import json
from utils import to_graph_mol, align_mol_to_frags

dataset = 'zinc'


def read_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    num_lines = len(lines)
    data = []
    for i, line in enumerate(lines):
        toks = line.strip().split(' ')
        if len(toks) == 3:
            smi_frags, abs_dist, angle = toks
            smi_mol = smi_frags
            smi_linker = ''
        elif len(toks) == 5:
            smi_mol, smi_linker, smi_frags, abs_dist, angle = toks
        else:
            print("Incorrect input format. Please check the README for useage.")
            exit()
        data.append({'smi_mol': smi_mol, 'smi_linker': smi_linker,
                     'smi_frags': smi_frags,
                     'abs_dist': [abs_dist, angle]})
        if i % 2000 == 0:
            print('Finished reading: %d / %d' % (i, num_lines), end='\r')
    print('Finished reading: %d / %d' % (num_lines, num_lines))
    return data


def preprocess(raw_data, dataset, name, test=False):
    print('Parsing smiles as graphs.')
    processed_data = []
    total = len(raw_data)
    for i, (smi_mol, smi_frags, smi_link, abs_dist) in enumerate([(mol['smi_mol'], mol['smi_frags'],
                                                                   mol['smi_linker'], mol['abs_dist']) for mol in
                                                                  raw_data]):
        if test:
            smi_mol = smi_frags
            smi_link = ''
        (mol_out, mol_in), nodes_to_keep, exit_points = align_mol_to_frags(smi_mol, smi_link, smi_frags)
        if mol_out == []:
            continue
        nodes_in, edges_in = to_graph_mol(mol_in, dataset)
        nodes_out, edges_out = to_graph_mol(mol_out, dataset)
        if min(len(edges_in), len(edges_out)) <= 0:
            continue
        processed_data.append({
            'graph_in': edges_in,
            'graph_out': edges_out,
            'node_features_in': nodes_in,
            'node_features_out': nodes_out,
            'smiles_out': smi_mol,
            'smiles_in': smi_frags,
            'v_to_keep': nodes_to_keep,
            'exit_points': exit_points,
            'abs_dist': abs_dist
        })
        if i % 500 == 0:
            print('Processed: %d / %d' % (i, total), end='\r')
    print('Processed: %d / %d' % (total, total))
    print('Saving data')
    with open('molecules_%s.json' % name, 'w') as f:
        json.dump(processed_data, f)
    print('Length raw data: \t%d' % total)
    print('Length processed data: \t%d' % len(processed_data))
