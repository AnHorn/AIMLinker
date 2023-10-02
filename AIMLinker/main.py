"""
Usage:
  main [options]

Options:
  -h, --help                    Show help messages
  -m, --min_atoms INT           Min atoms of linker
  -M, --max_atoms INT           Max atoms of linker
  -g, --generate_num INT        Number of generation per atom amounts
  -i, --input_file FILE Input   Mol file name
"""
import os
from collections import defaultdict

import pandas as pd
from docopt import docopt
from rdkit import Chem
from rdkit.Chem import PandasTools

import postprocess
from model import GGNNModel
from prepare_data import read_file, preprocess
from utils import compute_distance_and_angle


class FragmentLinking:
    def __init__(self, args):
        # Input parameters
        self.min_atoms = args["--min_atoms"]
        self.max_atoms = args["--max_atoms"]
        self.generate_num = args["--generate_num"]
        self.input_file = args["--input_file"]
        self.output_prefix = self.input_file.split(".")[0]

        # Basic settings
        os.environ["CUDA_VISIBLE_DEVICES"] = "0"

        # Get fragments from mol file
        self.mols_to_link = Chem.SDMolSupplier(self.input_file)
        for mol_to_link in self.mols_to_link:
            Chem.SanitizeMol(mol_to_link)

    def preprocess_data(self):
        data_path = "fragments_positioning_data.txt"
        if os.path.exists(data_path):
            os.remove(data_path)

        for mol_to_link in self.mols_to_link:
            # Get distance and angle between fragments
            dist, ang = compute_distance_and_angle(mol_to_link, "", Chem.MolToSmiles(mol_to_link))

            # Write data to file
            with open(data_path, 'a') as f:
                f.write("%s %s %s\n" % (Chem.MolToSmiles(mol_to_link), dist, ang))

            # Preprocess data
            raw_data = read_file(data_path)
            preprocess(raw_data, "zinc", "preprocess", True)

    def generate_linker(self):
        args = defaultdict(None)
        args["--dataset"] = "zinc"
        args["--config"] = '{"generation": true, \
                             "batch_size": 1, \
                             "number_of_generation_per_valid": ' + self.generate_num + ', \
                             "min_atoms": ' + str(self.min_atoms) + ', "max_atoms": ' + str(self.max_atoms) + ', \
                             "train_file": "molecules_preprocess.json", \
                             "valid_file": "molecules_preprocess.json", \
                             "output_name": "fragment_linking_result.smi"}'
        args["--freeze-graph-model"] = False
        args["--restore"] = "models/pretrained_model.pickle"

        # Setup model and generate molecules
        model = GGNNModel(args)
        model.train()

    def run(self):
        self.preprocess_data()
        self.generate_linker()


class PostProcess:
    def __init__(self, args):
        self.input_file = "fragment_linking_result.smi"
        self.output_prefix = args["--input_file"].split(".")[0]

        self.frag_smiles = ""
        self.gen_smiles_list = []

    def load_molecules(self):
        with open(self.input_file, "r") as f:
            for line in f:
                if line:
                    tmp = line.strip().split(" ")[0:3]
                    self.gen_smiles_list.append(tmp[2])

        self.frag_smiles = tmp[0]

    def filter(self):
        self.gen_smiles_list = postprocess.remove_duplicate(self.gen_smiles_list)
        self.gen_smiles_list = postprocess.remove_non_protac(self.gen_smiles_list, self.frag_smiles)
        self.gen_smiles_list = postprocess.remove_undruglike(self.gen_smiles_list)
        self.gen_smiles_list = postprocess.remove_disobey_bredt_rule(self.gen_smiles_list)

    def write_sdf_pass(self):
        df_linker = pd.DataFrame(columns=["Compound ID", "Smiles"])
        for i, smiles in enumerate(self.gen_smiles_list):
            # if self.pass_filter_list[i]:
            compound_name = "linker_%06d" % (i + 1)
            df_linker = df_linker.append({"Compound ID": compound_name, "Smiles": smiles}, ignore_index=True)

        PandasTools.AddMoleculeColumnToFrame(df_linker, "Smiles", "Molecule")
        PandasTools.WriteSDF(df_linker, f"{self.output_prefix}_result.sdf", molColName="Molecule",
                             properties=list(df_linker.columns))

    def run(self):
        self.load_molecules()
        self.filter()
        self.write_sdf_pass()


if __name__ == "__main__":
    arguments = docopt(__doc__)

    linker_prediction = FragmentLinking(arguments)
    linker_prediction.run()

    post_process = PostProcess(arguments)
    post_process.run()
