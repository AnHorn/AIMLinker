from rdkit import Chem
from rdkit.Chem import AllChem


def remove_duplicate(smiles_list):
    """
    Remove duplicates.
    """
    keep_set = set()
    for smiles in smiles_list:
        keep_set.add(smiles)

    return list(keep_set)


def remove_non_protac(smiles_list, fragment_smiles):
    """
    Remove non-PROTAC structure.
    
    fragment_smiles: String type, SMILES of two fragments, e.g. 
    "NC(=O)CC1N=C(c2ccc(Cl)cc2)c2c(sc(C)c2C)-n2c(C)nnc21.Oc1cccc2c1C(=O)N(C1CCC(=O)NC1=O)C2=O"
    """
    fragment_mol = Chem.MolFromSmiles(fragment_smiles)
    du = Chem.MolFromSmiles("*")
    hydrogen = Chem.MolFromSmiles("[H]")

    keep_list = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        clean_frag = Chem.RemoveHs(AllChem.ReplaceSubstructs(fragment_mol, du, hydrogen, True)[0])
        if mol.HasSubstructMatch(clean_frag):
            keep_list.append(smiles)

    return keep_list


def remove_undruglike(smiles_list):
    """
    Remove structures containing unfavorible substructure.
    """
    keep_list = []
    for smiles in smiles_list:
        if is_druglike(smiles):
            keep_list.append(smiles)

    return keep_list


def remove_disobey_bredt_rule(smiles_list):
    """
    Remove structures disobeying Bredt's Rule.
    """
    keep_list = []
    for smiles in smiles_list:
        if is_obey_bredt_rule(smiles):
            keep_list.append(smiles)

    return keep_list


def is_druglike(smiles):
    """
    Return true if the structure does not 
    contain unfavorible substructure. 
    """
    mol = Chem.MolFromSmiles(smiles)

    exclusion_mol_list = []
    with open("unfavorible_substructure.smi", "r") as f:
        lines = f.readlines()
    for smiles in lines:
        exclusion_mol_list.append(Chem.MolFromSmiles(smiles))

    for exclusion_mol in exclusion_mol_list:
        if mol.HasSubstructMatch(exclusion_mol):
            return False

    return True


def is_obey_bredt_rule(smiles):
    """
    Return true if the structure obeys Bredt's Rule. 
    Bredt's Rule describes that double bonds cannot be 
    placed at the bridgehead of a bridged ring system.
    """
    mol = Chem.MolFromSmiles(smiles)

    ssr = Chem.GetSymmSSSR(mol)
    for i in range(len(ssr)):
        for j in range(i + 1, len(ssr)):
            ssr1 = set(ssr[i])
            ssr2 = set(ssr[j])
            ssr1_ssr2 = ssr1.intersection(ssr2)
            # Two small rings form a bridge
            if len(ssr1_ssr2) > 2:
                bridge = ssr1_ssr2
                for atom_id in bridge:
                    neighbors = {x.GetIdx() for x in mol.GetAtomWithIdx(atom_id).GetNeighbors()}
                    neighbor_bridge = neighbors.intersection(bridge)
                    neighbor_ssr1 = neighbors.intersection(ssr1)
                    neighbor_ssr2 = neighbors.intersection(ssr2)

                    # Atom is bridge head
                    if len(neighbor_bridge) == 1 and len(neighbor_ssr1) == 2 and len(neighbor_ssr2) == 2:
                        for bond_id in mol.GetAtomWithIdx(atom_id).GetBonds():
                            # Bonds is double/aromatic bond
                            if bond_id.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.AROMATIC]:
                                ssr_union = ssr1.union(ssr2)
                                # The double/aromatic bond in both rings (ssr1 and ssr2)
                                if bond_id.GetBeginAtomIdx() in ssr_union and bond_id.GetEndAtomIdx() in ssr_union:
                                    return False

    return True
