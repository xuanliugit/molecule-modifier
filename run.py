from rdkit import Chem
from accfg import (set_atom_idx,
                   remove_fg_list_from_mol)
import argparse
import sys

def attach_fg_to_mol(mol, fg_mol, in_atom, out_atom):
    """
    Attach a functional group to a molecule.

    Args:
        mol (rdkit.Chem.Mol): The target molecule.
        fg_mol (rdkit.Chem.Mol): The functional group molecule.
        in_atom (int): The atom index in the functional group to attach from.
        out_atom (int): The atom index in the target molecule to attach to.
    Returns:
        new_mol (rdkit.Chem.Mol): The new molecule with the functional group attached.
    """
    mol = Chem.RWMol(mol)
    fg_mol = Chem.RWMol(fg_mol)
    fg_atom = fg_mol.GetAtom(in_atom)
    fg_bonds = fg_atom.GetBonds()
    if len(fg_bonds) != 1:
        raise ValueError("The in_atom should have exactly one bond.")
    fg_bond = fg_bonds[0]
    fg_neighbor = fg_bond.GetOtherAtom(fg_atom)
    fg_neighbor_idx = fg_neighbor.GetIdx()
    
    # Remove the bond and the in_atom from the functional group
    fg_mol.RemoveBond(in_atom, fg_neighbor_idx)
    fg_mol.RemoveAtom(in_atom)
    
    # Add the functional group atoms to the target molecule
    atom_map = {}
    for atom in fg_mol.GetAtoms():
        new_idx = mol.AddAtom(atom)
        atom_map[atom.GetIdx()] = new_idx
    
    # Add bonds from the functional group to the target molecule
    for bond in fg_mol.GetBonds():
        begin_idx = atom_map[bond.GetBeginAtomIdx()]
        end_idx = atom_map[bond.GetEndAtomIdx()]
        mol.AddBond(begin_idx, end_idx, bond.GetBondType())
    
    # Add a bond between the functional group and the target molecule
    mol.AddBond(atom_map[fg_neighbor_idx], out_atom, Chem.BondType.SINGLE)
    
    new_mol = mol.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol

def clear_atom_prop(mol):
    for atom in mol.GetAtoms():
        prop_names = list(atom.GetPropNames())
        for prop in prop_names:
            atom.ClearProp(prop)
    return mol

def convert_fg_name_smiles_to_smiles(text):
    smiles = text.split('(')[-1][:-1]
    return smiles

def clean_radical_in_fg(fg_mol):
    for atom in fg_mol.GetAtoms():
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + atom.GetNumRadicalElectrons())
        atom.SetNumRadicalElectrons(0)
    return fg_mol

def build_mol_from_modification(target_mapped_smiles, disconnect_list, connect_dict):
    """
    Build a molecule from modification information.

    Args:
        target_mapped_smiles (str): The SMILES representation of the target molecule.
        disconnect_list (list): A list of disconnected fragments. e.g., [['Arylchloride', [[0]]]]
        connect_dict (dict): A dictionary mapping fragment connections. e.g., {'C2 ([CH3:0][CH2:1])': [[1, 0, 'target molecule']]}
    Returns:
        ref_smiles (str): The SMILES representation of the rebuilt molecule.
    """
    target_mol = Chem.MolFromSmiles(target_mapped_smiles)
    #target_mapped_smiles = Chem.MolToSmiles(target_mol)
    #target_mol = remove_atom_mapping(target_mol)
    # Chem.SanitizeMol(target_mol)
    target_mol = set_atom_idx(target_mol,'atomNote')
    new_mol = remove_fg_list_from_mol(target_mol, disconnect_list)
    for fg_name_smiles in connect_dict:
        fg_smiles = convert_fg_name_smiles_to_smiles(fg_name_smiles)
        fg_mol = Chem.MolFromSmiles(fg_smiles)
        fg_mol = clean_radical_in_fg(fg_mol)
        for connect_info in connect_dict[fg_name_smiles]:
            # Check source info to target molecule
            in_atom = connect_info[0]
            out_atom = connect_info[1]
            target_info = connect_info[2]
            if target_info != 'target molecule':
                raise ValueError("Currently only support connecting to target molecule.")
            else:
                for atom in fg_mol.GetAtoms():
                    if atom.GetProp('molAtomMapNumber') == str(in_atom):
                        atom.SetNumExplicitHs(atom.GetNumExplicitHs()-1)
                        atom.SetProp('atomNote', 'RecoverPoint')
                for atom in new_mol.GetAtoms():
                    if atom.GetProp('molAtomMapNumber') == str(out_atom):
                        atom.SetNumExplicitHs(atom.GetNumExplicitHs()-1)
                        atom.SetProp('atomNote', 'RecoverPoint')
                new_mol = Chem.CombineMols(fg_mol, new_mol)
                recover_point = []
                for atom in new_mol.GetAtoms():
                    if atom.HasProp('atomNote') and atom.GetProp('atomNote') == 'RecoverPoint':
                        recover_point.append(atom.GetIdx())
                        atom.SetProp('atomNote', atom.GetProp('molAtomMapNumber'))
                new_mol_edit = Chem.RWMol(new_mol)
                new_mol_edit.AddBond(recover_point[0], recover_point[1], order=Chem.rdchem.BondType.SINGLE)
                new_mol = new_mol_edit.GetMol()
    
    return Chem.MolToSmiles(clear_atom_prop(new_mol))

def arg_parse():
    
    parser = argparse.ArgumentParser(description="Build molecule from modification information.")
    parser.add_argument('--target_mapped_smiles', type=str, required=True,
                        help='The SMILES representation of the target molecule with atom mapping.')
    parser.add_argument('--disconnect_list', type=str, required=True,
                        help='A list of disconnected fragments in string format. e.g., \'[[\"Arylchloride\", [[0]]]]\'')
    parser.add_argument('--connect_dict', type=str, required=True,
                        help='A dictionary mapping fragment connections in string format. e.g., \'{"C2 ([CH3:0][CH2:1])": [[1, 2, \"target molecule\"]]}\'')
    args = parser.parse_args()
    return args

def main():
    args = arg_parse()
    target_mapped_smiles = args.target_mapped_smiles
    disconnect_list = eval(args.disconnect_list)
    connect_dict = eval(args.connect_dict)
    # print("Target mapped SMILES:", target_mapped_smiles)
    # print("Disconnect list:", disconnect_list)
    # print("Connect dict:", connect_dict)
    rebuilt_smiles = build_mol_from_modification(target_mapped_smiles, disconnect_list, connect_dict)
    print(rebuilt_smiles)
    sys.stdout.flush()
    
    
if __name__ == "__main__":
    '''
    Example usage:
    python run.py --target_mapped_smiles "[Cl:0][c:1]1[cH:2][c:3]([Cl:4])[c:5](-[c:6]2[c:7]([Cl:8])[cH:9][cH:10][cH:11][c:12]2[Cl:13])[c:14]([Cl:15])[cH:16]1" --disconnect_list '[["Arylchloride", [[15]]]]' --connect_dict '{"C2 ([CH3:0][CH2:1])": [[1, 2, "target molecule"]]}'
    '''
    main()