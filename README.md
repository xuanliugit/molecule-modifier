# Molecule Modifier

Build a molecule from modification information.

    Args:
        target_mapped_smiles (str): The SMILES representation of the target molecule.
        disconnect_list (list): A list of disconnected fragments. e.g., [['Arylchloride', [[0]]]]
        connect_dict (dict): A dictionary mapping fragment connections. e.g., {'C2 ([CH3:0][CH2:1])': [[1, 0, 'target molecule']]}
    Returns:
        ref_smiles (str): The SMILES representation of the rebuilt molecule.

```python
python run.py --target_mapped_smiles "[Cl:0][c:1]1[cH:2][c:3]([Cl:4])[c:5](-[c:6]2[c:7]([Cl:8])[cH:9][cH:10][cH:11][c:12]2[Cl:13])[c:14]([Cl:15])[cH:16]1" --disconnect_list '[["Arylchloride", [[15]]]]' --connect_dict '{"C2 ([CH3:0][CH2:1])": [[1, 2, "target molecule"]]}'
```