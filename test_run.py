import pytest

from run import build_mol_from_modification


def test_build_mol_from_modification_recovers_expected_smiles():
    target_mapped_smiles = (
        "[Cl:0][c:1]1[cH:2][c:3]([Cl:4])[c:5](-[c:6]2[c:7]([Cl:8])[cH:9][cH:10][cH:11][c:12]2[Cl:13])[c:14]([Cl:15])[cH:16]1"
    )
    disconnect_list = [["Arylchloride", [[15]]]]
    connect_dict = {"C2 ([CH3:0][CH2:1])": [[1, 2, "target molecule"]]}

    result = build_mol_from_modification(target_mapped_smiles, disconnect_list, connect_dict)

    assert result == "CCc1c(Cl)ccc(-c2c(Cl)cccc2Cl)c1Cl"
