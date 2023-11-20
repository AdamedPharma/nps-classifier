#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Fragments


smiles = input("smiles: ")
# remove the radicals
smiles = max(smiles.split("."), key=len)
mol = Chem.MolFromSmiles(smiles)


def find_substructure(substructures: list):
    substructure = None
    for i, substructure in enumerate(substructures):

        if mol.HasSubstructMatch(substructure):
            matches = mol.GetSubstructMatches(substructure)
            return substructure, matches, i
    if substructure is None:
        return False


def main():
    structure_I = Chem.MolFromSmiles("NCCc1c[nH]c2ccccc12")
    structure_I_pyrrolidine = Chem.MolFromSmiles("C1CCCN1CCc1c[nH]c2ccccc12")
    structure_I_dihydrogen_phosphate = Chem.MolFromSmiles("NCCc1c[nH]c2cccc(OP(=O)(O)[O-])c12")
    structure_II = Chem.MolFromSmiles("[H][C@@]2(C(N)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34")
    structure_II_morpholine = Chem.MolFromSmiles("[H][C@@]2(C(N1CCOCC1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34")
    structure_II_pyrrolidine = Chem.MolFromSmiles("[H][C@@]2(C(N1CCCC1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34")
    structure_II_azetidine_I = Chem.MolFromSmiles("[H][C@@]2(C(N1C(C)C(C)C1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34")
    structure_II_azetidine_II = Chem.MolFromSmiles("[H][C@@]2(C(N1C(C)CC(C)1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34")
    substructures = [structure_II_azetidine_I, structure_II_azetidine_II, structure_II_pyrrolidine,
                     structure_II_morpholine, structure_II, structure_I_pyrrolidine, structure_I_dihydrogen_phosphate,
                     structure_I]
    condensed_systems = [structure_II_azetidine_I, structure_II_azetidine_II, structure_II_pyrrolidine,
                         structure_II_morpholine]

    if find_substructure(substructures) is not None:
        substructure, matches, i = find_substructure(substructures)
        # if substructure:
        mol2edit = Chem.RWMol(mol)
        if len(matches) == 1:
            for atom in reversed(sorted(matches[0])):
                mol2edit.RemoveAtom(atom)

        substituents = Chem.MolToSmiles(mol2edit, canonical=False).split(".")
        res = []
        ri = mol.GetRingInfo()
        for idx in matches[0]:
            atom = mol.GetAtomWithIdx(idx)  # get each atom from main structure
            neighbors = atom.GetNeighbors()
            if not all(n.GetIdx() in matches[0] for n in neighbors):
                add = [idx, atom.GetIsAromatic(), atom.GetSymbol(), ri.IsAtomInRingOfSize(atom.GetIdx(), 6)]
                if i > 4:
                    # nitrogen aliphatic atom can have one or two substituents
                    if len(atom.GetNeighbors()) == 2 and add[2] == "N" and add[1] is False: res.append(add)
                    if len(atom.GetNeighbors()) == 3 and add[2] == "N" and add[
                        1] is False and substructure != structure_I_pyrrolidine: res.extend([add, add])
                    # nitrogen aromatic atom can have only one substituent
                    if len(atom.GetNeighbors()) == 3 and add[2] == "N" and add[1] == True: res.append(add)

                    if not all(idx in matches[0] for idx in [atom.GetIdx() for atom in
                                                             atom.GetNeighbors()]):  # discard atoms from structure without substituents
                        # carbon aliphatic atom can have only one substituent
                        if len(atom.GetNeighbors()) == 3 and add[2] == "C" and add[1] is False: res.append(add)
                        # carbon aromatic atom can have only one substituent
                        if len(atom.GetNeighbors()) == 3 and add[2] == "C" and add[1] == True:
                            res.append(add)
                        # carbon atom from benzene ring can have multiple substituents
                        elif len(atom.GetNeighbors()) == 3 and add[2] == "C" and add[1] == True and add[3] == True:
                            res.append(add)
                else:  # II-structure
                    # nitrogen benzene atom can have one substituents
                    if len(atom.GetNeighbors()) == 3 and add[2] == "N" and ri.IsAtomInRingOfSize(atom.GetIdx(),
                                                                                                 5): res.append(add)
                    # nitrogen aromatic 5-ring atom can have one substituent
                    if len(atom.GetNeighbors()) == 3 and add[2] == "N" and add[1] is False and add[
                        3] is True: res.append(add)
                    # nitrogen aliphatic atom can have one or two substituents
                    if len(atom.GetNeighbors()) == 3 and add[2] == "N" and add[1] is False and add[
                        3] is False and substructure not in condensed_systems: res.extend([add, add])
                    if len(atom.GetNeighbors()) == 2 and add[2] == "N" and add[1] is False: res.append(add)

        def condensed(res: list):
            # if methylodioxy group is present and two aromatic carbons have substituent then one of them has to be removed
            to_remove = [None]
            for each in res:
                if each[3] is True and each[2] == "C":
                    to_remove = each
            res.remove(to_remove)
            return res

        if mol.HasSubstructMatch(Chem.MolFromSmiles("c2ccc1OCOc1c2")): res = condensed(res)

        def r12(s: Chem.rdchem.Mol) -> bool:
            # conditions for aliphatic nitrogen substituents
            if all(atom.GetAtomicNum() == 6 for atom in s.GetAtoms()) and len(
                    carbons) <= 6:  # alkyl group up to 6 carbons (allyl includeed)
                return True
            if mol.HasSubstructMatch(structure_I_pyrrolidine):  # check the pyrrolidine
                return True
            else:
                return False

        def r345(s: Chem.rdchem.Mol, cs: int) -> bool:
            # conditions for aromatic nitrogen and aliphatic carbon and aromatic carbon substituents
            if all(atom.GetAtomicNum() == 6 for atom in s.GetAtoms()):  # check if the substituent contains only carbon atoms
                if carbons <= cs:
                    return True
                else:
                    return False  # more than 3 carbons atoms
            else:
                return False  # atoms other than carbons

        def rn(s: Chem.rdchem.Mol) -> bool:

            if mol.HasSubstructMatch(
                    structure_I_dihydrogen_phosphate):  # check the dihydrogen phosphate, while True i cannot have another S
                return True
            if mol.HasSubstructMatch(
                    Chem.MolFromSmiles("c2ccc1OCOc1c2")) and num_heavy_atoms == 9:  # check methylodioxy group
                return True
            if s.HasSubstructMatch(
                    Chem.MolFromSmiles("OC(=O)C")) and num_heavy_atoms == 4:  # check the acethoxyl group
                return True
            if s.HasSubstructMatch(Chem.MolFromSmiles("OC")) and Fragments.fr_Al_OH(
                    s) == 0 and num_heavy_atoms == 2:  # check the methoxyl group
                return True
            if s.HasSubstructMatch(Chem.MolFromSmiles("O")) and num_heavy_atoms == 1:  # check the hydroxyl group
                return True
            if s.HasSubstructMatch(
                    Chem.MolFromSmiles("SC")) and num_heavy_atoms == 2:  # check the methylotiol group
                return True
            else:
                return False

        def r1_II(s: Chem.rdchem.Mol) -> bool:
            # conditions for 5-ring nitrogen
            if all(atom.GetAtomicNum() == 6 for atom in s.GetAtoms()) and carbons <= 3:  # check if the substituent contains only carbon atoms
                return True
            if num_heavy_atoms <= 4 and s.HasSubstructMatch(Chem.MolFromSmiles("COC")):
                return True
            else:
                return False

        def r2_II(s: Chem.rdchem.Mol) -> bool:
            # conditions for 6-ring nitrogen
            if all(atom.GetAtomicNum() == 6 for atom in s.GetAtoms()) and carbons <= 4:
                return True

        def r34_II(s: Chem.rdchem.Mol) -> bool:
            # conditions for aliphatic nitrogen
            if all(atom.GetAtomicNum() == 6 for atom in s.GetAtoms()) and carbons <= 5:
                return True
            if Fragments.fr_Al_OH(s) == 1 and num_heavy_atoms == 2:
                return True
            for system in condensed_systems:
                if s.HasSubstructMatch(system) and carbons <= 4:
                    return True
            else:
                return False

        result = []
        res.sort()

        for r, s in zip(res, substituents):
            num_heavy_atoms = s.GetNumHeavyAtoms()
            s = Chem.MolFromSmiles(s)
            carbons = [atom.GetAtomicNum() == 6 for atom in
                       s.GetAtoms()].count(True)
            if i > 4:
                if r[1] is False and r[2] == "N":
                    result.insert(0, r12(s))  # substituent conditions for aliphatic nitrogen
                if (r[1] is True and r[2] == "N") or (r[1] is False and r[2] == "C"):
                    result.insert(1, r345(s,
                                          3))  # substituent conditions for aromatic nitrogen and aliphatic carbon
                if r[1] is True and r[2] == "C" and r[3] is False:
                    result.insert(2, r345(s, 2))  # substituent conditions for aromatic carbon
                if r[1] is True and r[2] == "C" and r[3] is True:  # attach to benzene ring
                    result.insert(3, rn(s))  # attach to benzene ring
            else:
                if r[2] == "N":
                    result.insert(0, r1_II(s))  # substituent conditions for 5-ring nitrogen
                if r[1] is False and r[3] is True and r[2] == "N":
                    result.insert(1, r2_II(s))  # substituent conditions for 6-ring nitrogen
                if r[1] is False and r[3] is False and r[2] == "N":
                    result.insert(2, r34_II(s))  # substituent conditions for aliphatic nitrogen
        else:
            return all([val for val in result])
    # return False
    else:
        return False  # main structure was not found


if __name__ == '__main__':
    main()

