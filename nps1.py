#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import Descriptors
from collections import Counter

systems_map_I = {
    "O1C(CCN)=CC2=CC3=C(C=CO3)C=C12": "benzodifuranyl",
    "O1C=C(CCN)C2=CC3=C(C=CO3)C=C12": "benzodifuranyl",
    "O1C=CC2=C(CCN)C3=C(C=CO3)C=C12": "benzodifuranyl",
    "O1C(CCN)=CC2=CC3=C(OC=C3)C=C12": "benzodifuranyl",
    "O1C=C(CCN)C2=CC3=C(OC=C3)C=C12": "benzodifuranyl",
    "C1CC2=C(CCN)C3=C(CCO3)C=C2O1": "tetrahydrobenzofuranyl",
    "C1CC2=C(CCN)C3=C(OCC3)C=C2O1": "tetrahydrobenzofuranyl",
    "C1CC2=CC3=C(OCC3)C(CCN)=C2O1": "tetrahydrobenzofuranyl",
    "C1COC2=C(CCN)C3=C(OCCC3)C=C2C1": "tetrahydrobenzodipiranyl",
    "C1COC2=CC3=C(CCCO3)C(CCN)=C2C1": "tetrahydrobenzodipiranyl",
    "O1C(CCN)=CC2=C1C=CC=C2": "benzofuranyl",
    "O1C=C(CCN)C2=C1C=CC=C2": "benzofuranyl",
    "O1C=CC2=C1C(CCN)=CC=C2": "benzofuranyl",
    "O1C=CC2=C1C=CC(CCN)=C2": "benzofuranyl",
    "O1C=CC2=C1C=CC=C(CCN)2": "benzofuranyl",
    "O1C=CC2=C1C=C(CCN)C=C2": "benzofuranyl",
    "C1CC2=C(O1)C(CCN)=CC=C2": "dihydrobenzofuranyl",
    "C1CC2=C(O1)C=C(CCN)C=C2": "dihydrobenzofuranyl",
    "C1CC2=C(O1)C=CC(CCN)=C2": "dihydrobenzofuranyl",
    "C1CC2=C(O1)C=CC=C2(CCN)": "dihydrobenzofuranyl",
    "C1C=CC2=C1C=CC=C2CCN": "indenyl",
    "C1C=CC2=C1C=C(CCN)C=C2": "indenyl",
    "C1C(CCN)=CC2=C1C=CC=C2": "indenyl",
    "C1CC2=C(C1)C=CC=C2CCN": "indanyl",
    "C1CC2=C(C1)C=C(CCN)C=C2": "indanyl",
    "C1OC2=C(O1)C(CCN)=CC=C2": "metylodioksyfenyl",
    "C1OC2=C(O1)C=C(CCN)C=C2": "metylodioksyfenyl",
    "C1COC2=C(CCN)C=CC=C2O1": "etylodioksyfenyl",
    "C1COC2=CC(CCN)=CC=C2O1": "etylodioksyfenyl",
    "C1=C(CCN)C2=CC=CC=C2C=C1": "naftyl",
    "C1=CC2=CC(CCN)=CC=C2C=C1": "naftyl",
    "C1CCC2=C(CCN)C=CC=C2C1": "tetralinyl",
    "C1CCC2=CC(CCN)=CC=C2C1": "tetralinyl",
    "N1C=CC=C1CCN": "pirolil",
    "N1C=C(CCN)C=C1": "pirolil",
    "S1C(CCN)=CC=C1": "tiofuranyl",
    "S1C=CC(CCN)=C1": "tiofuranyl",
    "C1=CC(CCN)=NC=C1": "pirydyl",
    "C1=CC=NC=C1CCN": "pirydyl",
    "NCCC1=CC=NC=C1": "pirydyl",
    "O1C=CC=C1(CCN)": "furyl",
    "O1C=C(CCN)C=C1": "furyl",
    "C1CCCC1(CCN)": "cyklopentyl",
    "C1CCCCC1(CCN)": "cykloheksyl",
    "C1=CC=CC(CCN)=C1": "fenyl"
}


def find_substructure(systems_map_I: dict, mol: Chem.rdchem.Mol) -> bool | Chem.rdchem.Mol | tuple:
    substructure = None
    for system, name in systems_map_I.items():
        if mol.HasSubstructMatch(Chem.MolFromSmiles(system)):
            substructure = Chem.MolFromSmiles(system)
            matches = mol.GetSubstructMatches(substructure)
            return substructure, matches, name

    if substructure is None:
        return False


def get_central_atoms(mol: Chem.rdchem.Mol, suspected: tuple) -> list | bool:
    res = []
    rings = mol.GetRingInfo()  # information about the all rings in mol object
    atom_rings = rings.AtomRings()
    for idx in suspected:
        atom = mol.GetAtomWithIdx(idx)  # get each atom from main structure

        neighbors = atom.GetNeighbors()

        if not all(n.GetIdx() in suspected for n in neighbors): 

            ring = [list(idxs) for idxs in atom_rings if idx in idxs]  # list of idx of ring which the atom is part of
            is_suspected_part = any(set(lst).issubset(set(list(suspected))) for lst in ring)
            add = [idx, atom.GetIsAromatic(), atom.GetSymbol(), atom.IsInRing(), is_suspected_part]
            # aliphatic nitrogen
            if add[2] == "N" and add[3] is False and len(neighbors) == 3: res.extend([add, add])
            if add[2] == "N" and add[3] is False and len(neighbors) == 2: res.append(add)
            if add[2] == "N" and add[3] is True and len(neighbors) == 3:
                for s_idx in suspected:
                    atom = mol.GetAtomWithIdx(s_idx)
                    # if atom belongs to suspected cyclic system
                    if any(set(lst).issubset(set(list(suspected))) for lst in
                           [list(idxs) for idxs in atom_rings if s_idx in idxs]) and atom.GetSymbol() == "C":
                        if all([n.IsInRing() for n in atom.GetNeighbors()]):
                            pass  # nitrogen atom is in cyclic system condensed with suspected cyclic system
                    elif add[4] is False:
                        res.append(add) 

            # system carbon
            if add[2] == "C" and add[3] is True and add[4] == True and len(neighbors) == 3:
                if len(ring) == 1: res.append(add)
                if len(ring) > 1:
                    for s_idx in suspected:
                        atom = mol.GetAtomWithIdx(s_idx)
                        if atom.GetSymbol() == "N" and atom.IsInRing() is True:
                            pass  # nitrogen atom is in cyclic system condensed with system ring
                        if any(set(lst).issubset(set(list(suspected))) for lst in
                               [list(idxs) for idxs in atom_rings if s_idx in idxs]) and atom.GetSymbol() == "C":
                            if all([n.IsInRing() for n in atom.GetNeighbors()]):
                                pass

            # aliphatic carbon
            if add[2] == "C" and add[3] is False:

                if len(neighbors) == 3:
                    res.append(add)
                elif len(neighbors) == 4:
                    res.extend([add, add])

            if add[2] == "C" and add[3] is True and add[4] is False:
                add.append(True)  # True for being condensed
                res.append(add)  # aliphatic carbon condensed with

    return res


def if_in_ring(s: Chem.rdchem.Mol, in_ring: bool) -> tuple:
    is_in_ring = []
    for i in range(s.GetNumAtoms()):
        atom = s.GetAtomWithIdx(i)
        is_in_ring.append([i, atom.IsInRing()])

    in_ring_info = []
    i = 0
    while i < len(is_in_ring):
        if is_in_ring[i][1] is in_ring:
            in_ring_info.append(is_in_ring[i])
        elif is_in_ring[i][1] is not in_ring:
            break
        i += 1

    return in_ring_info, is_in_ring


def until(s, stop_a, in_add):
    s2edit = Chem.RWMol(s)
    at = min([atom.GetIdx() for atom in s.GetAtoms() if atom.GetSymbol() == stop_a]) + in_add
    if at + 1 < s.GetNumHeavyAtoms():
        for atom in range(0, at + 1):
            if atom == at:
                s2edit.RemoveAtom(atom)
                break
    else:
        return s

    return s2edit


# rs(s, num_heavy_atoms, carbons, permitted_atoms)
def rs(s: Chem.rdchem.Mol, num_heavy_atoms: int, carbons: int, permitted_atoms: list) -> bool:
    if num_heavy_atoms <= 8:
        if Fragments.fr_halogen(s) == 1 and num_heavy_atoms == 1:
            return True
        if all(atom.GetAtomicNum() == 6 for atom in s.GetAtoms()) and carbons <= 6:
            return True
        if 2 <= num_heavy_atoms <= 6 and sum(1 for atom in s.GetAtoms() if atom.GetAtomicNum() == 8) == 1:
            return True
        if Fragments.fr_COO(s) == 1:
            return True
        if Fragments.fr_sulfone(s) == 1 and carbons <= 6:
            return True
        if Fragments.fr_nitro(s) == 1:  # cannot find nitro group directly attach to benzene ring
            return True
        elif 6 < num_heavy_atoms <= 8 and all(atom.GetAtomicNum() in permitted_atoms for atom in s.GetAtoms()):
            return True
        else:
            return False
    else:
        return False  # too many atoms


# r12(s, part_s, ring_part_s, is_in_ring, carbons, num_heavy_atoms, permitted_atoms)
def r12(s: Chem.rdchem.Mol, part_s: Chem.rdchem.Mol, ring_part_s: Chem.rdchem.Mol, is_in_ring: int, carbons: int, num_heavy_atoms: int, permitted_atoms: list) -> bool:
    part_s_carbons = [atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()].count(True)
    if num_heavy_atoms - is_in_ring <= 10:
        if ring_part_s is not None:
            if all(atom.GetAtomicNum() in [6, 8, 16, 7] for atom in
                   ring_part_s.GetAtoms()) and 5 <= ring_part_s.GetNumHeavyAtoms() <= 7:
                return True
        else:
            if all(atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()) and part_s_carbons <= 6:
                return True
            if sum(1 for atom in part_s.GetAtoms() if atom.GetAtomicNum() == 8) == 1 and part_s_carbons <= 6:
                return True
            if Fragments.fr_Al_OH(part_s) == 1:
                return True
            if Fragments.fr_NH2(part_s) == 1 or Fragments.fr_NH1(part_s) or Fragments.fr_NH0(part_s) == 1:
                return True

            elif 6 < num_heavy_atoms <= 10 and all(atom.GetAtomicNum() in permitted_atoms for atom in s.GetAtoms()):
                return True
            else:
                return False
    else:
        return False  # to many atoms


# r3456(s, part_s, ring_part_s, condense, is_in_ring, carbons, num_heavy_atoms, permitted_atoms)
def r3456(s: Chem.rdchem.Mol, part_s: Chem.rdchem.Mol, ring_part_s: Chem.rdchem.Mol, condense: bool, is_in_ring: int,
          carbons: int, num_heavy_atoms: int, permitted_atoms: list) -> bool:
    part_s_carbons = [atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()].count(True)
    if num_heavy_atoms - is_in_ring <= 10:
        if ring_part_s is not None:

            if all(atom.GetAtomicNum() == 6 for atom in ring_part_s.GetAtoms()) and part_s_carbons <= 10:
                return True
        else:

            if condense:
                if all(atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()) and 4 <= part_s_carbons <= 6:
                    return True
            if Fragments.fr_halogen(s) == 1 and num_heavy_atoms == 1:
                return True
            if all(atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()) and part_s_carbons <= 10:
                return True
            if s.HasSubstructMatch(Chem.MolFromSmiles("O")):
                # if Fragments.fr_Al_OH(part_s) == 1:
                tocheck = until(s, "O", 0)
                if Fragments.fr_Al_OH(tocheck) == 1:
                    return True
                # if part_s.HasSubstructMatch(Chem.MolFromSmiles("COC")):
                tocheck = until(s, "O", 1)
                # print(Chem.MolToSmiles(tocheck))
                if tocheck.HasSubstructMatch(Chem.MolFromSmiles("COC")) and [atom.GetAtomicNum() == 6 for atom in tocheck.GetAtoms()].count(True) <= 10:
                    return True
                tocheck = until(s, "O", 1)
                if Fragments.fr_C_O_noCOO(tocheck) == 1 and [atom.GetAtomicNum() == 6 for atom in tocheck.GetAtoms()].count(True) <= 10:
                    return True
            if s.HasSubstructMatch(Chem.MolFromSmiles("S")):
            # if Fragments.fr_sulfone(part_s) == 1 :
                tocheck = until(s, "S", 1)
                if Fragments.fr_sulfone(tocheck) == 1 and [atom.GetAtomicNum() == 6 for atom in tocheck.GetAtoms()].count(True) <= 10:
                    return True
            # if Fragments.fr_C_O_noCOO(part_s) == 1:  # Number of carbonyl O, excluding COOH

            elif 6 < num_heavy_atoms <= 10 and all(atom.GetAtomicNum() in permitted_atoms for atom in s.GetAtoms()):
                return True
            else:
                return False
    else:
        return False


def main(smiles):

    smiles = max(smiles.split("."), key=len)  # remove the radicals
    mol = Chem.MolFromSmiles(smiles)

    if round(Descriptors.ExactMolWt(mol), 2) <= 500:  

        if find_substructure(systems_map_I, mol) is not None:
            substructure, matches, name = find_substructure(systems_map_I, mol)
            for suspected in matches:
                res = get_central_atoms(mol, suspected)
                if res is not False:

                    mol2edit = Chem.RWMol(mol)
                    for atom in reversed(sorted(suspected)):
                        mol2edit.RemoveAtom(atom)
                    s4condense = Chem.MolToSmiles(mol2edit, canonical=True).split(".")

                    substituents = []
                    mol2substituents = Chem.RWMol(mol)
                    for c_idx in res:
                        mol2substituents.RemoveAtom(c_idx[0])
                        substituent = Chem.MolToSmiles(mol2substituents, canonical=True).split(".")[0]

                        substituents.append(substituent)
                        mol2substituents = Chem.RWMol(mol)

                    s4condense = list((Counter(s4condense) - Counter(substituents)).elements())  # overwrite
                    result = []
                    if len(res) != len(substituents):  # validation for s order would be beneficial
                        print("Do weryfikacji")
                        return False
                    for r, s in zip(res, substituents):

                        if len(r) == 6:  # if condense is True
                            s = s4condense[0]
                        permitted_atoms = [6, 8, 7, 9, 17, 35, 53, 16]
                        s = Chem.MolFromSmiles(s)

                        non_ring_s, is_in_ring = if_in_ring(s, False)  # list of [bool, atom_idx]
                        ring_s, _ = if_in_ring(s, True)
                        is_in_ring = [i[1] for i in is_in_ring if i[1] is True].count(True)

                        carbons = [atom.GetAtomicNum() == 6 for atom in s.GetAtoms()].count(True)
                        num_heavy_atoms = s.GetNumHeavyAtoms()
                        atoms = s.GetAtoms()
                        part_s = None
                        ring_part_s = None

                        if len(non_ring_s) > 0:  # s starts with aliphatic atom
                            idx_to_remove = [i for i in range(len(non_ring_s), num_heavy_atoms)]
                            part_s = Chem.RWMol(s)
                            for atom in reversed(sorted(idx_to_remove)):
                                part_s.RemoveAtom(atom)

                        else:
                            idx_to_remove = [i for i in range(len(ring_s), num_heavy_atoms)]
                            ring_part_s = Chem.RWMol(s)
                            for atom in reversed(sorted(idx_to_remove)):
                                ring_part_s.RemoveAtom(atom)

                        if r[2] == "N":
                            result.append(r12(s, part_s, ring_part_s, is_in_ring, carbons, num_heavy_atoms, permitted_atoms))

                        if r[2] == "C" and r[4] is True:
                            result.append(rs(s, num_heavy_atoms, carbons, permitted_atoms))

                        if r[2] == "C" and r[4] is False:
                            condense = False
                            if len(r) == 6:  # if condense is True
                                condense = True
                            result.append(r3456(s, part_s, ring_part_s, condense, is_in_ring, carbons, num_heavy_atoms, permitted_atoms))

                    return all(i for i in result)

        else:
            return False  # mw above 500


if __name__ == '__main__':
  main()
