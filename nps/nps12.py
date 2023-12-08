from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import Descriptors
from collections import Counter
from typing import Union

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
systems_map_II = {
    "O1C(C(=O)CN)=CC2=CC3=C(C=CO3)C=C12": "benzodifuranyl",
    "O1C=C(C(=O)CN)C2=CC3=C(C=CO3)C=C12": "benzodifuranyl",
    "O1C=CC2=C(C(=O)CN)C3=C(C=CO3)C=C12": "benzodifuranyl",
    "O1C(C(=O)CN)=CC2=CC3=C(OC=C3)C=C12": "benzodifuranyl",
    "O1C=C(C(=O)CN)C2=CC3=C(OC=C3)C=C12": "benzodifuranyl",
    "C1CC2=C(C(=O)CN)C3=C(CCO3)C=C2O1": "tetrahydrobenzofuranyl",
    "C1CC2=C(C(=O)CN)C3=C(OCC3)C=C2O1": "tetrahydrobenzofuranyl",
    "C1CC2=CC3=C(OCC3)C(C(=O)CN)=C2O1": "tetrahydrobenzofuranyl",
    "C1COC2=C(C(=O)CN)C3=C(OCCC3)C=C2C1": "tetrahydrobenzodipiranyl",
    "C1COC2=CC3=C(CCCO3)C(C(=O)CN)=C2C1": "tetrahydrobenzodipiranyl",
    "O1C(C(=O)CN)=CC2=C1C=CC=C2": "benzofuranyl",
    "O1C=C(C(=O)CN)C2=C1C=CC=C2": "benzofuranyl",
    "O1C=CC2=C1C(C(=O)CN)=CC=C2": "benzofuranyl",
    "O1C=CC2=C1C=CC(C(=O)CN)=C2": "benzofuranyl",
    "O1C=CC2=C1C=CC=C(C(=O)CN)2": "benzofuranyl",
    "O1C=CC2=C1C=C(C(=O)CN)C=C2": "benzofuranyl",
    "C1CC2=C(O1)C(C(=O)CN)=CC=C2": "dihydrobenzofuranyl",
    "C1CC2=C(O1)C=C(C(=O)CN)C=C2": "dihydrobenzofuranyl",
    "C1CC2=C(O1)C=CC(C(=O)CN)=C2": "dihydrobenzofuranyl",
    "C1CC2=C(O1)C=CC=C2(C(=O)CN)": "dihydrobenzofuranyl",
    "C1C=CC2=C1C=CC=C2C(=O)CN": "indenyl",
    "C1C=CC2=C1C=C(C(=O)CN)C=C2": "indenyl",
    "C1C(C(=O)CN)=CC2=C1C=CC=C2": "indenyl",
    "C1CC2=C(C1)C=CC=C2C(=O)CN": "indanyl",
    "C1CC2=C(C1)C=C(C(=O)CN)C=C2": "indanyl",
    "C1OC2=C(O1)C(C(=O)CN)=CC=C2": "metylodioksyfenyl",
    "C1OC2=C(O1)C=C(C(=O)CN)C=C2": "metylodioksyfenyl",
    "C1COC2=C(C(=O)CN)C=CC=C2O1": "etylodioksyfenyl",
    "C1COC2=CC(C(=O)CN)=CC=C2O1": "etylodioksyfenyl",
    "C1=C(C(=O)CN)C2=CC=CC=C2C=C1": "naftyl",
    "C1=CC2=CC(C(=O)CN)=CC=C2C=C1": "naftyl",
    "C1CCC2=C(C(=O)CN)C=CC=C2C1": "tetralinyl",
    "C1CCC2=CC(C(=O)CN)=CC=C2C1": "tetralinyl",
    "N1C=CC=C1C(=O)CN": "pirolil",
    "N1C=C(C(=O)CN)C=C1": "pirolil",
    "S1C(C(=O)CN)=CC=C1": "tiofuranyl",
    "S1C=CC(C(=O)CN)=C1": "tiofuranyl",
    "C1=CC(C(=O)CN)=NC=C1": "pirydyl",
    "C1=CC=NC=C1C(=O)CN": "pirydyl",
    "NC(=O)CC1=CC=NC=C1": "pirydyl",
    "O1C=CC=C1(C(=O)CN)": "furyl",
    "O1C=C(C(=O)CN)C=C1": "furyl",
    "C1CCCC1(C(=O)CN)": "cyklopentyl",
    "C1CCCCC1(C(=O)CN)": "cykloheksyl",
    "C1=CC=CC(C(=O)CN)=C1": "fenyl"
}


def find_substructure(systems_map: dict, mol: Chem.rdchem.Mol) -> Union[bool, Chem.rdchem.Mol, tuple]:
    substructure = None
    for system, name in systems_map.items():
        if mol.HasSubstructMatch(Chem.MolFromSmiles(system)):
            substructure = Chem.MolFromSmiles(system)
            matches = mol.GetSubstructMatches(substructure)
            return substructure, matches, name

    if substructure is None:
        return False


def get_central_atoms(mol: Chem.rdchem.Mol, suspected: tuple, desc: list) -> Union[list, bool]:
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
            if add[2] == "N" and add[3] is False and len(neighbors) == 3:
                res.extend([add, add])
            if add[2] == "N" and add[3] is False and len(neighbors) == 2:
                res.append(add)
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
            if add[2] == "C" and add[3] is True and add[4] is True and len(neighbors) == 3:
                if len(ring) == 1:
                    res.append(add)
                if len(ring) > 1:
                    for s_idx in suspected:
                        atom = mol.GetAtomWithIdx(s_idx)
                        if atom.GetSymbol() == "N" and atom.IsInRing() is True:
                            pass  # nitrogen atom is in cyclic system condensed with system ring
                        if any(set(lst).issubset(set(list(suspected))) for lst in
                               [list(idxs) for idxs in atom_rings if s_idx in idxs]) and atom.GetSymbol() == "C":
                            if all([n.IsInRing() for n in atom.GetNeighbors()]):
                                pass  # system condensed with main system

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

        s = Chem.MolToSmiles(s)  # object must be re-converted
        s = Chem.MolFromSmiles(s)

        return s

    s2edit = Chem.MolToSmiles(s2edit)
    s2edit = Chem.MolFromSmiles(s2edit)

    return s2edit


def rs(s: Chem.rdchem.Mol, part_s: Chem.rdchem.Mol, num_heavy_atoms: int,
       carbons: int, permitted_atoms: list, desc: list) -> bool:
    second_num_s = s.GetNumHeavyAtoms() - part_s.GetNumHeavyAtoms()
    s_ring_atoms = [atom.IsInRing() for atom in s.GetAtoms()].count(True)
    s_part_ring_atoms = [atom.IsInRing() for atom in part_s.GetAtoms()].count(True)
    part_s_carbons = [atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()].count(True)

    if part_s is None:
        desc.append("Zawiera niedozwolony podstawnik cykliczny.")
        return False
    else:
        if num_heavy_atoms <= 14:
            if all(atom.GetAtomicNum() in permitted_atoms for atom in s.GetAtoms()):

                if all(atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()) and part_s_carbons <= 6:
                    desc.append(f"Zawiera łańcuch węglowy; {carbons} atomów węgla.")
                    return True

                if Fragments.fr_halogen(s) == 1 and num_heavy_atoms == 1:
                    desc.append(f"Zawiera atom {Chem.MolToSmiles(s)}.")
                    return True

                if part_s.HasSubstructMatch(Chem.MolFromSmiles("O")):
                    tocheck = until(part_s, "O", 1)
                    if [atom.GetAtomicNum() == 6 for atom in tocheck.GetAtoms()].count(True) <= 6:

                        if Fragments.fr_COO(tocheck) == 1:

                            desc.append(f"Zawiera grupę karboksylową; {part_s_carbons} atomów węgla. "
                                        f"Druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")
                            return True

                        else:
                            desc.append(f"Zawiera grupę alkoksylowa; {part_s_carbons} atomów węgla. "
                                        f"Druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")
                            return True
                if part_s_carbons <= 6:
                    if Fragments.fr_sulfone(part_s) == 1:
                        desc.append(f"Zawiera grupę alkilosulfonową; {part_s_carbons} atomów węgla.")
                        return True
                    if Fragments.fr_nitro(part_s) == 1:  # cannot find nitro group directly attach to benzene ring
                        desc.append(f"Zawiera grupę nitrową; {part_s_carbons} atomów węgla.")
                        return True

                if s.GetAtomWithIdx(0).GetAtomicNum() == 6:
                    for p_at in [Chem.Atom(at).GetSymbol() for at in permitted_atoms[1:]]:
                        if s.HasSubstructMatch(Chem.MolFromSmiles(p_at)):
                            smi = Chem.MolToSmiles(until(s, p_at, -1))
                            smi = Chem.MolFromSmiles(smi)
                            if all(atom.GetAtomicNum() == 6 for atom in smi.GetAtoms()) and smi.GetNumHeavyAtoms() <= 6:
                                desc.append(f"Zawiera łańcuch węglowy; {smi.GetNumHeavyAtoms()} atomów węgla; "
                                            f"druga część podstawnika zawiera {num_heavy_atoms - smi.GetNumHeavyAtoms()} atomów.")
                                return True

                else:
                    desc.append("Nie spełnia warunków. Zawiera niedozwolony atom lub grupę atomów.")
                    return False
        else:
            desc.append(f"Zawiera niedopuszczalną liczbę atomów: {num_heavy_atoms}.")
            return False  # too many atoms


def r12(s: Chem.rdchem.Mol, part_s: Chem.rdchem.Mol, ring_part_s: Chem.rdchem.Mol, is_in_ring: int,
        num_heavy_atoms: int, permitted_atoms: list, desc: list) -> bool:
    second_num_s = s.GetNumHeavyAtoms() - part_s.GetNumHeavyAtoms()
    s_ring_atoms = [atom.IsInRing() for atom in s.GetAtoms()].count(True)
    s_part_ring_atoms = [atom.IsInRing() for atom in part_s.GetAtoms()].count(True)
    part_s_carbons = [atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()].count(True)

    if num_heavy_atoms - is_in_ring <= 16:
        if ring_part_s is not None:
            if all(atom.GetAtomicNum() in [6, 8, 16, 7] for atom in
                   ring_part_s.GetAtoms()) and 5 <= ring_part_s.GetNumHeavyAtoms() <= 7:
                ring_at = ring_part_s.GetNumHeavyAtoms()
                desc.append(f"Zawiera pierścień, składający się z {ring_at} atomów.")
                return True
        else:
            if all(atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()) and part_s_carbons <= 6:
                desc.append(f"Zawiera łańcuch węglowy; {part_s_carbons} atomów węgla; "
                            f"druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")
                return True

            if Fragments.fr_Al_OH(part_s) == 1:
                desc.append(f"Zawiera grupę hydroksylową; {part_s_carbons} atomów węgla; "
                            f"druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")

                return True

            if sum(1 for atom in part_s.GetAtoms() if atom.GetAtomicNum() == 8) == 1 and part_s_carbons <= 6:
                desc.append(f"Zawiera grupę alkilokarbonylową; {part_s_carbons} atomów węgla; "
                            f"druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")

                return True

            if Fragments.fr_NH2(part_s) == 1 or Fragments.fr_NH1(part_s) or Fragments.fr_NH0(part_s) == 1:
                desc.append("Zawiera grupę aminową.")
                return True

            elif all(atom.GetAtomicNum() in permitted_atoms for atom in s.GetAtoms()):
                desc.append(f"Druga część podstawnika zawiera {num_heavy_atoms} dozwolonych atomów.")
                return True
            else:
                desc.append("Nie spełnia warunków. Zawiera niedozwolony atom lub grupę atomów.")
                return False
    else:
        desc.append(
            f"Zawiera niedopuszczalną liczbę atomów, które nie są częścią pierścienia: {num_heavy_atoms - is_in_ring}.")
        return False  # to many atoms


def r3456(s: Chem.rdchem.Mol, part_s: Chem.rdchem.Mol, ring_part_s: Chem.rdchem.Mol, condense: bool, is_in_ring: int,
          num_heavy_atoms: int, permitted_atoms: list, desc: list) -> bool:
    second_num_s = s.GetNumHeavyAtoms() - part_s.GetNumHeavyAtoms()
    s_ring_atoms = [atom.IsInRing() for atom in s.GetAtoms()].count(True)
    s_part_ring_atoms = [atom.IsInRing() for atom in part_s.GetAtoms()].count(True)
    part_s_carbons = [atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()].count(True)

    if num_heavy_atoms - is_in_ring <= 20:
        if ring_part_s is not None:

            if all(atom.GetAtomicNum() == 6 for atom in ring_part_s.GetAtoms()) and part_s_carbons <= 10:
                ring_at = ring_part_s.GetNumHeavyAtoms()
                desc.append(f"Zawiera pierścień, składający się z {ring_at} atomów.")
                return True
        else:

            if condense:
                if all(atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()) and 4 <= part_s_carbons <= 6:
                    desc.append(f"Zawiera pierścień skondensowany ze strukturą główną; {part_s_carbons} atomów węgla.")
                    return True

            if Fragments.fr_halogen(s) == 1 and num_heavy_atoms == 1:
                desc.append(f"Zawiera atom {Chem.MolToSmiles(s)}.")
                return True
            if all(atom.GetAtomicNum() == 6 for atom in part_s.GetAtoms()) and part_s_carbons <= 6:
                desc.append(f"Zawiera łańcuch węglowy; {part_s_carbons} atomów węgla; "
                            f"druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")
                return True

            if part_s.HasSubstructMatch(Chem.MolFromSmiles("O")):
                tocheck = until(part_s, "O", 0)
                if Fragments.fr_Al_OH(tocheck) == 1 or num_heavy_atoms == 1:
                    desc.append(f"Zawiera grupę hydroksylową; {part_s_carbons} atomów węgla.")
                    return True

                tocheck = until(part_s, "O", 1)
                if Fragments.fr_C_O_noCOO(tocheck) == 1 and [atom.GetAtomicNum() == 6 for atom in
                                                             tocheck.GetAtoms()].count(True) <= 10:
                    desc.append(f"Zawiera grupę alkilokarbonylową; {part_s_carbons} atomów węgla. "
                                f"Druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")
                    return True

                else:
                    desc.append(f"Zawiera grupę alkoksylowa; {part_s_carbons} atomów węgla. "
                                f"Druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")
                    return True

            if part_s.HasSubstructMatch(Chem.MolFromSmiles("S")):
                tocheck = until(part_s, "S", 1)
                if Fragments.fr_sulfone(tocheck) == 1 and [atom.GetAtomicNum() == 6 for atom in
                                                           tocheck.GetAtoms()].count(True) <= 10:
                    desc.append(f"Zawiera grupę alkilosulfonylową; {part_s_carbons} atomów węgla. "
                                f"Druga część podstawnika zawiera {second_num_s} atomów, w tym {s_ring_atoms - s_part_ring_atoms} atomów w pierścieniu.")
                    return True

            elif all(atom.GetAtomicNum() in permitted_atoms for atom in s.GetAtoms()):
                desc.append(f"Druga część podstawnika zawiera {num_heavy_atoms} dozwolonych atomów.")
                return True
            else:
                desc.append("Nie spełnia warunków. Zawiera niedozwolony atom lub grupę atomów.")
                return False
    else:
        desc.append(
            f"Zawiera niedopuszczalną liczbę atomów, które nie są częścią pierścienia: {num_heavy_atoms - is_in_ring}.")
        return False


def classifier(smiles: str, systems_map: dict) -> tuple:
    try:
        smiles = max(smiles.split("."), key=len)  # remove the radicals
        mol = Chem.MolFromSmiles(smiles)
        desc = []
        mol2move = mol

        mw = round(Descriptors.ExactMolWt(mol), 2)
        if mw <= 500:
            desc.append(f"Masa molowa: {mw}.")
            if find_substructure(systems_map, mol) is not False:
                substructure, matches, name = find_substructure(systems_map, mol)
                desc.append(f"Znaleziono układ cykliczny: {name}.")
                for suspected in matches:
                    i = 0
                    res = get_central_atoms(mol, suspected, desc)
                    if res is not False:
                        mol2edit = Chem.RWMol(mol)
                        for atom in reversed(sorted(suspected)):
                            mol2edit.RemoveAtom(atom)
                        substituents = Chem.MolToSmiles(mol2edit, canonical=True).split(".")

                        to_m = []
                        substituents4condense = []
                        mol2substituents = Chem.RWMol(mol)
                        res = sorted(res)
                        for c_idx in res:
                            mol2substituents.RemoveAtom(c_idx[0])
                            substituent = Chem.MolToSmiles(mol2substituents, canonical=True).split(".")[
                                0]  # move list to the next step instead od element
                            to_m.append(Chem.MolToSmiles(mol2substituents, canonical=True).split("."))

                            substituents4condense.append(substituent)
                            mol2substituents = Chem.RWMol(mol)

                        s4condense = list(
                            (Counter(substituents) - Counter(substituents4condense)).elements())  # overwrite
                        result = []

                        for r, s in zip(res, to_m):
                            s = [i for i in s if i in substituents]  # filtering
                            if len(r) == 6:  # if condense is True
                                s = s4condense[0]
                            permitted_atoms = [6, 8, 7, 9, 17, 35, 53, 16]

                            s = s[0]
                            s = Chem.MolFromSmiles(s)

                            non_ring_s, is_in_ring = if_in_ring(s, False)  # list of [bool, atom_idx]
                            ring_s, _ = if_in_ring(s, True)
                            is_in_ring = [i[1] for i in is_in_ring if i[1] is True].count(True)

                            carbons = [atom.GetAtomicNum() == 6 for atom in s.GetAtoms()].count(True)
                            num_heavy_atoms = s.GetNumHeavyAtoms()

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
                                desc.append(f"Podstawnik R1-2:")
                                aliphatic_nitrogen = r12(s, part_s, ring_part_s, is_in_ring, num_heavy_atoms,
                                                         permitted_atoms, desc)
                                result.append(aliphatic_nitrogen)

                            if r[2] == "C" and r[4] is True:
                                desc.append(f"Podstawnik R:")
                                ring_carbon = rs(s, part_s, num_heavy_atoms, carbons, permitted_atoms, desc)
                                result.append(ring_carbon)

                            if r[2] == "C" and r[4] is False:
                                condense = False
                                if len(r) == 6:  # if condense is True
                                    condense = True
                                desc.append(f"Podstawnik R3-6:")
                                aliphatic_carbon = r3456(s, part_s, ring_part_s, condense, is_in_ring, num_heavy_atoms,
                                                         permitted_atoms, desc)
                                result.append(aliphatic_carbon)

                        i += 1
                        desc = " ".join(desc)
                        if not all(i for i in result) and i < len(matches):
                            desc = []
                            continue
                        else:
                            return all(i for i in result), desc, suspected, mol2move

            else:
                desc.append("Struktura główna nie została znaleziona.")
                desc = " ".join(desc)
                return False, desc, None, mol2move

        else:
            desc.append(f"Dopuszczalna masa molowa została przekroczona: {mw}.")
            desc = " ".join(desc)
            return False, desc, None, mol2move  # mw above 500

    except Exception:
        return False, "Do weryfikacji", None, mol2move


# res, desc, suspected, mol2move = classifier(smiles, systems_map_I)
