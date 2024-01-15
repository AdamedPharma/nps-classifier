from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Fragments
from typing import Union

import warnings

warnings.filterwarnings("ignore")

systems_map_III = {
    "c1cc2c(cc1)ncc2": "indol-1,3-diyl",  # a
    "c1cc2c(cc1)n(nc2)": "indazol-1,3-diyl",  # b  # [nH0] - cannot have s
    "c3cc4ncnc4cc3": "benzimidazol-1,2-diyl izomer I",  # c
    "c1cc2ncn(c2cc1)": "benzimidazol-1,2-diyl izomer II",  # d # łącznik to **
    "c1c2c(ncc1)c(cn2)": "4-azaindol-1,3-diyl",  # e
    "c1c2c(ncc1)c(nn2)": "4-azaindazol-1,3-diyl",  # f
    "c1c2c(cnc1)c(cn2)": "5-azaindol-1,3-diyl",  # g
    "c1c2c(cnc1)c(nn2)": "5-azaindazol-1,3-diyl",  # h
    "c1cc2c(cn1)ncc2": "6-azaindol-1,3-diyl",  # i
    "c1cc2c(cn1)n(nc2)": "6-azaindazol-1,3-diyl",  # j
    "c1cc2c(nc1)ncc2": "7-azaindol-1,3-diyl",  # k
    "c1cc2c(nc1)n(nc2)": "7-azaindazol-1,3-diyl",  # l
    "c3cc2c1cc(ccc1nc2cc3)": "karbazol-1,4-diyl",  # m
    "*n3nc(cc3)c4ccccc4": "pirazol-1,5-diyl",  # m
    "*n2nc(cc2c1ccccc1)": "pirazol-1,3-diyl",  # o
    "O=C2c1ccc(cc1N(C=C2))": "4-chinolon-1,3-diyl"  # p
}


def find_smarts_substructure(systems_map: dict, mol: Chem.rdchem.Mol) -> Union[bool, Chem.rdchem.Mol, tuple]:
    substructure = None
    for system, name in systems_map.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(system)):
            substructure = Chem.MolFromSmarts(system)
            matches = mol.GetSubstructMatches(substructure)
            return substructure, matches, name

    if substructure is None:
        return False


def side_chain_linker(mol: Chem.rdchem.Mol, atom_rings: tuple, suspected: list):
    # get central N for side chain, then match the rignt chain
    # get central C for linker, then match the rignt S

    # ring_info_mol_reconvert = mol.GetRingInfo()
    # atom_rings = ring_info_mol_reconvert.AtomRings()  # number of rings in mol

    nitro_sc = []  # lst of nitro in ring, with attach side chain
    nitro_l = []  # lst of nitro in ring, with attach linker
    carbon_l = []  # lst of carbon in ring, with attach linker
    carbon_sc = []  # lst of carbon in ring, with attach side chain

    for ring in atom_rings:
        ring_l = len(ring)  # len of ring in mol
        for idx in list(ring):
            atom = mol.GetAtomWithIdx(idx)

            # symbol, ring size, idx, type, degree
            degree = atom.GetDegree()
            symbol = atom.GetSymbol()

            if symbol == "N" and atom.GetIdx() == idx and atom.GetDegree() == 3 and ring_l == 5:
                nitro_sc.append([symbol, ring_l, idx, "side", degree])

            elif symbol == "C" and atom.GetIdx() == idx and atom.GetDegree() == 3 and idx in suspected and ring_l == 5 and [
                idx in tup for tup in atom_rings].count(True) == 1:
                # [idx in tup for tup in atom_rings].count(True) == 1 - check if atom with idx is only in one ring
                carbon_l.append([symbol, ring_l, idx, "linker", degree])

            elif symbol == "N" and atom.GetIdx() == idx and idx in suspected and ring_l == 5:  # add degree
                nitro_l.append([symbol, ring_l, idx, "linker", degree])

            elif atom.GetSymbol() == "C" and atom.GetIdx() == idx and idx in suspected and ring_l == 5:  # add degree
                if not mol.HasSubstructMatch(Chem.MolFromSmiles("nn", sanitize=False)):
                    carbon_sc.append([symbol, ring_l, idx, "side", degree])

    return nitro_sc, carbon_l, nitro_l, carbon_sc


def side_chain_limit(s: list, desc: list):
    n_permitted_at = [6, 8, 16, 9, 17, 35, 53, 7]  # can be written once
    # ["C", "O", "S", "F", "Cl", "Br", "J", "N"]
    for side in s:
        side = Chem.MolFromSmiles(side, sanitize=False)
        if side is not None:
            if all(atom.GetAtomicNum() in n_permitted_at for atom in
                   side.GetAtoms()):  # check if all in permitted atoms
                if 3 <= side.GetNumHeavyAtoms() <= 11:
                    desc.append(f"Zawiera {side.GetNumHeavyAtoms()} dozwolonych atomów.")
                    return True
                else:
                    desc.append(f"Zawiera niedozwoloną liczbę atomów: {side.GetNumHeavyAtoms()}")
                    return False
            else:
                desc.append("Zawiera niedozwolone atomy.")
                return False
        else:
            desc.append("Zaweira niedozwolony skondensowany łańcuch boczny.")
            return False


def linker_limit(s: list, desc: list):
    # limitation for linker + attched
    n_permitted_at = [6, 8, 16, 9, 17, 35, 53, 7]
    for link in s:
        link = Chem.MolFromSmiles(link, sanitize=False)
        if link is not None:
            if all(atom.GetAtomicNum() in n_permitted_at for atom in
                   link.GetAtoms()):  # check if all in permitted atoms
                if link.HasSubstructMatch(Chem.MolFromSmiles("C=O")) or link.HasSubstructMatch(
                        Chem.MolFromSmiles("C=N")):
                    desc.append(
                        f"Łącznik zawiera wiązanie podwójne do atomu tlenu albo azotu. Grupa przyłączona zawiera {link.GetNumHeavyAtoms()} dozwolonych atomów.")
                    return True
                # handled as none after
            else:
                desc.append("Zawiera niedozwolone atomy.")
                return False

        else:
            desc.append("Brak łącznika oraz grupy przyłączonej.")
            return False


def classifier(smiles: str, systems_map: dict):

    try:
        smiles = max(smiles.split("."), key=len)  # remove the radicals
        mol = Chem.MolFromSmiles(smiles)
        desc = []
        mol2move = mol

        if find_smarts_substructure(systems_map_V, mol) is not False:
            substructure, matches, name = find_smarts_substructure(systems_map_V, mol)
            desc.append(f"Układ cykliczny grupy podstawowej: {name}.")
            mw = round(Descriptors.ExactMolWt(mol), 2)
            desc.append(f"Masa cząsteczkowa: {mw}.")

            suspected = matches[0]

            ring_info_mol_reconvert = mol.GetRingInfo()
            atom_rings = ring_info_mol_reconvert.AtomRings()  # number of rings in mol

            nitro_sc, carbon_l, nitro_l, carbon_sc = side_chain_linker(mol, atom_rings, suspected)
            idxs = []
            res = []
            if nitro_sc:
                idxs.append([idx[2] for idx in nitro_sc][0])
                res.append(nitro_sc[0])
            if carbon_l:
                idxs.append([idx[2] for idx in carbon_l][0])
                res.append(carbon_l[0])
            if nitro_l:
                idxs.append([idx[2] for idx in nitro_l][0])
                res.append(nitro_l[0])
            if carbon_sc:
                idxs.append([idx[2] for idx in carbon_sc][0])
                res.append(carbon_sc[0])

            mol2edit = Chem.RWMol(mol)
            for atom in reversed(sorted(suspected)):
                mol2edit.RemoveAtom(atom)
            substituents = Chem.MolToSmiles(mol2edit, canonical=True).split(".")

            to_m = []
            mol2substituents = Chem.RWMol(mol)
            for id in idxs:
                mol2substituents.RemoveAtom(id)
                substituent = Chem.MolToSmiles(mol2substituents, canonical=True).split(".")[0]
                # # move list to the next step instead of element
                to_m.append(Chem.MolToSmiles(mol2substituents, canonical=True).split("."))

            result = []
            for r, s in zip(res, to_m):
                if r[0] == "N" and r[3] == "side":
                    s = [i for i in s if i in substituents]
                    desc.append("Łańcuch boczny:")
                    scl = side_chain_limit(s, desc)
                    result.append(scl)

                elif r[0] == "C" and r[3] == "linker":
                    desc.append("Łącznik oraz grupa przyłączona:")
                    ll = linker_limit(s, desc)
                    if ll is not None:
                        result.append(ll)
                    else:
                        desc.append("Brak łącznika. Do weryfikacji")
                        desc = " ".join(desc)
                        return True, desc, None, mol2move

                elif r[0] == "N" and r[3] == "linker":
                    desc.append("Łącznik oraz grupa przyłączona:")
                    ll = linker_limit(s, desc)
                    if ll is not None:
                        result.append(ll)
                    else:
                        desc.append("Brak łącznika. Do weryfikacji")
                        desc = " ".join(desc)
                        return True, desc, None, mol2move

                if r[0] == "C" and r[3] == "side":
                    s = [i for i in s if i in substituents]
                    desc.append("Łańcuch boczny:")
                    scl = side_chain_limit(s, desc)
                    result.append(scl)

            desc.append("Zweryfikuj masę cząsteczkową grupy przyłączonej.")
            # linker can be only one, so remove duplicate
            desc = " ".join(desc)
            # if all(i for i in result)
            return all(i for i in result), desc, suspected, mol2move

        else:
            desc.append("Struktura główna nie została znaleziona.")
            desc = " ".join(desc)
            return False, desc, None, mol2move

    except Exception:
        return False, "Do weryfikacji", None


# cns = ["CCCCCN1C=C(C2=CC=CC=C21)C(=O)C3=CC=CC4=CC=CC=C43", "CC(C)[C@@H](C(=O)OC)NC(=O)C1=CN(C2=CC=CC=C21)CCCCCF", "CCCCCN1C=C(C2=CC=CC=C21)C(=O)C3=CC=C(C=C3)OC", "C1=CC=C2C(=C1)C(=CN2CCCCCF)C(=O)OC3=CC=CC4=C3N=CC=C4 ",
#        "C1=CC=C2C(=C1)C(=CN2CCCCCF)C(=O)C3=CC=CC=C3I", "C1=CC=C2C(=C1)C(=CN2CCCCCO)C(=O)C3=CC=CC=C3I", "CCCCCN1C=C(C2=CC=CC=C21)C(=O)OC3=CC=CC4=C3N=CC=C4", "CC(C)C(C(=O)N)NC(=O)C1=NN(C2=CC=CC=C21)CC3CCCCC3",
#        "CC(C)[C@@H](C(=O)OC)NC(=O)C1=NN(C2=CC=CC=C21)CC3CCCCC3", "CC1=C(C2=C3N1C(COC3=CC=C2)CN4CCOCC4)C(=O)C5=CC=CC6=CC=CC=C65", "CCCCCN1C=C(C2=CC=CC=C21)C(=O)C3C(C3(C)C)(C)C", "CC1(C(C1(C)C)C(=O)C2=CN(C3=CC=CC=C32)CCCCCF)C",
#        "CC1(C(C1(C)C)C(=O)C2=CN(C3=CC=CC=C32)CC4=CC=C(C=C4)F)C", "CC(C)C(C(=O)OC)NC(=O)C1=NN(C2=CC=CC=C21)CC3=CC=C(C=C3)F", "CC(C)(C)[C@@H](C(=O)OC)NC(=O)C1=NN(C2=CC=CC=C21)CC3=CC=C(C=C3)F",
#        "CCCCCN1C2=CC=CC=C2C(=N1)C(=O)NC34CC5CC(C3)CC(C5)C4", "Cc1cc2c(c(C)c1C)n(CCC)nc2C(=O)c3cc4ccccc4cc3"]

# smiles = "C1=CC=C2C(=C1)C(=CN2CCCCCF)C(=O)OC3=CC=CC4=C3N=CC=C4"
# res, desc, suspected, mol2move = classifier(smiles, systems_map_III)
# print(res, desc, suspected, mol2move)
