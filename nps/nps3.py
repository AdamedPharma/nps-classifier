from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Fragments
from typing import Union

systems_map_V = {
    # grupy podstawowe, gdzie * to jeden atom łańcucha bocznego
    "c1cc2c(cc1)n(*)cc2": "indol-1,3-diyl",  # a
    "c1cc2c(cc1)n(nc2)*": "indazol-1,3-diyl",  # b  # [nH0] - cannot have s
    "c3cc4ncn(*)c4cc3": "benzimidazol-1,2-diyl izomer I",  # c
    "c1cc2nc(*)n(c2cc1)**": "benzimidazol-1,2-diyl izomer II",  # d # łącznik to **
    "c1c2c(ncc1)c(cn2*)": "4-azaindol-1,3-diyl",  # e
    "c1c2c(ncc1)c(nn2*)": "4-azaindazol-1,3-diyl",  # f
    "c1c2c(cnc1)c(cn2*)": "5-azaindol-1,3-diyl",  # g
    "c1c2c(cnc1)c(nn2*)": "5-azaindazol-1,3-diyl",  # h
    "c1cc2c(cn1)n(*)cc2": "6-azaindol-1,3-diyl",  # i
    "c1cc2c(cn1)n(nc2)*": "6-azaindazol-1,3-diyl",  # j
    "c1cc2c(nc1)n(*)cc2": "7-azaindol-1,3-diyl",  # k
    "c1cc2c(nc1)n(nc2)*": "7-azaindazol-1,3-diyl",  # l
    "c3cc2c1cc(ccc1n(*)c2cc3)": "karbazol-1,4-diyl",  # m
    "*n3nc(cc3)c4ccccc4": "pirazol-1,5-diyl",  # m
    "*n2nc(cc2c1ccccc1)": "pirazol-1,3-diyl",  # o
    "O=C2c1ccc(cc1N(C=C2)*)": "4-chinolon-1,3-diyl"  # p
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


# def linker(mol: Chem.rdchem.Mol, atom_rings: tuple):
#     # get central C for linker, then match the rignt S
#     carbon_s = []  # lst of carbon in ring, with attach linker

#     for ring in atom_rings:
#         # print(ring)  # tuple with idx of atom in ring
#         ring_l = len(ring)  # len of ring in mol
#         for idx in list(ring):
#             atom = mol.GetAtomWithIdx(idx)


def side_chain(mol: Chem.rdchem.Mol, atom_rings: tuple, suspected: list):
    # get central N for side chain, then match the rignt chain
    # get central C for linker, then match the rignt S

    # ring_info_mol_reconvert = mol.GetRingInfo()
    # atom_rings = ring_info_mol_reconvert.AtomRings()  # number of rings in mol

    nitro_sc = []  # lst of nitro in ring, with attach side chain
    carbon_sc = []  # lst of carbon in ring, with attach linker

    for ring in atom_rings:
        # print(atom_rings)
        # print(ring)  # tuple with idx of atom in ring
        ring_l = len(ring)  # len of ring in mol
        for idx in list(ring):
            atom = mol.GetAtomWithIdx(idx)
            # print(atom.GetDegree())
            # print(atom.GetIdx(), atom.GetSymbol())
            if atom.GetSymbol() == "N" and atom.GetIdx() == idx and atom.GetDegree() == atom.GetTotalValence() and ring_l == 5:
                # atom.GetDegree() == atom.GetTotalValence() - discard atoms which must have no S
                nitro_sc.append([atom.GetSymbol(), ring_l, idx])  # information lst, number of nitro and ring size

            elif atom.GetSymbol() == "C" and atom.GetIdx() == idx and atom.GetDegree() == 3 and idx in suspected and ring_l == 5 and [
                idx in tup for tup in atom_rings].count(True) == 1:
                # [idx in tup for tup in atom_rings].count(True) == 1 - check if atom with idx is only in one ring
                carbon_sc.append([atom.GetSymbol(), ring_l, idx])

    return nitro_sc, carbon_sc


def side_chain_limit():
    n_permitted_at = ["C", "O", "S", "F", "Cl", "Br", "J", "N"]
    # if all()


def classifier(smiles: str, systems_map: dict):
    smiles = max(smiles.split("."), key=len)  # remove the radicals
    mol = Chem.MolFromSmiles(smiles)
    desc = []
    mol2move = mol

    if find_smarts_substructure(systems_map_V, mol) is not False:
        substructure, matches, name = find_smarts_substructure(systems_map_V, mol)
        desc.append(f"Układ cykliczny grupy podstawowej: {name}.")

        suspected = matches[0]

        ring_info_mol_reconvert = mol.GetRingInfo()
        atom_rings = ring_info_mol_reconvert.AtomRings()  # number of rings in mol
        nitro_sc, carbon_sc = side_chain(mol, atom_rings, suspected)

        mol2edit = Chem.RWMol(mol)
        for atom in reversed(sorted(suspected)):
            mol2edit.RemoveAtom(atom)
        substituents = Chem.MolToSmiles(mol2edit, canonical=True).split(".")
        # print(substituents)
        
        to_m = []
        mol2substituents = Chem.RWMol(mol)
        idxs_n = [idx[2] for idx in nitro_sc]
        
        # idxs = [idxs_n[0], idxs_c[0]]
        
        for n_id in idxs_n:
            print(n_id)
            mol2substituents.RemoveAtom(n_id)
            substituents = Chem.MolToSmiles(mol2substituents, canonical=True).split(".")
            substituent = Chem.MolToSmiles(mol2substituents, canonical=True).split(".")[0]
            # # move list to the next step instead of element
            to_m.append(Chem.MolToSmiles(mol2substituents, canonical=True).split("."))


        for s, m, res in zip(substituents, to_m, nitro_sc):
            if s in m:
                print(res, s)
            
        
        desc = " ".join(desc)
        return True, desc, suspected, mol2move

    else:
        desc.append("Struktura główna nie została znaleziona.")
        desc = " ".join(desc)
        return False, desc, None, mol2move


smiles = "Cc1cc2c(c(C)c1C)n(CCC)nc2C(=O)c3cc4ccccc4cc3"
res, desc, suspected, mol2move = classifier(smiles, systems_map_V)
res, desc, suspected, mol2move
