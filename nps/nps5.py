from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Fragments
from typing import Union


systems_map_V = {

    # _benzen R1 can be benzene or pyridine
    # _tiofuran R1 is tiofuran
    
    # loprazolam
    "CN1CCN(CC1)\C=C1/N=C2N(C1=O)C1=C([C,N]=[C,N][C,N]=[C,N]1)C=NC2": "lop_benzen loprazolamu",
    "CN1CCN(CC1)\C=C1/N=C2N(C1=O)C1=C([C,S]=[C,S]=[C,S]1)C=NC2": "lop_tiofuran loprazolam",

    # ketazolam
    "C1=NC=C2N1C1=C([C,N]=[C,N][C,N]=[C,N]1)C1OC=CC(=O)N1C2": "3keta_benzen ketazolamu",
    "C1=NC=C2N1C1=C([C,S]=[C,S]=[C,S]1)C1OC=CC(=O)N1C2": "3keta_tiofuran ketazolamu",
    "C1=NN=C2N1C1=C([C,N]=[C,N][C,N]=[C,N]1)C1OC=CC(=O)N1C2": "2keta_benzen ketazolamu",
    "C1=NN=C2N1C1=C([C,S]=[C,S]=[C,S]1)C1OC=CC(=O)N1C2": "2keta_tiofuran ketazolamu",
    "N1C2=C([C,N]=[C,N][C,N]=[C,N]2)C2OC=CC(=O)N2CC1": "1keta_benzen ketazolamu",
    "N1C2=C([C,S]=[C,S]=[C,S]2)C2OC=CC(=O)N2CC1": "1keta_tiofuran ketazolamu",

    # oksazolam
    "C1CN2C(O1)C1=C([C,N]=[C,N][C,N]=[C,N]1)N1C=NC=C1C2": "3oksa_benzen oksazolamu",
    "C1CN2C(O1)C1=C([C,S]=[C,S]=[C,S]1)N1C=NC=C1C2": "3oksa_tiofuran oksazolamu",
    "C1CN2C(O1)C1=C([C,N]=[C,N][C,N]=[C,N]1)N1C=NN=C1C2": "2oksa_benzen oksazolamu", 
    "C1CN2C(O1)C1=C([C,S]=[C,S]=[C,S]1)N1C=NN=C1C2": "2oksa_tiofuran oksazolamu",
    "C1CN2C(O1)C1=C([C,N]=[C,N][C,N]=[C,N]1)NCC2": "1oksa_benzen oksazolamu",
    "C1CN2C(O1)C1=C([C,S]=[C,S]=[C,S]1)NCC2": "1osa_tiofuran oksazolamu",

    # chlorodiazepoksyd
    "[O-][N+]1=CC2=C([C,N]=[C,N][C,N]=[C,N]2)N2C=NC=C2C1": "3chlo_benzen chlorodiazepoksydu",
    "[O-][N+]1=CC2=C([C,S]=[C,S]=[C,S]2)N2C=NC=C2C1": "3chlo_tiofuran chlorodiazepoksydu",
    "[O-][N+]1=CC2=C([C,N]=[C,N][C,N]=[C,N]2)N2C=NN=C2C1": "2chlo_benzen chlorodiazepoksydu",
    "[O-][N+]1=CC2=C([C,S]=[C,S]=[C,S]2)N2C=NN=C2C1": "2chlo_tiofuran chlorodiazepoksydu",
    "[O-][N+]1=CC2=C([C,N]=[C,N][C,N]=[C,N]2)NCC1": "1chlo_benzen chlorodiazepoksydu",
    "[O-][N+]1=CC2=C([C,S]=[C,S]=[C,S]2)NCC1": "1chlo_tiofuran chlorodiazepoksydu",
    
    
    # 1,5-benzodiazepine
    "N1C2=C([C,N]=[C,N][C,N]=[C,N]2)N2C=NC=C2CC1": "benzo_153_benzen 1,5-benzodiazepin",
    "N1C2=C([C,S]=[C,S]=[C,S]2)N2C=NC=C2CC1": "benzo_153_tiofuran 1,5-benzodiazepin",
    "N1C2=C([C,N]=[C,N][C,N]=[C,N]2)N2C=NN=C2CC1": "benzo_152_benzen 1,5-benzodiazepin",
    "N1C2=C([C,S]=[C,S]=[C,S]2)N2C=NN=C2CC1": "benzo_152_tiofuran 1,5-benzodiazepin",
    "N1C2=C([C,N]=[C,N][C,N]=[C,N]2)NCCC1": "benzo_151_benzen 1,5-benzodiazepin",
    "N1C2=C([C,S]=[C,S]=[C,S]2)NCCC1": "benzo_151_tiofuran 1,5-benzodiazepin",
    
    # 1,4-benzodiazepine
    "C1=NC=C2N1C1=C([C,S]=[C,S]=[C,S]1)C=NC2": "benzo_143_tiofuran 1,4-benzodiazepin",
    "C1=NC=C2N1C1=C([C,N]=[C,N][C,N]=[C,N]1)C=NC2": "benzo_143_benzen 1,4-benzodiazepin",
    "C1=NN=C2N1C1=C([C,S]=[C,S]=[C,S]1)C=NC2": "benzo_142_tiofuran 1,4-benzodiazepin",
    "C1=NN=C2N1C1=C([C,N]=[C,N][C,N]=[C,N]1)C=NC2": "benzo_142_benzen 1,4-benzodiazepin", 
    "N1C2=C([S,C]~[S,C]~[S,C]2)C=NCC1": "benzo_141_tiofuran 1,4-benzodiazepin",  
    "N1C2=C([C,N]=[C,N][C,N]=[C,N]2)C=NCC1=*": "benzo_141_benzen 1,4-benzodiazepin" 
    
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


def functional(s, desc: list):
    
    num_heavy_atoms = s.GetNumHeavyAtoms()
    if Fragments.fr_halogen(s) == 1 and num_heavy_atoms == 1:
        desc.append("Atom halogenu.")
        
    elif Fragments.fr_nitro(s) == 1:
        desc.append("Grupa nitrowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("N")) and num_heavy_atoms == 1:
        desc.append("Grupa aminowa.")
        
    elif s.HasSubstructMatch(Chem.MolFromSmiles("C1=CC=CC=C1")):
        desc.append("Grupa fenylowa.")

    elif Fragments.fr_pyridine(s):
        desc.append("Pirydyl.")

    elif s.HasSubstructMatch(Chem.MolFromSmarts("C1C~C~C(C~C1)")): # can be double bond inside
        desc.append("Cykloheksyl.")
        
    elif s.HasSubstructMatch(Chem.MolFromSmiles("O")) and num_heavy_atoms == 1:
        desc.append("Grupa hydroksylowa albo karbonylowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("C=O")) and num_heavy_atoms == 3:
        desc.append("Karboksyl.")
        
    elif num_heavy_atoms <= 2 and all(atom.GetAtomicNum() == 6 for atom in s.GetAtoms()):
        desc.append("Metyl albo etyl.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("CCOC(C)=O")):
        desc.append("Grupa etoksykarbonylowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("CN(C)C(C)=O")):
        desc.append("Grupa (N,N-dimetylo) karbamoilowa.")
        
    elif s.HasSubstructMatch(Chem.MolFromSmiles("CN(C)C")):
        desc.append("Grupa (N,N-dimetyloamino)metylowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("CCN(C)CC")):
        desc.append("Grupa (N,N-dietyloamino)metylowa.")
        
    elif s.HasSubstructMatch(Chem.MolFromSmiles("CCN(C)C")):
        desc.append("Grupa (N,N-dimetyloamino)etylowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("CCN(CC)CC")):
        desc.append("Grupa (N,N-dietyloamino)etylowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("FC(F)F")):
        desc.append("Grupa (trifluoro)metylowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("CC#C")):
        desc.append("Grupa prop-2-yn-1-ylowa.")

    elif s.HasSubstructMatch(Chem.MolFromSmiles("N=C")):
        desc.append("Grupa iminowa lub  N-metyloiminowa.")
        
    else:
        desc.append(f"Podstawnik {Chem.MolToSmiles(s)} nie spełnia warunków. Zawiera niedozwolony atom lub grupę atomów.")

    return desc


def classifier(smiles: str, systems_map: dict):
    
    try:
        smiles = max(smiles.split("."), key=len)  # remove the radicals
        mol = Chem.MolFromSmarts(smiles)
        desc = []
        mol2move = mol
        mw = round(Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles)), 2)
        permitted_atoms = [17, 9, 35, 53, 6, 8, 7, 16]
        # carabon, oxygen, nitrogen, sulfur, fluorine, chlorine, bromine, iodine
        if mw <= 600:
            if all(atom.GetAtomicNum() in permitted_atoms for atom in mol.GetAtoms()):
                if find_smarts_substructure(systems_map_V, mol) is not False:
                    substructure, matches, name = find_smarts_substructure(systems_map_V, mol)
                    
                    name = name.split(" ")[1]
                    desc.append(f"Znaleziono pochodną {name}.")
                    
                    suspected = matches[0]
                    mol2edit = Chem.RWMol(mol)
                    for atom in reversed(sorted(suspected)):
                        mol2edit.RemoveAtom(atom)
                    substituents = Chem.MolToSmiles(mol2edit, canonical=True).split(".")
                    substituents = [Chem.MolFromSmiles(s) for s in substituents]
                    
                    if substituents:
                        desc.append("Podstawnik/i: ")
                        
                    for s in substituents:
                        desc = functional(s, desc)
                    
                    desc = " ".join(desc)
                    return True, desc, mol2move
    
                else:
                    desc.append("Struktura główna nie została znaleziona.")
                    desc = " ".join(desc)
                    return False, desc, None
                    
            else:
                desc.append("Substancja zawiera niedozwolone atomy")
                desc = " ".join(desc)
                return False, desc, None
    
        else:
            desc.append(f"Dopuszczalna masa molowa została przekroczona: {mw}.")
            desc = " ".join(desc)
            return False, desc, None  # mw above 600
            
    except Exception:
        return False, "Do weryfikacji", None


# res, desc, mol2move = classifier("CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)CN(C)C(C)=O", systems_map_V)
