from rdkit import Chem
from rdkit.Chem import Descriptors


systems_map = {
    "NC(=O)[C@]1C=C2c3cccc4c3c(cn4)C[C@@]2NC1": "structure_II_II",
    "[H][C@@]2(C(N1C(C)CC(C)1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34": "structure_II_azetidine_II",
    "[H][C@@]2(C(N1C(C)C(C)C1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34": "structure_II_azetidine_I",
    "[H][C@@]2(C(N1CCCC1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34": "structure_II_pyrrolidine",
    "[H][C@@]2(C(N1CCOCC1)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34": "structure_II_morpholine",
    "[H][C@@]2(C(N)=O)C=C1c3cccc4[nH]cc(C[C@@]1([H])NC2)c34": "structure_II",
    "NCCc1c[nH]c2ccccc12": "structure_I",
    "NCCc1cnc2ccccc12": "structure_I_I",
    "C1CCCN1CCc1c[nH]c2ccccc12": "structure_I_pyrrolidine",
    "NCCc1c[nH]c2cccc(OP(=O)(O)[O-])c12": "structure_I_dihydrogen_phosphate"
}


def find_smarts_substructure(systems_map: dict, mol: Chem.rdchem.Mol) -> bool | Chem.rdchem.Mol | tuple: 
    substructure = None
    for system, name in systems_map.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(system)):
            substructure = Chem.MolFromSmarts(system)
            matches = mol.GetSubstructMatches(substructure)
            return substructure, matches, name

    if substructure is None:
        return False
        

def classifier(smiles: str) -> tuple:
    
    try:
        smiles = max(smiles.split("."), key=len)  # remove the radicals
        mol = Chem.MolFromSmiles(smiles)
        desc = []
        permitted_at = []
        mw = round(Descriptors.ExactMolWt(mol), 2)
        if mw <= 500:    
            desc.append(f"Masa molowa: {mw}.")
            if find_smarts_substructure(systems_map, mol) is not False:
                    substructure, matches, name = find_smarts_substructure(systems_map, mol)
                    suspected = list(matches[0])
                    if "structure_II" in name:
                        desc.append(f"Znaleziono strukturę II.")
                        desc = " ".join(desc)
                        return True, desc, suspected
                    else:
                        desc = " ".join(desc)
                        desc.append(f"Znaleziono strukturę I.")
                        return True, desc, suspected
        
            else:
                desc.append("Struktura główna nie została znaleziona.")
                return False, desc, None
        
        else:
            desc.append(f"Dopuszczalna masa molowa została przekroczona: {mw}.")
            desc = " ".join(desc)
            return False, desc, None  # mw above 500
            
    except Exception:
        return False, ["Do weryfikacji"], None


# smiles = "CC(C)CC1C(=O)N2CCCC2C3(N1C(=O)C(O3)(C(C)C)NC(=O)C4CN(C5CC6=CNC7=CC=CC(=C67)C5=C4)C)O"
# classifier(smiles)
