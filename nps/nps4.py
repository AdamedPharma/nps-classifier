from rdkit import Chem
from rdkit.Chem import Descriptors
from typing import Union


systems_map_IV = {
    "structure_I": "C(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1",
    "structure_III": "C(c1ccccc1)CN1CCC(Nc2ccccc2)CC1",
    "structure_II": "C(=O)N(c1ccccc1)C1CCNCC1"
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


def classifier(smiles: str) -> tuple:
    try:

        smiles = max(smiles.split("."), key=len)  # remove the radicals
        mol = Chem.MolFromSmiles(smiles)
        mw = round(Descriptors.ExactMolWt(mol), 2)
        permitted_atoms = [17, 9, 35, 53, 6, 8, 7, 16]
        # fluorine, chlorine, bromine, iodine, carbon, oxygen, nitrogen, sulfur
        desc = []
        if mw <= 500:
            if all(atom.GetAtomicNum() in permitted_atoms for atom in mol.GetAtoms()):

                desc.append(f"Masa molowa: {mw}.")
                if find_smarts_substructure(systems_map_IV, mol) is not False:
                    substructure, matches, name = find_smarts_substructure(systems_map_IV, mol)
                    suspected = list(matches[0])

                    if name == "structure_I":
                        desc.append(f"Znaleziono strukturę I.")
                        desc = " ".join(desc)
                        return True, desc, suspected, mol2move
                        
                    if name == "structure_II":
                        desc.append(f"Znaleziono strukturę II.")
                        desc = " ".join(desc)
                        return True, desc, suspected, mol2move

                    if name == "structure_III":
                        desc.append(f"Znaleziono strukturę III.")
                        desc = " ".join(desc)
                        return True, desc, suspected, mol2move

                else:
                    desc.append("Struktura główna nie została znaleziona.")
                    desc = " ".join(desc)
                    return False, desc, None, mol2move

            else:
                desc.append("Substancja zawiera niedozwolone atomy")
                desc = " ".join(desc)
                return False, desc, None, mol2move
                
        else:
            desc.append(f"Dopuszczalna masa molowa została przekroczona: {mw}.")
            desc = " ".join(desc)
            return False, desc, None, mol2move  # mw above 500

    except Exception:
        return False, "Do weryfikacji", None, mol2move


# smiles = "CCC(=O)N(C1=CC=CC=C1)C2(CCN(CC2)CCC3=CC=CS3)COC"
# classifier(smiles)
