from nps import nps12
from nps.nps12 import systems_map_I, systems_map_II
from rdkit.Chem import Draw
from rdkit import Chem
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO


def main(smiles: str):

    _, _, suspected_I, mol = nps12.classifier(smiles, systems_map_I)
    _, _, suspected_II, mol = nps12.classifier(smiles, systems_map_II)

    suspected = ()

    if suspected_II:
        suspected = suspected_II
    if suspected_I:
        suspected = suspected_I

    img = Draw.MolsToGridImage([mol], molsPerRow=1,
                               highlightAtomLists=[list(suspected)], subImgSize=(1200, 1200))
    main_mol = np.array(img)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.axis("off")
    plt.imshow(main_mol)
    imagefile = BytesIO()
    img.save(imagefile, format="PNG")
    imagedata = imagefile.getvalue()

    return imagedata

# main("O=C(CN[C@@H](Cc1ccc(OC)cc1)C)c2ccc(O)c(c2)N")
