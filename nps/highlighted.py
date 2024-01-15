from nps import nps12
from nps import nps3
from nps import nps4
from nps import nps5
from nps import nps6

from nps.nps12 import systems_map_I, systems_map_II
from nps.nps3 import systems_map_III
from nps.nps4 import systems_map_IV
from nps.nps5 import systems_map_V
from nps.nps6 import systems_map_VI

from rdkit.Chem import Draw
from rdkit import Chem
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO


def main(smiles: str):
    
    _, _, suspected_I, mol = nps12.classifier(smiles, systems_map_I)
    _, _, suspected_II, mol = nps12.classifier(smiles, systems_map_II)
    # _, _, suspected_III, mol = nps3.classifier(smiles, systems_map_III)
    # _, _, suspected_IV, mol = nps4.classifier(smiles, systems_map_IV)
    # _, _, suspected_V, mol = nps5.classifier(smiles, systems_map_V)
    # _, _, suspected_VI, mol = nps6.classifier(smiles, systems_map_VI)
    
    suspected = ()
    
    # if suspected_VI:
    #     suspected = suspected_VI
    # elif suspected_V:
    #     suspected = suspected_V
    # elif suspected_IV:
    #     suspected = suspected_IV
    # elif suspected_II:
        # suspected = suspected_II
    if suspected_I:
        suspected = suspected_I
    elif suspected_II:
        suspected = suspected_II

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
