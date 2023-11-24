
from nps12 import classifier
from nps.nps12 import systems_map_I, systems_map_II
# from nps12 import systems_map_I, systems_map_II
import io

smiles = "O1C(CCN)=CC2=CC3=C(C=CO3)C=C12"


def highlighted_mol(smiles: str):

    from nps12 import classifier
    # from nps.nps12 import systems_map_I, systems_map_II
    from nps12 import systems_map_I, systems_map_II
    import io
    _, _, suspected = classifier(smiles, systems_map_I)

    from rdkit.Chem import Draw
    from rdkit import Chem
    import matplotlib.pyplot as plt
    import numpy as np

    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolsToGridImage([mol], molsPerRow=1,
                               highlightAtomLists=[list(suspected)], subImgSize=(1200, 1200))
    main_mol = np.array(img)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.axis("off")
    plt.imshow(main_mol)
    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    buf.seek(0)

    with open("output.png", "wb") as f:
        f.write(buf.read())

    buf.close()


highlighted_mol(smiles)
