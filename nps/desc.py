from nps import nps12
from nps.nps12 import systems_map_I, systems_map_II


def main(smiles: str):

    res_I, desc_I, _ = nps12.classifier(smiles, systems_map_I)
    res_II, desc_II, _ = nps12.classifier(smiles, systems_map_II)

    if res_II:
        description = desc_II
    elif res_I:
        description = desc_I
    else:
        description = ""
        
    return description

      