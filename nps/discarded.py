from nps import nps12
from nps.nps12 import systems_map_I

def main(smiles: str):

    res_I, desc_I, _, _ = nps12.classifier(smiles, systems_map_I)
  
    if res_I is False:
        description = desc_I
  
    return description
