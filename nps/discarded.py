from nps import nps12
from nps.nps12 import systems_map_I

def main(smiles: str):

    description = ""
    if nps12.classifier(smiles, systems_map_I):
        res_I, desc_I, _, _ = nps12.classifier(smiles, systems_map_I)
    else:
        description = "Do weryfikacji"
  
    if res_I is False:
        description = desc_I
    else:
        description = ""
  
    return description
