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


def main(smiles: str):

    res_I, desc_I, _, _ = nps12.classifier(smiles, systems_map_I)
    res_II, desc_II, _, _ = nps12.classifier(smiles, systems_map_II)
    res_III, desc_III, _, _ = nps3.classifier(smiles, systems_map_III)
    res_IV, desc_IV, _, _ = nps4.classifier(smiles, systems_map_IV)
    res_V, desc_V, _, _ = nps5.classifier(smiles, systems_map_V)
    res_VI, desc_VI, _, _ = nps6.classifier(smiles, systems_map_VI)

    if res_VI:
        description = desc_VI
    elif res_V:
        description = desc_V
    elif res_IV:
        description = desc_IV
    elif res_II:
        description = desc_II
    elif res_I:
        description = desc_I
    elif res_III:
        description = desc_III
    else:
        description = ""
        
    return description
    
