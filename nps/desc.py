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

    description = ""
    if nps6.classifier(smiles, systems_map_VI):
        res_VI, desc_VI, _, _ = nps6.classifier(smiles, systems_map_VI)
        if res_VI:
            description = desc_VI

    if nps5.classifier(smiles, systems_map_V):
        res_V, desc_V, _, _ = nps5.classifier(smiles, systems_map_V)
        if res_V is True and res_VI is False:
            description = desc_V

    if nps4.classifier(smiles, systems_map_IV):
        res_IV, desc_IV, _, _ = nps4.classifier(smiles, systems_map_IV)
        if res_IV is True and res_VI is False and res_V is False:
            description = desc_IV

    if nps12.classifier(smiles, systems_map_II):
        res_II, desc_II, _, _ = nps12.classifier(smiles, systems_map_II)
        if res_II is True and res_VI is False and res_V is False and res_IV is False:
            description = desc_II

    if nps12.classifier(smiles, systems_map_I):
        res_I, desc_I, _, _ = nps12.classifier(smiles, systems_map_I)
        if res_I is True and res_VI is False and res_V is False and res_IV is False and res_II is False:
            description = desc_I
        
    if nps3.classifier(smiles, systems_map_III):
        res_III, desc_III, _, _ = nps3.classifier(smiles, systems_map_III)
        if res_III is True and res_VI is False and res_V is False and res_IV is False and res_II is False and res_I is False:
            description = desc_III
    else:
        description = "" 

    return description
