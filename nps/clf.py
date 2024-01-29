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

    predict = ""
    if nps6.classifier(smiles, systems_map_VI):
        res_VI, _, _, _ = nps6.classifier(smiles, systems_map_VI)
        if res_VI:
            predict = "NPS_VI"

    elif nps5.classifier(smiles, systems_map_V):
        res_V, _, _, _ = nps5.classifier(smiles, systems_map_V)
        if res_V:
            predict = "NPS_V"

    elif nps4.classifier(smiles, systems_map_IV):
        res_IV, _, _, _ = nps4.classifier(smiles, systems_map_IV)
        if res_IV:
            predict = "NPS_IV"

    elif nps12.classifier(smiles, systems_map_II):
        res_II, _, _, _ = nps12.classifier(smiles, systems_map_II)
        if res_II:
            predict = "NPS_II"

    elif nps12.classifier(smiles, systems_map_I):
        res_I, _, _, _ = nps12.classifier(smiles, systems_map_I)
        if res_I:
            predict = "NPS_I"
        
    elif nps3.classifier(smiles, systems_map_III):
        res_III, _, _, _ = nps3.classifier(smiles, systems_map_III)
        if res_III:
            predict = "NPS_III"
    else:
        predict = "nie podlega" 

    return predict
