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

    predict = "nie podlega" 
    if nps6.classifier(smiles, systems_map_VI):
        res_VI, _, _, _ = nps6.classifier(smiles, systems_map_VI)
        if res_VI:
            predict = "NPS_VI"

    if nps5.classifier(smiles, systems_map_V):
        res_V, _, _, _ = nps5.classifier(smiles, systems_map_V)
        if res_V and predict == "nie podlega":
            predict = "NPS_V"

    if nps4.classifier(smiles, systems_map_IV):
        res_IV, _, _, _ = nps4.classifier(smiles, systems_map_IV)
        if res_IV and predict == "nie podlega":
            predict = "NPS_IV"

    if nps12.classifier(smiles, systems_map_II):
        res_II, _, _, _ = nps12.classifier(smiles, systems_map_II)
        if res_II and predict == "nie podlega":
            predict = "NPS_II"

    if nps12.classifier(smiles, systems_map_I):
        res_I, _, _, _ = nps12.classifier(smiles, systems_map_I)
        if res_I and predict == "nie podlega":
            predict = "NPS_I"
        
    if nps3.classifier(smiles, systems_map_III):
        res_III, _, _, _ = nps3.classifier(smiles, systems_map_III)
        if res_III and predict == "nie podlega":
            predict = "NPS_III" 

    return predict
