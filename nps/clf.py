from nps import nps12
from nps.nps12 import systems_map_I, systems_map_II


def main(smiles: str):

    res_I, _, _, _ = nps12.classifier(smiles, systems_map_I)
    res_II, _, _, _ = nps12.classifier(smiles, systems_map_II)

    if res_II:
        predict = "NPS_II"
    elif res_I:
        predict = "NPS_I"
    else:
        predict = "nie podlega"

    return predict
