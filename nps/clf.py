from nps import nps12
from nps.nps12 import systems_map_I


def main(smiles: str):

    res_I, _, _ = nps12.classifier("CCC", systems_map_I)
    res_II, _, _ = nps12.classifier("CCC", systems_map_I)

    if res_I:
        predict = "NPS_I"
    elif res_II:
        predict = "NPS_II"
    else:
        predict = "nie podlega"

    return predict
