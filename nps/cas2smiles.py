#!/usr/bin/env python

import json
import time
import urllib.error
import urllib.request


# Function for resolving given CAS number into CID
def cas_to_cid(cas):
    path_prolog = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    path_compound = '/compound/'
    path_name = 'name/'
    path_cas = cas
    path_cas_rest = '/cids/JSON'

    url = path_prolog + path_compound + path_name + path_cas + path_cas_rest

    try:
        request = urllib.request.urlopen(url)
    except urllib.error.HTTPError:
        print('HTTPError while requesting cas', cas)
        return ''

    if request is not None:
        reply = request.read()
        if reply is not None and len(reply) > 0:
            json_out = json.loads(reply)
            cid = json_out['IdentifierList']['CID'][0]
            return cid
    return ''


# Function for searching and extracting SMILES code with entering CID
def cid_to_smiles(cid):
    path_prolog = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
    path_compound = '/compound/'
    path_name = 'cid/'
    path_cid = str(cid)
    path_cid_rest = '/property/IsomericSMILES/JSON'  # get isomeric smiles

    url = path_prolog + path_compound + path_name + path_cid + path_cid_rest

    try:
        request = urllib.request.urlopen(url)
    except urllib.error.HTTPError:
        print('HTTPError while requesting cid', cid)
        return ''

    if request is not None:
        reply = request.read()
        if reply is not None and len(reply) > 0:
            json_out = json.loads(reply)
            smiles = json_out['PropertyTable']['Properties'][0]['IsomericSMILES']
            return smiles
    return ''


# Both functions described above are now called in the third function map_cas_list_to_csv
def map_cas_list_to_csv(list_cas):
    output = ''
    for cas in list_cas:
        cid = cas_to_cid(cas)
        if len(str(cid)) > 0:
            smiles = cid_to_smiles(cid)
            if len(smiles) > 0:
                line = cas + '|' + str(cid) + '|' + smiles
                output = output + line + '\n'

                time.sleep(0.8)
    return output


def main():
    import sys

    # Load list with CAS numbers where SMILES code is to be requested
    list_cas = ["152237-40-6"]  # list of CAS's ex. "152237-40-6"
    output = map_cas_list_to_csv(list_cas)  # 'CAS|CID|SMILES\n'
    smiles = output.split('|')[-1]

    return smiles


if __name__ == '__main__':
    main()


