# This script optimizes the get ancestors function for every entity in the corpora
# It outputs a file with the ancestors of each entity in json format.

import ssmpy as ssm
import xmltodict
from multiprocessing import Pool
from collections import ChainMap
import json
import sys
import requests
import traceback


def get_ids():
    infile = open(sys.argv[1], 'r', encoding='utf-8')
    data = infile.readlines()
    infile.close()

    s = set()
    for line in data:
        line = line.rstrip()
        line = line.split(' ')
        identifier = line[-1]
        identifier = identifier.replace(':', '_')
        s.add(identifier)
    return s


def get_path_to_root(entity_id):
    if entity_id.startswith('CHEBI'):
        ssm.semantic_base('bin/DiShIn/chebi.db')
    if entity_id.startswith('GO'):
        ssm.semantic_base('bin/DiShIn/go.db')

    e1 = ssm.get_id(entity_id.replace(':', '_'))

    a = ssm.common_ancestors(e1, e1)
    a = [ssm.get_name(x) for x in a]

    return a


def get_path_api(entity_id):
    if entity_id.startswith('GO'):
        go_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms//ancestors"
        go_index = go_url.find('/ancestors')
        go_url = go_url[:go_index] + entity_id.replace('_',
                                                       ':') + go_url[go_index:]
        response = requests.get(go_url)
        json_data = json.loads(response.text)
        if json_data['results'][0]['isObsolete']:
            path = []
        else:
            path = json_data['results'][0]['ancestors']

    elif entity_id.startswith('CHEBI'):
        chebi_url = "https://www.ebi.ac.uk/webservices/chebi/2.0/test/getOntologyParents?chebiId="
        chebi_url = chebi_url + entity_id.replace('_', ':')
        response = requests.get(chebi_url)
        obj_data = xmltodict.parse(response.text)
        json_data = json.dumps(obj_data)
        chebi_data = json.loads(json_data)
        chebi_data = chebi_data['S:Envelope']['S:Body'][
            'getOntologyParentsResponse']['return']['ListElement']
        if type(chebi_data) is list:
            path = [e['chebiId'] for e in chebi_data]
        elif type(chebi_data) is dict:
            path = [chebi_data['chebiId']]

    return path


def get_ancestors(id):
    """Obtain the path to lowest common ancestor of each entity of each pair and path from LCA to root
    :param sentence_labels: list of (e1, e2)
    :param sentence_entities: dictionary mapping entity ID to ((e_start, e_end), text, paths_to_root)
    :return: left and right paths to LCA
    """
    root_path = get_path_to_root(id)
    api_path = get_path_api(id)
    set_path = list(set(root_path + api_path))
    path = [
        i.replace(':', '_') for i in set_path
        if i.startswith('CHEBI') or i.startswith('GO')
    ]
    print('get_ancestors', id)
    path = {id: path}
    return path


if __name__ == "__main__":
    pool = Pool(20)
    ids = get_ids()
    outfile = open(sys.argv[2], 'w')
    try:
        results = pool.map(get_ancestors, ids, 427)
    except Exception:
        print(traceback.format_exc())
    data = dict(ChainMap(*results))
    # Serializing json
    json_object = json.dumps(data, indent=4)
    # Writing to sample.json
    outfile.write(json_object)