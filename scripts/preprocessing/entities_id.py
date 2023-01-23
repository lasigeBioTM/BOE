import ontology_preprocessing
import sys, time
from multiprocessing import Pool
from fuzzywuzzy import process
from fuzzywuzzy import fuzz

is_a_graph, name_to_id_g, synonym_to_id_g, id_to_name, id_to_index = ontology_preprocessing.load_go(
)

is_a_graph, name_to_id_c, synonym_to_id_c, id_to_name, id_to_index = ontology_preprocessing.load_chebi(
)


def map_to_ontology(text, name_to_id, synonym_to_id):
    s = time.time()
    entity_id = None
    if text in name_to_id:
        entity_id = name_to_id[text]
    elif text in synonym_to_id:
        entity_id = synonym_to_id[text][0]
    else:  #add syn
        e = process.extractOne(text,
                               name_to_id.keys(),
                               scorer=fuzz.token_sort_ratio)
        if e[1] < 70:
            entity_syn = process.extract(text,
                                         synonym_to_id.keys(),
                                         limit=10,
                                         scorer=fuzz.token_sort_ratio)
            if entity_syn[0][1] > e[1]:
                e = entity_syn[0]
        if e[0] in name_to_id:
            entity_id = name_to_id[e[0]]
        elif e[0] in synonym_to_id:
            entity_id = synonym_to_id[e[0]][0]
    e = time.time()
    print('map_to_ontology')
    print(e - s)
    return (text, entity_id)


def get_genes(text):
    return map_to_ontology(text, name_to_id_g, synonym_to_id_g)


def get_chemicals(text):
    return map_to_ontology(text, name_to_id_c, synonym_to_id_c)


def get_entities():
    infile = open(sys.argv[1], 'r', encoding='utf-8')
    data = infile.readlines()
    infile.close()

    entities_info = {'GENE': set(), 'CHEMICAL': set()}
    for line in data:
        line = line.strip()
        line = tuple(line.split('\t'))
        if 'GENE' in line[2]:
            entities_info['GENE'].add(line[5])
        else:
            entities_info['CHEMICAL'].add(line[5])
    print('Entities file processed')
    return entities_info


if __name__ == "__main__":
    pool = Pool(20)
    names = get_entities()
    gfile = open(sys.argv[2], 'w', encoding='utf-8')
    cfile = open(sys.argv[3], 'w', encoding='utf-8')
    genes = []
    chemicals = []
    try:
        genes = pool.map(get_genes, names['GENE'])
        chemicals = pool.map(get_chemicals, names['CHEMICAL'])
    except Exception as e:
        print(e)
    gfile.write('\n'.join('%s %s' % x for x in genes))
    cfile.write('\n'.join('%s %s' % x for x in chemicals))