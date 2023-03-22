import obonet
import networkx


def load_chebi(path='data/chebi.obo'):
    print('\nLoading the ChEBI Ontology from {}...'.format(path))
    graph = obonet.read_obo(path)
    graph.add_node('BE:00000', name='ROOT')  # biomedical entity general root concept
    graph.add_edge('CHEBI:24431', 'BE:00000', edgetype='is_a')  # chemical_entity
    graph = graph.to_directed()
    is_a_graph = networkx.MultiDiGraph([(u, v, d) for u, v, d in graph.edges(data=True) if d['edgetype'] == 'is_a'])

    id_to_name = {}
    name_to_id = {}
    for id_, data in graph.nodes(data=True):
        try:
            id_to_name[id_] = data['name']
            name_to_id[data['name']] = id_
        except KeyError:
            pass

    id_to_index = {e: i + 1
                   for i, e in enumerate(graph.nodes())
                   }  # ids should start on 1 and not 0
    id_to_index[''] = 0
    
    synonym_to_id = {}
    print('Synonyms to IDs...')
    for n in graph.nodes(data=True):
        for syn in n[1].get('synonym', []):
            syn_name = syn.split('"')
            if len(syn_name) > 2:
                syn_name = syn.split('"')[1]
                synonym_to_id.setdefault(syn_name, []).append(n[0])

    print('Done:', len(name_to_id), 'IDs with', len(synonym_to_id), 'synonyms.')

    return is_a_graph, name_to_id, synonym_to_id, id_to_name, id_to_index


def load_genego(path='data/genego.txt'):
    gene_go = {}
    with open(path, mode='r') as in_file:
        for line in in_file:
            line = line.rstrip()
            line = line.split(' ')
            gene_go[line[0]] = line[1]
    return gene_go


def load_go(path='data/go.obo'):
    print('\nLoading the Gene Ontology from {}...'.format(path))
    graph = obonet.read_obo(path)
    graph.add_node('BE:00000', name='ROOT')  # biomedical entity general root concept
    graph.add_edge('GO:0008150', 'BE:00000', edgetype='is_a')  # biological_process
    graph = graph.to_directed()

    is_a_graph = networkx.MultiDiGraph([(u, v, d)
                                        for u, v, d in graph.edges(data=True)
                                        if d['edgetype'] == 'is_a'])

    id_to_name = {
        id_: data['name']
        for id_, data in graph.nodes(data=True) if 'namespace' in data
        if data['namespace'] == 'biological_process'
    }  # only biological_process
    
    name_to_id = {
        data['name']: id_
        for id_, data in graph.nodes(data=True) if 'namespace' in data
        if data['namespace'] == 'biological_process'
    }  # only biological_process
    name_to_id.update(load_genego())

    id_to_index = {e: i + 1
                   for i, e in enumerate(graph.nodes())
                   }  # ids should start on 1 and not 0
    id_to_index[''] = 0
    
    synonym_to_id = {}
    print('Synonyms to IDs...')
    for n in graph.nodes(data=True):
        for syn in n[1].get('synonym', []):
            syn_name = syn.split('"')
            if len(syn_name) > 2:
                syn_name = syn.split('"')[1]
                synonym_to_id.setdefault(syn_name, []).append(n[0])

    print('Done:', len(name_to_id), 'IDs with', len(synonym_to_id), 'synonyms.')

    return is_a_graph, name_to_id, synonym_to_id, id_to_name, id_to_index


def load_data(path):
    all_data = {}
    with open(path, mode='r', encoding='utf-8') as in_file:
        for line in in_file:
            line = line.rstrip()
            line = line.split(' ')
            identifier = line[-1]
            text = ' '.join(line[:-1])
            all_data[text] = identifier
    return all_data
