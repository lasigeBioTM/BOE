from __future__ import unicode_literals, print_function
from itertools import combinations, product
from subprocess import PIPE, Popen
from multiprocessing import Process
from ssmpy import ssm

import networkx as nx
import pysbd
import en_core_web_sm
import numpy as np
import time
import sys
import logging
import json
import multiprocessing

import ontology_preprocessing

# logger configuration
logging.basicConfig(filename="biont.log",
                    format='%(asctime)s %(message)s',
                    filemode='w')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

sys.path.append('bin/DiShIn/')  # access bin/DiShIn/

# Input Parameters
sst_light_directory = 'bin/sst-light-0.4/'
temporary_directory = sys.argv[-1]
biomedical_entity_to_ontology = {'gene': 'G', 'drug': 'C'}

neg_gv_list = {}

label_to_pair_type = {
    'NO_RELATION': 0,
    'INHIBITOR': 1,
    'PART-OF': 2,
    'SUBSTRATE': 3,
    'ACTIVATOR': 4,
    'INDIRECT-DOWNREGULATOR': 5,
    'ANTAGONIST': 6,
    'INDIRECT-UPREGULATOR': 7,
    'AGONIST': 8,
    'DIRECT-REGULATOR': 9,
    'PRODUCT-OF': 10,
    'AGONIST-ACTIVATOR': 11,
    'AGONIST-INHIBITOR': 12,
    'SUBSTRATE_PRODUCT-OF': 13
}
pair_type_to_label = {v: k for k, v in label_to_pair_type.items()}

#import debugpy

# Allow other computers to attach to debugpy at this IP address and port.
#debugpy.listen(("0.0.0.0", 5678))

# Pause the program until a remote debugger is attached
#debugpy.wait_for_client()

# --------------------------------------------------------------
#                 DIVIDED BY SENTENCES ABSTRACTS
# --------------------------------------------------------------


def divided_by_sentences(abstract):
    """Divides abstracts by sentences

    :param abstract:load
    :return: list of sentences of the abstract
    """

    seg = pysbd.Segmenter(language="en", clean=False)
    sentences = seg.segment(abstract)

    return sentences


def get_abstracts_info(abstract_file):
    abstracts_info_dict = {}
    with open(abstract_file, encoding='utf-8') as af:
        for line in af:
            line = line.strip()
            line = line.split('\t')
            # k: abstract_id v: text
            abstracts_info_dict[line[0]] = [line[1]] + divided_by_sentences(
                line[2])
        print('Abstract file processed')

    return abstracts_info_dict


def get_entities_info(entities_file):
    entities_info_dict = {}
    with open(entities_file, encoding='utf-8') as ef:
        for line in ef:
            line = line.strip()
            line = tuple(line.split('\t'))
            if line[0] in entities_info_dict:
                entities_info_dict[line[0]].append(line[1:])
            else:
                entities_info_dict[line[0]] = [line[1:]]
    for k, v in entities_info_dict.items():
        entities_info_dict[k] = sorted(v, key=lambda x: int(x[2]))
    print('Entities file processed')

    return entities_info_dict


def get_relations_info(relations_file):
    relations_info_dict = {}
    if relations_file:
        with open(relations_file, encoding='utf-8') as rf:
            for line in rf:
                line = line.strip()
                line = line.split('\t')
                if line[0] in relations_info_dict:
                    relations_info_dict[line[0]].append(
                        (line[1], line[2].split(':')[1],
                         line[3].split(':')[1]))
                else:
                    relations_info_dict[line[0]] = [
                        (line[1], line[2].split(':')[1], line[3].split(':')[1])
                    ]
        print('Relations file processed')

    return relations_info_dict


if temporary_directory == 'temp/':
    abstracts_info = get_abstracts_info(
        'corpora/training/drugprot_training_abstracs.tsv')
    entities_info = get_entities_info(
        'corpora/training/drugprot_training_entities.tsv')
    relations_info = get_relations_info(
        'corpora/training/drugprot_training_relations.tsv')
    g_data = ontology_preprocessing.load_data('corpora/training/t_g_ents.txt')
    c_data = ontology_preprocessing.load_data('corpora/training/t_c_ents.txt')
    g_json = open('../corpora/json/t_g_paths.json')
    g_paths = json.load(g_json)
    c_json = open('../corpora/json/t_c_paths.json')
    c_paths = json.load(c_json)

elif temporary_directory == 'temp_test/':
    abstracts_info = get_abstracts_info(
        'corpora/test-background/test_background_abstracts.tsv')
    entities_info = get_entities_info(
        'corpora/test-background/test_background_entities.tsv')
    relations_info = get_relations_info('')
    all_data = ontology_preprocessing.load_data(
        'corpora/test-background/tb_ent.txt')
    paths = {}

elif temporary_directory == 'temp_sample/':
    abstracts_info = get_abstracts_info(
        'corpora/sample_corpora/sample_training_abstracts.tsv')
    entities_info = get_entities_info(
        'corpora/sample_corpora/sample_training_entities.tsv')
    relations_info = get_relations_info(
        'corpora/sample_corpora/sample_training_relations.tsv')
    # Same as dev for testing
    g_data = ontology_preprocessing.load_data(
        'corpora/development/dev_g_ents.txt')
    c_data = ontology_preprocessing.load_data(
        'corpora/development/dev_c_ents.txt')
    g_json = open('../corpora/json/dev_g_paths.json')
    g_paths = json.load(g_json)
    c_json = open('../corpora/json/dev_c_paths.json')
    c_paths = json.load(c_json)

elif temporary_directory == 'temp_dev/':
    abstracts_info = get_abstracts_info(
        'corpora/development/drugprot_development_abstracs.tsv')
    entities_info = get_entities_info(
        'corpora/development/drugprot_development_entities.tsv')
    relations_info = get_relations_info(
        'corpora/development/drugprot_development_relations.tsv')
    g_data = ontology_preprocessing.load_data(
        'corpora/development/dev_g_ents.txt')
    c_data = ontology_preprocessing.load_data(
        'corpora/development/dev_c_ents.txt')
    g_json = open('../corpora/json/dev_g_paths.json')
    g_paths = json.load(g_json)
    c_json = open('../corpora/json/dev_c_paths.json')
    c_paths = json.load(c_json)

# --------------------------------------------------------------
#                    UPDATE SENTENCES OFFSETS
# --------------------------------------------------------------


def get_overlapped_ents(entities_per_sentence):
    overlapped = set()
    overlapped_offset = 0
    ent = ()

    for previous, current in zip(entities_per_sentence,
                                 entities_per_sentence[1:]):
        if current[2] <= overlapped_offset and current != ent:
            overlapped.add(current)
        if previous[2] <= overlapped_offset and previous != ent:
            overlapped.add(previous)
        if current[1] <= previous[1] <= current[2] and current[
                2] > overlapped_offset:
            overlapped.add(previous)
            overlapped_offset = current[2]
            ent = current
        elif previous[1] <= current[1] <= previous[2] and previous[
                2] > overlapped_offset:
            overlapped.add(current)
            overlapped_offset = previous[2]
            ent = previous
        if ent:
            overlapped.add(ent)

    return overlapped


def get_new_offsets_sentences(abstract_id):
    entities_per_sentence = {}
    abstract = abstracts_info[abstract_id].copy()
    annotation_lines = entities_info[abstract_id].copy()

    limit_1 = 0
    limit_2 = 0
    sentence_id = 0

    for sentence in abstract:
        limit_2 += len(sentence)
        k = ('a' + abstract_id + '.s' + str(sentence_id), sentence)
        entities_per_sentence[k] = []

        for annotation in annotation_lines[:]:
            offset_1 = int(annotation[2]) + 1
            offset_2 = int(annotation[3]) + 1

            if limit_1 <= int(offset_1 - 1) <= limit_2 and limit_1 <= int(
                    offset_2 - 1) <= limit_2:
                updated_offset_1 = int(offset_1) - limit_1
                updated_offset_2 = int(offset_2) - limit_1

                entities_per_sentence[k].append(
                    (annotation[0], updated_offset_1, updated_offset_2,
                     annotation[4], annotation[1]))

                annotation_lines.remove(annotation)
            else:
                break

        overlapped = list(get_overlapped_ents(entities_per_sentence[k]))
        if overlapped:
            '''
            for ent in overlapped:
                entities_per_sentence[k].remove(ent)
            '''
            # Create extra sentences if overlapped entities exist
            overlapped += entities_per_sentence[k].copy()
            overlapped = list(set(overlapped))
            del entities_per_sentence[k]
            extra_entities_per_sentence = []

            for c in list(combinations(overlapped, 2)):
                # GENE-DRUG or DRUG-GENE and not overlapped
                if c[0][4][0] != c[1][4][0] and not get_overlapped_ents(c):
                    extra_entities_per_sentence.append(list(c))

            for i in range(len(extra_entities_per_sentence)):
                if len(extra_entities_per_sentence[i]) > 1:
                    x = ('a' + abstract_id + '.s' + str(sentence_id) + '.x' +
                         str(i), sentence)
                    entities_per_sentence[x] = []
                    for e in extra_entities_per_sentence[i]:
                        entities_per_sentence[x].append(e)
        else:
            ents = entities_per_sentence[k].copy()
            del entities_per_sentence[k]
            tuple_per_sentence = []

            for c in list(combinations(ents, 2)):
                # GENE-DRUG or DRUG-GENE and not overlapped
                if c[0][4][0] != c[1][4][0]:
                    tuple_per_sentence.append(list(c))

            for i in range(len(tuple_per_sentence)):
                if len(tuple_per_sentence[i]) > 1:
                    x = ('a' + abstract_id + '.s' + str(sentence_id) + '.x' +
                         str(i), sentence)
                    entities_per_sentence[x] = []
                    for e in tuple_per_sentence[i]:
                        entities_per_sentence[x].append(e)

        sentence_id += 1
        limit_1 += len(sentence) + 1

    return entities_per_sentence


entities_per_abstract = {
    k: get_new_offsets_sentences(k)
    for k in entities_info.keys()
}

# --------------------------------------------------------------
#                     GET SENTENCE ENTITIES
# --------------------------------------------------------------


def get_sentence_entities():
    """

    :return:
    """

    entities = {}
    positive_entities = set()

    for k in entities_info.keys():

        for sentence, entities_sentence in entities_per_abstract[k].items():

            sentence_entities = {}

            entity_number = 1
            for entity in entities_sentence:
                if 'GENE' in entity[4]:
                    entity_id = g_data[entity[3]]
                elif 'CHEMICAL' in entity[4]:
                    entity_id = c_data[entity[3]]

                if k in relations_info:
                    for pair in relations_info[k]:

                        if entity[0] == pair[1]:
                            positive_entities.add(sentence[0] + '.u' +
                                                  entity[0])

                        elif entity[0] == pair[2]:
                            positive_entities.add(sentence[0] + '.u' +
                                                  entity[0])

                sentence_entities[sentence[0] + '.u' +
                                  entity[0]] = (eval('[' + str(entity[1]) +
                                                     ', ' + str(entity[2]) +
                                                     ']'), entity[3],
                                                entity_id)
                entity_number += 1

            entities[sentence[0]] = sentence_entities

    return entities, positive_entities


# --------------------------------------------------------------
#                     GET ENTITIES ANCESTORS
# --------------------------------------------------------------


def get_path_to_root(entity_id):
    """

    :param entity_id:
    :return:
    """

    if entity_id.startswith('CHEBI'):
        ssm.semantic_base('bin/DiShIn/chebi.db')
    if entity_id.startswith('GO'):
        ssm.semantic_base('bin/DiShIn/go.db')

    e1 = ssm.get_id(entity_id.replace(':', '_'))

    a = ssm.common_ancestors(e1, e1)
    a = [ssm.get_name(x) for x in a]

    return a


def get_ancestors(sentence_labels, sentence_entities):
    """Obtain the path to lowest common ancestor of each entity of each pair and path from LCA to root

    :param sentence_labels: list of (e1, e2)
    :param sentence_entities: dictionary mapping entity ID to ((e_start, e_end), text, paths_to_root)
    :return: left and right paths to LCA
    """
    right_paths = []
    left_paths = []
    common_ancestors = []

    for p in sentence_labels:
        left_ent = sentence_entities[p[0]][2]
        left_ent = left_ent.replace(':', '_')
        right_ent = sentence_entities[p[1]][2]
        right_ent = right_ent.replace(':', '_')

        if left_ent in g_paths:
            left_path = g_paths[left_ent]
        elif left_ent in c_paths:
            left_path = c_paths[left_ent]
        else:
            left_path = get_path_to_root(left_ent)
            paths[left_ent] = left_path

        if right_ent in g_paths:
            right_path = g_paths[right_ent]
        elif right_ent in c_paths:
            right_path = c_paths[right_ent]
        else:
            right_path = get_path_to_root(right_ent)
            paths[right_ent] = right_path

        left_paths.append(left_path)
        right_paths.append(right_path)

    return common_ancestors, (left_paths, right_paths)


# --------------------------------------------------------------
#               PARSE CORPUS SENTENCES USING SPACY
# --------------------------------------------------------------


def prevent_sentence_segmentation(doc):
    """

    :param doc:
    :return:
    """

    for token in doc:
        # This will entirely disable spaCy's sentence detection
        token.is_sent_start = False

    return doc


nlp = en_core_web_sm.load(disable=['ner'])
nlp.add_pipe(prevent_sentence_segmentation,
             name='prevent-sbd',
             before='parser')


def parse_sentence_spacy(sentence_text, sentence_entities):
    """

    :param sentence_text:
    :param sentence_entities:
    :return:
    """

    # Clean text to make tokenization easier
    # Whitespace before and after entity to ease tokenization
    for e in sentence_entities:
        idx = sentence_entities[e][0]
        if idx[0] != 1:
            sentence_text = sentence_text[:idx[0] -
                                          2] + ' ' + sentence_text[idx[0] - 1:]
        if idx[1] != len(sentence_text):
            sentence_text = sentence_text[:idx[1] -
                                          1] + ' ' + sentence_text[idx[1]:]
        else:
            sentence_text = sentence_text[:idx[1] - 1] + ' .'
        sentence_text = sentence_text[:idx[0] - 1] + 'a' * len(
            sentence_entities[e][1]) + sentence_text[idx[1] - 1:]

    # Clean text to make tokenization easier
    sentence_text = sentence_text.replace('[', ' ').replace('\'', ' ').replace('/', ' ') \
        .replace(';', ',').replace('*', ' ').replace(':', ',').replace(']', ' ')

    # Replace entity in text
    for e in sentence_entities:
        idx = sentence_entities[e][0]
        if idx[0] == 1:
            sentence_text = sentence_entities[e][1] + sentence_text[idx[1] -
                                                                    1:]
        else:
            sentence_text = sentence_text[:idx[0] - 1] + sentence_entities[e][
                1] + sentence_text[idx[1] - 1:]
        sentence_text = sentence_text[:idx[0] - 1] + sentence_text[idx[0] - 1:idx[1] - 1].replace(' ', '_')\
            .replace(':', '_') + sentence_text[idx[1] - 1:]

    # Use spacy to parse a sentence
    parsed = nlp(sentence_text)

    return parsed


def run_sst(base_dir, token_seq):
    """

    :param base_dir:
    :param token_seq:
    :return:
    """

    chunk_size = 500
    wordnet_tags = {}
    sent_ids = list(token_seq.keys())

    chunks = [
        sent_ids[i:i + chunk_size] for i in range(0, len(sent_ids), chunk_size)
    ]

    for i, chunk in enumerate(chunks):
        sentence_file = open('{}/sentences_{}.txt'.format(
            temporary_directory + base_dir.split('/')[1], i),
                             'w',
                             encoding='utf-8')

        for sent in chunk:
            sentence_file.write("{}\t{}\t.\n".format(
                sent, '\t'.join(token_seq[sent])))

        sentence_file.close()
        sst_args = [
            sst_light_directory + 'sst', 'bitag',
            '{}/MODELS/WSJPOSc_base_20'.format(sst_light_directory),
            '{}/DATA/WSJPOSc.TAGSET'.format(sst_light_directory),
            '{}/MODELS/SEM07_base_12'.format(sst_light_directory),
            '{}/DATA/WNSS_07.TAGSET'.format(sst_light_directory),
            '{}/sentences_{}.txt'.format(
                temporary_directory + base_dir.split('/')[1], i), '0', '0'
        ]

        p = Popen(sst_args, stdout=PIPE)
        p.communicate()

        with open('{}/sentences_{}.txt.tags'.format(
                temporary_directory + base_dir.split('/')[1], i),
                  encoding='utf-8') as f:
            output = f.read()

        sstoutput = parse_sst_results(output)
        wordnet_tags.update(sstoutput)

    return wordnet_tags


def parse_sst_results(results):
    """

    :param results:
    :return:
    """

    sentences = {}
    lines = results.strip().split('\n')

    for l in lines:
        values = l.split('\t')
        wntags = [x.split(' ')[-1].split('-')[-1] for x in values[1:]]
        sentences[values[0]] = wntags

    return sentences


def parse_sentences_spacy(base_dir, entities):
    """

    :param base_dir:
    :param entities:
    :return:
    """

    # First iterate all documents, and preprocess all sentences
    parsed_sentences = {}
    token_seq = {}

    for k in entities_info.keys():

        for sentence, entities_sentence in entities_per_abstract[k].items():

            parsed_sentence = parse_sentence_spacy(sentence[1],
                                                   entities[sentence[0]])
            parsed_sentences[sentence[0]] = parsed_sentence
            tokens = []

            for t in parsed_sentence:
                tokens.append(
                    t.text.replace(' ', '_').replace('\t',
                                                     '_').replace('\n', '_'))
            token_seq[sentence[0]] = tokens

    wordnet_tags = run_sst(base_dir, token_seq)

    return parsed_sentences, wordnet_tags


def get_network_graph_spacy(document):
    """Convert the dependencies of the spacy document object to a networkX graph

    :param document: spacy parsed document object
    :return: networkX graph object and nodes list
    """

    edges = []
    nodes = []

    # Ensure that every token is connected
    for s in document.sents:
        edges.append(('ROOT', '{0}-{1}'.format(s.root.lower_, s.root.i)))

    for token in document:
        nodes.append('{0}-{1}'.format(token.lower_, token.i))

        for child in token.children:
            edges.append(('{0}-{1}'.format(token.lower_, token.i),
                          '{0}-{1}'.format(child.lower_, child.i)))

    return nx.Graph(edges), nodes


def get_head_tokens_spacy(entities, sentence, positive_entities):
    """

    :param entities: dictionary mapping entity IDs to (offset, text)
    :param sentence: sentence parsed by spacy
    :param positive_entities:
    :return: dictionary mapping head tokens word-idx to entity IDs
    """

    sentence_head_tokens_type_1 = {}
    sentence_head_tokens_type_2 = {}
    pos_gv = set()
    neg_gv = set()

    for eid in entities:

        offset = (entities[eid][0][0], entities[eid][0][1])
        entity_tokens = sentence.char_span(offset[0] - 1, offset[1] - 1)

        if not entity_tokens:
            logging.warning(('No tokens found:', entities[eid], sentence.text,
                             '|'.join([t.text for t in sentence])))

        else:
            head_token = '{0}-{1}'.format(
                entity_tokens.root.lower_.replace(' ', '_'),
                entity_tokens.root.i)

            if eid in positive_entities:
                pos_gv.add(entity_tokens.root.head.lower_)
            else:
                neg_gv.add(entity_tokens.root.head.lower_)

            if head_token in sentence_head_tokens_type_1:
                logging.warning(
                    ('Head token conflict:',
                     sentence_head_tokens_type_1[head_token], entities[eid]))
            elif head_token in sentence_head_tokens_type_2:
                logging.warning(
                    ('Head token conflict:',
                     sentence_head_tokens_type_2[head_token], entities[eid]))

            if biomedical_entity_to_ontology['drug'] == entities[eid][2][0]:
                sentence_head_tokens_type_1[head_token] = eid
            elif biomedical_entity_to_ontology['gene'] == entities[eid][2][0]:
                sentence_head_tokens_type_2[head_token] = eid

    return sentence_head_tokens_type_1, sentence_head_tokens_type_2, pos_gv, neg_gv


def process_sentence_spacy(base_dir,
                           sentence,
                           sentence_entities,
                           sentence_pairs,
                           positive_entities,
                           wordnet_tags=None,
                           mask_entities=True,
                           min_sdp_len=0,
                           max_sdp_len=15):
    """Process sentence to obtain labels, instances and classes for a ML classifier

    :param base_dir:
    :param sentence: sentence processed by spacy
    :param sentence_entities: dictionary mapping entity ID to ((e_start, e_end), text, paths_to_root)
    :param sentence_pairs: dictionary mapping pairs of known entities in this sentence to pair types
    :param positive_entities:
    :param wordnet_tags:
    :param mask_entities:
    :param min_sdp_len:
    :param max_sdp_len:
    :return: labels of each pair (according to sentence_entities, word vectors and classes (pair types according to sentence_pairs)
    """
    left_word_vectors = []
    right_word_vectors = []
    left_wordnets = []
    right_wordnets = []
    classes = []
    labels = []

    graph, nodes_list = get_network_graph_spacy(sentence)
    sentence_head_tokens_type_1, sentence_head_tokens_type_2, pos_gv, neg_gv = get_head_tokens_spacy(
        sentence_entities, sentence, positive_entities)

    entity_offsets = [sentence_entities[x][0][0] for x in sentence_entities]

    for e in product(sentence_head_tokens_type_1.keys(),
                     sentence_head_tokens_type_2.keys()):
        try:
            e1_text = sentence_entities[sentence_head_tokens_type_1[e[0]]]
            e2_text = sentence_entities[sentence_head_tokens_type_2[e[1]]]

        except KeyError:
            print(e[0], e[1])
            print(sentence)
            print(sentence_head_tokens_type_1)
            print(sentence_head_tokens_type_2)

        if e1_text[1].lower() == e2_text[1].lower():
            continue

        if e1_text[0][0] > e2_text[0][0]:
            middle_text = sentence.text[e2_text[0][1]:e1_text[0][0]]
            e1, e2 = e[1], e[0]
        else:
            middle_text = sentence.text[e1_text[0][1]:e2_text[0][0]]
            e2, e1 = e[1], e[0]

        if any(punctuation in middle_text.strip() for punctuation in '!?.'):
            continue

        try:
            sdp = nx.shortest_path(graph, source=e1, target=e2)

            if len(sdp) < min_sdp_len or len(sdp) > max_sdp_len:
                continue

            neg = False
            is_neg_gv = False
            for i, element in enumerate(sdp):
                token_text = element.split('-')[0]
                if (i == 1 or i == len(sdp) - 2) and token_text in neg_gv_list:
                    pass
                    print('Skipped gv {} {}:'.format(token_text, str(sdp)))

            if neg or is_neg_gv:
                continue

            vector = []
            wordnet_vector = []
            negations = 0
            head_token_position = None

            for i, element in enumerate(sdp):
                if element != 'ROOT':
                    token_idx = int(
                        element.split('-')[-1])  # get the index of the token
                    sdp_token = sentence[token_idx]  # get the token object

                    if mask_entities and sdp_token.idx in entity_offsets:
                        vector.append('entity')
                    else:
                        vector.append(sdp_token.text)
                    if wordnet_tags:
                        wordnet_vector.append(wordnet_tags[token_idx])

                    head_token = '{}-{}'.format(
                        sdp_token.head.lower_,
                        sdp_token.head.i)  # get the key of head token

                    # Head token must not have its head in the path, otherwise that would be the head token
                    # In some cases the token is its own head
                    if head_token not in sdp or head_token == element:
                        head_token_position = i + negations

            if head_token_position is None:
                print('Head token not found:', e1_text, e2_text, sdp)
                sys.exit()
            else:
                left_vector = vector[:head_token_position + 1]
                right_vector = vector[head_token_position:]
                left_wordnet = wordnet_vector[:head_token_position + 1]
                right_wordnet = wordnet_vector[head_token_position:]

            left_word_vectors.append(left_vector)
            right_word_vectors.append(right_vector)
            left_wordnets.append(left_wordnet)
            right_wordnets.append(right_wordnet)

        except nx.exception.NetworkXNoPath:
            logging.warning('No path:', e1_text, e2_text, graph.nodes())
            left_word_vectors.append([])
            right_word_vectors.append([])
            left_wordnets.append([])
            right_wordnets.append([])

        except nx.NodeNotFound:
            logging.warning(('Node not found:', e1_text, e2_text, e1, e2,
                             list(sentence), graph.nodes()))
            left_word_vectors.append([])
            right_word_vectors.append([])
            left_wordnets.append([])
            right_wordnets.append([])

        if sentence_head_tokens_type_1.get(e1):
            try:
                labels.append((sentence_head_tokens_type_1[e1],
                               sentence_head_tokens_type_2[e2]))
            except KeyError:
                print(e1, e2)
                print(sentence)
                print(sentence_head_tokens_type_1)
                print(sentence_head_tokens_type_2)

            if (sentence_head_tokens_type_1[e1],
                    sentence_head_tokens_type_2[e2]) in sentence_pairs:
                classes.append(
                    sentence_pairs[(sentence_head_tokens_type_1[e1],
                                    sentence_head_tokens_type_2[e2])])
            else:
                classes.append(0)
        else:
            try:
                labels.append((sentence_head_tokens_type_1[e2],
                               sentence_head_tokens_type_2[e1]))
            except KeyError:
                print(e1, e2)
                print(sentence)
                print(sentence_head_tokens_type_1)
                print(sentence_head_tokens_type_2)

            if (sentence_head_tokens_type_1[e2],
                    sentence_head_tokens_type_2[e1]) in sentence_pairs:
                classes.append(
                    sentence_pairs[(sentence_head_tokens_type_1[e2],
                                    sentence_head_tokens_type_2[e1])])
            else:
                classes.append(0)

    return labels, (left_word_vectors, right_word_vectors), (
        left_wordnets, right_wordnets), classes, pos_gv, neg_gv


# --------------------------------------------------------------
#    PARSE CORPUS WITH SDP VECTORS FOR EACH RELATION INSTANCE
# --------------------------------------------------------------


def get_sdp_instances(base_dir, parser='spacy'):
    """Parse corpus, return vectors of SDP of each relation instance

    :param base_dir: directory containing the sentences
    :param parser:
    :return: labels (eid1, eid2), instances (vectors), classes (0/1), common ancestors, l/r ancestors, l/r wordnet
    """

    entities, positive_entities = get_sentence_entities()
    if parser == 'spacy':
        parsed_sentences, wordnet_sentences = parse_sentences_spacy(
            base_dir, entities)

    else:
        parsed_sentences = wordnet_sentences = None

    left_instances = []
    right_instances = []
    left_ancestors = []
    right_ancestors = []
    common_ancestors = []
    left_wordnet = []
    right_wordnet = []
    classes = []
    labels = []
    all_pos_gv = set()
    all_neg_gv = set()

    s = time.time()
    counter = 0

    print(len(entities_info.keys()))
    for k in entities_info.keys():  # 3500 train 10750 test

        q = multiprocessing.Queue()
        processes = []

        for sentence, entities_sentence in entities_per_abstract[k].items():
            p = Process(target=add_helper,
                        args=(
                            q,
                            base_dir,
                            parser,
                            entities,
                            positive_entities,
                            parsed_sentences,
                            wordnet_sentences,
                            k,
                            sentence,
                            entities_sentence,
                        ))
            logging.info((p.name, sentence[0]))
            processes.append(p)

        for proc in processes:
            proc.start()

        for proc in processes:
            try:
                out = q.get(True, 0.1)
                if out[0]:
                    labels += out[0]
                else:
                    labels += []

                if out[1]:
                    left_instances += out[1]
                else:
                    left_instances += []

                if out[2]:
                    right_instances += out[2]
                else:
                    right_instances += []

                if out[3]:
                    classes += out[3]
                else:
                    classes += []

                if out[4]:
                    common_ancestors += out[4]
                else:
                    common_ancestors += []

                if out[5]:
                    left_ancestors += out[5]
                else:
                    left_ancestors += []

                if out[6]:
                    right_ancestors += out[6]
                else:
                    right_ancestors += []

                if out[7]:
                    left_wordnet += out[7]
                else:
                    left_wordnet += []

                if out[8]:
                    right_wordnet += out[8]
                else:
                    right_wordnet += []

                if out[9]:
                    all_pos_gv.update(out[9])
                else:
                    all_pos_gv.update()

                if out[10]:
                    all_neg_gv.update(out[10])
                else:
                    all_neg_gv.update()

            except Exception as e:
                print(e)

        for proc in processes:
            proc.join()

        # TO BE SAFE
        counter += 1
        if temporary_directory == 'temp/':
            divider = 500
        elif temporary_directory == 'temp_dev/':
            divider = 250
        else:
            divider = 1075

        if counter % divider == 0:
            print('\n\n')
            print(counter)
            print('\n\n')
            train_labels, x_train, y_train, x_train_ancestors, x_train_subpaths, x_train_wordnet = labels, (
                left_instances, right_instances), classes, common_ancestors, (
                    left_ancestors, right_ancestors), (left_wordnet,
                                                       right_wordnet)

            if temporary_directory == 'temp/':
                preprocess_what = 'train'
            elif temporary_directory == 'temp_test/':
                preprocess_what = 'test'
            elif temporary_directory == 'temp_sample/':
                preprocess_what = 'sample'
            else:
                preprocess_what = 'dev'

            np.save(
                temporary_directory + 'drug_gene/' + preprocess_what +
                str(counter) + '_y_ck.npy', y_train)
            np.save(
                temporary_directory + 'drug_gene/' + preprocess_what +
                str(counter) + '_labels_ck.npy', train_labels)
            np.save(
                temporary_directory + 'drug_gene/' + preprocess_what +
                str(counter) + '_x_words_ck.npy', x_train)
            np.save(
                temporary_directory + 'drug_gene/' + preprocess_what +
                str(counter) + '_x_wordnet_ck.npy', x_train_wordnet)
            np.save(
                temporary_directory + 'drug_gene/' + preprocess_what +
                str(counter) + '_x_subpaths_ck.npy', x_train_subpaths)
            np.save(
                temporary_directory + 'drug_gene/' + preprocess_what +
                str(counter) + "_x_ancestors_ck.npy", x_train_ancestors)

    e = time.time()
    print(e - s)

    return labels, (left_instances,
                    right_instances), classes, common_ancestors, (
                        left_ancestors, right_ancestors), (
                            left_wordnet,
                            right_wordnet), all_neg_gv, all_pos_gv


def add_helper(q, base_dir, parser, entities, positive_entities,
               parsed_sentences, wordnet_sentences, k, sentence,
               entities_sentence):

    sentence_pairs = {}
    sentence_entities = entities_sentence
    # entity id; offset1; offset2; entity text; entity type.

    all_pairs = {}
    if k in relations_info:
        for pair in relations_info[k]:
            if pair[1] in all_pairs:
                all_pairs[pair[1]].append((pair[0], pair[2]))
            else:
                all_pairs[pair[1]] = []
                all_pairs[pair[1]].append((pair[0], pair[2]))

    comb = combinations([e[0] for e in sentence_entities], 2)

    for pair in comb:
        sentence_pairs[(sentence[0] + '.u' + pair[0], sentence[0] + '.u' +
                        pair[1])] = label_to_pair_type['NO_RELATION']

    for entity in sentence_entities:
        if entity[0] in all_pairs:
            for matching_entities in sentence_entities:
                for pair in all_pairs[entity[0]]:
                    if pair[1] == matching_entities[0]:
                        sentence_pairs[(sentence[0] + '.u' + entity[0],
                                        sentence[0] + '.u' + matching_entities[0])] = \
                            label_to_pair_type[pair[0]]

    if len(sentence_pairs) > 0:  # skip sentences without pairs
        sentence_entities = entities[sentence[0]]
        parsed_sentence = parsed_sentences[sentence[0]]
        wordnet_sentence = wordnet_sentences[sentence[0]]

        if parser == 'spacy':
            sentence_labels, sentence_we_instances, sentence_wn_instances, sentence_classes, pos_gv, neg_gv = \
                process_sentence_spacy(base_dir, parsed_sentence, sentence_entities, sentence_pairs,
                                       positive_entities, wordnet_sentence)
        else:
            sentence_labels = sentence_we_instances = sentence_classes = sentence_wn_instances = pos_gv = neg_gv = None
        sentence_ancestors, sentence_subpaths = get_ancestors(
            sentence_labels, sentence_entities)

        labels = sentence_labels
        left_instances = sentence_we_instances[0]
        right_instances = sentence_we_instances[1]
        classes = sentence_classes
        common_ancestors = sentence_ancestors
        left_ancestors = sentence_subpaths[0]
        right_ancestors = sentence_subpaths[1]
        left_wordnet = sentence_wn_instances[0]
        right_wordnet = sentence_wn_instances[1]

        q.put([
            labels, left_instances, right_instances, classes, common_ancestors,
            left_ancestors, right_ancestors, left_wordnet, right_wordnet,
            pos_gv, neg_gv
        ])

    return
