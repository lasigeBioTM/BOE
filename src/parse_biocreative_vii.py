from __future__ import unicode_literals, print_function
import sys
import json
import string
import time
import logging
# logger configuration
logging.basicConfig(filename="biont.log",
                    format='%(asctime)s %(message)s',
                    filemode='w')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
from multiprocessing import Process, JoinableQueue
from itertools import combinations, product
from subprocess import PIPE, Popen

import ontology_preprocessing

import networkx as nx
import pysbd
import en_core_web_sm
from spacy import Language
import numpy as np
#import debugpy
# Allow other computers to attach to debugpy at this IP address and port.
#debugpy.listen(("0.0.0.0", 5678))
# Pause the program until a remote debugger is attached
#debugpy.wait_for_client()
# Input Parameters
neg_gv_list = {}
sst_light_directory = 'bin/sst-light-0.4/'
temporary_directory = sys.argv[-1]
biomedical_entity_to_ontology = {'gene': 'G', 'drug': 'C'}
preprocessing_method = sys.argv[2]  # main or extra
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


def divided_by_sentences(abstract):
    seg = pysbd.Segmenter(language="en", clean=False)
    sentences = seg.segment(abstract)

    return sentences


def get_abstracts_info(abstract_file):
    abstracts_info_dict = {}
    with open(abstract_file, encoding='utf-8') as af:
        for line in af:
            line = line.strip()
            line = line.split('\t')
            abstracts_info_dict[line[0]] = [line[1]] + \
            divided_by_sentences(line[2])
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
    abstracts_info = get_abstracts_info('corpora/training/drugprot_training_abstracs.tsv')
    entities_info = get_entities_info('corpora/training/drugprot_training_entities.tsv')
    relations_info = get_relations_info('corpora/training/drugprot_training_relations.tsv')

    g_data = ontology_preprocessing.load_data('corpora/training/t_g_ents.txt')
    c_data = ontology_preprocessing.load_data('corpora/training/t_c_ents.txt')
    
    g_json = open('corpora/json/t_g_paths.json')
    g_paths = json.load(g_json)
    c_json = open('corpora/json/t_c_paths.json')
    c_paths = json.load(c_json)

elif temporary_directory == 'temp_dev/':
    abstracts_info = get_abstracts_info('corpora/development/drugprot_development_abstracs.tsv')
    entities_info = get_entities_info('corpora/development/drugprot_development_entities.tsv')
    relations_info = get_relations_info('corpora/development/drugprot_development_relations.tsv')

    g_data = ontology_preprocessing.load_data('corpora/development/dev_g_ents.txt')
    c_data = ontology_preprocessing.load_data('corpora/development/dev_c_ents.txt')
    
    g_json = open('corpora/json/dev_g_paths.json')
    g_paths = json.load(g_json)
    c_json = open('corpora/json/dev_c_paths.json')
    c_paths = json.load(c_json)


def get_overlapped_ents(entities_per_sentence):
    overlapped = set()
    overlapped_offset = 0
    ent = ()

    for previous, current in zip(entities_per_sentence, entities_per_sentence[1:]):
        if current[2] <= overlapped_offset and current != ent:
            overlapped.add(current)
        if previous[2] <= overlapped_offset and previous != ent:
            overlapped.add(previous)
        
        if current[1] <= previous[1] <= current[2] and \
            current[2] > overlapped_offset:
            overlapped.add(previous)
            overlapped_offset = current[2]
            ent = current
        
        elif previous[1] <= current[1] <= previous[2] and \
            previous[2] > overlapped_offset:
            overlapped.add(current)
            overlapped_offset = previous[2]
            ent = previous

    return overlapped


def get_new_offsets_sentences(abstract_id):
    entities_per_sentence = {}
    abstract = abstracts_info[abstract_id].copy()
    annotation_lines = entities_info[abstract_id].copy()

    sentence_id = 0
    sentence_limit_1 = 0
    sentence_limit_2 = 0
    for sentence in abstract:
        k = ('a' + abstract_id + '.s' + str(sentence_id), sentence)
        entities_per_sentence[k] = []
        sentence_limit_1 = sentence_limit_2
        sentence_limit_2 += len(sentence)
        for annotation in annotation_lines[:]:
            if sentence_limit_1 <= int(annotation[2]) and int(annotation[3]) <= sentence_limit_2:
                try:
                    i_ent = int(annotation[2]) - sentence_limit_1
                    i_ent =  i_ent - 1 if i_ent > 0 else 0
                    s_ent = sentence.index(annotation[4], i_ent, len(sentence))
                    e_ent = s_ent + int(annotation[3]) - int(annotation[2])
                    entities_per_sentence[k].append((annotation[0], s_ent, e_ent, annotation[4], annotation[1]))
                    annotation_lines.remove(annotation)
                
                except ValueError:
                    print(annotation, sentence, i_ent)
            else:
                break
        
        overlapped = list(get_overlapped_ents(entities_per_sentence[k]))
        if preprocessing_method == 'main':
            if overlapped:
                for ent in overlapped:
                    entities_per_sentence[k].remove(ent)

        elif preprocessing_method == 'extra':
            if overlapped:
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

    return entities_per_sentence

entities_per_abstract = {
    k: get_new_offsets_sentences(k)
    for k in entities_info.keys()
}


def get_sentence_entities():
    entities = {}
    positive_entities = set()

    for k in entities_info.keys():
        for sentence, entities_sentence in entities_per_abstract[k].items():
            sentence_entities = {}
            for entity in entities_sentence:
                if 'GENE' in entity[4]:
                    entity_id = g_data[entity[3]]
                elif 'CHEMICAL' in entity[4]:
                    entity_id = c_data[entity[3]]
                
                if k in relations_info:
                    for pair in relations_info[k]:
                        if entity[0] == pair[1]:
                            positive_entities.add(sentence[0] + '.u' + entity[0])
                        elif entity[0] == pair[2]:
                            positive_entities.add(sentence[0] + '.u' + entity[0])
                
                sentence_entities[sentence[0] + '.u' + entity[0]] = \
                    (eval('[' + str(entity[1]) + ', ' + str(entity[2]) + ']'), entity[3], entity_id)
                
            entities[sentence[0]] = sentence_entities

    return entities, positive_entities


def get_ancestors(sentence_labels, sentence_entities):
    right_paths = []
    left_paths = []

    for p in sentence_labels:
        left_ent = sentence_entities[p[0]][2]
        left_ent = left_ent.replace(':', '_')
        right_ent = sentence_entities[p[1]][2]
        right_ent = right_ent.replace(':', '_')

        if left_ent in g_paths:
            left_path = g_paths[left_ent]
            right_path = c_paths[right_ent]
        elif left_ent in c_paths:
            left_path = c_paths[left_ent]
            right_path = g_paths[right_ent]

        left_paths.append(left_path)
        right_paths.append(right_path)

    return (left_paths, right_paths)


@Language.component("psbd")
def prevent_sentence_segmentation(doc):
    # This will entirely disable spaCy's sentence detection
    for token in doc:
        token.is_sent_start = False

    return doc

nlp = en_core_web_sm.load(disable=['ner'])
nlp.add_pipe("psbd", name='prevent-sbd', before='parser')


def parse_sentence_spacy(sentence_text, sentence_entities):
    # Replace entities with string of a
    for e in sentence_entities:
        s_ent = sentence_entities[e][0][0]
        e_ent = sentence_entities[e][0][1]
        name_ent = sentence_entities[e][1]
        a_str = 'a'*len(name_ent)
        if s_ent != 0:
            sentence_text = sentence_text[:s_ent-1] + \
                            ' ' + a_str + ' ' + \
                            sentence_text[e_ent+1:]
        else:
            sentence_text = a_str + ' ' + \
                            sentence_text[e_ent+1:]
    # Clean text to make tokenization easier
    sentence_text = sentence_text.translate(str.maketrans(string.punctuation, ' '*len(string.punctuation)))
    # Replace entities in text
    for e in sentence_entities:
        s_ent = sentence_entities[e][0][0]
        e_ent = sentence_entities[e][0][1]
        name_ent = sentence_entities[e][1].translate(str.maketrans(string.punctuation, ' '*len(string.punctuation))).replace(' ', '_')
        if s_ent != 0:
            sentence_text = sentence_text[:s_ent] + name_ent + sentence_text[e_ent:]
        else:
            sentence_text = name_ent + sentence_text[e_ent:]
    # Use spacy to parse a sentence
    parsed = nlp(sentence_text)
    return parsed


def run_sst(token_seq):
    chunk_size = 500
    wordnet_tags = {}
    sent_ids = list(token_seq.keys())
    chunks = [sent_ids[i:i + chunk_size] for i in range(0, len(sent_ids), chunk_size)]
    for i, chunk in enumerate(chunks):
        with open('{}sentences_{}.txt'.format(temporary_directory, i), 'w', encoding='utf-8') as f:
            for sent in chunk:
                f.write("{}\t{}\t.\n".format(sent, '\t'.join(token_seq[sent])))

        sst_args = [
            sst_light_directory + 'sst', 
            'bitag',
            '{}MODELS/WSJPOSc_base_20'.format(sst_light_directory),
            '{}DATA/WSJPOSc.TAGSET'.format(sst_light_directory),
            '{}MODELS/SEM07_base_12'.format(sst_light_directory),
            '{}DATA/WNSS_07.TAGSET'.format(sst_light_directory),
            '{}sentences_{}.txt'.format(temporary_directory, i), 
            '0', 
            '0']

        p = Popen(sst_args, stdout=PIPE)
        p.communicate()

        with open('{}sentences_{}.txt.tags'.format(temporary_directory, i), encoding='utf-8') as f:
            output = f.read()

        sstoutput = parse_sst_results(output)
        wordnet_tags.update(sstoutput)

    return wordnet_tags


def parse_sst_results(results):
    sentences = {}
    lines = results.strip().split('\n')

    for l in lines:
        values = l.split('\t')
        wntags = [x.split(' ')[-1].split('-')[-1] for x in values[1:]]
        sentences[values[0]] = wntags

    return sentences


def parse_sentences_spacy(entities):
    # First iterate all documents, and preprocess all sentences
    parsed_sentences = {}
    token_seq = {}
    for k in entities_info.keys():
        for sentence, entities_sentence in entities_per_abstract[k].items():
            parsed_sentence = parse_sentence_spacy(sentence[1], entities[sentence[0]])
            parsed_sentences[sentence[0]] = parsed_sentence
            tokens = []
            for t in parsed_sentence:
                tokens.append(t.text.replace(' ', '_').replace('\t', '_').replace('\n', '_'))
            token_seq[sentence[0]] = tokens

    wordnet_tags = run_sst(token_seq)

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
        entity_tokens = sentence.char_span(offset[0], offset[1])

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


def process_sentence_spacy(sentence,
                           sentence_entities,
                           sentence_pairs,
                           positive_entities,
                           wordnet_tags=None,
                           mask_entities=True,
                           min_sdp_len=0,
                           max_sdp_len=15):
    """Process sentence to obtain labels, instances and classes for a ML classifier
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
    sentence_head_tokens_type_1, sentence_head_tokens_type_2, pos_gv, neg_gv = \
        get_head_tokens_spacy(sentence_entities, sentence, positive_entities)

    entity_offsets = [sentence_entities[x][0][0] for x in sentence_entities]

    for e in product(sentence_head_tokens_type_1.keys(), sentence_head_tokens_type_2.keys()):
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


def get_sdp_instances():
    s = time.time()
    entities, positive_entities = get_sentence_entities()
    parsed_sentences, wordnet_sentences = parse_sentences_spacy(entities)

    left_instances = []
    right_instances = []
    left_ancestors = []
    right_ancestors = []
    left_wordnet = []
    right_wordnet = []
    classes = []
    labels = []
    all_pos_gv = set()
    all_neg_gv = set()

    counter = 0
    print(len(entities_info.keys()))
    for k in entities_info.keys():
        q = JoinableQueue()
        processes = []
        for sentence, entities_sentence in entities_per_abstract[k].items():
            p = Process(target = add_helper,
                        args = (
                            q,
                            entities,
                            positive_entities,
                            parsed_sentences,
                            wordnet_sentences,
                            k,
                            sentence,
                            entities_sentence,
                        ))
            processes.append(p)

        for proc in processes:
            proc.start()

        for proc in processes:
            out = q.get()
            labels += out[0]
            left_instances += out[1]
            right_instances += out[2]
            classes += out[3]
            left_ancestors += out[4]
            right_ancestors += out[5]
            left_wordnet += out[6]
            right_wordnet += out[7]
            all_pos_gv.update(out[8])
            all_neg_gv.update(out[9])
            q.task_done()

        for proc in processes:
            proc.join()

        q.join()
        # TO BE SAFE
        counter += 1
        if temporary_directory == 'temp/':
            divider = 500
        elif temporary_directory == 'temp_dev/':
            divider = 250
        if counter % divider == 0:
            print('\n\n')
            print(counter)
            print('\n\n')
            train_labels, x_train, y_train, x_train_subpaths, x_train_wordnet = \
            labels, (left_instances, right_instances), classes, (left_ancestors, right_ancestors), (left_wordnet, right_wordnet)

            np.save(temporary_directory + str(counter) + '_y_ck.npy', y_train)
            np.save(temporary_directory + str(counter) + '_labels_ck.npy', train_labels)
            np.save(temporary_directory + str(counter) + '_x_words_ck.npy', x_train)
            np.save(temporary_directory + str(counter) + '_x_wordnet_ck.npy', x_train_wordnet)
            np.save(temporary_directory + str(counter) + '_x_subpaths_ck.npy', x_train_subpaths)

    e = time.time()
    print(e - s)

    return labels, (left_instances, right_instances), classes, (left_ancestors, right_ancestors), (left_wordnet, right_wordnet), all_neg_gv, all_pos_gv


def add_helper(q, entities, positive_entities, parsed_sentences, wordnet_sentences, k, sentence, sentence_entities):
    left_instances = []
    right_instances = []
    left_ancestors = []
    right_ancestors = []
    left_wordnet = []
    right_wordnet = []
    classes = []
    labels = []
    pos_gv = set()
    neg_gv = set()

    sentence_pairs = {}
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
        sentence_pairs[(sentence[0] + '.u' + pair[0], sentence[0] + '.u' + pair[1])] = label_to_pair_type['NO_RELATION']

    for entity in sentence_entities:
        if entity[0] in all_pairs:
            for matching_entities in sentence_entities:
                for pair in all_pairs[entity[0]]:
                    if pair[1] == matching_entities[0]:
                        sentence_pairs[(sentence[0] + '.u' + entity[0], sentence[0] + '.u' + matching_entities[0])] = label_to_pair_type[pair[0]]

    if len(sentence_pairs) > 0:  # skip sentences without pairs
        sentence_entities = entities[sentence[0]]
        parsed_sentence = parsed_sentences[sentence[0]]
        wordnet_sentence = wordnet_sentences[sentence[0]]
        
        sentence_labels, sentence_we_instances, sentence_wn_instances, sentence_classes, pos_gv, neg_gv = \
        process_sentence_spacy(parsed_sentence, sentence_entities, sentence_pairs, positive_entities, wordnet_sentence)
        
        sentence_subpaths = get_ancestors(sentence_labels, sentence_entities)

        labels = sentence_labels
        left_instances = sentence_we_instances[0]
        right_instances = sentence_we_instances[1]
        classes = sentence_classes
        left_ancestors = sentence_subpaths[0]
        right_ancestors = sentence_subpaths[1]
        left_wordnet = sentence_wn_instances[0]
        right_wordnet = sentence_wn_instances[1]

    q.put([
        labels, left_instances, right_instances, classes,
        left_ancestors, right_ancestors, left_wordnet, right_wordnet,
        pos_gv, neg_gv
    ])

    return
