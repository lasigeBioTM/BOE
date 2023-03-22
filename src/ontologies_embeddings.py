import sys
import random
import os
# Just disables the warning, doesn't enable AVX/FMA
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['CUDA_DEVICE_ORDER'] = 'PCI_BUS_ID'
os.environ['CUDA_VISIBLE_DEVICES'] = '0'

import ontology_preprocessing
from parse_biocreative_vii import get_sdp_instances
from models import get_model, n_classes, max_ancestors_length, max_sentence_length

import numpy as np
import matplotlib.pyplot as plt
from gensim.models.keyedvectors import KeyedVectors
from keras import utils
from keras.models import model_from_json
from keras.utils import pad_sequences
from keras.callbacks import ModelCheckpoint
# To run on GPU
import tensorflow as tf
from tensorflow.python.keras import backend as K
config = tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True  # dynamically grow the memory used on the GPU
K.set_session(tf.compat.v1.Session(config=config))
tf.compat.v1.disable_v2_behavior()
#from tensorflow.python.client import device_lib
#print(device_lib.list_local_devices())
# Input Parameters
n_epochs = 40
batch_size = 64
validation_split = 0.2
dict_type_entity_ontology = {
    'DRUG': 'ontology_preprocessing.load_chebi()',
    'GENE': 'ontology_preprocessing.load_go()'
}
temporary_directory = sys.argv[-1]
sst_light_directory = 'bin/sst-light-0.4/'
data_directory = 'data/'
models_directory = 'models/'
results_directory = 'results/'


def get_data():
    labels = []
    left_instances = []  # word indexes
    right_instances = []
    left_ancestors = []
    right_ancestors = []
    left_wordnet = []  # wordnet IDs
    right_wordnet = []
    all_pos_gv = set()  # anti positive governors
    all_neg_gv = set()
    classes = np.empty((0, ))

    dir_labels, dir_instances, dir_classes, dir_ancestors, dir_wordnet, \
    neg_gv, pos_gv = get_sdp_instances()

    dir_classes = np.array(dir_classes)
    labels += dir_labels
    left_instances += dir_instances[0]
    right_instances += dir_instances[1]
    left_ancestors += dir_ancestors[0]
    right_ancestors += dir_ancestors[1]
    left_wordnet += dir_wordnet[0]
    right_wordnet += dir_wordnet[1]
    classes = np.concatenate((classes, dir_classes), axis=0)

    all_pos_gv.update(pos_gv)
    all_neg_gv.update(neg_gv)

    return labels, (left_instances, right_instances), classes, (left_ancestors, right_ancestors), (left_wordnet, right_wordnet)


def preprocess_ids(x_data, id_to_index, max_len):
    data = []
    for i, seq in enumerate(x_data):
        idxs = []
        for d in seq:
            if d and d.startswith('CHEBI') or d.startswith('GO'):
                if d.replace('_', ':') not in id_to_index:
                    pass
                else:
                    idxs.append(id_to_index[d.replace('_', ':')])
        data.append(idxs)
    data = pad_sequences(data, maxlen=max_len)

    return data


def get_w2v(file_name='{}/PubMed-w2v.bin'.format(data_directory)):
    word_vectors = KeyedVectors.load_word2vec_format(file_name, binary=True)

    return word_vectors


def preprocess_sequences(x_data, embeddings_index):
    data = []
    for i, seq in enumerate(x_data):
        idxs = []
        for w in seq:
            if w.lower() in embeddings_index.vocab:
                idxs.append(embeddings_index.vocab[w.lower()].index)
        if None in idxs:
            print(seq, idxs)
        data.append(idxs)
    data = pad_sequences(data, maxlen=max_sentence_length, padding='post')

    return data


def get_wordnet_indexes():
    embedding_indexes = {}
    with open('{}/DATA/WNSS_07.TAGSET'.format(sst_light_directory), 'r') as f:
        lines = f.readlines()
        i = 0
        for l in lines:
            if l.startswith('I-'):
                continue
            embedding_indexes[l.strip().split('-')[-1]] = i
            i += 1

    return embedding_indexes


def preprocess_sequences_glove(x_data, embeddings_index):
    data = []
    for i, seq in enumerate(x_data):
        idxs = [embeddings_index.get(w) for w in seq if w in embeddings_index]
        if None in idxs:
            print(seq, idxs)
        data.append(idxs)
    data = pad_sequences(data, maxlen=max_sentence_length)

    return data


def concatenation_ancestors(x_subpaths_train, list_order):
    is_a_graph, name_to_id, synonym_to_id, id_to_name, id_to_index_entity_left = \
            eval(dict_type_entity_ontology['DRUG'])
    is_a_graph, name_to_id, synonym_to_id, id_to_name, id_to_index_entity_right = \
            eval(dict_type_entity_ontology['GENE'])
    id_to_index = {**id_to_index_entity_left, **id_to_index_entity_right}

    x_ids_left = preprocess_ids(x_subpaths_train[0], id_to_index, max_ancestors_length)
    x_ids_right = preprocess_ids(x_subpaths_train[1], id_to_index, max_ancestors_length)

    return x_ids_left[list_order], x_ids_right[list_order], id_to_index


def write_plots(history, model_name):
    plt.figure()
    plt.plot(history.history['f1'])
    plt.plot(history.history['val_f1'])
    plt.title('model eval')
    plt.ylabel('score')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    plt.savefig(results_directory + '/' + '{}_acc.png'.format(model_name))

    plt.figure()
    plt.plot(history.history['loss'])
    plt.plot(history.history['val_loss'])
    plt.title('model loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
    plt.savefig(results_directory + '/' + '{}_loss.png'.format(model_name))


def join_channels(model_name, channels, y_train, train_labels, n_classes, 
                  x_words_train, x_wordnet_train, x_subpaths_train, id_to_index):
    # Remove and replace previous model files
    if os.path.isfile('{}/{}.json'.format(models_directory, model_name)):
        os.remove('{}/{}.json'.format(models_directory, model_name))

    if os.path.isfile('{}/{}.h5'.format(models_directory, model_name)):
        os.remove('{}/{}.h5'.format(models_directory, model_name))

    y_train = utils.to_categorical(y_train, num_classes=n_classes)
    list_order = np.arange(len(y_train))
    print()
    print(list_order)

    random.seed(1)
    random.shuffle(list_order)
    y_train = y_train[list_order]

    train_labels = train_labels[list_order]
    print('\nTraining order:', list_order)

    inputs = {}

    if 'words' in channels:  ##### 1ST CHANNEL
        word_vectors = get_w2v()
        w2v_layer = word_vectors.get_keras_embedding(train_embeddings=False)
        x_words_left = preprocess_sequences(x_words_train[0], word_vectors)
        x_words_right = preprocess_sequences(x_words_train[1], word_vectors)
        del word_vectors
        inputs['left_words'] = x_words_left[list_order]
        inputs['right_words'] = x_words_right[list_order]
    else:
        w2v_layer = None

    if 'wordnet' in channels:  ##### 2ND CHANNEL
        wn_index = get_wordnet_indexes()
        x_wn_left = preprocess_sequences_glove(x_wordnet_train[0], wn_index)
        x_wn_right = preprocess_sequences_glove(x_wordnet_train[1], wn_index)
        inputs['left_wordnet'] = x_wn_left[list_order]
        inputs['right_wordnet'] = x_wn_right[list_order]
    else:
        wn_index = None

    if 'concatenation_ancestors' in channels:  ##### 3RD CHANNEL
        x_left, x_right, id_to_index = concatenation_ancestors(x_subpaths_train, list_order)
        inputs['left_ancestors'] = x_left
        inputs['right_ancestors'] = x_right
    else:
        id_to_index = None
    # Model
    model = get_model(w2v_layer, channels, wn_index, id_to_index)
    del w2v_layer
    del wn_index
    del id_to_index
    # Serialize model to JSON
    model_json = model.to_json()
    with open('{}/{}.json'.format(models_directory, model_name), 'w') as json_file:
        json_file.write(model_json)

    checkpoint = ModelCheckpoint(filepath='{}/{}.h5'.format(models_directory, model_name),
                                 verbose=1,
                                 save_best_only=True)
    class_weight = {
        0: 0.5,
        1: 1.0,
        2: 1.7670781432557048,
        3: 1.0,
        4: 1.2890649159309184,
        5: 1.3608608726449865,
        6: 1.6744392894003468,
        7: 1.3246812160674841,
        8: 2.0630715593582787,
        9: 1.0,
        10: 1.728335057605479,
        11: 5.186499263874312,
        12: 5.9888457363992496,
        13: 5.334919268992585
    }
    
    history = model.fit(inputs, {'output': y_train},
                        validation_split=validation_split,
                        epochs=n_epochs,
                        batch_size=batch_size,
                        verbose=2,
                        callbacks=[checkpoint],
                        class_weight=class_weight)

    write_plots(history, model_name)
    print('Saved model to disk.')


def predict(model_name, channels, test_labels, x_words_test, x_wn_test, x_subpaths_test, id_to_index):
    inputs = {}
    if 'words' in channels:  ##### 1ST CHANNEL
        word_vectors = get_w2v()
        x_words_test_left = preprocess_sequences([['entity'] + x[1:] for x in x_words_test[0]], word_vectors)
        x_words_test_right = preprocess_sequences([x[:-1] + ['entity'] for x in x_words_test[1]], word_vectors)
        inputs['left_words'] = x_words_test_left
        inputs['right_words'] = x_words_test_right
        del word_vectors

    if 'wordnet' in channels:  ##### 2ND CHANNEL
        wn_index = get_wordnet_indexes()
        x_wordnet_test_left = preprocess_sequences_glove(x_wn_test[0], wn_index)
        x_wordnet_test_right = preprocess_sequences_glove(x_wn_test[1], wn_index)
        inputs['left_wordnet'] = x_wordnet_test_left
        inputs['right_wordnet'] = x_wordnet_test_right
        del wn_index

    if 'concatenation_ancestors' in channels:  ##### 3RD CHANNEL
        x_ids_left = preprocess_ids(x_subpaths_test[0], id_to_index, max_ancestors_length)
        x_ids_right = preprocess_ids(x_subpaths_test[1], id_to_index, max_ancestors_length)
        inputs['left_ancestors'] = x_ids_left
        inputs['right_ancestors'] = x_ids_right
        del id_to_index
    # Load JSON
    json_file = open('{}/{}.json'.format(models_directory, model_name), 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(loaded_model_json)
    # Load weights
    loaded_model.load_weights('{}/{}.h5'.format(models_directory, model_name))

    print('Loaded model {}/{} from disk.'.format(models_directory, model_name))

    scores = loaded_model.predict(inputs)
    # Write results to file
    from parse_biocreative_vii import pair_type_to_label

    with open(results_directory + '/{}_results.txt'.format(model_name), 'w') as f:
        f.write('\t'.join(['Entity_1', 'Entity_2', 'Predicted_class\n']))
        for i, pair in enumerate(test_labels):
            f.write('\t'.join((pair[0], pair[1], pair_type_to_label[(np.argmax(scores[i]))] + '\n')))


def main():
    type_of_action = sys.argv[1]  # preprocess, train or test
    ##### PREPROCESS #####
    if type_of_action == 'preprocess':
        train_labels, x_train, y_train, x_train_subpaths, x_train_wordnet = get_data()
        np.save(temporary_directory + 'y.npy', y_train)
        np.save(temporary_directory + 'labels.npy', train_labels)
        np.save(temporary_directory + 'x_words.npy', x_train)
        np.save(temporary_directory + 'x_wordnet.npy', x_train_wordnet)
        np.save(temporary_directory + 'x_subpaths.npy', x_train_subpaths)

    ##### TRAIN #####
    elif type_of_action == 'train':
        model_name = sys.argv[2]  # model_1, model_2 etc.
        channels = sys.argv[3:-1]  # channels to use or string with only one channel
        y_train = np.load(temporary_directory + 'y.npy', allow_pickle=True)
        train_labels = np.load(temporary_directory + 'labels.npy', allow_pickle=True)
        x_words_train = np.load(temporary_directory + 'x_words.npy', allow_pickle=True)
        x_wordnet_train = np.load(temporary_directory + 'x_wordnet.npy', allow_pickle=True)
        x_subpaths_train = np.load(temporary_directory + 'x_subpaths.npy', allow_pickle=True)

        join_channels(model_name, channels, y_train, train_labels, n_classes, x_words_train, 
                      x_wordnet_train, x_subpaths_train, max_ancestors_length)

    ##### TEST / PREDICT #####
    elif type_of_action == 'test':
        model_name = sys.argv[2]  # model name
        channels = sys.argv[3:-1]  # channels to use

        id_to_index = {}
        is_a_graph_chebi, name_to_id_chebi, synonym_to_id_chebi, id_to_name_chebi, id_to_index_chebi = \
        ontology_preprocessing.load_chebi()
        id_to_index = {**id_to_index, **id_to_index_chebi}

        is_a_graph_go, name_to_id_go, synonym_to_id_go, id_to_name_go, id_to_index_go = \
        ontology_preprocessing.load_go()
        id_to_index = {**id_to_index, **id_to_index_go}

        test_labels = np.load(temporary_directory + 'labels.npy', allow_pickle=True)
        x_words_test = np.load(temporary_directory + 'x_words.npy', allow_pickle=True)
        x_wordnet_test = np.load(temporary_directory + 'x_wordnet.npy', allow_pickle=True)
        x_subpaths_test = np.load(temporary_directory + 'x_subpaths.npy', allow_pickle=True)
        predict(model_name, channels, test_labels, x_words_test, x_wordnet_test, x_subpaths_test, id_to_index)

        del id_to_index

    else:
        print('The type of action was not properly defined, it has to be preprocess, train or test.')


if __name__ == '__main__':
    main()
