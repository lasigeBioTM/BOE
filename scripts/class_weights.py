# This script optimizes the weights of each label according with its size in the dataset
# It outputs the score for each label.

import numpy as np
import math
import sys


# labels_dict : {ind_label: count_label}
# mu : parameter to tune
def create_class_weight(labels_dict, mu=0.15):
    total = np.sum(list(labels_dict.values()))
    keys = labels_dict.keys()
    class_weight = dict()

    for key in keys:
        score = math.log(mu * total / float(labels_dict[key]))
        class_weight[key] = score if score > 1.0 else 1.0

    return class_weight


relations = open(sys.argv[1], 'r', encoding='utf-8')
rel = relations.readlines()
labels_dict = {}

for line in rel:
    line = line.strip()
    line = line.split('\t')
    if line[1] not in labels_dict:
        labels_dict[line[1]] = 1
    else:
        labels_dict[line[1]] += 1

print(create_class_weight(labels_dict))