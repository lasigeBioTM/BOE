# This script removes the extra relations in drugprot_training_relations.tsv and drugprot_development_relations.tsv
# The first relation found for a pair is the one selected.

import sys

relations = open(sys.argv[1], 'r', encoding='utf-8')
rel = relations.readlines()
rel_type = {}

for line in rel:
    line = line.strip()
    line = line.split('\t')
    if line[0] not in rel_type:
        rel_type[line[0]] = {}
        rel_type[line[0]][(line[2][5:], line[3][5:])] = line[1]
    else:
        if (line[2][5:], line[3][5:]) not in rel_type[line[0]]:
            rel_type[line[0]][(line[2][5:], line[3][5:])] = line[1]

new_relations = open(sys.argv[2], 'w', encoding='utf-8')

for key, value in rel_type.items():
    for k, v in value.items():
        new_relations.write('%s\t%s\t%s\t%s\n' %
                            (key, v, 'Arg1:' + k[0], 'Arg2:' + k[1]))

new_relations.close()
relations.close()