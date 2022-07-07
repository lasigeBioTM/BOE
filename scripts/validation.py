# This script converts BiOnt's result in a friendly format for evaluation.

import sys

results = open(sys.argv[1], 'r', encoding='utf-8')
results.readline()  # skip header
res = results.readlines()
results.close()

res_rel = {}
for line in res:
    line = line.strip()
    line = line.split('\t')

    arg1 = line[0].split('.')
    arg2 = line[1].split('.')
    rel_type = line[2]

    if arg1[0][1:] not in res_rel:
        res_rel[arg1[0][1:]] = {}

    for i in range(len(arg1)):
        if 'uT' in arg1[i]:
            ent1 = arg1[i]

    for i in range(len(arg2)):
        if 'uT' in arg2[i]:
            ent2 = arg2[i]

    res_rel[arg1[0][1:]][(ent1[1:], ent2[1:])] = rel_type
print('Results processed')

outfile = open(sys.argv[2], 'w', encoding='utf-8')

for key, value in res_rel.items():
    for k, v in value.items():
        outfile.write('%s\t%s\t%s\t%s\n' %
                      (key, v, 'Arg1:' + k[0], 'Arg2:' + k[1]))
outfile.close()