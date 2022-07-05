#This script outputs a file with the PubMed ids.

import sys

pmids = open(sys.argv[1], 'r', encoding='utf-8')
lines = pmids.readlines()
pmids.close()
out_file = open(sys.argv[2], 'w')
for line in lines:
    line = line.strip()
    line = line.split('\t')
    out_file.write('%s\n' % (line[0]))
out_file.close()