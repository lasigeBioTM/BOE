# This script takes the gold standard relations and logs the overlapped relations.
# Logs the number of overlapped entities, the number of relations with these entities
# and the number of predicted relations.
import sys
import itertools
from collections import Counter


def get_relations(relations_file):
    gs_relations_d = dict()
    gs_lines = relations_file.readlines()
    relations_file.close()

    for line in gs_lines:
        line = line.strip()
        line = line.split('\t')
        abstract = line[0]
        relation_type = line[1]
        arg1 = line[2][5:]
        arg2 = line[3][5:]
        if abstract in gs_relations_d and arg1 in gs_relations_d[abstract] \
            and arg2 in gs_relations_d[abstract][arg1]:
            gs_relations_d[abstract][arg1][arg2].append(relation_type)
        elif abstract in gs_relations_d and arg1 in gs_relations_d[abstract]:
            gs_relations_d[abstract][arg1][arg2] = [relation_type]
        elif abstract in gs_relations_d:
            gs_relations_d[abstract][arg1] = {arg2: [relation_type]}
        else:
            gs_relations_d[abstract] = {arg1: {arg2: [relation_type]}}

    return gs_relations_d

def get_entities(entities_file):
    gs_entities_d = dict()
    gs_lines = entities_file.readlines()
    entities_file.close()

    for line in gs_lines:
        line = line.strip()
        line = line.split('\t')
        abstract = line[0]
        arg = line[1]
        arg_type = line[2]
        pos_s = line[3]
        pos_e = line[4]
        if abstract in gs_entities_d and (pos_s, pos_e) in gs_entities_d[abstract]:
            gs_entities_d[abstract][(pos_s, pos_e)].append((arg_type, arg))
        elif abstract in gs_entities_d:
            gs_entities_d[abstract][(pos_s, pos_e)] = [(arg_type, arg)]
        else:
            gs_entities_d[abstract] = {(pos_s, pos_e): [(arg_type, arg)]}
    
    return gs_entities_d

def get_overlapped(entities_abstract):
    overlapped = set()
    overlapped_offset = 0
    entity = ()

    for previous, current in zip(entities_abstract, entities_abstract[1:]):
        if current[1] <= overlapped_offset and current != entity:
            overlapped.add(current)
        if previous[1] <= overlapped_offset and previous != entity:
            overlapped.add(previous)

        if current[0] <= previous[0] <= current[1] and \
        current[1] > overlapped_offset:
            overlapped.add(previous)
            overlapped_offset = current[1]
            entity = current
        elif previous[0] <= current[0] <= previous[1] and \
        previous[1] > overlapped_offset:
            overlapped.add(current)
            overlapped_offset = previous[1]
            entity = previous

        if entity:
            overlapped.add(entity)

    return overlapped

def get_overlapped_entities_relations(gs_relations_d, gs_entities_d, abstract, overlapped_positions):
    overlapped_entities = list(itertools.chain(*[gs_entities_d[abstract][str(p[0]), str(p[1])] for p in overlapped_positions]))
    abstract_relations = [(arg1, arg2, relations) for arg1, arg2_d in gs_relations_d[abstract].items() for arg2, relations in arg2_d.items()]
    overlapped_entities_relations_d = dict()
    relations_type_d = dict()
    multiple_relations = 0

    for ent in overlapped_entities:
        arg = ent[1]
        for r in abstract_relations[:]:
            arg1 = r[0]
            arg2 = r[1]
            relations = r[2]
            if arg in r:
                if abstract in overlapped_entities_relations_d \
                    and arg1 in overlapped_entities_relations_d[abstract]:
                    overlapped_entities_relations_d[abstract][arg1][arg2] = relations
                elif abstract in overlapped_entities_relations_d:
                    overlapped_entities_relations_d[abstract][arg1] = {arg2: relations}
                else:
                    overlapped_entities_relations_d[abstract] = {arg1: {arg2: relations}}

                if len(relations) > 1:
                    multiple_relations += 1
                
                for t in relations:
                    if t in relations_type_d:
                        relations_type_d[t] += 1
                    else:
                        relations_type_d[t] = 1
                abstract_relations.remove(r)

    return overlapped_entities_relations_d, relations_type_d, multiple_relations

def get_overlap_stats(pd_relations_d, overlapped_entities_relations_d):
    predictions_type_d = dict()
    for abstract in overlapped_entities_relations_d:
        if abstract in pd_relations_d:
            for arg1, arg2_d in overlapped_entities_relations_d[abstract].items():
                if arg1 in pd_relations_d[abstract]:
                    for arg2, relations in arg2_d.items():
                        if arg2 in pd_relations_d[abstract][arg1]:
                            for rel in relations:
                                if rel in pd_relations_d[abstract][arg1][arg2]:
                                    print(abstract, arg1, arg2, relations)
                                    if rel in predictions_type_d:
                                        predictions_type_d[rel] += 1
                                    else:
                                        predictions_type_d[rel] = 1
    return predictions_type_d

def write_overlap_stats(relations_file, predictions_file, entities_file, outfile):
    gs_relations_d = get_relations(relations_file)
    gs_entities_d = get_entities(entities_file)
    abstracts = 0
    overlapped_entities = 0
    overlapped_entities_relations_d = dict()
    relations_type_d = dict()
    multiple_relations = 0

    for abstract in gs_relations_d:
        entities_position = sorted([(int(p[0]), int(p[1])) for p in list(gs_entities_d[abstract].keys())])
        overlapped = get_overlapped(entities_position)
        if overlapped:
            abstracts += 1
            overlapped_entities += (len(overlapped))
            d = get_overlapped_entities_relations(gs_relations_d, gs_entities_d, abstract, overlapped)
            overlapped_entities_relations_d.update(d[0])
            relations_type_d = Counter(relations_type_d) + Counter(d[1])
            multiple_relations += d[2]

    pd_relations_d = get_relations(predictions_file)
    predictions_type_d = dict.fromkeys(relations_type_d.keys(), 0)
    predictions_type_d.update(get_overlap_stats(pd_relations_d, overlapped_entities_relations_d))

    outfile.write('%d Abstracts w/ Relations & Overlapped Entities\n' % (abstracts))
    outfile.write('%d Abstracts w/ Overlapped Entities w/ Relations\n\n' % (len(overlapped_entities_relations_d.keys())))
    outfile.write('%d Overlapped Entities in Abstracts w/ Relations\n\n' % (overlapped_entities))
    outfile.write('%d Overlapped Relations\n' % (sum(relations_type_d.values())))
    outfile.write('%d Multiple Overlapped Relations w/ Overlapped Entities\n' % (multiple_relations))

    for relation_type, number in relations_type_d.items():
        outfile.write('%s\t%d\t%d\t\n' % (relation_type, number, predictions_type_d[relation_type]))

def write_multiple_relations(relations_file, outfile):
    gs_relations_d = get_relations(relations_file)

    for abstract, arg1 in gs_relations_d.items():
        for arg1, arg2 in arg1.items():
            for arg2, relations in arg2.items():
                if len(relations) > 1:
                    outfile.write('%s\t%s\t%s\t%s\n' %
                                    (abstract, arg1, arg2, relations))


relations_file = open(sys.argv[1], 'r', encoding='utf-8')
predictions_file = open(sys.argv[2], 'r', encoding='utf-8')
entities_file = open(sys.argv[3], 'r', encoding='utf-8')
outfile_stats = open(sys.argv[4], 'w', encoding='utf-8')
write_overlap_stats(relations_file, predictions_file, entities_file, outfile_stats)

#outfile_relations = open(sys.argv[2], 'w', encoding='utf-8')
#write_multiple_relations(relations_file, outfile_relations)
