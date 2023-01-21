# This script outputs, for each relation type, 
# the number of false positives by relation type in the predictions. 
# Then, for each relation type, 
# the number of entities present in both relation types by name and ontology id.
import sys

def get_entities(entities_file):
    gs_entities_d = dict()
    gs_lines = entities_file.readlines()
    entities_file.close()

    for line in gs_lines:
        line = line.strip()
        line = line.split('\t')
        abstract = line[0]
        arg = line[1]
        arg_name = line[5]
        if abstract in gs_entities_d:
            gs_entities_d[abstract][arg] = arg_name
        else:
            gs_entities_d[abstract] = {arg: arg_name}
    
    return gs_entities_d

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

def get_ontology(ontology_file):
    ontology_d = dict()
    file_lines = ontology_file.readlines()
    ontology_file.close()

    for line in file_lines:
        line = line.strip()
        line = line.split(' ')
        arg_name = ' '.join(line[:-1])
        id = line[-1].replace(':', '_')
        ontology_d[arg_name] = id
    
    return ontology_d

def get_fp_stats(entities_file, relations_file, predictions_file, chemicals_file, genes_file):
    gs_entities_d = get_entities(entities_file)
    gs_relations_d = get_relations(relations_file)
    chemicals_d = get_ontology(chemicals_file)
    genes_d = get_ontology(genes_file)
    pd_relations_d = get_relations(predictions_file)
    entities_stats_d = dict()
    relations_stats_d = dict()
    for abstract, arg1_d in gs_relations_d.items():
        for arg1, arg2_d in arg1_d.items():
            for arg2, relations in arg2_d.items():
                for t in relations:
                    c_name = gs_entities_d[abstract][arg1]
                    c_id = chemicals_d[c_name]
                    g_name = gs_entities_d[abstract][arg2]
                    g_id = genes_d[g_name]
                    if t not in entities_stats_d:
                        entities_stats_d[t]={'chemical_names': set(),
                                              'chemical_ids': set(),
                                              'gene_names': set(),
                                              'gene_ids': set()}
                    entities_stats_d[t]['chemical_names'].add(c_name)
                    entities_stats_d[t]['chemical_ids'].add(c_id)
                    entities_stats_d[t]['gene_names'].add(g_name)
                    entities_stats_d[t]['gene_ids'].add(g_id)

                    if abstract in pd_relations_d and arg1 in pd_relations_d[abstract] \
                        and arg2 in pd_relations_d[abstract][arg1] and t not in pd_relations_d[abstract][arg1][arg2]:
                        pd_t = pd_relations_d[abstract][arg1][arg2][0]
                        if t in relations_stats_d and pd_t in relations_stats_d[t]:
                            relations_stats_d[t][pd_t] += 1
                        elif t in relations_stats_d:
                            relations_stats_d[t][pd_t] = 1
                        else:
                            relations_stats_d[t] = {pd_t: 1}
    return entities_stats_d, relations_stats_d

def write_fp_stats(entities_file, relations_file, predictions_file, chemicals_file, genes_file, outfile):
    d = get_fp_stats(entities_file, relations_file, predictions_file, chemicals_file, genes_file)
    case_study = ["PRODUCT-OF", "ACTIVATOR", "INDIRECT-DOWNREGULATOR", "AGONIST", "INDIRECT-UPREGULATOR"]
    entities_stats_d = d[0]
    relations_stats_d = d[1]
    for t, pd_t_d in relations_stats_d.items():
        if t in case_study:
            t_c = (entities_stats_d[t]['chemical_names'], entities_stats_d[t]['chemical_ids'])
            t_g = (entities_stats_d[t]['gene_names'], entities_stats_d[t]['gene_ids'])
            outfile.write('%s\tchem_names\t%d\tchem_ids\t%d\tgene_names\t%d\tgene_ids\t%d\t\n' % 
                    (t, len(t_c[0]), len(t_c[1]), len(t_g[0]), len(t_g[1])))
            
            for t_fp, t_fp_count in pd_t_d.items():
                if t_fp != 'NO_RELATION':
                    t_fp_c = (entities_stats_d[t_fp]['chemical_names'], entities_stats_d[t_fp]['chemical_ids'])
                    t_fp_g = (entities_stats_d[t_fp]['gene_names'], entities_stats_d[t_fp]['gene_ids'])
                    outfile.write('%d\t%s\tchem_names\t%d\tchem_ids\t%d\tgene_names\t%d\tgene_ids\t%d\n' % 
                    (t_fp_count, t_fp, len(t_c[0]&t_fp_c[0]), len(t_c[1]&t_fp_c[1]), len(t_g[0]&t_fp_g[0]), len(t_g[1]&t_fp_g[1])))
            outfile.write('\n')

ent_file = open(sys.argv[1], 'r', encoding='utf-8')
rel_file = open(sys.argv[2], 'r', encoding='utf-8')
pred_file = open(sys.argv[3], 'r', encoding='utf-8')
c_file = open(sys.argv[4], 'r', encoding='utf-8')
g_file = open(sys.argv[5], 'r', encoding='utf-8')
outfile = open(sys.argv[6], 'w', encoding='utf-8')

write_fp_stats(ent_file, rel_file, pred_file, c_file, g_file, outfile)