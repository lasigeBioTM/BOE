# This script outputs, for each relation type, 
# the number of false negatives by relation type in the predictions. 
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

def get_false_negatives(predictions_file):
    pd_relations_d = get_relations(predictions_file)
    fn_relations_d = dict()
    for abstract, arg1_d in dev_relations_d.items():
        for arg1, arg2_d in arg1_d.items():
            for arg2, relations in arg2_d.items():
                for t in relations:
                    if abstract in pd_relations_d and arg1 in pd_relations_d[abstract] \
                        and arg2 in pd_relations_d[abstract][arg1] and t not in pd_relations_d[abstract][arg1][arg2]:
                        pd_t = pd_relations_d[abstract][arg1][arg2][0]
                        if t in fn_relations_d and pd_t in fn_relations_d[t]:
                            fn_relations_d[t][pd_t] += 1
                        elif t in fn_relations_d:
                            fn_relations_d[t][pd_t] = 1
                        else:
                            fn_relations_d[t] = {pd_t: 1}
    
    return fn_relations_d

def get_train_entities():
    entities_d = dict()
    for abstract, arg1_d in train_relations_d.items():
        for arg1, arg2_d in arg1_d.items():
            for arg2, relations in arg2_d.items():
                for t in relations:
                    if t not in entities_d:
                        entities_d[t] = {'chemical_names': set(),
                                        'chemical_ids': set(),
                                        'gene_names': set(),
                                        'gene_ids': set()}
                    c_name = train_entities_d[abstract][arg1]
                    c_id = chemicals_d[c_name]
                    g_name = train_entities_d[abstract][arg2]
                    g_id = genes_d[g_name]
                    entities_d[t]['chemical_names'].add(c_name)
                    entities_d[t]['chemical_ids'].add(c_id)
                    entities_d[t]['gene_names'].add(g_name)
                    entities_d[t]['gene_ids'].add(g_id)
    return entities_d

def get_train_relations(c_entities, g_entities, type, false_negative_type):
    relations_d = dict()
    for abstract, arg1_d in train_relations_d.items():
        for arg1, arg2_d in arg1_d.items():
            for arg2, relations in arg2_d.items():
                for t in relations:
                    if t in [type, false_negative_type]:
                        if t not in relations_d:
                            relations_d[t] = set()
                        c_name = train_entities_d[abstract][arg1]
                        c_id = chemicals_d[c_name]
                        g_name = train_entities_d[abstract][arg2]
                        g_id = genes_d[g_name]
                        condition = c_name in c_entities[0] or \
                                    g_name in g_entities[0] or \
                                    c_id in c_entities[1] or \
                                    g_id in g_entities[1]                                    
                        if condition:
                            relations_d[t].add(str(abstract + arg1 + arg2))
    return relations_d

def write_fn_stats(predictions, outfile):
    fn_relations_d = get_false_negatives(predictions)
    case_study = ['PART-OF', 'INDIRECT-DOWNREGULATOR', 'INDIRECT-UPREGULATOR', 'ACTIVATOR', 'AGONIST', 'PRODUCT-OF', 'AGONIST-ACTIVATOR', 'AGONIST-INHIBITOR', 'SUBSTRATE_PRODUCT-OF']
    train_entities_d = get_train_entities()
    for t, pd_t_d in fn_relations_d.items():
        if t in case_study:
            train_t_c = (train_entities_d[t]['chemical_names'], train_entities_d[t]['chemical_ids'])
            train_t_g = (train_entities_d[t]['gene_names'], train_entities_d[t]['gene_ids'])
            
            outfile.write('%s\tchem_names\t%d\tchem_ids\t%d\tgene_names\t%d\tgene_ids\t%d\t\n' % 
                    (t, len(train_t_c[0]), len(train_t_c[1]), len(train_t_g[0]), len(train_t_g[1])))
            
            share_c = (set(), set())
            share_g = (set(), set())
            common_relations_train_t = set()
            for t_fn, t_fn_count in pd_t_d.items():
                if t_fn != 'NO_RELATION':
                    train_t_fn_c = (train_entities_d[t_fn]['chemical_names'], train_entities_d[t_fn]['chemical_ids'])
                    train_t_fn_g = (train_entities_d[t_fn]['gene_names'], train_entities_d[t_fn]['gene_ids'])
                    train_common_chemical = (train_t_c[0]&train_t_fn_c[0], train_t_c[1]&train_t_fn_c[1])
                    train_common_gene = (train_t_g[0]&train_t_fn_g[0], train_t_g[1]&train_t_fn_g[1])
                    unique_chemical = (train_common_chemical[0]-share_c[0], train_common_chemical[1]-share_c[1])
                    unique_gene = (train_common_gene[0]-share_g[0], train_common_gene[1]-share_g[1])
                    share_c[0].update(train_common_chemical[0])
                    share_c[1].update(train_common_chemical[1])
                    share_g[0].update(train_common_gene[0])
                    share_g[1].update(train_common_gene[1])

                    outfile.write('%d\t%s\tchem_names\t%d\tchem_ids\t%d\tgene_names\t%d\tgene_ids\t%d\n' % 
                    (t_fn_count, t_fn, len(unique_chemical[0]), len(unique_chemical[1]), len(unique_gene[0]), len(unique_gene[1])))
                    
                    relations_train = get_train_relations(unique_chemical, unique_gene, t, t_fn)
                    unique_relations_train_t = relations_train[t] - common_relations_train_t
                    common_relations_train_t.update(relations_train[t])
                    
                    outfile.write('%s\trelations\t%d\t%s\trelations\t%d\n' % (t, len(unique_relations_train_t), t_fn, len(relations_train[t_fn])))
                    
            outfile.write('%s\tshared_chem_names\t%d\tshared_chem_ids\t%d\tshared_gene_names\t%d\tshared_gene_ids\t%d\t\n' % 
                    (t, len(share_c[0]), len(share_c[1]), len(share_g[0]), len(share_g[1])))
            outfile.write('\n')

dev_relations_d = get_relations(open('./corpora/development/drugprot_development_relations.tsv', 'r', encoding='utf-8'))
train_entities_d = get_entities(open('./corpora/training/drugprot_training_entities.tsv', 'r', encoding='utf-8'))
train_relations_d = get_relations(open('./corpora/training/drugprot_training_relations.tsv', 'r', encoding='utf-8'))

chemicals_d = dict(get_ontology(open('./corpora/training/t_c_ents.txt', 'r', encoding='utf-8')),
            **get_ontology(open('./corpora/development/dev_c_ents.txt', 'r', encoding='utf-8')))
genes_d = dict(get_ontology(open('./corpora/training/t_g_ents.txt', 'r', encoding='utf-8')),
            **get_ontology(open('./corpora/development/dev_g_ents.txt', 'r', encoding='utf-8')))

write_fn_stats(open(sys.argv[1], 'r', encoding='utf-8'), open(sys.argv[2], 'w', encoding='utf-8'))