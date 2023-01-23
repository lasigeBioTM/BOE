'''
import requests
import xmltodict, json
import json

query = "http://mygene.info/v3/query?q=symbol:"
d = {}
with open('universe.txt', mode='r') as in_file, \
     open('universeid.txt', mode='w') as out_file:
    for g in in_file:
        response = requests.get(query+g)
        json_data = json.loads(response.text)
        if json_data['hits'] and json_data['hits'][0]['_id'].isdecimal():
            d[g.rstrip()] = json_data['hits'][0]['_id']
        else: 
            response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}[sym]+AND+human[ORGN]".format(g))
            json_data = xmltodict.parse(response.text)
            print(g.rstrip(), json_data['eSearchResult']['IdList']['Id'])
            if json_data['eSearchResult']['IdList']['Id']:
                d[g.rstrip()] = json_data['eSearchResult']['IdList']['Id']
            else:
                d[g.rstrip()] = 'NULL'
    for key, value in d.items(): 
        out_file.write('%s %s\n' % (key, value))
    in_file.close()
    out_file.close()
'''
'''
query = "http://mygene.info/v3/query?q=symbol:"
d = {}
with open('universeid.txt', mode='r') as in_file, \
     open('gene_symbol_thesaurus.txt', mode='r') as genes_file, \
     open('saurus.txt', mode='w') as out_file:
    for g in in_file:
        g = g.rstrip()
        g = g.split(' ')
        d[g[0]] = g[1]
    next(genes_file)
    for g in genes_file:
        g = g.rstrip()
        g = g.split('\t')
        if g[0] not in d:
            response = requests.get(query+g[0])
            print(response.status_code)
            while response.status_code != 200:
                response = requests.get(query+g[0])
            json_data = json.loads(response.text)
            if json_data['hits'] and json_data['hits'][0]['_id'].isdecimal():
                d[g[0]] = json_data['hits'][0]['_id']
            else: 
                response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}[sym]+AND+human[ORGN]".format(g[0]))
                while response.status_code != 200:
                    response = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}[sym]+AND+human[ORGN]".format(g[0]))
                json_data = xmltodict.parse(response.text)
                if json_data['eSearchResult']['IdList']:
                    d[g[0]] = json_data['eSearchResult']['IdList']['Id']
                    print(g[0], json_data['eSearchResult']['IdList']['Id'])
                else:
                    d[g[0]] = 'NULL'
    for key, value in d.items(): 
        out_file.write('%s %s\n' % (key, value))
    in_file.close()
    genes_file.close()
    out_file.close()
'''
'''
d = {}
with open('saurus.txt', mode='r') as in_file, \
     open('thesaurus.txt', mode='w') as out_file:
    for g in in_file:
        g = g.rstrip()
        g = g.split(' ')
        if '[' in g[1]:
            g[1] = g[1][1:-1].split(', ')
            g[1] = g[1][0][1:-1]
        d[g[0]] = g[1]
    for key, value in d.items(): 
        out_file.write('%s %s\n' % (key, value))
    in_file.close()
    out_file.close()
'''
'''
d = {}
with open('thesaurus.txt', mode='r') as in_file, \
     open('gene_symbol_thesaurus.txt', mode='r') as genes_file, \
     open('genethesaurus.txt', mode='w') as out_file:
    for g in in_file:
        g = g.rstrip()
        g = g.split(' ')
        d[g[0]] = g[1]
    next(genes_file)
    for g in genes_file:
        g = g.rstrip()
        g = g.split('\t')
        if g[1] != 'NA':
            g[1] = g[1].split(',')
            for i in g[1]:
                d[i] = d.get(g[0])
    for key, value in d.items(): 
        out_file.write('%s %s\n' % (key, value))
    in_file.close()
    genes_file.close()
    out_file.close()
'''
import os
import collections

def dict_g2go(file_g2go):
    """Creates a dictionary of type {gene1:[(GO_ID, Evidence, GO_name, category),
    (GO_ID, Evidence, GO_name, category), ...], }
    :param file_g2go: file with relations gene to GO
    :return: dict of type {gene1:[(GO_ID, Evidence, GO_name, category),
             (GO_ID, Evidence, GO_name, category), ...], }
    """
    gene2go = open(file_g2go, 'r', encoding = 'utf-8')
    gene2go.readline()  # skip header
    relations_g2go = gene2go.readlines()
    gene2go.close()
    relations_g2go.pop()
    dict_gene_go = {}

    for line in relations_g2go:
        line = line.split('\t')
        gene_id = line[1]
        go = line[2]
        evidence = line[3]
        name = line[5]
        category = line[7][:-1]

        if gene_id not in dict_gene_go and category == 'Process':
            dict_gene_go[gene_id] = []
            dict_gene_go[gene_id].append((go, evidence, name, category))
        else:
            if category == 'Process':
                dict_gene_go[gene_id].append((go, evidence, name, category))

    return dict_gene_go




"""Generates a file for each abstract with the correspondent phenotype and GO annotations,
    creates a dictionary of type {gene_id:go_id, gene_id2:go_id, } and
    a dictionary of type {gene_name:go_name, gene_name:go_name, }
:param annotations_path: divided by sentences annotations path
:param file_g2go: file with relations gene to GO
:param destination_path: destination path
:return: file for each abstract with the correspondent phenotype and GO annotations,
            creates a dictionary of type {gene_id:go_id, gene_id2:go_id, } and
            a dictionary of type {gene_name:go_name, gene_name:go_name, }
            annotation file example:
            26 29  negative regulation of cell proliferation   GO_0008285
            279	288	bilateral	HP_0012832
            313	323	unilateral	HP_0012833
"""
dict_gene_id_go = dict_g2go('C:\\Users\\rcass\\Documents\\C_FCUL\\scripts\\gene2go')
dict_gene_go_id = {}
with open('C:\\Users\\rcass\\Documents\\C_FCUL\\scripts\\genethesaurus.txt', mode='r') as in_file, \
open('C:\\Users\\rcass\\Documents\\C_FCUL\\scripts\\genego.txt', mode='w') as out_file:
    for line in in_file:
        line = line.rstrip()
        line = line.split(' ')
        if line[1] in dict_gene_id_go:
            #for value_tup in dict_gene_id_go[line[1]]:
            list_evidence = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP',
                                'HDA', 'HMP', 'HGI', 'HEP', 'ISS', 'ISO', 'ISA',
                                'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA',
                                'TAS', 'NAS', 'IC', 'ND', 'IEA']  # order criteria

            d = {}
            d[line[1]] = {}
            for value in dict_gene_id_go[line[1]]:
                if value[1] not in d[line[1]]:
                    d[line[1]][value[1]] = [value]
                else:
                    d[line[1]][value[1]].append(value)

            for v in list_evidence:
                if v in d[line[1]]:
                    dict_gene_go_id[line[0]] = d[line[1]][v][0][0].replace(':', '_')
                    break

    for key, value in dict_gene_go_id.items(): 
        out_file.write('%s %s\n' % (key, value))
    in_file.close()
    out_file.close()
