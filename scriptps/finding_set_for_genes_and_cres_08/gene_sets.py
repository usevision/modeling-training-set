#!/usr/bin/env python3

import argparse as ap
import numpy as np

'''Usage: ./gene_sets.py --genes ~/mm10_genome/gencode.vM4.genes.gtf --train train_genes.bed --test test_genes.bed --ref referenceLoci_genes.bed'''

parser = ap.ArgumentParser(description='go through list of genes in train/test/reference set and make sure that each gene is in a unique set')
parser.add_argument('--genes', action='store', nargs=1, type=str, required=True, help='gtf file gencode.vM4.genes.gtf with only genes')
parser.add_argument('--train', action='store', nargs=1, type=str, required = True, help='bed file with train.bed intersected with gencode.vM4.genes.gtf')
parser.add_argument('--test', action='store', nargs=1, type=str, required = True, help='bed file with test.bed intersected with gencode.vM4.genes.gtf')
parser.add_argument('--ref', action='store', nargs=1, type=str, required=True, help='bed file with referenceLoci_fullRegion.bed intersected with gencode.vM4.genes.gtf')
args=parser.parse_args()
genes_file = args.genes[0]
train_file = args.train[0]
test_file = args.test[0]
ref_file = args.ref[0]

set_dict = {'train':0,
            'test':1,
            'ref':2}

'''add genes to a dictionary to track what sets it overlaps, the genomic location, and its common name'''
geneIDs = {}
for line in open(genes_file):
    fields=line.strip('\r\n').split('\t')
    gene_ID = fields[8].split(";")[0].split(" ")[1].replace('"','')
    common_name = fields[8].split(";")[4].split(" ")[-1].replace('"','')
    chr = fields[0]
    start = fields[3]
    end = fields[4]
    if chr != 'chrM':
        if gene_ID not in geneIDs:
            geneIDs[gene_ID]={'sets':[], 'loc':(chr, start, end), 'common_name':common_name} #checked manually that the common name for any gene that overlapped multiple sets where one of those sets was the reference set was not a reference loci of interest

'''genes that overlap training set'''
setOI = 'train'
for line in open(train_file):
    fields=line.strip('\r\n').split('\t')
    gene_ID_tr = fields[11].split(";")[0].split(" ")[1].replace('"','')
    geneIDs[gene_ID_tr]['sets'].append(set_dict[setOI])

'''genes that overlap test set'''
setOI = 'test'
for line in open(test_file):
    fields=line.strip('\r\n').split('\t')
    gene_ID_te = fields[11].split(";")[0].split(" ")[1].replace('"','')
    geneIDs[gene_ID_te]['sets'].append(set_dict[setOI])

'''genes that overlap the ref set'''
setOI = 'ref'
for line in open(ref_file):
    fields=line.strip('\r\n').split('\t')
    gene_ID_ref = fields[12].split(";")[0].split(" ")[1].replace('"','')
    geneIDs[gene_ID_ref]['sets'].append(set_dict[setOI])

'''if a gene overlaps more than one set'''
for key in geneIDs:
    uniq = set(geneIDs[key]['sets']) #gene could be in multiple windows within the same set; this together with line 66-67 will reduce those
    uniq_len = len(uniq)
    if uniq_len > 1:
        uniq = list(uniq)
        uniq_index = np.random.randint(uniq_len) #get a random index less than the length of the list of sets
        assigned_set = uniq[uniq_index] #assign a set based on this random index
        geneIDs[key]['sets']=[assigned_set]
    else:
        geneIDs[key]['sets'] = list(uniq)

test_gene_file = open('test_genes_uniq.bed','w+')
train_gene_file = open('train_genes_uniq.bed','w+')
ref_gene_file = open('ref_genes_uniq.bed','w+')

'''go through the dictionary and write the gene to whatever set it belongs to'''
for key in geneIDs:
    chr, start, stop = geneIDs[key]['loc']
    if geneIDs[key]['sets'][0] == 0:
        train_gene_file.write(chr + '\t' + start + '\t' + stop + '\t' + key + '\t' + geneIDs[key]['common_name'] + '\n')
    elif geneIDs[key]['sets'][0] == 1:
        test_gene_file.write(chr + '\t' + start + '\t' + stop + '\t' + key + '\t' + geneIDs[key]['common_name'] + '\n')
    elif geneIDs[key]['sets'][0] == 2:
        ref_gene_file.write(chr + '\t' + start + '\t' + stop + '\t' + key + '\t' + geneIDs[key]['common_name'] + '\n')
