#!/usr/bin/env python3

import argparse as ap
import numpy as np

'''Usage: ./cre_sets.py --cres VISIONmusHem_ccREs_locs --train ccREs_int_train.bed --test ccREs_int_test.bed --ref ccREs_int_ref.bed'''

parser = ap.ArgumentParser(description='go through list of cres in train/test/reference set and make sure that each cre is in a unique set')
parser.add_argument('--cres', action='store', nargs=1, type=str, required=True, help='bed file with cre_locs')
parser.add_argument('--train', action='store', nargs=1, type=str, required = True, help='bed file with train.bed intersected with cre_locs')
parser.add_argument('--test', action='store', nargs=1, type=str, required = True, help='bed file with test.bed intersected with cre_locs')
parser.add_argument('--ref', action='store', nargs=1, type=str, required=True, help='bed file with referenceLoci_fullRegion.bed intersected with cre_locs')
args=parser.parse_args()
cre_file = args.cres[0]
train_file = args.train[0]
test_file = args.test[0]
ref_file = args.ref[0]

set_dict = {'train':0,
            'test':1,
            'ref':2}

'''add ccREs to a dictionary to track what sets it overlaps & the genomic location'''
ccREs = {}
for line in open(cre_file):
    fields=line.strip('\r\n').split('\t')
    chr = fields[0]
    start = fields[1]
    end = fields[2]
    cre_ID = "_".join((chr,start,end))
    if chr != 'chrM':
        if cre_ID not in ccREs:
            ccREs[cre_ID]={'sets':[], 'loc':(chr, start, end)} #should check that the known ref CREs remain in the reference region. Way utilized assigned anything (in this case 1) that overlapped multiple including reference to reference

'''cres that overlap training set'''
setOI = 'train'
for line in open(train_file):
    fields=line.strip('\r\n').split('\t')
    cre_ID = "_".join((fields[0], fields[1], fields[2]))
    ccREs[cre_ID]['sets'].append(set_dict[setOI])

'''cres that overlap test set'''
setOI = 'test'
for line in open(test_file):
    fields=line.strip('\r\n').split('\t')
    cre_ID = "_".join((fields[0], fields[1], fields[2]))
    ccREs[cre_ID]['sets'].append(set_dict[setOI])

'''cres that overlap the ref set'''
setOI = 'ref'
for line in open(ref_file):
    fields=line.strip('\r\n').split('\t')
    cre_ID = "_".join((fields[0], fields[1], fields[2]))
    ccREs[cre_ID]['sets'].append(set_dict[setOI])

'''if a cres overlaps more than one set'''
for key in ccREs:
    uniq = set(ccREs[key]['sets']) #cres could be in multiple windows within the same set; this together with line 67-68 will reduce those
    uniq_len = len(uniq)
    if uniq_len > 1:
        uniq = list(uniq)
        if 2 in uniq:
            ccREs[key]['sets'] = 2
        else:
            uniq_index = np.random.randint(uniq_len) #get a random index less than the length of the list of sets
            assigned_set = uniq[uniq_index] #assign a set based on this random index
            ccREs[key]['sets']= assigned_set
    else:
        ccREs[key]['sets'] = list(uniq)[0]

test_cre_file = open('test_cre_uniq.bed','w+')
train_cre_file = open('train_cre_uniq.bed','w+')
ref_cre_file = open('ref_cre_uniq.bed','w+')

'''go through the dictionary and write the cre to whatever set it belongs to'''
for key in ccREs:
    chr, start, stop = ccREs[key]['loc']
    if ccREs[key]['sets'] == 0:
        train_cre_file.write(chr + '\t' + start + '\t' + stop + '\n')
    elif ccREs[key]['sets'] == 1:
        test_cre_file.write(chr + '\t' + start + '\t' + stop + '\n')
    elif ccREs[key]['sets'] == 2:
        ref_cre_file.write(chr + '\t' + start + '\t' + stop + '\n')

test_cre_file.close()
train_cre_file.close()
ref_cre_file.close()
