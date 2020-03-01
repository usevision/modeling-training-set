#!/usr/bin/env python3

import argparse as ap
import logging
import datetime
import numpy as np
import subprocess
import matplotlib.pyplot as plt

logging.basicConfig(filename='furtherAsses.log', level=logging.DEBUG)
logging.info(datetime.datetime.now())

parser = ap.ArgumentParser(description='further assess training and test split')
parser.add_argument('--trainSet', action = 'store', nargs=1, type=str, required = True, help='file in bed format with train set')
parser.add_argument('--testSet', action = 'store', nargs=1, type=str, required=True, help='file in bed format with test set')
parser.add_argument('--LADs', action='store', nargs='+', type=str, required=True, help='files with LAD information') #~/taylorLab/VISION/data/cscore/*_10K_cscore.bg
parser.add_argument('--IDEAS', action='store', nargs='+', type=str, required=True, help='files with IDEAS information for quiescent only') #~/taylorLab/VISION/data/IDEAS/ideasVisionV20p8Seg*.bed.o2.bed
parser.add_argument('--RNAseq_both', action='store', nargs=1, type=str, required=True, help='file with RNAseq expression information for 16 cell types (including Amit data)') #~/taylorLab/VISION/data/RNAseq/both_RNAseq.tab
args = parser.parse_args()
trainSet = args.trainSet[0]
testSet = args.testSet[0]
LAD_files = args.LADs
IDEAS_files = args.IDEAS
RNAseq_both = args.RNAseq_both[0]

def bedtools_intersect(fileA, fileB):
    bedtools_intersect = subprocess.check_output("bedtools intersect -a {} -b {} -wo".format(fileA, fileB), shell=True).decode("utf-8").splitlines()
    bedtools_dont_intersect = subprocess.check_output("bedtools intersect -a {} -b {} -v".format(fileA, fileB), shell=True).decode("utf-8").splitlines()
    return(bedtools_intersect, bedtools_dont_intersect)

'''Using recorded split, further assess for ....'''
'''....LADs''' #proportion of windows (which have a cscore suggesting LADs) that overlap regions
cant_determine_train = []
train_values = {}
cant_determine_test = []
test_values = {}
for i,LAD_file in enumerate(LAD_files):
    bedtools_intersect_train, bedtools_complement_train = bedtools_intersect(trainSet, LAD_file)
    for line in bedtools_complement_train:
        fields = line.strip('\r\n').split('\t')
        valuesOI = (fields[0], fields[1], fields[2])
        cant_determine_train.append(valuesOI)
    for line in bedtools_intersect_train:
        fields=line.strip('\r\n').split('\t')
        if float(fields[6]) < 0:
            LAD = False
        else:
            LAD = True
        valuesOI = (fields[0], fields[1], fields[2])
        if valuesOI not in train_values:
            train_values[valuesOI]={0:[],
                                    1:[],
                                    2:[]}

        train_values[valuesOI][i].append(LAD)
    bedtools_intersect_test, bedtools_complement_test = bedtools_intersect(testSet, LAD_file)
    for line in bedtools_complement_test:
        fields = line.strip('\r\n').split('\t')
        valuesOI = (fields[0], fields[1], fields[2])
        cant_determine_test.append(valuesOI)

    for line in bedtools_intersect_test:
        fields=line.strip('\r\n').split('\t')
        valuesOI=(fields[0], fields[1], fields[3])
        if float(fields[6]) < 0:
            LAD = False
        else:
            LAD = True
        if valuesOI not in test_values:
            test_values[valuesOI]={0:[],
                                   1:[],
                                   2:[]}
        test_values[valuesOI][i].append(LAD)

LADs_train_first = []
LADs_test_first = []
LADs_train_second = []
LADs_test_second = []
LADs_train_third = []
LADs_test_third = []

for key in train_values:
    first = np.sum(train_values[key][0])/len(train_values[key][0])
    LADs_train_first.append(first)
    second = np.sum(train_values[key][1])/len(train_values[key][1])
    LADs_train_second.append(second)
    third = np.sum(train_values[key][2])/len(train_values[key][2])
    LADs_train_third.append(third)

for key in test_values:
    first = np.sum(test_values[key][0])/len(test_values[key][0])
    LADs_test_first.append(first)
    second = np.sum(test_values[key][1])/len(test_values[key][1])
    LADs_test_second.append(second)
    third = np.sum(test_values[key][2])/len(test_values[key][2])
    LADs_test_third.append(third)

fig, (ax1, ax2) = plt.subplots(1,2)
fig.suptitle('proportion of LADs within region')
ax1.set_ylabel('number of regions')
ax1.set_xlabel('proportion LADs')
y1tr, bins1tr, patches = ax1.hist(LADs_train_first, bins=10, color='dodgerblue')
bincenters1tr = 0.5*(bins1tr[1:]+bins1tr[:-1])
y2tr, bins2tr, patches =ax1.hist(LADs_train_second, bins=10, color='darkorange')
bincenters2tr = 0.5*(bins2tr[1:]+bins2tr[:-1])
y3tr, bins3tr, patches = ax1.hist(LADs_train_third, bins=10, color='saddlebrown')
bincenters3tr = 0.5*(bins3tr[1:]+bins3tr[:-1])
ax1.set_title('Train Set')
ax2.set_xlabel('proportion LADs')
ax2.set_title('Test Set')
y1t, bins1t, patches = ax2.hist(LADs_test_first, bins=10, color='dodgerblue')
bincenters1t = 0.5*(bins1t[1:]+bins1t[:-1])
y2t, bins2t, patches = ax2.hist(LADs_test_second, bins=10, color='darkorange')
bincenters2t = 0.5*(bins2t[1:]+bins2t[:-1])
y3t, bins3t, patches = ax2.hist(LADs_test_third, bins=10, color='saddlebrown')
bincenters3t = 0.5*(bins3t[1:]+bins3t[:-1])
plt.legend(['GSM2514768', 'HPC7_Rep1', 'HPC7_Rep2'])
#ax2.legend(['GSM2514768', 'HPC7_Rep1', 'HPC7_Rep2'])
fig.savefig("LAD_split_hist.png")
plt.close(fig)

fig, (ax1, ax2) = plt.subplots(1,2)
fig.suptitle('proportion of LADs within region')
ax1.set_ylabel('number of regions')
ax1.set_xlabel('proportion LADs')
ax1.plot(bincenters1tr, y1tr, color='dodgerblue')
ax1.plot(bincenters2tr, y2tr, color='darkorange')
ax1.plot(bincenters3tr, y3tr, color='saddlebrown')
ax1.set_title('Train Set')
ax2.set_xlabel('proportion LADs')
ax2.set_title('Test Set')
ax2.plot(bincenters1t, y1t, color='dodgerblue')
ax2.plot(bincenters2t, y2t, color='darkorange')
ax2.plot(bincenters3t, y3t, color='saddlebrown')
plt.legend(['GSM2514768', 'HPC7_Rep1', 'HPC7_Rep2'])
#ax2.legend(['GSM2514768', 'HPC7_Rep1', 'HPC7_Rep2'])
fig.savefig("LAD_split_line.png")
plt.close(fig)

'''....genes expressed across 0-16 cell types''' #data not available for HPC7 or Eryfl
bedtools_intersect_all_train, bedtools_complement_all_train = bedtools_intersect(trainSet, RNAseq_both)
bedtools_intersect_all_test, bedtools_complement_all_test = bedtools_intersect(testSet, RNAseq_both)

def get_list_of_expressed(bedtools_intersect):
    expressed = []
    for line in bedtools_intersect:
        counter_expressed = 16
        fields = line.strip('\r\n').split('\t')
        expression = fields[6].split('=')[1:]
        for value in expression:
            decimal, bleh = value.split(';')
            valueOI = float(decimal)
            if valueOI == 0.00:
                counter_expressed += -1
        expressed.append(counter_expressed)
    return (expressed)

train_expressed = get_list_of_expressed(bedtools_intersect_all_train)
print(max(train_expressed))
print(min(train_expressed))
print(np.average(train_expressed))
test_expressed = get_list_of_expressed(bedtools_intersect_all_test)
print(max(test_expressed))
print(min(test_expressed))
print(np.average(test_expressed))

fig, (ax1, ax2) = plt.subplots(1,2)
plt.suptitle("Genes expressed in n-cell types")
ax1.hist(train_expressed, bins=16)
ax1.set_xlabel("n-cell types")
ax1.set_ylabel("number of genes")
ax1.set_title("Train Set")
ax1.xaxis.set_tick_params(rotation=45)
ax2.hist(test_expressed, bins=16)
ax2.set_xlabel("n-cell types")
ax2.set_title("Test Set")
ax2.xaxis.set_tick_params(rotation=45)
#plt.tight_layout()
fig.savefig("expression_n_cellTypes.png")
plt.close(fig)
