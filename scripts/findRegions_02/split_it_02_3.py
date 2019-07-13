#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

#This script fully filters CHB regions to produce proposed_splits or regions that are 4-7 Mbp in length

furtherPartCounter = 0

def furtherPartitionIt(furtherPartitionDict, furtherPartCounter, distance, FTWT, chromosome, start, end):
    furtherPartCounter+=1
    theLen = len(furtherPartitionDict[distance])
    for i,value in enumerate(furtherPartitionDict[distance]):
        transformed_value = value*1000000
        if i == 0:
            FTWT.write(chromosome + '\t' + str(start)+'\t'+str(start+transformed_value)+'\t'+'together_'+str(furtherPartCounter)+'\n')
            new_start = start+transformed_value+1
        elif i != 0 and i != theLen-1:
            new_end = new_start+transformed_value
            FTWT.write(chromosome + '\t' + str(new_start)+'\t'+str(new_end)+'\t'+'together_'+str(furtherPartCounter)+'\n')
            new_start = new_end+1
        elif i == theLen-1:
            FTWT.write(chromosome + '\t' + str(new_start)+'\t'+str(end)+'\t'+'together_'+str(furtherPartCounter)+'\n')
    return (furtherPartCounter)

furtherPartition={4:[4],
                  5:[5],
                  6:[6],
                  7:[7],
                  8:[4,4],
                  9:[4,5],
                  10:[5,5],
                  11:[5,6],
                  12:[6,6],
                  13:[6,7],
                  14:[4,5,5],
                  15:[4,5,6],
                  16:[5,5,6],
                  17:[5,6,6],
                  18:[6,6,6],
                  19:[4,4,5,6],
                  20:[4,5,5,6],
                  21:[5,5,5,6],
                  22:[5,5,6,6]}


filePath = '/Users/kateweaver/taylorLab/VISION/train_test_split/labelRepresentationMethod/correct_resultsAfterBashScript/higher_threshold/mergedFiles/filter4Mbp/candidate_regions_75K_5K_coverage_yes_TA_18_3.0_7500.bed' #the cat >> and self annotated file from filtered_{}_regions_

proposed_splits = open('proposed_splits.bed', 'w+')
for i, line in enumerate(open(filePath)):
    fields=line.strip('\r\n').split('\t')
    chromosome = fields[0]
    start=int(fields[1])
    end=int(fields[2])
    distance = round(float(fields[3]))
    chrLoc = fields[4] #I manually annotated this final field originally in the candidate_regions_..._.bed file
    if chrLoc == 'nextToLast': #don't want to write these coordinates until I know whether the final region of the chromosome is stand alone or not
        nextToLast_start = start
        nextToLast_end = end
        previousDistance = distance ##manually had checked that none of these were by themselves greater than 7

    if chrLoc == 'last':
        if distance >= 4:
            proposed_splits.write(chromosome+'\t'+str(nextToLast_start)+'\t'+str(nextToLast_end)+'\t'+'by_itself'+'\n') #go ahead and print the next to last
            if len(furtherPartition[distance])>1: #does the last need split further first?
                furtherPartCounter = furtherPartitionIt(furtherPartition, furtherPartCounter, distance, proposed_splits, chromosome, start, end)
            else:
                proposed_splits.write(chromosome+'\t'+str(start)+'\t'+str(end)+'\t'+'by_itself' + '\n')
        else: #combine the next to the last with the last
            if len(furtherPartition[distance+previousDistance])==1: #is the combined still under 7?
                proposed_splits.write(chromosome+'\t'+str(nextToLast_start)+'\t'+str(end)+'\t'+'by_itself'+'\n')
            else:
                furtherPartCounter = furtherPartitionIt(furtherPartition, furtherPartCounter, distance+previousDistance, proposed_splits, chromosome, nextToLast_start, end)
    elif chrLoc == 'notLast':
        previousDistance = distance
        if len(furtherPartition[distance])>1:
            furtherPartCounter = furtherPartitionIt(furtherPartition, furtherPartCounter, distance, proposed_splits, chromosome, start, end)
        else:
            proposed_splits.write(chromosome + '\t' + str(start) + '\t' + str(end) + '\t' + 'by_itself' + '\n')
