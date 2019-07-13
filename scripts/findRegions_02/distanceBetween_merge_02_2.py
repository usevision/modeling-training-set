#!/usr/bin/env python3

'''Usage ./distanceBetween_merge.py
This file filters out CHB regions which are most definitely less than 4Mbp (set by the filter variable line 46) from the CHB before it. Reports CHB boundaries in filtered_{}_split_loc_ files and then mostly filtered region boundaries in filtered_{}_regions_'''

import matplotlib.pyplot as plt

def happyWithSite(axes, fileTWT2, fileTWT4, fields, chromosome, num_filtered, previousStart, region_start, distances, rounded):
    axes.scatter(rel_start, num_filtered, marker='>')
    end = int(fields[2])
    fileTWT4.write(chromosome + '\t' + str(region_start) + '\t' + str(end) + '\n')
    rel_end = end/chrSizes[chromosome]
    axes.scatter(rel_end,num_filtered, marker='<')
    num_filtered += 1
    fileTWT2.write(chromosome + '\t' + str(previousStart) + '\t' + str(int((region_start + end)/2)) + '\t' + str(rounded) + '\n')
    previousStart = int(((region_start + end)/2) + 1)
    distances.append(rounded)
    return (end, num_filtered, previousStart, distances)

chrSizes = {'chr1':195471971,
            'chr2':182113224,
            'chrX':171031299,
            'chr3':160039680,
            'chr4':156508116,
            'chr5':151834684,
            'chr6':149736546,
            'chr7':145441459,
            'chr10':130694993,
            'chr8':129401213,
            'chr14':124902244,
            'chr9':124595110,
            'chr11':122082543,
            'chr13':120421639,
            'chr12':120129022,
            'chr15':104043685,
            'chr16':98207768,
            'chr17':94987271,
            'chrY':91744698,
            'chr18':90702639,
            'chr19':61431566}

chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
                 'chr18', 'chr19', 'chrX', 'chrY']

filter = 4

#filePath = '/Users/kateweaver/taylorLab/VISION/train_test_split/labelRepresentationMethod/correct_resultsAfterBashScript/higher_threshold/mergedFiles/chrAwk/merged_potentialWindows_{}.{}.{}.bed'
filePath = '/Users/kateweaver/taylorLab/VISION/train_test_split/labelRepresentationMethod/correct_resultsAfterBashScript/higher_threshold/mergedFiles/chrAwk/merged_potentialWindows_{}.{}.bed'
Sep='_'
Windows = '75K'
Slide = '5K'
COC = 'coverage'
Affirms = ['yes', 'no']
TAprefs = ['T', 'TA']
Nums = '18'
Thresholds = ['3.0','4.25','5.0']
Top = '7500'
#merge_num = 'multiple'
merge_num = 'both' #had toyed around with only using CHBs which were supported by multiple overlapping regions to produce the final ROI (>1 in final column of merged_potentialWindows.... files). 'both' here signifies that I'm including those regions which had nothing overlapping it as well
for affirm in Affirms:
    for TApref in TAprefs:
        for threshold in Thresholds:
            restOfFile = Windows+Sep+Slide+Sep+COC+Sep+affirm+Sep+TApref+Sep+Nums+Sep+threshold+Sep+Top
            fileToWriteTo1 = open('filtered_{}_num_{}.txt'.format(merge_num, restOfFile), 'w+') #this file will report on the number of CHBs that meet the basic filter requirement
            fileToWriteTo3 = open('beginningOfChr_{}_{}.txt'.format(merge_num, restOfFile),'w+') #this file will report on the first CHB distance (distance between beginning of boundary and beginning of chromosome)
            fileToWriteTo5 = open('endOfChr_{}_{}.txt'.format(merge_num, restOfFile), 'w+') #this file will report on the last CHB distance (distance between end of boundary and end of chromosome)
            for chromosome in chromosomes:
                fig1,ax1=plt.subplots()
                fig1.suptitle("{} Location of Regions OI".format(chromosome))
                ax1.set_xlabel('Relative Chromosome Location')
                ax1.set_ylabel('Region number')
                ax1.axvline(0,linestyle='--')
                ax1.axvline(1, linestyle='--')
                num_filtered = 0
                #fileOI = filePath.format(restOfFile, chromosome, merge_num)
                fileOI = filePath.format(restOfFile, chromosome)
                fileToWriteTo2 = open('filtered_{}_regions_{}_{}.bed'.format(merge_num, restOfFile, chromosome), 'w+') #this file will report on regions for partitioning in bed format, specifically with distance between in final column. In this last column, the only row that should be below the filter set in this script is the final row. #Went through and annotated the files by hand aftewards to make a single file of all the chromosomes (cat >>) and then had a final final column that said if it was first, next to last, or last region
                fileToWriteTo4 = open('filtered_{}_splitLoc_{}_{}.bed'.format(merge_num, restOfFile, chromosome), 'w+') #this file will report on CHB regions in bed format
                end = 0
                previousStart = 0
                distances = []
                for i, line in enumerate(open(fileOI)):
                    if i == 0:
                        firstValue = True
                    fields=line.strip('\r\n').split('\t')
                    region_start = int(fields[1])
                    rel_start = region_start/chrSizes[chromosome]
                    distance = region_start - end
                    rounded = float("{0:.2f}".format(distance/1000000))
                    if firstValue:
                        if rounded/filter > 1 or (rounded%filter == 0 and rounded >= filter):
                            fileToWriteTo3.write(chromosome + '\t' + str(rounded) + '\n')
                            end, num_filtered, previousStart, distances = happyWithSite(ax1, fileToWriteTo2, fileToWriteTo4, fields, chromosome, num_filtered, previousStart, region_start, distances, rounded)
                            firstValue=False
                    else:
                        if rounded/filter > 1 or (rounded%filter == 0 and rounded >= filter):
                            end, num_filtered, previousStart, distances = happyWithSite(ax1, fileToWriteTo2, fileToWriteTo4, fields, chromosome, num_filtered, previousStart, region_start, distances, rounded)
                fileToWriteTo2.write(chromosome+'\t'+str(previousStart)+'\t'+str(chrSizes[chromosome])+ '\t' + str((chrSizes[chromosome]-previousStart)/1000000)+'\n')
                fileToWriteTo5.write(chromosome + '\t' + str((chrSizes[chromosome]-previousStart)/1000000)+'\n')
                fileToWriteTo1.write(chromosome + '\t' + str(num_filtered) + '\n')
                fig2,ax2=plt.subplots()
                fig2.suptitle("{} Distances/1000000 between Regions OI".format(chromosome))
                ax2.axvline(1, linestyle='--', color='black')
                ax2.axvline(3, linestyle='--', color='indigo')
                ax2.axvline(5, linestyle='--', color='darkslategray')
                ax2.axvline(10, linestyle='--', color='darkgoldenrod')
                ax2.set_xlabel('Distance/1000000 bp')
                ax2.set_ylabel('Count')
                ax2.hist(distances, bins=len(list(set(distances))))

                fig1.savefig('filtered_location_{}_{}.png'.format(restOfFile, chromosome))
                plt.close(fig1)

                fig2.savefig('distances_{}_filtered_{}_{}.png'.format(merge_num, restOfFile, chromosome))
                plt.close(fig2)
