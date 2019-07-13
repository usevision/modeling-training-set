#!/usr/bin/env python3

import fasta
import sys
import matplotlib.pyplot as plt
import numpy as np

'''Usage: ./computeMeanGC_byChr.py chr{}
ran for each chromosome individually (as concurrent but separate jobs); lines were appended to gcMean_proposed_splits.bed'''
#chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5','chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11','chr12', 'chr13', 'chr14', 'chr15', 'chr16','chr17','chr18', 'chr19', 'chrX', 'chrY']
#filePath = '/Users/kateweaver/taylorLab/VISION/train_test_split/labelRepresentationMethod/correct_resultsAfterBashScript/higher_threshold/mergedFiles/proposed_splits_{}.bed'
filePath = '/home/kweave23/VISION_train_test_split/proposed_splits_{}.bed'

def computeGC(sequence):
    Gs = sequence.casefold().count('G'.casefold())
    Cs = sequence.casefold().count('C'.casefold())
    if Gs==0 and Cs==0:
        GCcontent = 0
    else:
        length = len(sequence)
        GCcontent = (Gs+Cs)/length
    return GCcontent

'''sliding windows''' #1 Mbp windows; slide by 500bp; average a proposed split region
fileToWriteTo = open('gcMean_proposed_splits.bed','a+')

chromosome = sys.argv[1]
#file = '/Users/kateweaver/mm10_genome/{}.fa'.format(chromosome)
file='/project/vision/Data/mm10_genome/{}.fa'.format(chromosome)
reader = fasta.FASTAReader(open(file))
gcMeans = []
starts = []
region_starts=[]
gcList = []
for ident, sequence in reader:
    for i, line in enumerate(open(filePath.format(chromosome))):
        if i%25 == 0 and i !=0:
            print("line 38", i, flush=True)
        fields=line.strip('\r\n').split('\t')
        chromosome = fields[0]
        start = int(fields[1])
        region_starts.append(start)
        end = int(fields[2])
        together = fields[3]

        window=0
        slides=0
        subsequence = sequence[start:end+1]
        seqLen = len(subsequence)
        for i in range(0, seqLen-1000000, 500):
            gc = computeGC(subsequence[i:i+1000000])
            gcList.append(gc)
            starts.append(start+i)
            window += gc
            slides += 1
        finalStart = list(range(0,seqLen-1000000,500))[-1]
        gc=computeGC(subsequence[finalStart:seqLen])
        gcList.append(gc)
        starts.append(finalStart)
        window += gc
        slides += 1
        gcMean = window/slides
        gcMeans.append(gcMean)
        fileToWriteTo.write(ident + '\t' + str(start) + '\t' + str(end) + '\t' + str(gcMean) + '\t'+together+'\t' + str(slides) + '\t' + str(seqLen) + '\n')
        fileToWriteTo.flush()
    #plot the sliding window for each region on the overall chromosome
    sort_index = np.argsort(starts)
    starts=np.array(starts)
    gcList=np.array(gcList)
    start_sorted = starts[sort_index]
    gcList_sorted = gcList[sort_index]
    fig, ax = plt.subplots()
    plt.plot(start_sorted, gcList_sorted)
    plt.scatter(region_starts, gcMeans)
    plt.suptitle('1 Mbp Windows (sliding 500 bp and taking mean at proposed splits) {}'.format(chromosome))
    ax.set_ylabel('GC content')
    ax.set_xlabel('Nucleotide Start Position')
    fig.savefig('slidingWindow_{}.png'.format(chromosome))
    plt.close(fig)
    print('finished with ', ident, flush=True)
