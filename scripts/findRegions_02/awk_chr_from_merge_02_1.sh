#!/bin/bash

#this script separates chromosomes from the merged CHB files into separate files for easier use in distanceBetween_merge.py and split_it.py

#PATH1P1=/Users/kateweaver/taylorLab/VISION/train_test_split/labelRepresentationMethod/resultsAfterBashScript/mergedFiles/merged_potentialWindows_
PATH1P1=/Users/kateweaver/taylorLab/VISION/train_test_split/labelRepresentationMethod/correct_resultsAfterBashScript/higher_threshold/mergedFiles/merged_potentialWindows_
PATH1P2=.bed

function awk_it {
  for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
  do
    awk '{if ($1 == "'$CHR'") {print}}' $1 > $2.$CHR$3
  done
}

function call_it {
  SEP=_
  TOP=7500
  SLIDE=5K
  COC=coverage
  for WINDOW in 55K 75K
  do
    for AFFIRM in yes no
    do
      for TAPREF in T TA
      do
        for NUM in 15 18
        do
          for THRESHOLD in 3.0 4.25 5.0
          do
            RESTOFFILE=$WINDOW$SEP$SLIDE$SEP$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$THRESHOLD$SEP$TOP
            FILEOI=$1$RESTOFFILE$2
            awk_it $FILEOI $1$RESTOFFILE $2
          done
        done
      done
    done
  done
}

call_it $PATH1P1 $PATH1P2
