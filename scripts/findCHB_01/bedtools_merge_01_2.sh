#!/bin/bash

#This script is merging potential CHB windows for many combinations of parameters from original labelRepresentation.py script
SLIDE=5K
EXPANDS=5000
TOP=7500
SEP=_
DIREND1=coverage/
DIREND2=7500/
COC=coverage
END1=.txt
END2=.bed
PATH1=/Users/kateweaver/taylorLab/VISION/train_test_split/labelRepresentationMethod/correct_resultsAfterBashScript/higher_threshold/
FILLER1=potentialWindows_params_
FILLER2=merged_potentialWindows_


for WINDOW in 55K 75K
do
  PATH2=$PATH1$WINDOW$SEP$SLIDE$SEP$DIREND1
  for AFFIRM in yes no
  do
    for TAPREF in T TA
    do
      for NUM in 15 18
      do
        for THRESHOLD in 3.0 4.25 5.0
        do
          PATH3=$PATH2$WINDOW$SEP$SLIDE$SEP$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$THRESHOLD$SEP$DIREND2
          cd $PATH3
          #echo $PATH3
          if [ $WINDOW = 55K ]
          then
            #echo yes1
            EXPANDW=55000
          elif [ $WINDOW = 75K ]
          then
            #echo yes2
            EXPANDW=75000
          else
            echo yikes
          fi
          F1=$PATH3$FILLER1$EXPANDW$SEP$EXPANDS$SEP$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$THRESHOLD$SEP$TOP
          F2=$PATH3$FILLER2$WINDOW$SEP$SLIDE$SEP$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$THRESHOLD$SEP$TOP
          tail -n +2 $F1$END1 | sort -k1,1 -k2,2n > $F1$END2
          bedtools merge -i $F1$END2 -c 1 -o count > $F2$END2
        done
      done
    done
  done
done
