#!/bin/bash

#this script actually runs labelRepresentation.py and specifies parameters for argparse()
#this script only includes 2 of 48+ iterations/combinations of thresholds/parameters I used. Later scripts often merge, etc for all/many combinations that I hadn't already ruled out
#COC - coverage or count
#TAPREF - TA or T
#AFFIRM - yes or no
#NUM - 15 or 18
#RATIO lower_thresholds(1.35, 1.5, 2.0) & higher_thresholds(3.0, 4.25, 5.0)
#TOP 7500
#WINDOW 75K 55K
#SLIDE 5K

SEP=_
END=.txt
COC=coverage
TAPREF=TA
AFFIRM=yes
NUM=18
#RATIO=1.35
RATIO=3.0
TOP=7500
WINDOW=55000
SLIDE=5000
NEW=/home/kweave23/VISION_train_test_split/labelRepresentation/55K_5K_$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$RATIO$SEP$TOP

mkdir $NEW
cd $NEW
python ./../labelRepresentation.py --window $WINDOW --slide $SLIDE --COC $COC --countQ $AFFIRM --TApref $TAPREF --numCellTypes $NUM --thresholdRatio $RATIO --thresholdTop $TOP
awk '{print $1}' potentialWindows_params_$WINDOW$SEP$SLIDE$SEP$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$RATIO$SEP$TOP$END | sort | uniq -c > chrom.txt
echo finished with {$COC}_{$AFFIRM}_{$TAPREF}_{$NUM}_{$RATIO}
cd ..

WINDOW=75000
SLIDE=5000
NEW=/home/kweave23/VISION_train_test_split/labelRepresentation/75K_5K_$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$RATIO$SEP$TOP

mkdir $NEW
cd $NEW
python ./../labelRepresentation.py --window $WINDOW --slide $SLIDE --COC $COC --countQ $AFFIRM --TApref $TAPREF --numCellTypes $NUM --thresholdRatio $RATIO --thresholdTop $TOP
awk '{print $1}' potentialWindows_params_$WINDOW$SEP$SLIDE$SEP$COC$SEP$AFFIRM$SEP$TAPREF$SEP$NUM$SEP$RATIO$SEP$TOP$END | sort | uniq -c > chrom.txt
echo finished with {$COC}_{$AFFIRM}_{$TAPREF}_{$NUM}_{$RATIO}
