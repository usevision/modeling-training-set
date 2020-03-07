#!/bin/bash

#Usage; ./run_cre_prop_in_state.sh train_file test_file ref_file

source activate basic

date; time bedtools intersect -a $1 -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames >train_region_int_all.bed
date; time bedtools intersect -a $1 -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > train_region_v_all.bed
wc -l train_region_v_all.bed

date; time bedtools intersect -a $2 -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > test_region_int_all.bed
date; time bedtools intersect -a $2 -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > test_region_v_all.bed
wc -l test_region_v_all.bed

date; time bedtools intersect -a $3  -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -wa -wb -filenames > ref_region_int_all.bed
date; time bedtools intersect -a $3  -b /project/vision/Data/IDEAS/ideasVisionV20p8Seg*.bed -v -filenames > ref_region_v_all.bed
wc -l ref_region_v_all.bed

IDEAS=/project/vision/Data/IDEAS/ideasVisionV20p8Seg
BED=.bed
for CT in B Cd4 Cd8 Cfue Cfum Clp Cmp Er4 Eryad Eryfl G1e Gmp Hpc7 Imk Mk Lsk Mep Mon Neu Nk
do
  for type in ref train test
  do
    date; time awk -v var=$IDEAS$CT$BED '{if ($4 == var) {print}}' ${type}_region_int_all.bed > ${type}_region_int_$CT$BED
    date; time awk '{print $1, $2, $3, $6"_"$7"_"$8}' OFS='\t' ${type}_region_int_$CT$BED > ${type}_region_istate_$CT$BED
    date; time bedtools groupby -i ${type}_region_istate_$CT$BED -g 1,2,3 -c 4 -o collapse > ${type}_region_collapse_$CT$BED
  done
done

date; time python ccre_prop_in_state.py train_region_collapse_{}.bed 469 train_region_state_prop.npz
date; time python ccre_prop_in_state.py test_region_collapse_{}.bed 53 test_region_state_prop.npz
date; time python ccre_prop_in_state.py ref_region_collapse_{}.bed 12 ref_region_state_prop.npz

conda deactivate
