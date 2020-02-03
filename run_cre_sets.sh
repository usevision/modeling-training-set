#!/bin/bash

source activate basic

date; time cut -f1,2,3 VISIONmusHem_ccREs.bed > VISIONmusHem_ccREs_locs.bed
date; time bedtools intersect -a VISIONmusHem_ccREs_locs.bed -b train.bed -wa -wb > ccREs_int_train.bed
date; time bedtools intersect -a VISIONmusHem_ccREs_locs.bed -b test.bed -wa -wb > ccREs_int_test.bed
date; time bedtools intersect -a VISIONmusHem_ccREs_locs.bed -b referenceLoci_fullRegion.bed -wa -wb > ccREs_int_ref.bed

wc -l VISIONmusHem_ccREs_locs.bed
wc -l ccREs_int_train.bed ccREs_int_test.bed ccREs_int_ref.bed

date; time python cre_sets.py --cres VISIONmusHem_ccREs_locs.bed --train ccREs_int_train.bed --test ccREs_int_test.bed --ref ccREs_int_ref.bed

wc -l test_cre_uniq.bed train_cre_uniq.bed ref_cre_uniq.bed

conda deactivate
