# modeling-training-set
Training sets for VISION modeling

`.bed` files containing a list of full partitioned regions (genomic coordinates) are in the directory [full_region_bed](https://github.com/usevision/modeling-training-set/tree/master/full_region_bed)<br />
`.bed` files containing a list of partitioned genes (based on the full partitioned regions and [this python script](https://github.com/usevision/modeling-training-set/blob/master/scripts/finding_set_for_genes_and_cres_08/gene_sets.py)) are in the directory [specific_gene_bed](https://github.com/usevision/modeling-training-set/tree/master/specific_gene_bed)<br />
`.bed` files containing a list of partitioned cCREs (based on the full partitioned regions and [this python script](https://github.com/usevision/modeling-training-set/blob/master/scripts/finding_set_for_genes_and_cres_08/cre_sets.py) which was run by [this bash script](https://github.com/usevision/modeling-training-set/blob/master/scripts/finding_set_for_genes_and_cres_08/run_cre_sets.sh)) are in the directory [specific_cre_bed](https://github.com/usevision/modeling-training-set/tree/master/specific_cre_bed)<br /> 
