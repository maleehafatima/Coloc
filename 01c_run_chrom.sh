#!/bin/bash

#./01c_run_chrom.sh /homes/data/COLOC_input_data/uniq_pred_db_hg38.chrchrom.maf0.01.R20.8.anno.txt.gz /homes/data/COLOC_input_data/annotation_all_aptamers_ENSG.txt /homes/data/COLOC_input_data/uniq_pred_db_hg38.chrchrom.maf0.01.R20.8.geno.txt.gz /homes/data/COLOC_input_data/annotation_all_aptamers_ENSG.txt ASW AFA CEU
#For file paths in place of chromosome or population, use chrom and pop
snp_annot=$1 #first argument is path to SNP annotation file
gene_annot=$2 #path to gene annotation file
genotype_file=$3 #path to genotype file
expression_file=$4 #path to the expression file
pops=${@:5:99}

for chrom in {1..3}
do 
	nohup Rscript 01b_run_pull_snps_driving.R $chrom $snp_annot $gene_annot $genotype_file $expression_file $pops > output/LD_matrix/nohup_1Mb_chrom${chrom}_2.out &
done
