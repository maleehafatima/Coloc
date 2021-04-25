#!/bin/bash
# Author: Elyse Geoffroy
# This file makes the LD matrices for colocalization using plink with the dosage bed/bim/fam files and a list of sig snps to include inthe matrix for each gene.
pop=$1 #first argument is the population abbreviation
chr=$2 #second argument is the chromosome
out=$3 #output directory

for file in $"{out}/LD_matrix"/${pop}/${pop}_chr_${chr}_*
do
	filename="${file##*/}"
	substring=('_1Mb_of_gene.txt') #file endings
	filename=${filename%"${substring}"} #% operations extracts the substring from the filename

	plink --bfile ${out}/LD_matrix/${pop}/${pop}_chr${chr}_dose --r square gz yes-really --extract ${file} --write-snplist --out ${out}/LD_matrix/${pop}/${pop}_1Mb_coords_LDMatrix/${filename}_1Mb_LD
	#--bfile specifies the input bfiles
	#--r calculates the inter-variant allele count correlations
	#square specifies the shape of the output matrix
	#gz causes the outputs to be gzipped
	#yes-really requests an unfiltered, non-distributed all pairs computation on more than 400k variants, needed because outputs will be extremely large
	#extract pulls out the significant SNPs from the bfiles specified in the significant SNP list
	#write--snplist writes the list IDs of significant SNPs used to an output file
done
