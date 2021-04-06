#!/bin/bash

pop=$1 #stores first argument in variable pop
vcf=$2 #second argument is the path to the directory containing directories for each population that contain vcf files

mkdir ${PWD}/output/LD_matrix/${pop}/old_bims

#for chr in {1..22} #Loop through chromosomes 1-22
for chr in {2..3} #for testing purpose
do
	#echo ${vcf}/${pop}.chr${chr}.phase3.vcf.gz
	#echo ${PWD}/output/LD_matrix/${pop}/${pop}_chr${chr}_dose
	plink --vcf ${vcf}/${pop}.chr${chr}.phase3.vcf.gz --make-bed --out ${PWD}/output/LD_matrix/${pop}/${pop}_chr${chr}_dose
	#--vcf flag specifies the vcf input files #1000 genomes genotype dosage data
	#--make-bed reformats dosage files into .bed/bim/fam files #PLINK binary formats

	#Reformat rsids in the outputted bim files to match chr:position:allele1:allele2
	mv ${PWD}/output/LD_matrix/${pop}/${pop}_chr${chr}_dose.bim ${PWD}/output/LD_matrix/${pop}/old_bims/${pop}_chr${chr}_dose.bim
	Rscript 02b_reformat_rsids.R ${PWD}/output/LD_matrix/${pop}/ ${chr} ${pop}
done
