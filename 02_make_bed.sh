#!/bin/bash

pop=$1 #stores first argument in variable pop
for chr in {1..22} #Loop through chromosomes 1-22
do
	./plink --vcf /home/egeoffroy/MESA_TOPMED_Imputation/${pop}/chr${chr}.dose.vcf.gz --make-bed --out /home/egeoffroy/LD_matrix/${pop}/${pop}_chr${chr}_dose
	#--vcf flag specifies the vcf input files #1000 genomes genotype dosage data
	#--make-bed reformats dosage files into .bed/bim/fam files #PLINK binary formats
done
