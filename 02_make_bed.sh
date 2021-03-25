#!/bin/bash

pop=$1 #stores first argument in variable pop
vcf=$2 #second argument is the path to the directory containing directories for each population that contain vcf files
for chr in {1..22} #Loop through chromosomes 1-22
do
	./plink --vcf ${vcf}/${pop}/chr${chr}.dose.vcf.gz --make-bed --out output/LD_matrix/${pop}/${pop}_chr${chr}_dose
	#--vcf flag specifies the vcf input files #1000 genomes genotype dosage data
	#--make-bed reformats dosage files into .bed/bim/fam files #PLINK binary formats
done
