#!/bin/bash

pop=$1 #stores first argument in variable pop
vcf=$2 #second argument is the path to the directory directories for each population that contain vcf files
out=$3 #third argument is the path to the directory for the output bed/bim/fam files
for chr in {1..22} #Loop through chromosomes 1-22
do
	./plink --vcf ${vcf}/${pop}/chr${chr}.dose.vcf.gz --make-bed --out ${out}/${pop}/${pop}_chr${chr}_dose
	#--vcf flag specifies the vcf input files #1000 genomes genotype dosage data
	#--make-bed reformats dosage files into .bed/bim/fam files #PLINK binary formats
done
