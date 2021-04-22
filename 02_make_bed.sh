#!/bin/bash

pop=$1 #stores first argument in variable pop
chr=$2 #chromosome number
vcf=$3 #third argument is the path to the directory containing directories for each population that contain vcf files
out=$4 #output directory

#Run plink command on correct vcf file
plink --vcf ${vcf} --make-bed --mind 0.01 --maf 0.01 --out ${out}/LD_matrix/${pop}/${pop}_chr${chr}_dose
#--vcf flag specifies the vcf input files #1000 genomes genotype dosage data
#--make-bed reformats dosage files into .bed/bim/fam files #PLINK binary formats

