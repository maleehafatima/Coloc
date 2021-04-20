# GWAS-eQTL Colocalization Pipeline
---
## Languages and Packages
---
**R**
```
library(dplyr)
library(reshape2)
library(glmnet)
library(methods)
library(doMC)
library(doRNG)
library(tidyr)
library(tibble)
```
**bash**

**Python**
```
import os
import os.path
import argparse
```
## Software ##
---
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- [COLOC 3.2-1](https://github.com/chr1swallace/coloc)

## Running the pipeline ##

## Workflow ##
---
### Scripts 1a-b ###
These scripts use imputed transcript levels from the eQTL population to generate lists of significant SNPs that are within 1Mb of a predicted gene. The goal is to pull all significant cis-acting variants. This step can be executed by calling script 1c, which runs script 1b for every chromosome:

#### Input Files: ####
- SNP annotation file
- Gene annotation file
- Genotype file
- Expression file
#### Output Files: ####
- output/LD_matrix/{pop}/{pop}\_chr\_{chrom}\_{gene}\_1Mb\_of\_gene.txt

### Script 2 ###
This script takes in genotype dosage files in vcf file format and converts them to PLINK binary formats, .bed, .bim, and .fam. This step can be executed by calling script 2.
```
chmod u+x 02_make_bed.sh 
./02_make_bed.sh {pop} {path/to/directory/containing/population_folders_containing_dosage_files/}
```
#### Input Files: ####
- Genotype dosage files
#### Output Files: ####
- output/LD_matrix/{pop}/{pop}\_chr{chrom}\_dose.bed
- output/LD_matrix/{pop}/{pop}\_chr{chrom}\_dose.bim
- output/LD_matrix/{pop}/{pop}\_chr{chrom}\_dose.fam

### Script 3 ###
This script uses the significant SNP lists produced by script 1, along with the bfiles produced by script 2, to calculate LD matrices between every significant SNP around every predicted gene. This step can be executed by calling script 3.
```
chmod u+x 03_make_LD_matrix.sh
./03_make_LD_matrix.sh {pop}
```
#### Input Files: ####
- output/LD_matrix/{pop}/{pop}\_chr\_{chrom}\_{gene}\_1Mb\_of\_gene.txt
- output/LD_matrix/{pop}/{pop}\_chr{chrom}\_dose.bed
- output/LD_matrix/{pop}/{pop}\_chr{chrom}\_dose.bim
- output/LD_matrix/{pop}/{pop}\_chr{chrom}\_dose.fam
#### Output Files: ####
- output/LD_matrix/{pop}/{pop}\_1Mb\_coords\_LDMatrix/{pop}\_chr\_{chrom}\_{gene}\_1Mb\_LD.ld.gz
- output/LD_matrix/{pop}/{pop}\_1Mb\_coords\_LDMatrix/{pop}\_chr\_{chrom}\_{gene}\_1Mb\_LD.snplist.gz

### Script 4 ###
This script takes the GWAS summary statistics, eQTL, and .frq files input by the user and formats them properly in preparation of COLOC analysis.
#### Input Files: ####
- input/{pop}\_prot\_hg38.frq
- input/GWAS_SS/WojcikG_{phenotype}.txt.gz
- input/cis_eQTLs\_{pop}\_WG\_all\_cis.txt.gz
#### Output Files: ####
- output/pQTL/{pop}/pQTL\_{pop}\_{phenotype}.txt.gz
- output/GWAS_TOPMED/{pop}/GWAS\_TOPMED\_{pop}\_{phenotype}.txt.gz

### Scripts 5a-b ###
These scripts take the input arguments from the user and run COLOC analysis for each population and phenotype specified by the user using the outputs from scripts 3 and 4. A summary table of the COLOC reults are written to an outfile.
#### Input Files: ####
- output/pQTL/{pop}/pQTL\_{pop}\_{phenotype}.txt.gz
- output/GWAS_TOPMED/{pop}/GWAS\_TOPMED\_{pop}\_{phenotype}.txt.gz
- output/LD_matrix/{pop}/{pop}\_1Mb\_coords\_LDMatrix/{pop}\_chr\_{chrom}\_{gene}\_1Mb\_LD.ld.gz
- output/LD_matrix/{pop}/{pop}\_1Mb\_coords\_LDMatrix/{pop}\_chr\_{chrom}\_{gene}\_1Mb\_LD.snplist.gz
#### Output Files: ####
- output/Coloc_output/{pop}\_{phenotype}.txt
