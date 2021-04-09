library(dplyr)
library(data.table)

#For running script 1: check if they have the necessary directories otherwise make output -> LD_matrix -> folders for every population

## Get & set working directory
#Tell user to cd into Coloc before running
wd <- getwd()
wd <- paste(wd,'/',sep='')


## Import necessary scripts
source(paste(wd,'05_coloc.R',sep=''))
source(paste(wd, '05b_run_coloc.R', sep=''))


## Import arguments 
args <- commandArgs(trailingOnly=TRUE)
print(args)
input_dir <- args[1]
pops_1 <- args[2]
pops_4 <- args[3]
phenos <- args[4]
pop_sizes <- as.numeric(args[5])
#gene_id<-args[6]  #NOTE: probs don't need this argument

genes<-readRDS(gene_id) #NOTE: need to figure out how this works with 5b


## MAKE DIRECTORIES
#Create output dir
dir.create(paste(wd, 'output', sep = ''), showWarnings = F)
#List of output dirs
dirs <- c('LD_matrix', 'gene_lists', 'GWAS_TOPMED', 'pQTL', 'Coloc_output')
#Create output subdirs
for(dir in dirs){
  dir.create(paste(wd,'output/', dir, sep = ''), showWarnings = F)
}
#Create pop dirs in subdirs
for(pop in pops_1){
  dir.create(paste(wd,'output/LD_matrix/', pop, sep = ''), showWarnings = F)
  dir.create(paste(wd,'output/LD_matrix/', pop, '/', pop, '_1Mb_coords_LDMatrix', sep = ''), showWarnings = F)
  dir.create(paste(wd,'output/gene_lists/', pop, sep = ''), showWarnings = F)
}
for(pop in pops_4){
  dir.create(paste(wd,'output/GWAS_TOPMED/', pop, sep = ''), showWarnings = F)
  dir.create(paste(wd,'output/pQTL/', pop, sep = ''), showWarnings = F)
}






#Call 5b as no hup in this script

## Run COLOC analysis  
for(i in 1:length(pop_sizes)){
  for(pheno in phenos){
    coloc_analysis(pop_1=as.character(pops_1[i]), pop_4=as.character(pops_4[i]), pop_size=pop_sizes[i], phenotype=pheno, gene_id=genes, ld_matrix=ld_matrix)
  }
}

