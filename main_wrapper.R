library(dplyr)
library(data.table)

#Get & set working directory
wd <- getwd()
wd <- paste(wd,'/',sep='')

#Import necessary scripts
source(paste(wd,'05_coloc.R',sep=''))
source(paste(wd, '05b_run_coloc.R', sep=''))

#make output folder 

#Import arguments 
args <- commandArgs(trailingOnly=TRUE)
print(args)
pops <- args[1]
phenos <- args[2]
eqtl <- args[3] #NOTE: probs don't need this argument
gwas<- args[4] #NOTE: probs don't need this argument
LD_dir<- args[5] #NOTE: probs don't need this argument
gene_id<-args[6]
pop_sizes<-as.numeric(args[7])
ld_matrix<-args[8]
gwas_size<-as.numeric(args[9]) #NOTE: probs don't need this argument

genes<-readRDS(gene_id) #NOTE: need to figure out how this works with 5b





#Call 5b as no hup in this script

#(From script 5c) just becomes part of end of wrapper 
for(i in 1:length(pops)){
        for(pheno in phenos){
                coloc_analysis(pop=as.character(pops[i]), pop_size=pop_sizes[i], phenotype=pheno, gene_id=genes, ld_matrix=ld_matrix)
        }
}

