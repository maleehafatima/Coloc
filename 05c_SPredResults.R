# Author: Elyse Geoffroy
# This script calls function in 05b to run coloc. Parameter required is population id 

library(dplyr)
library(stringr)
library(data.table)
source('/home/rschubert1/scratch/test_coloc_for_elyse/05b_run_coloc.R')

#Note: population sample sizes are hard-coded for the MESA populations. This will have to become a new parameter in COMP BIO project pipeline
#Note: will have to change this script to use an input list of genes instead of SPrediXcan results

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
print(args)
pop <- args[1]
print(pop)


files <- list.files('/home/igregga/topmed/proteome/PAGE', pattern = '_bonferroni_all_pheno.csv', recursive = F, full.names=T)
files2 <- list.files('/home/egeoffroy/topmed/SPrediXcan/Wojcik/baseline', pattern = '_all_pheno.csv', recursive = F, full.names=T)
files <- c(files, files2)
SPrediXcan_Results <- data.frame()

for(file in files){
	file <- fread(file, header= T, sep = ',', stringsAsFactors=F)
	SPrediXcan_Results <- rbind(SPrediXcan_Results, file)
}
#print(head(SPrediXcan_Results))

#pops <- c('AFA')
#pops <- c('AFA', 'ALL', 'CAU', 'CHN', 'HIS')

if(pop=='AFA'){sample_size<-183}
if(pop=='ALL'){sample_size<-971}
if(pop=='CAU'){sample_size<-416}
if(pop=='CHN'){sample_size<-71}
if(pop=='HIS'){sample_size<-301}

SPred_pop <- SPrediXcan_Results %>% filter(Model == pop) %>% select(gene, Phenotype)
SPred_pop$Phenotype <- str_replace_all(SPred_pop$Phenotype, 'LDL_choleseterol', 'LDL_cholesterol')
# print(SPred_pop)
phenos <- unique(SPred_pop$Phenotype)

#Note: this is the important part of this script, above not really needed for our pipeline construction 
for(pheno in phenos){
	data <- SPred_pop %>% filter(Phenotype == pheno)
	genes <- unique(data$gene)
	coloc_analysis(pop=as.character(pop), pop_size=sample_size, phenotype=pheno, gene_id=genes)
  
}
