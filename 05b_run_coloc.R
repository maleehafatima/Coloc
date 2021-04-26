library(dplyr)
library(data.table)

#Get input arguments from command line
args <- commandArgs(trailingOnly = TRUE)
gwas <- args[1] #specify exact file path of formatted gwas file in wrapper for this pop & phenotype
eqtl <- args[2] #specify exact file path of formatted eqtl file in wrapper for this pop & phenotype
ld <- args[3] #specify directory for ld matrices for this pop (aka pop_1Mb_coords_LDMatrix)
out <- args [4] #output directory
pop <- args[5]
pop_size <- as.numeric(args[6])
pheno <- args[7]
genes <- tail(args,-7) #if user input gene ids

#Get & set working directory (tell user to cd into Coloc repo)
wd <- getwd()
wd <- paste(wd,'/',sep='')

#Uses 05_coloc script to run coloc analysis
source(paste(wd,'05_coloc.R',sep=''))

"%&%" = function(a,b) paste(a,b,sep="")

#If user input a list of genes
if(length(genes) != 0){
  #RDS to save gene list
  rds<- out %&% "/gene_lists/gene_list_" %&% pop %&% "_" %&% pheno %&% ".RDS"
  saveRDS(genes,rds)
}else{rds <- NULL}


## Check if eQTL/GWAS files exist
if(!file.exists(eqtl)){
  cat(eqtl, "not exists. Exiting.\n")
  return(NA)
}
if(!file.exists(gwas)){
  cat(gwas, "not exists. Exiting.\n")
  return(NA)
}
  
    
## Get more parameters for main function
F_gwas <- fread(gwas, header = T, stringsAsFactors=F)
# print(F_gwas)
gwas_size <- F_gwas$`sample_size`[1]
# print(gwas_size)
    
## Run main function from 05_coloc
main(eqtl=eqtl, 
     gwas=gwas, 
     mode = 'bse', 
     gene_list=rds,
     directory=ld,
     eqtlGeneCol='gene_id', 
     eqtlSNPCol='variant_id', 
     eqtlMAFCol='maf', 
     eqtlSeCol=6, 
     eqtlBetaCol=5, 
     eqtlSampleSize=pop_size, 
     gwasSNPCol=1,
     gwasBetaCol=2, 
     gwasSeCol=3, 
     gwasSampleSize=gwas_size,
     method="cond", 
     outFile=paste(out,'/Coloc_output/', pop, '_', pheno, '.txt', sep = ''), ld_header = 'T')






