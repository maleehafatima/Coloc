library(dplyr)
library(data.table)

#Get & set working directory
#wd <- getwd()
#wd <- paste(wd,'/',sep='')
wd <- '/homes/anovak9/Coloc/'
print(wd)

#Uses 05_coloc script to run coloc analysis
source(paste(wd,'05_coloc.R',sep=''))

"%&%" = function(a,b) paste(a,b,sep="")

coloc_analysis <- function(pop_1=NULL, pop_4=NULL, pop_size=NULL, phenotype=NULL, gene_id=NULL){
  
  #RDS to save gene list
  #rds<- wd %&% "output/gene_lists/" %&% pop_1 %&% "/gene_list_" %&% pop_1 %&% "_" %&% phenotype %&% ".RDS"
  #saveRDS(gene_id,rds)
  
  ## Get pQTL/GWAS directories 
  eqtl <- paste(wd, "output/pQTL/", pop_4, "/pQTL_", pop_4, '_', phenotype, '.txt.gz', sep = '')
  gwas <- paste(wd, "output/GWAS_TOPMED/", pop_4, "/GWAS_TOPMED_", pop_4, '_', phenotype, '.txt.gz', sep = '')
  
  ## Check if pQTL/GWAS files exist
  if(!file.exists(eqtl)){
    cat(eqtl, "not exists. Exiting.\n")
    return(NA)
  }
  if(!file.exists(gwas)){
    cat(gwas, "not exists. Exiting.\n")
    return(NA)
  }
  
  
  #if(!is.null(ld_matrix)){
    
    ## Get more parameters for main function
    F_gwas <- fread(gwas, header = T, stringsAsFactors=F)
    # print(F_gwas)
    gwas_size <- F_gwas$`sample_size`[1]
    # print(gwas_size)
    LD_dir<-paste(wd, 'output/LD_matrix/', pop_1, '/', pop_1, '_1Mb_coords_LDMatrix', sep ='')
    
    ## Run main function from 05_coloc
    main(eqtl=eqtl, 
         gwas=gwas, 
         mode = 'bse', 
         gene_list=gene_id,
         directory=LD_dir,
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
         #LD=ld_matrix, 
         method="cond", 
         outFile=paste(wd,'output/Coloc_output/', pop_1, '_', pheno, '.txt', sep = ''), ld_header = 'T')
 # }
}

pop_1 <- 'ASW'
pop_4 <- 'AFA'
pop_sizes <- 183
phenos <- c('BMI', 'C-reactive', 'HDL_cholesterol', 'Total_cholesterol')
ld_matrix <- 'T'

for(i in 1:length(pop_1)){
  for(pheno in phenos){
    coloc_analysis(pop_1=as.character(pop_1[i]), pop_4=as.character(pop_4[i]), pop_size=pop_sizes[i], phenotype=pheno)
  }
}

eqtl <- paste(wd, "output/pQTL/", pop_4, "/pQTL_", pop_4, '_', pheno, '.txt.gz', sep = '')
eqtldf<-as.data.frame(check.fread(eqtl,header=T,stringsAsFactors = F)) %>% distinct()
genes <- eqtldf[,1] %>% unlist %>% unname %>% unique()
length(ngenes)
gene <- genes[2]
LD_dir<-paste(wd, 'output/LD_matrix/', pop_1, '/', pop_1, '_1Mb_coords_LDMatrix', sep ='')
if(!is.null(list.files(LD_dir, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1])){
  ld_matrix1 <- list.files(LD_dir, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1]
} 
!is.null(list.files(LD_dir, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1])

ld_files <- list.files(LD_dir, pattern = paste('_1Mb_LD.ld.gz', sep = ''))
if(length(grep(pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), ld_files)) > 0){
  ld_matrix1 <- grep(pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), ld_files)
} else {
  print('Does not exist')
}



