library(dplyr)
library(data.table)

#Get & set working directory
wd <- getwd()
wd <- paste(wd,'/',sep='')

#Uses 05_coloc script to run coloc analysis
source(paste(wd,'05_coloc.R',sep=''))

"%&%" = function(a,b) paste(a,b,sep="")

#Note: should have an input list of the genes to test in which populations and which chrom they are on?

coloc_analysis <- function(pop=NULL, pop_size=NULL, phenotype=NULL, gene_id=NULL, ld_matrix =NULL){
  
  ## Reformat gene id
  gene <- gsub("\\..*","",gene_id)
  print(gene)
  #RDS is for gene list?
  rds<- wd %&% "gene_lists/" %&% pop %&% "/gene_list_" %&% pop %&% "_" %&% phenotype %&% ".RDS"
  saveRDS(gene_id,rds)
  
  ## Get pQTL/GWAS directories 
  pqtl <- paste(wd, "output/pQTL/", pop, "/pQTL_", pop, '_', phenotype, '.txt.gz', sep = '')
  gwas <- paste(wd, "output/GWAS_TOPMED/", pop, "/GWAS_TOPMED_", pop, '_', phenotype, '.txt.gz', sep = '')
  
  ## Check if pQTL/GWAS files exist
  if(!file.exists(eqtl)){
    cat(eqtl, "not exists. Exiting.\n")
    return(NA)
  }
  if(!file.exists(gwas)){
    cat(eqtl, "not exists. Exiting.\n")
    return(NA)
  }
  
  
  if(!is.null(ld_matrix)){
    
    ## Get more parameters for main function
    F_gwas <- fread(gwas, header = T, stringsAsFactors=F)
    # print(F_gwas)
    gwas_size <- F_gwas$`sample_size`[1]
    # print(gwas_size)
    LD_dir<-paste(wd, 'LD_matrix/', pop, '/', pop, '_1Mb_coords_LDMatrix', sep ='')
    
    ## Run main function from 05_coloc
    main(eqtl=eqtl, 
         gwas=gwas, 
         mode = 'bse', 
         gene_list=gene, #NOTE: need to figure out this parameter w/RDS stuff
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
         LD=ld_matrix, 
         method="cond", 
         outFile=paste(wd,'output/Coloc_output/', pop, '_', pheno, '.txt', sep = ''), ld_header = 'T')
  }
}
