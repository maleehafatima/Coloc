# This is where the pQTL and GWAS SS files will be formatted for coloc v2 
# Author: Elyse Geoffroy

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))

#Get input arguments from command line
args <- commandArgs(trailingOnly = TRUE)
gwas <- args[1] #specify exact file path of gwas sum stat for this pheno in wrapper command
pqtl <- args[2] #specify exact file path in wrapper command for this pop
frq <- args[3] #not sure if this will work, need to specify exact file path in wrapper command
out <- args [4] #output directory
pop <- args[5]
pop_size <- as.numeric(args[6])
pheno <- args[7]
chrs <- as.numeric(args[8])


#Output: 'ab' with no space
"%&%" <- function(a,b) paste(a,b,sep="") 


## If formatted GWAS & pQTL files don't exist for this pheno, create them
if(!file.exists(out %&% "/GWAS_TOPMED/" %&% pop %&% "/GWAS_TOPMED_" %&% pop %&% "_" %&% pheno %&% ".txt") 
   & !file.exists(out %&% "/pQTL/" %&% pop %&% "/pQTL_" %&% pop %&% "_" %&% pheno %&% ".txt")){
  
  ## Store pop's MAF & SNP data from frq file
  frq <- fread(frq)
  frq <- frq %>% dplyr::select(SNP, MAF)
        
  ## MAKE GWAS_write DF
  GWAS_result <- fread(gwas, header = T)
  #reformat chr_pos column
  GWAS_result$chr_pos <- paste(gsub("chr", "", GWAS_result$chromosome), GWAS_result$base_pair_location, sep = ":")
  #Filter & rename necessary columns for formatted GWAS DF
  GWAS_for_COLOC <- GWAS_result %>% dplyr::select(SNP_hg38, Beta, SE, `Effect-allele-frequency`, `Sample-size`) 
  colnames(GWAS_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")
  #Filter out rows with missing values (COLOC doesn't like missing data)
  GWAS_for_COLOC <- GWAS_for_COLOC[complete.cases(GWAS_for_COLOC),] 
  GWAS_write <- GWAS_for_COLOC
        
  ## MAKE pQTL_write DF
  pQTL_write <- data.frame(gene_id = character(), variant_id = character(), maf = numeric(), pval_nominal = numeric(), slope = numeric(), slope_se = numeric(), stringsAsFactors = F)
  for(chr in chrs){ 
    #read in matrix pQTL results
    mpqtl <- fread(pqtl, nThread = 40) 
    #make your own standard error since it's not in the mpQTL output
    mpqtl$se <- mpqtl$beta / mpqtl$statistic 
    #add n_samples column
    mpqtl$n_samples <- pop_sample_size[pop]
    #add MAF to pQTL DF (left_join combines mpqtl & frq dfs)
    mpQTL_for_COLOC <- left_join(mpqtl, frq, by = c("snps" = "SNP")) 
    #subset & rename cols needed for COLOC input
    mpQTL_for_COLOC <- mpQTL_for_COLOC %>% dplyr::select(gene, snps, MAF, pvalue, beta, se) 
    colnames(mpQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
    #Filter out rows with missing values
    mpQTL_for_COLOC <- mpQTL_for_COLOC[complete.cases(mpQTL_for_COLOC),]
    pQTL_write <- rbind(pQTL_write, mpQTL_for_COLOC)
    }
        
  ## Remove SNPs not found in both DFs
  snps_in_both <- intersect(GWAS_write$panel_variant_id, pQTL_write$variant_id)  
  GWAS_write <- subset(GWAS_write, panel_variant_id %in% snps_in_both)
  #print(head(GWAS_write))
  pQTL_write <- subset(pQTL_write, variant_id %in% snps_in_both)
  #order rows by gene id
  pQTL_write <- pQTL_write[order(pQTL_write$gene_id),]
        
  ## Write out & gzip formatted pQTL file
  fwrite(unique(pQTL_write), out %&% "/pQTL/" %&% pop %&% "/pQTL_" %&% pop %&% "_" %&% pheno %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
  gzip(out %&% "/pQTL/" %&% pop %&% "/pQTL_" %&% pop %&% "_" %&% pheno %&% ".txt", destname = out %&% "/pQTL/" %&% pop %&% "/pQTL_" %&% pop %&% "_" %&% pheno %&% ".txt.gz")         
        
  ## Write out & gzip formatted GWAS file
  fwrite(unique(GWAS_write), out %&% "/GWAS_TOPMED/" %&% pop %&% "/GWAS_TOPMED_" %&% pop %&% "_" %&% pheno %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
  gzip(out %&% "/GWAS_TOPMED/" %&% pop %&% "/GWAS_TOPMED_" %&% pops %&% "_" %&% pheno %&% ".txt", destname = out %&% "/GWAS_TOPMED/" %&% pop %&% "/GWAS_TOPMED_" %&% pop %&% "_" %&% pheno %&% ".txt.gz")
  print("Completed with " %&% pop %&% ", for " %&% pheno %&% ".")
        
} else{
  #files do exist
  next
  }

