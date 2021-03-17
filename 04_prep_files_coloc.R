# This is where the pQTL and GWAS SS files will be formatted for coloc v2 
# Author: Elyse Geoffroy

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))

#Output: 'ab' with no space
"%&%" <- function(a,b) paste(a,b,sep="") #'WBC', 'Platelet'

#phenos come from GWAS output file
phenos <- c('BMI')
chrs <- c(1:22)

#pops <- c('AFA', 'ALL', 'CAU', 'CHN', 'HIS')
#pop_sample_size <- c(183, 971, 416, 71, 301)

pops <- c( 'CAU')
pop_sample_size <-  c( 416)

#Maybe can put .frq in if block

for(pop in pops){ 
  #store pop's MAF & SNP data from frq file
  frq <- fread(paste("/home/egeoffroy/LD_matrix/", pops[pop], "_prot_hg38.frq", sep = ''))
  frq <- frq %>% dplyr::select(SNP, MAF)

  for(pheno in phenos){ 
    #if formatted GWAS & pQTL files don't exist, create them
    if(!file.exists("/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt") & !file.exists("/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt")){
      GWAS_result <- fread(paste("/home/egeoffroy/Wojcik/Wojcik_build38/gzipped_versions/WojcikG_", pheno, ".txt.gz", sep = '') , header = T)
    
      GWAS_result$chr_pos <- paste(gsub("chr", "", GWAS_result$chromosome), GWAS_result$base_pair_location, sep = ":")
    
      #Make formatted GWAS df
      GWAS_for_COLOC <- GWAS_result %>% dplyr::select(SNP_hg38, Beta, SE, `Effect-allele-frequency`, `Sample-size`) 
      colnames(GWAS_for_COLOC) <- c("panel_variant_id", "effect_size", "standard_error", "frequency", "sample_size")
    
      #Filter out rows with missing values (COLOC doesn't like missing data)
      GWAS_for_COLOC <- GWAS_for_COLOC[complete.cases(GWAS_for_COLOC),] 
    
      GWAS_write <- GWAS_for_COLOC
    
      #Make formatted pQTL df
      pQTL_write <- data.frame(gene_id = character(), variant_id = character(), maf = numeric(), pval_nominal = numeric(), slope = numeric(), slope_se = numeric(), stringsAsFactors = F)

      #Note: it may be more useful to also have chromosome as a parameter for the overall pipeline
      for(chr in chrs){ 
        #read in matrix eQTL results
        mpqtl <- fread("/home/egeoffroy/LD_matrix/cis_eQTLs_" %&% pops[pop] %&% "_WG_all_cis.txt.gz", nThread = 40) 
        #make your own standard error since it's not in the meQTL output
        mpqtl$se <- mpqtl$beta / mpqtl$statistic 
        #add n_samples column
        mpqtl$n_samples <- pop_sample_size[pop]
        #add freq to COLOC input
        #left_join combines mpqtl & frq dfs
        mpQTL_for_COLOC <- left_join(mpqtl, frq, by = c("snps" = "SNP")) 
        #subset to COLOC input
        mpQTL_for_COLOC <- mpQTL_for_COLOC %>% dplyr::select(gene, snps, MAF, pvalue, beta, se) 
        colnames(mpQTL_for_COLOC) <- c("gene_id", "variant_id", "maf", "pval_nominal", "slope", "slope_se")
        mpQTL_for_COLOC <- mpQTL_for_COLOC[complete.cases(mpQTL_for_COLOC),]

        pQTL_write <- rbind(pQTL_write, mpQTL_for_COLOC)
      }
    
      #Filter dfs
      snps_in_both <- intersect(GWAS_write$panel_variant_id, pQTL_write$variant_id)  
      GWAS_write <- subset(GWAS_write, panel_variant_id %in% snps_in_both)
      #print(head(GWAS_write))
      pQTL_write <- subset(pQTL_write, variant_id %in% snps_in_both)
      #order rows by gene id
      pQTL_write <- pQTL_write[order(pQTL_write$gene_id),]
      
      #write out the files that will be used to run coloc
    	fwrite(unique(pQTL_write), "/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", quote = F, sep = "\t", na = "NA", row.names = F, col.names = T)
    	gzip("/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", destname = "/home/egeoffroy/LD_matrix/coloc/pQTL_" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz") #script may only take .gz values so can't hurt to be too careful
    	fwrite(unique(GWAS_write), "/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", row.names = F, col.names = T, sep = "\t", quote = F, na = "NA")
    	gzip("/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt", "/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED" %&% pops[pop] %&% "_" %&% pheno %&% ".txt.gz")
    	print("Completed with " %&% pops[pop] %&% ", for " %&% pheno %&% ".")
    	
    } else{
        #files do exist
        next
      }
  }
  #End for loop for phenos
}
#End for loop for pops