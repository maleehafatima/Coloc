
suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
suppressMessages(library(doMC))
suppressMessages(library(doRNG))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))

"%&%" <- function(a,b) paste(a,b, sep = "") 
#Shortened notation for concatenating strings

## Couldn't use these line in pipeline due to confidentiality issues 
# get_gene_expression <- function(gene_expression_file_name, gene_annot) { #row are obs, columns are features
#   expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = NULL)))
#   expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df)))) #%>% mutate(id=gsub("\\.[0-9]+","",id))
#   expr_df
# }

get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- read.table(snp_annot_file_name, header = T, stringsAsFactors = F) %>%
    filter(!((refAllele == 'A' & effectAllele == 'T') |
               (refAllele == 'T' & effectAllele == 'A') |
               (refAllele == 'C' & effectAllele == 'G') |
               (refAllele == 'G' & effectAllele == 'C')) &
             !(is.na(rsid))) %>%
    distinct(varID, .keep_all = TRUE)
  snp_annot
}

get_maf_filtered_genotype <- function(genotype_file_name) {
  gt_df <- read.table(genotype_file_name, header = T, stringsAsFactors = F) %>% distinct(snp_ID,.keep_all=T) %>% column_to_rownames(var="snp_ID")
  gt_df
} 

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding',"aptamer", 'pseudogene', 'lincRNA',"aptamer","VALUE")){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types) 
  gene_df
}

## Couldn't use these line in pipeline due to confidentiality issues 
# get_gene_type <- function(gene_annot, gene) {
#   filter(gene_annot, gene_id == gene)$gene_type
# }

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window) & !is.na(rsid)) & (pos <= (coords[2] + cis_window)))
  print(c('number of snps: ', nrow(snp_info)))
  if (nrow(snp_info) == 0)
    return(NA)
  
   snp_info$SNP <- paste(snp_info$rsid, snp_info$refAllele, snp_info$effectAllele, sep = ':')
   cis_gt <- snp_info$SNP
   names(cis_gt) <- NULL
   print(cis_gt)
}

main <- function(snp_annot_file, gene_annot_file, genotype_file, expression_file,
                 covariates_file=NULL, chrom, pop, maf=0.01, n_folds=10, n_train_test_folds=5,
                 seed=NA, cis_window=1000000, alpha=0.5, null_testing=FALSE) {
  

  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  
  ## Couldn't use these line in pipeline due to confidentiality issues 
  #expr_df <- get_gene_expression(expression_file, gene_annot)
  #genes <- colnames(expr_df)
  #samples <- rownames(expr_df)
  
  genes <- list(gene_annot$gene_id)
  n_genes <- length(genes)
  snp_annot <- get_filtered_snp_annot(snp_annot_file)
  gt_df <- get_maf_filtered_genotype(genotype_file)
  for (i in 1:n_genes) {
    cat(i, "/", n_genes, "\n")
    gene <- unlist(genes[i])
    print(gene)
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    coords <- get_gene_coords(gene_annot, gene)
    
    cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)
    str(gene)
    str(gene_type)
    str(coords)
    str(cis_gt)
    
    print(cis_gt)
    if(length(cis_gt) > 2){
    	write.table(cis_gt, paste('output/LD_matrix/',pop,'/', pop, '_chr_', chrom, '_', gene, '_1Mb_of_gene.txt', sep = ''), quote = F, row.names=F, col.names=F)
    }
  }
}
    
