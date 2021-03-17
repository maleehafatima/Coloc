setwd("/home/egeoffroy/LD_matrix/")
source('01_pull_snps_driving.R') # this script pulls snps within 1Mb of the transcriptome start and end sites
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
pops = 'ASW'
pip <- '0.001'

for(pop1 in pops){
snp_annot_file <- "/home/rschubert1/data/TOPMED_Proteome/" %&% pop1 %&% "/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.anno.txt.gz"
gene_annot_file <- "/data/rschubert1/TOPMED_Proteome/annotation_all_aptamers_ENSG.txt"
genotype_file <- "/home/rschubert1/data/TOPMED_Proteome/" %&% pop1 %&% "/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.geno.txt.gz"
expression_file <- "/data/rschubert1/TOPMED_Proteome/annotation_all_aptamers_ENSG.txt"

main(snp_annot_file=snp_annot_file,
        gene_annot_file,
        genotype_file=genotype_file,
        expression_file=expression_file,
        pop=pop1,
        chrom=as.numeric(chrom),
        null_testing=FALSE)
}

