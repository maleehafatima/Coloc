#Test command
#Rscript 01b_run_pull_snps_driving.R 3 /homes/data/COLOC_input_data/uniq_pred_db_hg38.chrchrom.maf0.01.R20.8.anno.txt.gz /homes/data/COLOC_input_data/annotation_all_aptamers_ENSG.txt /homes/data/COLOC_input_data/uniq_pred_db_hg38.chrchrom.maf0.01.R20.8.geno.txt.gz /homes/data/COLOC_input_data/annotation_all_aptamers_ENSG.txt ASW AFA CEU
#Running 01c
#./01c_run_chrom.sh /homes/data/COLOC_input_data/uniq_pred_db_hg38.chrchrom.maf0.01.R20.8.anno.txt.gz /homes/data/COLOC_input_data/annotation_all_aptamers_ENSG.txt /homes/data/COLOC_input_data/uniq_pred_db_hg38.chrchrom.maf0.01.R20.8.geno.txt.gz /homes/data/COLOC_input_data/annotation_all_aptamers_ENSG.txt ASW AFA CEU

source('01_pull_snps_driving.R') # this script pulls snps within 1Mb of the transcriptome start and end sites
"%&%" <- function(a,b) paste(a,b, sep='')
argv <- commandArgs(trailingOnly = TRUE)
chrom <- argv[1]
snp_annot_file <- argv[2]
gene_annot_file <- argv[3] 
genotype_file <- argv[4]
expression_file <- argv[3]
out <- argv[5]
pop <- argv[6] 
pip <- '0.001' #what to do with pip

main(snp_annot_file=snp_annot_file,
        gene_annot_file,
        genotype_file=genotype_file,
        expression_file=expression_file,
        eQTL_bim_file=out%&%"/LD_matrix/"%&%pop%&%"/"%&%pop%&%"_chr"%&%chrom%&%"_dose.bim",
        out_dir = out,
        pop=pop,
        chrom=as.numeric(chrom),
        null_testing=FALSE)

