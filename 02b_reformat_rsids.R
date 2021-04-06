suppressMessages(library(dplyr))
"%&%" <- function(a,b) paste(a,b, sep = "") 

argv <- commandArgs(trailingOnly = TRUE)
bim <- argv[1] #first argument in the path to directory with the bim file
chr <- argv[2] #second argument is the chromosome
pop <- argv[3] #third argument is the population abbreviation

SNP_info <- read.table(bim%&%"old_bims/"%&%pop%&%"_chr"%&%chr%&%"_dose.bim") %>% 
  rename(chr=V1,rsid=V2,pos=V3,bp=V4,allele1=V5,allele2=V6) %>%
  mutate(rsid = "chr"%&%chr%&%":"%&%bp%&%":"%&%allele1%&%":"%&%allele2,.keep="used",.after=chr)

write.table(SNP_info,file=bim%&%pop%&%"_chr"%&%chr%&%"_dose.bim",quote=F,row.names=F,col.names=F,sep="\t")
