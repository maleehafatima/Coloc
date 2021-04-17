import os
import os.path
import argparse

'''
Set up argparser
'''
parser = argparse.ArgumentParser(description = 'Run GWAS-eQTL colocalization pipeline.')
#Add arguments
parser.add_argument('--gwas', required=True, help = 'GWAS summary statistics directory path')
parser.add_argument('--eqtl', required=True, help = 'eQTL data directory path')
parser.add_argument('--ld', required=True, help = 'LD matrices directory path')
#Add frq file dir?
parser.add_argument('--out', required=True, help = 'Specify directory for output files')
parser.add_argument('--pops1', action = 'extend', required=True, help = 'Populations used for vcf files')
parser.add_argument('--pops4', action = 'extend', required=True, help = 'Populations used for frq and eQTL files')
parser.add_argument('--pop_sizes', type = int, action = 'extend', required=True, help = 'Size of populations, in respective order to pop list')
parser.add_argument('--phenotypes', action = 'extend', required=True, help = 'Phenotypes to test')
parser.add_argument('--chrs', type = int, action = 'extend', required=True, help = 'Indicate what chromosomes to run')
#Flag to specify whether to run all genes or genes specified by user
#Maybe make default -1?
parser.add_argument('-gene_id', type = str, action = 'extend', nargs='?', default= False, help = 'Input specific genes ids to run colocalization on, default is set to False')
args = parser.parse_args()

'''
Store arguments
'''
gwas = args.gwas
eqtl = args.eqtl
ld = args.ld
out = args.out
pops1 = args.pops1
pops4 = args.pops4
pop_sizes = args.pop_sizes
phenos = args.phenotypes
chrs = args.chrs
if args.gene_id != False:
    gene_ids = args.gene_id

'''
Check/create necessary output directories
''' 
path = args.out + '/'
#Make main subdirectories
dirs = ['LD_matrix', 'gene_lists', 'GWAS_TOPMED', 'pQTL', 'Coloc_output']
for dir in dirs:
    temp_path = path + dir
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
#Make population directories in subdirectories 
for pop in pops1:
    temp_path = path + 'LD_matrix/' + pop
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
    temp_path = temp_path + '/' + pop + '_1Mb_coords_LDMatrix'
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
    temp_path = path + 'gene_lists/' + pop
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
#Make population directories in subdirectories 
for pop in pops4:
    temp_path = path + 'GWAS_TOPMED/' + pop
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
    temp_path = path + 'pQTL/' + pop
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
print('Directories created.')

'''
Run Scripts
'''
## Run scripts w/specified gene list
if args.gene_id != False:
    for i in len(pops1):
        for pheno in phenos:
            x = 1 #random line to get rid of error for now

## Run all genes in chromosomes
else:
    for i in len(pops1):
        for pheno in phenos:
            ## Scripts 1
            
            script1cmd = 'nohup Rscript 01b_run_pull_snps_driving.R $chrom $snp_annot $gene_annot $genotype_file $expression_file $pops > output/LD_matrix/nohup_1Mb_chrom$' + chr + '_2.out &'

            os.system('snp_annot=$1 #first argument is path to SNP annotation file')
            os.system('gene_annot=$2 #path to gene annotation file')
            os.system('genotype_file=$3 #path to genotype file')
            os.system('expression_file=$4 #path to the expression file')
            os.system('pops=${@:5:99}')
            os.system('for chr in chrs:')
            os.system('do')
            os.system(script1cmd)
            os.system('done')
                
            print('Pulling SNPs completed.')


            ## Scripts 2


            ## Scripts 3

            print('LD matrices created.')


            ## Script 4
            '''
            Command line arguments needed:
            args <- commandArgs(trailingOnly = TRUE)
            gwas <- args[1] #specify exact file path of gwas sum stat for this pheno in wrapper command
            pqtl <- args[2] #specify exact file path in wrapper command for this pop
            frq <- args[3] #not sure if this will work, need to specify exact file path in wrapper command
            out <- args [4] #output directory
            pop <- args[5]
            pop_size <- args[6]
            pheno <- args[7]
            chrs <- args[8]
            '''

            print('Input files formatted.')


            ## Scripts 5
            #call Rscript as nohup
            '''
            Command line arguments needed: 
            args <- commandArgs(trailingOnly = TRUE)
            gwas <- args[1] #specify exact file path of formatted gwas file in wrapper for this pop & phenotype
            pqtl <- args[2] #specify exact file path of formatted pqtl file in wrapper for this pop & phenotype
            ld <- args[3] #specify directory for ld matrices for this pop (aka pop_1Mb_coords_LDMatrix)
            out <- args [4] #output directory
            pop <- args[5]
            pop_size <- args[6]
            pheno <- args[7]
            #Way to do conditional tail argument?
            genes <- tail(args,-7)
            '''

            print('Coloc analysis finished')

print('Pipeline completed running.')
