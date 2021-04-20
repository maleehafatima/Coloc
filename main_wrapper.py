import os
import os.path
import argparse

'''
Set up argparser
'''
parser = argparse.ArgumentParser(description = 'Run GWAS-eQTL colocalization pipeline.')
#Add arguments
parser.add_argument('--gwas', required=True, help = 'GWAS summary statistics directory path')
parser.add_argument('--eqtl', required=True, help='Exact eQTL txt file path')
parser.add_argument('--vcf', required=True, help = 'Directory containing the vcf files from the eQTL population')
parser.add_argument('--snp_annot', required=True, help = 'SNP annotation files directory path')
parser.add_argument('--gene_annot', required=True, help = 'Exact gene annotation file path')
parser.add_argument('--geno', required=True, help = 'genotype files directory path')
parser.add_argument('--expr', required=True, help = 'Exact expression file path')
parser.add_argument('--frq', required = True, help = 'Exact frq file path')
parser.add_argument('--out', required=True, help = 'Specify main output directory for all output file')
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
vcf = args.vcf
ld = args.ld
snp_annot = args.snp_annot
gene_annot = args.gene_annot
geno = args.geno
expr = args.expr
frq = args.frq
out = args.out
pops1 = args.pops1
pops4 = args.pops4
pop_sizes = args.pop_sizes
phenos = args.phenotypes
chrs = args.chrs
#Convert lists from int to str
chrs = map(str, chrs)
pop_sizes = map(str, pop_sizes)
#Determine if user input list of genes
if args.gene_id != False:
    gene_ids = args.gene_id
#Get Coloc repo working directory
wd = os.getcwd()

'''
Check/create necessary output directories
''' 
path = args.out + '/'
#Make main subdirectories
dirs = ['LD_matrix', 'gene_lists', 'GWAS_TOPMED', 'eQTL', 'Coloc_output']
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
    temp_path = path + 'eQTL/' + pop
    if os.path.isdir(temp_path) == False:
        os.mkdir(temp_path)
#Save LD ouput dir
ld = path + 'LD_matrix/'

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
        
        pop = pops1[i]
        pop4 = pops4[i]
        pop_size = pop_sizes[i]
        
        for pheno in phenos:
            
            print('Run ' + pop + ' for ' + pheno + '.')

            ## Script 2
            for chr in chrs:
                #maybe implement subprocess for parallelization
                os.system("chmod u+x 02_make_bed.sh") #Make the script executable
                cmd = "./02_make_bed.sh "+pop+" "+vcf
                os.system(cmd)
            print('Bfiles made.')

            ## Scripts 1

            for chr in chrs:
                script1cmd = 'nohup Rscript 01b_run_pull_snps_driving.R ' + chr + ' ' + snp_annot + ' ' + gene_annot + ' ' + geno \
                    + ' ' + expr + ' ' + pop + ' > output/LD_matrix/nohup_1Mb_chrom' + chr + '_2.out &'
                os.system(script1cmd)

            print('Pulling SNPs completed.')

            ## Script 3
            for chr in chrs:
                os.system("chmod u+x 03_make_LD_matrix.sh")
                os.system("./03_make_LD_matrix.sh "+pop+" "+chr)
		    
            print('LD matrices created.')


            ## Script 4
            '''
            Command line arguments needed for Rscript:
            gwas <- args[1] #specify exact file path of gwas sum stat for this pheno in wrapper command
            pqtl <- args[2] #specify exact file path in wrapper command for this pop
            frq <- args[3] #not sure if this will work, need to specify exact file path in wrapper command
            out <- args [4] #output directory
            pop <- args[5]
            pop_size <- args[6]
            pheno <- args[7]
            chrs <- tail(args, -7)
            '''
            #Get list of files in gwas directory
            gwas_files = os.listdir(gwas)
            #Get gwas file specific to phenotype
            for file in gwas_files:
                if pheno in file:
                    pheno_gwas_file = file
            #Run script 4 command
            chrs_unlist = ' '.join(chrs)
            cmd = 'Rscript '+ wd + '/04_prep_files_coloc.R '+ gwas + '/' + pheno_gwas_file + ' ' + eqtl + ' ' + frq + ' ' + out \
                    + ' ' + pop4 + ' ' + pop_size + ' ' + pheno + ' ' + chrs_unlist
            os.system(cmd)

            print('Input files formatted.')


            ## Scripts 5
            #call Rscript as nohup
            '''
            Command line arguments needed for Rscript: 
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
            #Get formatted gwas file specific to pop and phenotype
            out_gwas = out + '/GWAS_TOPMED/' + pop4 
            gwas_files = os.listdir(out_gwas)
            for file in gwas_files:
                if pheno in file:
                    if pop4 in file: 
                        pheno_pop_gwas = file
            #Get formatted eqtl file specific to pop & phenotype
            out_eqtl = out + '/eQTL/' + pop4
            eqtl_files = os.listdir(out_eqtl)
            for file in eqtl_files:
                if pheno in file:
                    if pop4 in file:
                        pheno_pop_eqtl = file
            #Get LD dir for this pop
            pop_ld = ld + pop + '/' + pop + '_1Mb_coords_LDMatrix'
            #Run script 5 command
            cmd = 'nohup Rscript ' + wd + '/05b_run_coloc.R' + pheno_pop_gwas + ' ' + pheno_pop_eqtl + ' ' + pop_ld + ' ' + out \
                    + ' ' + pop + ' ' + pop_size + ' ' + pheno
            os.system(cmd)

            print('Coloc analysis finished')


print('Pipeline completed running.')
