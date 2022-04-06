from recombination_tools import *
from pysam import VariantFile
import argparse, os, random, sys

#########################
#   USER PARAMETERS     #	
#########################

parser = argparse.ArgumentParser(description = 'To do')
parser.add_argument('-R', '--rec_maps', help = 'folder containing all the recombination maps', required=True)
parser.add_argument('-P', '--ped', help = 'ped (pedigree) to simulate', required=True)
parser.add_argument('-V', '--vcf', help = 'vcf containing genomes of the founders of the pedigree', required=True)
parser.add_argument('--meta', help = 'Set this option to "True" if you want the algorithm to output the simulated \
    recombination sites and the seeds used for haplotype choosing', action="store_true")
parser.add_argument('--coverages', help = 'Pass a list of coverages to this argument if you want to get \
    simulated Genotype Likelihoods for each offspring at each position', nargs='+') 
parser.add_argument('--error', help = 'error rate for the Genotype Likelihood computation', type=float, \
    default = 0.0001)
args = parser.parse_args()

####################################
#  SIMULATING RECOMBINATION SITES  #
####################################

map_files=[args.rec_maps+file for file in os.listdir(args.rec_maps)] #list the recombination map available
recombination_maps=read_all_read_maps(map_files) #read the recombination maps and store it in a dictionnary

child_parents, child_rec_sites, all_childs = ped_to_rec(args.ped, recombination_maps, args.meta)

###############
#  WRITE VCF  #
###############


#Changing the format of the simulated data
child_rec_sites_reversed=[[] for _ in range(len(child_rec_sites[0][0]))]
for chr in range(len(child_rec_sites[0][0])):
	for ind_rec_site in child_rec_sites:
		child_rec_sites_reversed[chr].append((ind_rec_site[0][chr], ind_rec_site[1][chr]))

bcfInput = VariantFile(args.vcf)

#Only read information from the samples we have interest in
rec_parents = {parents[0] for parents in child_parents.values()} | {parents[1] for parents in child_parents.values()} 
bcfInput_parents = list(set(bcfInput.header.samples) & rec_parents)
bcfInput.subset_samples(bcfInput_parents)

#Write the new header
original_header=str(bcfInput.header)[:-1]
original_header+='\t'+'\t'.join(all_childs)
print(original_header)

#Useful variables
output_samples = list(bcfInput.header.samples) + all_childs
original_alleles_length = len(bcfInput_parents) 
toss_a_coin_per_chr = tossing_coins(args.ped+".rec", all_childs, args.meta)
parent_position=[(output_samples.index(child_parents[child][0]), output_samples.index(child_parents[child][1])) for child in all_childs]

#If --coverages
if args.coverages:
    coverage_header = str(bcfInput.header).split('\n')
    coverage_objs = [Coverage(coverage, args.error, all_childs, coverage_header) for coverage in args.coverages]

for record in bcfInput:
    if record.chrom[:3] == "chr":
        chromosome = record.chrom.split("chr")[1]
    else:
        chromosome = record.chrom
    chromosome_pos =  int(chromosome) - 1
    position = record.pos
    alleles=[al['GT'] for al in record.samples.values()]

    if args.coverages:
        split_rec = str(record).split("\t")[:8]
        split_rec.append("DP:PL")

    #-------Take the allele resulting from the meiosis for each parent-------#
    for parent_rec_sites, pos_parent, tossing in zip(child_rec_sites_reversed[chromosome_pos], parent_position, toss_a_coin_per_chr[chromosome]):
            out_haplo=tuple(str(alleles[parent][haplo_chooser(rec_sites, toss, position)]) for rec_sites, parent, toss in zip(parent_rec_sites, pos_parent, tossing))
            alleles.append(out_haplo)

    if args.coverages:
        split_rec = str(record).split("\t")[:8]
        split_rec.append('DP:PL')
        for cov in coverage_objs:
            cov.simulate_new_line(split_rec.copy())
                    
    #-------Add alleles to the record-------#

    output_record=str(record)[:-1]
    output_record+='\t'+'\t'.join(['|'.join(all) for all in alleles[original_alleles_length:]])
    print(output_record)
bcfInput.close()

#Close coverages files
if args.coverages:
    for cov in coverage_objs:
        cov.close_file()