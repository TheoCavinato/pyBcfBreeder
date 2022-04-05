# pyBcfBreeder
Based on recombination maps, a pedigree and a phased vcf/bcf, pyBcfBreeder simulate new genomes.

## Quickstart
In this example, you have a pedigree **example_ped.ped**, a VCF containing the genomes of the founders of this pedigree **example_vcf.bcf** and the recombination maps you want to use to simulate recombination are located in **maps/**. You can run: \
`python3 pyBcfBreeder.py --rec_maps maps/ --ped example_pedigree.ped --vcf example_vcf.bcf > result.vcf`. \
to obtain a new vcf **result.vcf** containing simulated genomes of the offspring of the pedigree with genomes of the founders.
If you have bcftools installed on your machine, you can **bgzip** the output as follow: \
`python3 pyBcfBreeder.py --rec_maps maps/ --ped example_pedigree.ped --vcf example_vcf.bcf | bgzip -c  > result.vcf.gz`

## Input
The pedigree hould be in the .ped format defined by PLINK, the Variant Call Format file containg the founders genome can be a .vcf/.vcf.gz/.bcf/

## Dependencies
The following packages are needed for the algorithm to work:
`numpy pysam`


## How does it work?
First, pyBcfBreeder use a poisson point process to simulate recombination sites in Morgan, and convert these sites into base pair using recombination maps.
Then, for each member of the pedigree, and for each position in the bcf, alleles will be copy from parental haplotypes depending on the previously generated recombination sites.
