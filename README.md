# pyBcfBreeder
Based on recombination maps, a pedigree and a phased vcf/bcf, pyBcfBreeder simulate new genomes.

## Quickstart
`python3 pyBcfBreeder.py --rec_maps <path_to_recombination_maps> --ped <pedigree> --vcf <vcf_with_founders_genomes> > result.vcf` \
For instance, using the data available in this repository you can do `python3 pyBcfBreeder.py --rec_maps maps/ --ped example_pedigree.ped --vcf example_vcf.bcf > result.vcf`. \
If you have bcftools installed on your machine, you can **bgzip** the output as follow: `python3 pyBcfBreeder.py --rec_maps maps/ --ped example_pedigree.ped --vcf example_vcf.bcf | bzip-c  > result.vcf`

## Input
The pedigree hould be in the .ped format defined by PLINK, the Variant Call Format file containg the founders genome can be a .vcf/.vcf.gz/.bcf/

## Dependencies
The following packages are needed for the algorithm to work:


## How does it work?
First, pyBcfBreeder use a poisson point process to simulate recombination sites in Morgan, and convert these sites into base pair using recombination maps.
Then, for each member of the pedigree, and for each position in the bcf, alleles will be copy from parental haplotypes depending on the previously generated recombination sites.
