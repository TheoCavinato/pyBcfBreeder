# pyBcfBreeder
Based on recombination maps, a pedigree and a phased vcf/bcf, pyBcfBreeder simulate new genomes.

## Quickstart
Just run `python3 pyBcfBreeder.py --rec_maps maps/ --ped example_pedigree.ped --vcf example_vcf.bcf | bgzip -c > test.vcf.gz`. \
If you want to use other recombination maps, replace the recombination maps in **maps** by yours.

## How does it work?
First, pyBcfBreeder use a poisson point process to simulate recombination sites in Morgan, and convert these sites into base pair using recombination maps.
Then, for each member of the pedigree, and for each position in the bcf, alleles will be copy from parental haplotypes depending on the previously generated recombination sites.
