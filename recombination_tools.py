import random, math, numpy, sys

#########################################
#   TOOLS TO READ RECOMBINATION MAPS    #
#########################################
def read_one_rec_map(file_path, rec_maps):
    rec_map_file=open(file_path,'r')
    print(file_path)

    header=rec_map_file.readline().strip().split('\t')
    col_cM=header.index("cM")
    col_pB=header.index("pos")
    col_chr=header.index("chr")

    pos_M=[]
    pos_pB=[]

    rec_map_info={}

    while True:
        line=rec_map_file.readline()
        line_splt=line.strip().split("\t")
        if not line:
            rec_maps[str(chr)]=(pos_M,pos_pB)
            break

        pos_M.append(float(line_splt[col_cM])/100)
        pos_pB.append(float(line_splt[col_pB]))
        chr=line_splt[col_chr]
        if chr in ("chrX", 'X'):
            break
    
    return rec_map_info

# read the recombination maps and store it in a dictionnary, where
# each key is a chromosome
# each value is a tuple where first element is the positions in Morgan, second element is the positions in bp
def read_all_read_maps(map_files):
    rec_maps={}
    for rec_map in map_files:
        read_one_rec_map(rec_map, rec_maps)
    return(rec_maps)

############################################
#   TOOLS TO SIMULATE RECOMBINATION SITES  #
############################################

class Chr:
    def __init__(self, chr_size_M, rng):
        self.chr_size_M=chr_size_M
        self.rng = rng
        self.rec_pos_M=[]
        self.generate_rec_pos()
        self.rec_pos_bp=[]

    def inverse_CDF(self,x):
        return -(math.log(1-x)/1)

    def generate_rec_pos(self):
        rec_pos=0
        iterator=0
        if self.chr_size_M!=0:
            while True:
                if rec_pos==0:
                    rec_pos=self.inverse_CDF(self.rng.uniform(0,1))
                else:
                    rec_pos=self.rec_pos_M[(iterator-1)]+self.inverse_CDF(self.rng.uniform(0,1))
                if rec_pos > self.chr_size_M:
                    break
                self.rec_pos_M.append(rec_pos)
                iterator+=1
    
    def convert_cM_to_bp(self, list_bp, list_M):
        for rec_pos in self.rec_pos_M:
            pos = numpy.searchsorted(list_M, rec_pos)
            conversion=self.linear_conversion(rec_pos, list_M[pos-1], list_M[pos], list_bp[pos-1], list_bp[pos])
            self.rec_pos_bp.append(conversion)

    def linear_conversion(self, X, A_1, A_2, B_1, B_2):
        return ((X-A_1)/(A_2-A_1))*(B_2-B_1)+B_1


#LINK RECOMBINATION SITES TO INDIVIDUALS FROM THE PEDIGREE
def ped_to_rec(sorted_ped_path, rec_maps, show):
    sorted_ped=open(sorted_ped_path, 'r')
    if show:
        output_file=open(sorted_ped_path+".rec",'w')

    child_parents, child_rec_sites, all_childs = {}, [], []

    #put 0 for lacking recombination maps
    chromosomes = [str(x) for x in range(1,23)]
    for chromosome in chromosomes:
        if chromosome not in rec_maps.keys():
            rec_maps[chromosome] = ([0],[0])

    #retrieve chromosome sizes about the rec_maps:
    all_chr_sizes_M = {chr:max(rec_maps[chr][0]) for chr in rec_maps.keys()}
    
    for line in sorted_ped:

        new_rec_line=line.split()
        child_id, parent_1_ID, parent_2_ID = new_rec_line[1:4]

        if parent_1_ID!='NA' and parent_2_ID!='NA': # if the child is not a founder

            rand_seeds=[]
            for _ in range(2): #one recombination list per parent
                parental_rec_sites=[]

                # set a simulation seed for the individual
                rand_seed = random.randrange(sys.maxsize)
                rng_chr = random.Random(rand_seed)

                for chromosome in chromosomes:
                    chromosome_size_M=all_chr_sizes_M[chromosome]
                    curr_chr_rec_map_M=rec_maps[chromosome][0]
                    curr_chr_rec_map_bp=rec_maps[chromosome][1]
                    new_chromosome=Chr(chromosome_size_M, rng_chr)
                    new_chromosome.convert_cM_to_bp(curr_chr_rec_map_bp, curr_chr_rec_map_M)
                    parental_rec_sites.append(new_chromosome.rec_pos_bp)

                new_rec_line.append(parental_rec_sites)
                rand_seeds.append(str(rand_seed))

            new_rec_line.extend(rand_seeds)

            all_childs.append(child_id)
            child_parents[child_id]=(parent_1_ID, parent_2_ID)
            child_rec_sites.append((new_rec_line[6], new_rec_line[7]))

            if show:
                new_rec_line[6] = str(new_rec_line[6])
                new_rec_line[7] = str(new_rec_line[7])
                output_file.write("\t".join(new_rec_line))
                output_file.write("\n")
    sorted_ped.close()
    if show:
        output_file.close()
    return  child_parents, child_rec_sites, all_childs

#################################
#   TOOLS TO WRITE A NEW VCF    #
#################################

def haplo_chooser(meiosis,toss_result, position):
    haplo_chooser=toss_result
    if len(meiosis):
            for rec_site in meiosis:
                    if rec_site > position:
                            break
                    haplo_chooser+=1

    if haplo_chooser%2:
            return 0
    else:
            return 1

def tossing_coins(out_path, all_childs, show):
    if show:
        metadata_writer=open(out_path+".meta.txt", 'w')
    toss_a_coin_per_chr = {}
    for chr in range(1,23):
        rand_seed = random.randrange(sys.maxsize)
        rng_chr = random.Random(rand_seed)
        toss_a_coin=[( rng_chr.randint(0,1), rng_chr.randint(0,1) ) for _ in range(len(all_childs))]
        toss_a_coin_per_chr[str(chr)]=toss_a_coin
        if show:
            metadata_writer.write(str(rand_seed))
            metadata_writer.write("\n")
    if show:
        metadata_writer.close()
    return toss_a_coin_per_chr

#####################################################
#  TOOLS TO SIMULATE PL AND DP (--coverage option)  #
#####################################################

class Coverage:
    def __init__(self, p_coverage, p_error, p_samples, p_header):
        self.coverage = float(p_coverage)
        self.error = p_error
        self.samples = p_samples
        self.samples_nbr = len(p_samples)
        self.output_file = open("PL_"+p_coverage+"x.vcf", 'w')
        self.new_line = []
        self.write_header(p_header)
    
    def close_file(self):
        self.output_file.close()
    
    def write_header(self, split_header):
        format_checker = False
        for line in split_header[:-2]:
            if line[:8] == "##FORMAT":
                format_checker=True
            if line[:8] != "##FORMAT" and format_checker:
                self.output_file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">')
                self.output_file.write('\n')
                self.output_file.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt">')
                self.output_file.write('\n')
                format_checker=False

            self.output_file.write(line)
            self.output_file.write('\n')
        last_line = split_header[-2].split('\t')[:9]
        for sample in self.samples:
            last_line.append(sample)
        self.output_file.write('\t'.join(last_line))
        self.output_file.write('\n')
    
    def emission_probas(self, H1, H2, read_allele):
        if H1 == H2 == read_allele:
            return math.log10(1 - self.error)
        elif H1 != read_allele and H2 != read_allele:
            return math.log10(self.error/3)
        else:
            return math.log10(1/2 - self.error/3)

    def GL_to_PL(self, GLs):
        log_values = [-10*GL for GL in GLs]
        min_value = min(log_values)
        PLs = [str(int(round(log_value - min_value))) for log_value in log_values]
        return ','.join(PLs)

    def simulate_DP_PL(self, GT):
        #simulate a number of reads covering the position
        DP = numpy.random.poisson(self.coverage, 1)[0]
        if DP==0:
            PL = "0,0,0"
            GT_str = "./."

        else:
            GT_int = [int(hap) for hap in GT]
            GT_str = '/'.join(GT)

            H0_nbr = numpy.random.binomial(DP, 0.5) #X alleles for hap0
            H1_nbr = DP - H0_nbr #DP - X alleles for hap1

            #calculate the likelihoods of the genotypes 00, 01 and 11
            p_D_00 = self.emission_probas( 0, 0, GT_int[0])*H0_nbr + self.emission_probas( 0, 0, GT_int[1])*H1_nbr
            p_D_01 = self.emission_probas( 0, 1, GT_int[0])*H0_nbr + self.emission_probas( 0, 1, GT_int[1])*H1_nbr
            p_D_11 = self.emission_probas( 1, 1, GT_int[0])*H0_nbr + self.emission_probas( 1, 1, GT_int[1])*H1_nbr

            PL = self.GL_to_PL((p_D_00, p_D_01, p_D_11))

        return GT_str + ':' + str(DP)+':'+PL
    
    def start_new_line(self, line):
        self.new_line = line
    
    def add_genotype(self, out_haplo):
        self.new_line.append(self.simulate_DP_PL(out_haplo))
    
    def write_new_line(self):
        self.output_file.write('\t'.join(self.new_line))
        self.output_file.write('\n')