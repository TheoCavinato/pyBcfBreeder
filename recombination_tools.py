import random, math, numpy, sys

#########################################
#   TOOLS TO READ RECOMBINATION MAPS    #
#########################################
def read_one_rec_map(file_path, rec_maps):
    rec_map_file=open(file_path,'r')

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
            rec_maps[str(chr)]=[pos_M,pos_pB]
            break

        pos_M.append(float(line_splt[col_cM])/100)
        pos_pB.append(float(line_splt[col_pB]))
        chr=line_splt[col_chr]
        if chr in ("chrX", 'X'):
            break
    
    return rec_map_info

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


#GIVE RECOMBINATION SITES TO NEW INDIVIDUALS FROM THE PEDIGREE
def ped_to_rec(sorted_ped_path, rec_maps, show):
    sorted_ped=open(sorted_ped_path, 'r')
    if show:
        output_file=open(sorted_ped_path+".rec",'w')

    child_parents, child_rec_sites, all_childs = {}, [], []

    #retrieve chromosome sizes about the rec_maps:
    all_chr_sizes_M = {chr:max(rec_maps[chr][0]) for chr in rec_maps.keys()}
    
    for line in sorted_ped:
        new_rec_line=line.split()
        child_id, parent_1_ID, parent_2_ID = new_rec_line[1:4]
        if parent_1_ID!='NA' and parent_2_ID!='NA': # if we know the parents of the children
            rand_seeds=[]
            for _ in range(2): #one recombination list per parent
                parental_rec_sites=[]

                # sort chromosomes by assending order
                chromosomes = [int(chr) for chr in rec_maps.keys()]
                chromosomes.sort()
                
                #set a seed for the individual
                rand_seed = random.randrange(sys.maxsize)
                rng_chr = random.Random(rand_seed)

                for chromosome in chromosomes:
                    
                    chromosome=str(chromosome)
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

