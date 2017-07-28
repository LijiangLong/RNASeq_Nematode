import pysam, time, operator, sys, scipy, collections, shutil, os, pyfaidx, vcf
from functools import reduce
from statistics import median
from collections import defaultdict
from skbio.alignment import StripedSmithWaterman
from subprocess import call, Popen, PIPE
from operator import attrgetter, itemgetter
from Bio import pairwise2
from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM


strand_converter = {True: -1,
                    False: 1,
                    '+': 1,
                    '-': -1}

COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}                
    
def reverse_complement(sequence):
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])
   
def coverage(bam_obj):
    stats = pysam.idxstats(bam_obj.filename).rstrip().split('\n')
    tot_reads = sum([int(x.split('\t')[2]) for x in stats])

    tot_bases = sum(bam_obj.lengths)
    count = 0
    read_len = []
    for read in bam_obj.fetch():
        read_len.append(len(read.seq))
        count += 1
        if count > 100000:
            break
    return tot_reads/tot_bases*median(read_len)

def avg_read_len(bam_obj):
    count = 0
    read_len = []
    for read in bam_obj.fetch():
        read_len.append(len(read.seq))
        count += 1
        if count > 100000:
            break
    return median(read_len)


def avg_coverage(bam_obj, chrom, pos1, pos2):
    n_reads = 0
    if pos2 - pos1 > 10000:
        return bam_obj.count(chrom, pos1, pos2)*avg_read_len(bam_obj)/(pos2-pos1+1)
    else:
        for puc in bam_obj.pileup(chrom, pos1, pos2):
            if puc.pos >= pos1 and puc.pos <= pos2:
                for pur in puc.pileups:
                    if pur.is_del or pur.is_refskip:
                        pass
                    else:
                        n_reads += 1
    return n_reads/(pos2 - pos1 + 1)

def genotype_bam(bam_obj, chrom, pos):
    pu = bam_obj.pileup(chrom, pos, pos+1)
    for puc in pu:
        if puc.pos == pos:
            return genotype_pucpileup(puc.pileups)
    return []

def genotype_pucpileup(puc_pileup):
    base_g = []
    for pur in puc_pileup:
        if pur.is_del or pur.is_refskip:
            continue
        elif pur.indel < 0:
            qp = pur.query_position
            base_g.append(pur.alignment.seq[qp] + str(pur.indel) )
        else:
            qp = pur.query_position
            if pur.indel > 0:
                base_g.append(pur.alignment.seq[qp:qp+1+pur.indel])
            else:
                base_g.append(pur.alignment.query_sequence[qp])
    return base_g

def overlap_bam(bam_obj, chrom, pos):
    reads = bam_obj.fetch(chrom, pos, pos+1)
    out_dist = []
    clipped = 0
    for read in reads:
        if not read.is_unmapped and not read.is_secondary and read.mapq > 0:
            if 'S' in read.cigarstring:
                clipped += 1
            start = read.get_blocks()[0][0]
            stop = read.get_blocks()[-1][1]-1
            out_dist.append(min([pos - start, stop - pos]))
    return out_dist, clipped

    
class clipped_obj():
    def __init__(self, bam_obj, ref_obj):
        self.bam_obj = bam_obj
        self.ref_obj = ref_obj
        self.chimera_reads = defaultdict(list)
        self.single_reads = defaultdict(list)
        self.close_hit = {}
        self.chimera_locs = set()

        
        for read in bam_obj.fetch():
            if not read.is_secondary and read.mapq > 10:
                new_read = clipped_read(read, self.bam_obj, self.ref_obj)
                if new_read is None:
                    continue
                if new_read.position == (None,None,None,None):
                    continue
                if new_read.is_simple:
                    if new_read.is_chimera:
                        self.chimera_reads[new_read.position].append(new_read)
                        self.chimera_locs.add((new_read.original_pos[0],None,1,''))
                        self.chimera_locs.add((None,new_read.original_pos[1],1,''))
                        
                    #else:
                    #    self.single_reads[new_read.position].append(new_read)

    def ret_long_clip(self, position):
        lengths = []
        consensus = []
        for cl_read in self.single_reads[position]:
            lengths.append(len(cl_read.clipped_seq))
        if len(lengths) == 0:
            return ''
        max_len = max(lengths)
        for cl_read in self.single_reads[position]:
            if len(cl_read.clipped_seq) == max_len:
                return cl_read
                   
                        
    def ret_consensus(self, position):
                   
        lengths = []
        consensus = []
        if position[0] is not None:
            index = 1
        if position[1] is not None:
            index = 0
        for cl_read in self.single_reads[position]:
            lengths.append(len(cl_read.clipped_seq[index]))
        if len(lengths) == 0:
            return ''
        for i in range(0,max(lengths)):
            temp_dict = defaultdict(int)
            total = 0
            for cl_read, length in zip(self.single_reads[position], lengths):
                if i < length:
                    temp_dict[cl_read.clipped_seq[index][i*(2*index-1)]] += 1
                    total += 1
                    cons_char = max(iter(temp_dict.items()), key=operator.itemgetter(1))[0]
            if temp_dict[cons_char]/total >= .74 and total > 1:
                consensus.append(cons_char)
            elif total == 1:
                consensus.append('')
            else: 
                consensus.append('N')
        consensus = ''.join(consensus)
        return consensus

    def ret_nearby_cov(self, bam_obj, position, delta):
        if position[0] == None:
            chrom = bam_obj.references[position[1][0]]
            stop = position[1][1]
            return avg_coverage(bam_obj, chrom, stop-delta, stop)
        else:
            chrom = bam_obj.references[position[0][0]]
            start = position[0][1]
            return avg_coverage(bam_obj, chrom, start, start+delta)
                   
    def ret_position(self, position):
        if position[0] is not None:
            return position[0][0]*25000000 + position[0][1]
        else:
            return position[1][0]*25000000 + position[1][1]


class clipped_read():
    def __init__(self, read, s_bam, ref_obj):
        self.read = read
        self.position = (None,None,None,None)
        self.is_clipped = False
        self.is_simple = False
        self.is_chimera = False
        self.clipped_seq = ''
        self.is_del = False
        self.is_inv = False
        self.is_dup = False
        self.is_transloc = False
        self.original_pos = (None,None)
        flipped = 1
        #Is it clipped?
        if ('S' in read.cigarstring or 'H' in read.cigarstring) and 'M' in read.cigarstring:
            self.is_clipped = True
        if self.is_clipped and len(read.cigar) == 2:
            self.is_simple = True
        if self.is_simple:
            try:
                SA_tag = read.get_tag('SA').split(',')
                self.is_chimera = True
                chi_cigarstring = SA_tag[3]
                if chi_cigarstring.count('S') != 1 or chi_cigarstring.count('M') != 1:
                    self.is_chimera = False
            except KeyError:
                pass
            #Calculate location
            if read.cigar[0][0] == 4:
                # e.g. 34S67M
                right_chr = read.rname
                right_pos = read.pos # Already 0-coordinate system
                if self.is_chimera:
                    right_strand = strand_converter[read.is_reverse]
                    left_chr = s_bam.gettid(SA_tag[0])
                    left_strand = strand_converter[SA_tag[2]]
                    if left_strand == right_strand:
                        # e.g. 36M65S
                        chi_mapped_bases = int(chi_cigarstring.split('M')[0])
                        left_pos = int(SA_tag[1]) + chi_mapped_bases - 2 # Convert 1-coordinate to 0-coordinate
                    elif left_strand != right_strand:
                        # e.g. 65S36M
    
                        chi_mapped_bases = int(chi_cigarstring.split('S')[1].split('M')[0])
                        t_left_pos = int(SA_tag[1]) - 1 # Convert 1-coordinate to 0-coordinate
                        left_pos = min(right_pos, t_left_pos)
                        right_pos = max(right_pos, t_left_pos)-1
                        flipped = -1
                        
                    mapped = (chi_mapped_bases-1,read.cigar[0][1])

                else:
                    self.position = (None,(right_chr,right_pos), 1, '')
                    self.clipped_seq = read.seq[:read.cigar[0][1]]
        

            elif read.cigar[0][0] == 0:
                # e.g. 35M65S
                left_chr = read.rname
                left_pos = read.pos + read.cigar[0][1] - 1 # Already 0-coordinate
                if self.is_chimera:
                    left_strand = strand_converter[read.is_reverse]
                    right_chr = s_bam.gettid(SA_tag[0])
                    right_strand = strand_converter[SA_tag[2]]
                    if left_strand == right_strand:
                        # e.g. 34S67M
                        try:
                            chi_mapped_bases = int(chi_cigarstring.split('S')[1].split('M')[0])
                        except ValueError:
                            return None
                        right_pos = int(SA_tag[1]) - 1 # Convert 1-coordinate to 0-coordinate
                    elif left_strand != right_strand:
                        # e.g. 36M65S
                        try:
                            chi_mapped_bases = int(chi_cigarstring.split('M')[0])
                        except ValueError:
                            return None
                        t_right_pos = int(SA_tag[1]) + chi_mapped_bases - 2 # Convert 1-coordinate to 0-coordinate
                        right_pos = max(t_right_pos, left_pos)
                        left_pos = min(t_right_pos, left_pos)+1

                    mapped = (read.cigar[0][1]-1,len(read.seq) - chi_mapped_bases)
                else:
                    self.position = ((left_chr,left_pos), None, 1, '')
                    self.clipped_seq = read.seq[read.cigar[0][1]:]
            if self.is_chimera:
                self.original_pos = ((left_chr,left_pos), (right_chr, right_pos))
                overlap = mapped[1] - mapped[0] - 1
                if overlap > 0:
                    ins_bases = read.seq[mapped[0]+1:mapped[1]]
                else:
                    left_pos = left_pos + overlap*flipped
                    mapped = (mapped[0] + overlap, mapped[1])
                    ins_bases = ''
                self.position = ((left_chr,left_pos), (right_chr, right_pos), left_strand*right_strand, ins_bases)
                if left_chr == right_chr and left_strand == right_strand and left_pos < right_pos:
                    self.is_del = True
                    #print('Del')
                elif left_chr == right_chr and left_strand == right_strand and left_pos > right_pos:
                    self.is_dup = True
                   # print('Dup')
                elif left_chr != right_chr:
                    self.is_transloc = True
                elif left_strand != right_strand:
                    self.is_inv = True
               # print(str(mapped) + '\t' + str(overlap) + '\t' + read.cigarstring + '\t' + chi_cigarstring + '\t' + ins_bases)
               # print(read.seq[:mapped[0]+1] + ' ' + ref_obj.fetch(ref_obj.references[left_chr], left_pos - mapped[0], left_pos + 1))
               # print(read.seq[mapped[1]:] + ' ' + ref_obj.fetch(ref_obj.references[right_chr], right_pos, right_pos + mapped[1]))
                
                        
    def __sub__(self, other):
        if self.position[1] != None:
        # clip1 = (None, chr1, pos) (chr1, pos, None)
            if other.position[0] == None:
                return None
            if self.position[1][0] != other.position[0][0]:
                return None
            return self.position[1][1] - other.position[0][1]
        if self.position[0] != None:
        # clip1 = (chr1, pos, None) (None, chr1, pos)
            if other.position[1] == None:
                return None
            if self.position[0][0] != other.position[1][0]:
                return None
            return  other.position[1][1] - self.position[0][1]
        
class Polys():
    def __init__(self, poly_list, species, genome_version, strains):
        self.species = species
        self.genome_version = genome_version
        if type(strains) is str:
            strains = [strains]
        self.to_geno_strains = set(strains)
        self.geno_strains = set(strains)
        if type(poly_list) is dict:
            self.polys = list(poly_list.values())
        else:
            self.polys = poly_list
        #self.out_vcf = SO.get_vcf_file(strains[0], self.genome_version)
        self.ref_obj = pysam.FastaFile(DI.ret_genome_versions(self.species, self.genome_version).unzipped_ref_file)
        
    def genotypePolys(self):
        for strain in self.to_geno_strains:
            bam_obj = pysam.AlignmentFile(FM.ret_asb(strain, self.species, self.genome_version, 'all'))
            for poly in self.polys:
                poly.genotype_bam(bam_obj, strain, add = True)
                    
    def vcfLoader(self, infile, add_geno = False):
        with open(infile) as f:
            sample_order = []
            for line in f:
                if line[:2] == '##':
                    continue
                else:
                    tokens = line.rstrip().split('\t')
                    if len(tokens) < 2:
                        continue
                    if tokens[0] == '#CHROM':
                        if add_geno:
                            for strain in tokens[9:]:
                                self.geno_strains.add(strain.replace('Sample',''))
                                sample_order.append(strain.replace('Sample', ''))
                            continue
                        else:
                            continue
                    chrom = tokens[0]
                    pos = int(tokens[1]) - 1
                    ref = tokens[3]
                    alt = tokens[4]
                    if alt == '<DEL>':
                        stop = int(line.split('END=')[1].split(';')[0])
                        ref = self.ref_obj.fetch(chrom, pos, stop)
                        alt = self.ref_obj.fetch(chrom, pos, pos+1)

                    snp_type = tokens[7].split('SVTYPE=')[1].split(';')
                    if snp_type == 'DUP':
                        tpoly = poly(chrom, pos, ref, alt, self.ref_obj, is_dup = True)
                    elif snp_type == 'INV':
                        tpoly = poly(chrom, pos, ref, alt, self.ref_obj, is_inv = True)
                    else:
                        tpoly = poly(chrom, pos, ref, alt, self.ref_obj)

                    if add_geno:
                        for i,tgeno in enumerate(tokens[9:]):
                            genos = tgeno.split(':')
                            tpoly.geno[sample_order[i]] = (genos[0],(int(genos[1]),int(genos[2]),int(genos[3])))
                            
                    self.polys.append(tpoly)


    def create_VCFfile(self, filename, geno = True):
        self.polys = sorted(self.polys, key = attrgetter('chrom', 'pos'))
        strain_order = sorted(self.geno_strains)
        with open(filename, 'w') as f:
            print(self.create_VCFheader(geno), file = f)
            prev_poly = None
            for poly in self.polys:
                if poly.empty:
                    continue
                if not self.polys_overlap(poly, prev_poly):
                    print(poly.ret_VCFrecord(strain_order, geno), file = f)
                    prev_poly = poly
                else:
                    print('Overlapping polys (2nd poly not added):', file = sys.stderr)
                    print(prev_poly.ret_VCFrecord(strain_order, geno), file = sys.stderr)
                    print(poly.ret_VCFrecord(strain_order, geno), file = sys.stderr)
                    

    def create_VCFheader(self, geno):
        out_string = '##fileformat=VCFv4.2\n'
        out_string += '##fileDate=' + time.strftime('%Y%m%d') + '\n'
        out_string += '##source=BamTools.py\n'
        out_string += '##reference=' + self.ref_obj.filename.decode('utf-8') + '\n'
        for i in range(0, self.ref_obj.nreferences):
            name = self.ref_obj.references[i]
            length = self.ref_obj.lengths[i]
            out_string += '##contig=<ID=' + name + ',length=' +str(length) + ',species="Metraclima zebra">\n'
        out_string += '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n'
        out_string += '##INFO=<ID=FST,Number=1,Type=Float,Description="Fst">\n'
        out_string += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        out_string += '##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference Reads">\n'
        out_string += '##FORMAT=<ID=AR,Number=1,Type=Float,Description="Alt Reads">\n'
        out_string += '##FORMAT=<ID=BR,Number=1,Type=Float,Description="Ambiguous Reads">\n'
        if geno:
            out_string += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(['Sample' + x for x in sorted(self.geno_strains)])
        else:
            out_string += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'

        return out_string

    def annotate_SNPs(self, outfile, snpeff_database):
        print ('Annotating')
        call(['snpeff', 'eff', snpeff_database, self.out_vcf], stdout = open(outfile, 'w'))
        #call(['python3', '/Users/pmcgrath7/Dropbox/Projects/Sequencing/Scripts/modify_snpEff.py', tempfile, self.genome_version], stdout = open(SO.get_vcf_file(self.mut_strain, self.genome_version, annotation = True), 'w'))
        #call(['rm','-f',tempfile, self.out_vcf])

    def polys_overlap(self, poly1, poly2):
        if poly2 is None:
            return False
        if poly1.chrom != poly2.chrom:
            return False
        if poly1.is_inv:
            if not poly2.is_inv:
                return False
        if abs(poly1.start - poly2.start) > 6:
            return False
        if abs(poly1.stop - poly2.stop) > 6:
            return False
        return True
            
    def plot_cluster(self):
        num_polys = len(self.polys)
        num_samples = len(self.geno_strains)
        geno_data = numpy.empty(shape = (num_polys, num_samples))
        for i, poly in enumerate(self.polys):
            geno_data[i] = poly.ret_geno_vector()
        
        
class poly():
    def __init__(self, chrom, pos, ref, alt, ref_o, is_inv = False, is_dup = False):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.geno = {}
        self.ref_o = ref_o
        self.start = self.pos
        self.Fst = ''
        self.empty = True
        if is_inv:
            self.stop = ref
            self.alt = 'INV'
        else:
            self.stop = self.pos + len(self.ref)
        self.is_inv = is_inv
        self.is_dup = is_dup

    def __cmp__(self, other):
        if self.chrom != other.chrom:
            return self.chrom.__cmp__(other.chrom)
        else:
            return self.pos.__cmp__(other.pos)
        
    def ret_ref_allele(self, flanking):
        if self.is_inv:
            seq = self.ref_o.fetch(self.chrom, max(0,self.start-flanking), self.stop + flanking).upper()
            if len(seq) > 3*flanking:
                return seq[:2*flanking] + seq[-2*flanking:]
            else:
                return seq
        else:
            return self.ref_o.fetch(self.chrom, max(0,self.start - flanking), self.stop + flanking).upper()

    def ret_alt_allele(self, flanking):
        if self.is_inv:
            seq = self.ref_o.fetch(self.chrom, max(0,self.start - flanking), self.start) + reverse_complement(self.ref_o.fetch(self.chrom, self.start, self.stop+1).upper()) + self.ref_o.fetch(self.chrom, self.stop+2, self.stop + flanking)
            if len(seq) > 3*flanking:
                return seq[:2*flanking] + seq[-2*flanking:]
            else:
                return seq
        else:
            return self.ref_o.fetch(self.chrom, max(0,self.start - flanking), self.start).upper() + self.alt  + self.ref_o.fetch(self.chrom, self.stop, self.stop + flanking).upper() 
        
    def genotype_bam(self, bam_file, strain, read_length = 100, add = False, troubleshooting = False):
        unique_reads = {}
        ref_seq = self.ret_ref_allele(read_length)
        alt_seq = self.ret_alt_allele(read_length)
        if troubleshooting != False:
            print(ref_seq, file = troubleshooting)
            print(alt_seq, file = troubleshooting)
        reads1 = bam_file.fetch(self.chrom, max(0,self.start-10), self.start + 10)
        for read in reads1:
            if read.pos - 5 < self.start and read.pos + len(read.seq) > self.start:
                unique_reads[read.qname] = read
        reads2 = bam_file.fetch(self.chrom, max(0,self.stop-10),  self.stop + 10)
        for read in reads2:
            if read.pos -5 < self.stop and read.pos + len(read.seq) > self.stop:
                unique_reads[read.qname] = read
        geno = [0,0,0]
        ref_SSW = StripedSmithWaterman(ref_seq)
        alt_SSW = StripedSmithWaterman(alt_seq)
       
        for read in unique_reads.values():
            if not read.is_secondary and read.mapq > 10:
                ref_score = max(ref_SSW(read.seq).optimal_alignment_score, ref_SSW(reverse_complement(read.seq)).optimal_alignment_score)
                alt_score = max(alt_SSW(read.seq).optimal_alignment_score, alt_SSW(reverse_complement(read.seq)).optimal_alignment_score)
                if ref_score - alt_score > 10:
                    geno[0]+=1
                elif ref_score - alt_score < -10:
                    geno[1]+=1
                else:
                    geno[2]+=1
                if troubleshooting != False:
                    print(read.seq + '\t' + str(ref_SSW(read.seq).optimal_alignment_score) + '\t' + str(alt_SSW(read.seq).optimal_alignment_score), file = troubleshooting)
                    print(reverse_complement(read.seq) + '\t' + str(ref_SSW(reverse_complement(read.seq)).optimal_alignment_score) + '\t' + str(alt_SSW(reverse_complement(read.seq)).optimal_alignment_score), file = troubleshooting)
                    
        if len(unique_reads) == 0:
            geno_o = ('./.',tuple(geno))
        if geno[0] == 0 and geno[1] == 0:
            geno_o = ('./.',tuple(geno))
        elif geno[0] > 4*geno[1]:
            geno_o = ('0/0',tuple(geno))
        elif geno[0] < 1/4*geno[1]:
            geno_o = ('1/1',tuple(geno))
            self.empty = False
        else:
            geno_o = ('0/1',tuple(geno))
            self.empty = False
        if add == True:
            self.geno[strain] = geno_o
        return geno_o
    
    def ret_VCFrecord(self, strain_order, geno):

        if self.is_inv:
            if self.stop - self.start < 100000:
                out_string = self.chrom + '\t' + str(self.pos+1) + '\t.\t' + str(self.ref_o.fetch(self.chrom, self.start-1, self.stop).upper()) + '\t' + str(reverse_complement(self.ref_o.fetch(self.chrom, self.start-1, self.stop).upper()))  + '\t30\tPASS\tNS='  + str(len(self.geno)) + ';END=' + str(self.stop) + ';SVTYPE=INV;SVLEN=' + str(self.stop - self.start + 1) + '\tGT:DP:RF:MF'
            else:
                out_string = self.chrom + '\t' + str(self.pos+1) + '\t.\t' + self.ref_o.fetch(self.chrom, self.start-1, self.start).upper() + '\t<INV>\t30\tPASS\tNS=' + str(len(self.geno)) + ';END=' + str(self.stop) + ';SVTYPE=INV;SVLEN=' + str(self.stop - self.start + 1) + '\tGT:DP:RF:MF'
                
        elif self.is_dup:
            if self.stop - self.start < 100000:
                out_string = self.chrom + '\t' + str(self.pos+1) + '\t.\t' + str(self.ref) + '\t' + self.alt + '\t30\tPASS\tNS=' + str(len(self.geno)) + ';END=' + str(self.stop) + ';SVTYPE=DUP;SVLEN=' + str(self.stop - self.start + 1) + '\tGT:DP:RF:MF'
            else:
                out_string = self.chrom + '\t' + str(self.pos+1) + '\t.\t' + self.ref_o.fetch(self.chrom, self.start-1, self.start).upper() + '\t<DUP>\t30\tPASS\tNS=' + str(len(self.geno)) + ';END=' + str(self.stop) + ';SVTYPE=DUP;SVLEN=' + str(self.stop - self.start + 1) + '\tGT:DP:RF:MF'
        else:
            if self.stop - self.start < 100000:
                out_string = self.chrom + '\t' + str(self.pos+1) + '\t.\t' + str(self.ref) + '\t' + self.alt + '\t30\tPASS\tNS=' + str(len(self.geno)) + ';END=' + str(self.stop) + ';SVTYPE=DEL;SVLEN=' + str(self.stop - self.start + 1) + '\tGT:DP:RF:MF'
            else:
                out_string = self.chrom + '\t' + str(self.pos+1) + '\t.\t' + self.ref_o.fetch(self.chrom, self.start-1, self.start).upper() + '\t<DEL>\t30\tPASS\tNS=' + str(len(self.geno)) + ';END=' + str(self.stop) + ';SVTYPE=DEL;SVLEN=' + str(self.stop - self.start + 1) + '\tGT:DP:RF:MF'
        if geno:
            for strain in strain_order:
                out_string += '\t' + self.geno[strain][0] + ':' + str(self.geno[strain][1][0]) + ':' + str(self.geno[strain][1][1]) + ':' + str(self.geno[strain][1][2])
        return out_string

    def ret_geno_vector(self):
        out_genos = [0]*len(self.geno)
        for i, geno in enumerate(self.geno):
            if geno[0] == '0/0':
                out_genos[i] = 0
            if geno[0] == '0/1':
                out_genos[i] = 1
            if geno[0] == '1/1':
                out_genos[i] = 2
        return out_genos
                
    
