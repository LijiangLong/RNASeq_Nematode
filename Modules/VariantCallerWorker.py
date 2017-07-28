from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM
import Modules.BamTools as BT
from Modules.deBrujinTools import DBGraph
from Bio import pairwise2
import sys, os, pysam, pyfaidx, operator, vcf
from collections import defaultdict
from subprocess import call

class SNVCallerGenotyper():
    def __init__(self, species, Genome_version):
        self.species = species
        self.Genome_version = Genome_version
        # Create pysam objects for reference and bam files
        self.ref_file = DI.ret_genome_versions(self.species, self.Genome_version).unzipped_ref_file
        self.ref = pysam.FastaFile(DI.ret_genome_versions(self.species, self.Genome_version).ref_file)
  
    def identifyVariants(self, discovery_samples):
        if len(discovery_samples) == 1 and SI.valid_group(discovery_samples[0]) is not False:
            self.discovery_samples = set(SI.valid_group(discovery_samples[0]))
            self.discovery_group = discovery_samples[0]
        else:
            self.discovery_samples = set(discovery_samples)

        for sample in sorted(self.discovery_samples):
            sample_vcf = FM.ret_vcf_file('Polys:SNVs', self.species, self.Genome_version, sample, 'None', sample)
            print('Identifying SNVs from: ' + sample, file = sys.stderr)
            mut_bam = FM.ret_asb(sample, self.species, self.Genome_version, 'all')
            call(['freebayes', '-f', self.ref_file, mut_bam], stdout = open(sample_vcf, 'w'))
            call(['bgzip', sample_vcf])
            call(['tabix', '-p', 'vcf', sample_vcf + '.gz'])

    def analyzeVariants(self, included_samples, excluded_samples, genotyped_samples):
        if len(included_samples) == 1 and SI.valid_group(included_samples[0]) is not False:
            self.included_samples = set(SI.valid_group(included_samples[0]))
            self.included_group = included_samples[0]
        else:
            self.included_samples = set(included_samples)

        if len(excluded_samples) == 1 and SI.valid_group(excluded_samples[0]) is not False:
            self.excluded_samples = set(SI.valid_group(excluded_samples[0]))
            self.excluded_group = excluded_samples[0]
        else:
            self.excluded_samples = set(excluded_samples)

        if len(genotyped_samples) == 1 and SI.valid_group(genotyped_samples[0]) is not False:
            self.genotyped_samples = set(SI.valid_group(genotyped_samples[0]))
            self.genotyped_group = genotyped_samples[0]
        else:
            self.genotyped_samples = set(genotyped_samples)

        self.genotyped_samples.update(self.included_samples)
        self.genotyped_samples.update(self.excluded_samples)
            
        included_vcfs = []
        temp_vcf = 'TempIncludedSamples.vcf'
        for sample in sorted(self.included_samples):
            included_vcfs.append(FM.ret_vcf_file('Polys:SNVs', self.species, self.Genome_version, sample, 'None', sample) + '.gz')

        call(['vcf-merge', '-d'] + included_vcfs, stdout = open(temp_vcf, 'w'))

        genotyped_bams = []
        temp2_vcf = 'TempIncludedSamplesGenotyped.vcf'
        for sample in self.genotyped_samples:
            genotyped_bams.append(FM.ret_asb(sample, self.species, self.Genome_version, 'all'))

        call(['freebayes', '-@', temp_vcf, '-l', '-f', self.ref_file] + genotyped_bams, stdout = open(temp2_vcf, 'w'))
        geno_file = vcf.VCFReader(filename = temp2_vcf)
        self.genotyped_vcf = 'TempIncludedSamplesGenotyped.vcf'
        out_vcf = vcf.Writer(open('TempIncludedSamplesGenotypedFiltered.vcf', 'w'), geno_file)
        for r in geno_file:
            included = False
            excluded = False
            for s in r.samples:
                if s.sample in included_samples:
                    if s.gt_nums == '1/1':
                        included = True
                if s.sample in excluded_samples:
                    if s.gt_nums in ['1/1','0/1']:
                        excluded = True

            if included and not excluded:
                out_vcf.write_record(r)

        
class LargeVarCallerGenotyper():
    def __init__(self, species, Genome_version):
        # Save input parameters
        self.species = species
        self.Genome_version = Genome_version
        self.ref = pysam.FastaFile(DI.ret_genome_versions(self.species, self.Genome_version).ref_file)

        # Set globals
        self.max_len = 2000000
        self.bin_size = 100 # Reads are added up over this level
        self.overlap = 100
        self.min_reads = 5

        
    def identifyVariants(self, discovery_samples):
        if len(discovery_samples) == 1 and SI.valid_group(discovery_samples[0]) is not False:
            self.discovery_samples = set(SI.valid_group(discovery_samples[0]))
            self.discovery_group = discovery_samples[0]
        else:
            self.discovery_samples = set(discovery_samples)
        print(self.discovery_samples)

        for sample in sorted(self.discovery_samples):
            sample_vcf = FM.ret_vcf_file('Polys:LargeVar', self.species, self.Genome_version, sample, 'None', sample)
            print('Identifying large variants from: ' + sample, file = sys.stderr)
            mut_all = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.Genome_version, 'all'))
            mut_clipped = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.Genome_version, 'clipped'))
            mut_clipped_obj = BT.clipped_obj(mut_clipped, self.ref)
            self.tempPolys = {}
            # Calculate avg cov        
            self.mut_cov = BT.coverage(mut_all)
            for loc, reads in mut_clipped_obj.chimera_reads.items():
                chrom = mut_clipped.references[loc[0][0]]
                start = loc[0][1]
                stop = loc[1][1]

                if len(reads) > 0 and stop - start < self.max_len:
                    if reads[0].is_del:
                        mut_cov = BT.avg_coverage(mut_all, chrom, start, stop)
                        refbase = self.ref.fetch(chrom,start,stop).upper()
                        altbase = self.ref.fetch(chrom,start,start+1).upper() + loc[3]
                        poly = BT.poly(chrom, start, refbase, altbase, self.ref)
                        mut_geno = poly.genotype_bam(mut_all, sample)
                        if mut_geno[0] != '0/0':
                            try:
                                mut_cov_ratio = (mut_cov/self.mut_cov)
                            except ZeroDivisionError:
                                continue
                            if mut_geno[0] == '1/1' and  mut_cov_ratio < .8:
                                self.tempPolys[(chrom, start, refbase, altbase)] = poly
                            elif mut_geno[0] == '0/1' and mut_cov_ratio < .8:
                                self.tempPolys[(chrom, start, refbase, altbase)] = poly

                    if reads[0].is_dup:
                        start = max(0,loc[1][1])
                        stop = loc[0][1]
                        mut_cov = BT.avg_coverage(mut_all, chrom, start, stop)
                        refbase = self.ref.fetch(chrom,stop,stop+1).upper()
                        altbase = refbase + self.ref.fetch(chrom,start,stop+1).upper() + loc[3]
                        poly = BT.poly(chrom, stop, refbase, altbase, self.ref, is_dup = True)
                        mut_geno = poly.genotype_bam(mut_all, sample)

                        if mut_geno[0] != '0/0':
                            try:
                                mut_cov_ratio = (mut_cov/mut_cov)
                            except ZeroDivisionError:
                                continue
                            self.tempPolys[(chrom, start, refbase, altbase)] = poly

                    if reads[0].is_inv:
                        pos1 = loc[0][1]
                        pos2 = loc[1][1]
                        start = min(pos1,pos2)+1
                        stop = max(pos1,pos2)
                        if start < 0 or stop < 0:
                            continue

                        alt = ''
                        poly = BT.poly(chrom, start, stop, alt, self.ref, is_inv=True)
                        mut_geno = poly.genotype_bam(mut_all, sample)
                        if mut_geno[0] != '0/0':
                            self.tempPolys[(chrom, start, stop, alt)] = poly

            # Identify large insertions
            discordant = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.Genome_version, 'discordant'))
        
            for i in range(0,discordant.nreferences):
                contig = discordant.references[i]
                length = discordant.lengths[i]
                current_hit = False
                hit_bin = False
                for j in range(0, length-self.bin_size, self.overlap):
                    reads = self.coverage_calculator(discordant, contig, j, j+self.bin_size)
                    #Does bin have necessary number of forward and reverse discordant reads?
                    if reads[0] > self.min_reads and reads[0] > 5 * reads[1]:
                        if j + 4*self.bin_size < length:
                            for k in range(1,4):
                                reads2 = self.coverage_calculator(discordant, contig, j + k*self.bin_size, j + (k+1)*self.bin_size)
                                if reads2[1] > self.min_reads and reads2[1] > 5*reads2[0]:
                                    hit_bin = True

                    # Have we reached the end of the hit?
                    if current_hit == True and hit_bin == False:
                        mid_pos = int(j + self.bin_size/2)
                        #Try to identify sequence of 
                        g = DBGraph(50, contig + ':' + str(mid_pos) + '_' + str(mid_pos), self.ref)
                        all_reads = self.collect_reads(contig + ':' + str(mid_pos) + '_' + str(mid_pos), discordant)
                        for read in all_reads:
                            g.add_read(read)
                        a = g.identify_polymorphism()
                        if a is not None:
                            self.tempPolys.append(a)
                        current_hit = False
                    
                    if hit_bin == True:
                        current_hit = True
                    hit_bin = False
                            
            Polys = BT.Polys(self.tempPolys, self.species, self.Genome_version, sample)
            Polys.genotypePolys()
            Polys.create_VCFfile(sample_vcf)
                    
            mut_all.close()
            mut_clipped.close()
            discordant.close()
            
    def analyzeVariants(self, included_samples, excluded_samples, genotyped_samples):
        if len(included_samples) == 1 and SI.valid_group(included_samples[0]) is not False:
            self.included_samples = set(SI.valid_group(included_samples[0]))
            self.included_group = included_samples[0]
        else:
            self.included_samples = set(included_samples)

        if len(excluded_samples) == 1 and SI.valid_group(excluded_samples[0]) is not False:
            self.excluded_samples = set(SI.valid_group(excluded_samples[0]))
            self.excluded_group = excluded_samples[0]
        else:
            self.excluded_samples = set(excluded_samples)

        if len(genotyped_samples) == 1 and SI.valid_group(genotyped_samples[0]) is not False:
            self.genotyped_samples = set(SI.valid_group(genotyped_samples[0]))
            self.genotyped_group = genotyped_samples[0]
        else:
            self.genotyped_samples = set(genotyped_samples)

        self.genotyped_samples.update(self.included_samples)
        self.genotyped_samples.update(self.excluded_samples)
            
        included_vcfs = []
        temp_vcf = 'TempIncludedLargeSamples.vcf'
        for sample in sorted(self.included_samples):
            included_vcfs.append(FM.ret_vcf_file('Polys:LargeVars', self.species, self.Genome_version, sample, 'None', sample) + '.gz')

        call(['vcf-merge', '-d'] + included_vcfs, stdout = open(temp_vcf, 'w'))

        largePolys = BT.Polys([], self.species, self.Genome_version, self.genotyped_samples)
        largePolys.vcfLoader(temp_vcf)
        largePolys.genotypePolys()
        temp2_vcf = 'TempIncludedSamplesLargeGenotyped.vcf'
        largePolys.create_VCFfile(temp2_vcf)
        
        geno_file = vcf.VCFReader(filename = temp2_vcf)
        out_vcf = vcf.Writer(open('TempIncludedSamplesLargeGenotypedFiltered.vcf', 'w'), geno_file)
        for r in geno_file:
            included = False
            excluded = False
            for s in r.samples:
                if s.sample in included_samples:
                    if s.gt_nums == '1/1':
                        included = True
                if s.sample in excluded_samples:
                    if s.gt_nums in ['1/1','0/1']:
                        excluded = True

            if included and not excluded:
                out_vcf.write_record(r)
            
    def coverage_calculator(self, bam_file, contig, start, stop):
        f_reads = 0
        r_reads = 0
        reads = bam_file.fetch(contig, start, stop)
        for read in reads:
            if read.is_unmapped:
                continue
            if 'S' not in read.cigarstring and 'H' not in read.cigarstring and read.mapq > 40:
                if read.is_reverse:
                    r_reads += 1
                else:
                    f_reads += 1
        return (f_reads, r_reads)

    def collect_reads(self, transposon, bam_file):
        all_reads = []
        mate_pairs = set()
        contig = transposon.split(':')[0]
        start = int(transposon.split(':')[1].split('_')[0])
        stop = int(transposon.split(':')[1].split('_')[1])
        
        reads_l = bam_file.fetch(contig, start - 250, start)
        for read in reads_l:
            if not read.is_reverse:
                all_reads.append(read.seq)
                mate_pairs.add((read.query_name, read.next_reference_name, read.mpos))

        reads_r = bam_file.fetch(contig, stop, stop+250)
        for read in reads_r:
            if read.is_reverse:
                all_reads.append(read.seq)
                mate_pairs.add((read.query_name, read.next_reference_name, read.mpos))

        for mp in mate_pairs:
            mp_reads = bam_file.fetch(mp[1], max(0,mp[2] - 1), mp[2] + 1)
            for mp_read in mp_reads:
                if mp_read.query_name == mp[0]:
                    if mp_read.is_reverse:
                        all_reads.append(mp_read.seq)
                    else:
                        all_reads.append(BT.reverse_complement(mp_read.seq))
        return(all_reads)

class DBGraph:
    def __init__(self, k_size, name, ref_o):
        self.k = k_size
        self.ref_o = ref_o
        self.nodes = {}
        self.name = name
                
    def add_read(self, read):
        for i in range(0,len(read) - self.k):
            seq1 = read[i:self.k+i]
            seq2 = read[i+1:self.k+1+i]
            if seq1 not in self.nodes:
                self.nodes[seq1] = db_node(seq1)
            else:
                self.nodes[seq1].num += 1
            if seq2 not in self.nodes:
                self.nodes[seq2] = db_node(seq2)
            self.nodes[seq1].f_edges[seq2] += 1
            self.nodes[seq2].r_edges[seq1] += 1

    def clean_graph(self):
        del_keys = []
        for node in self.nodes.values():
            if sum(node.f_edges.values()) <= 2 and sum(node.r_edges.values()) <= 2:
                for key, value in node.f_edges.items():
                    self.nodes[key].r_edges.pop(node.sequence)
                for key, value in node.r_edges.items():
                    self.nodes[key].f_edges.pop(node.sequence)
                del_keys.append(node.sequence)
        for key in del_keys:
            self.nodes.pop(key)
            
    def follow_node(self, seq):
        out_seq = seq
        out_count = 0
        count = 0
        while(True and count < 2000):
            try:
                new_seq = max(self.nodes[seq].f_edges.items(), key=operator.itemgetter(1))[0]
                out_seq += new_seq[-1]
                seq = new_seq
                count += 1
            except ValueError:
                break
        if count != 2000:
            return out_seq
        else:
            return ''
            
    def identify_polymorphism(self):
        self.clean_graph()
        long_seqs = []
        for seq, node in self.nodes.items():
            if sum(node.r_edges.values()) == 0:
                created_seq = self.follow_node(node.sequence)
                if len(created_seq) > 300:
                    long_seqs.append(created_seq)
        good_seqs = [True]*len(long_seqs)
        for i in range(0,len(long_seqs)):
            for j in range(0,len(long_seqs)):
                if i != j:
                    if long_seqs[i][-50:] == long_seqs[j][-50:]:
                        if len(long_seqs[i]) > len(long_seqs[j]):
                            good_seqs[j] = False
                        else:
                            good_seqs[i] = False
        out_seqs = []
        for i in range(0,len(long_seqs)):
            if good_seqs[i]:
                out_seqs.append(long_seqs[i])
        #print(out_seqs)
            
        if len(out_seqs) in [1,2]:
            contig = self.name.split(':')[0]
            pos1 = int(self.name.split(':')[1].split('_')[0])
            pos2 = int(self.name.split(':')[1].split('_')[1])            

            ref_seq = str(self.ref_o.fetch(contig, pos1-250, pos2+250)).replace('N','')
            self.out_seqs = out_seqs
            self.ref_seq = ref_seq
            if len(out_seqs) == 1:                
                c = out_seqs[0]
                out_aln = pairwise2.align.globalms(c, ref_seq, 2, -1, -2, -.1)
            else:
                a = out_seqs[0]
                b = out_seqs[1]
                a1 = pairwise2.align.globalms(a+'NNNNN' + b, ref_seq, 2, -1, -2, -.1)
                a2 = pairwise2.align.globalms(b+'NNNNN' + a, ref_seq, 2, -1, -2, -.1)
                if a1[0][2] > a2[0][2]:
                    out_aln = a1
                else:
                    out_aln = a2

            out_str = ''
            out_strs = []
            for c in out_aln[0][1]:
                if c is not '-':
                    out_str += c
                else:
                    if len(out_str) > 50:
                        out_strs.append(out_str)
                    out_str = ''
            if len(out_str) > 50:
                out_strs.append(out_str)
                
            if len(out_strs) == 2:
                i1 = out_aln[0][1].find(out_strs[0]) + len(out_strs[0])
                i2 = out_aln[0][1].find(out_strs[1])
                count = 0
                for i,c in enumerate(out_aln[0][1]):
                    if c is not '-':
                        count += 1
                    if i == i1:
                        break
                inserted_bases = out_aln[0][0][i1-1:i2].replace('-','')
                deleted_bases = out_aln[0][1][i1:i2].replace('-','')
                position = pos1-250+count
                reference_bases = str(self.ref_o.fetch(contig, position-1, position + len(out_aln[0][1][i1:i2].replace('-',''))))
                vcf_line = '\t'.join([contig, str(position), '.', reference_bases, inserted_bases, '30', 'PASS', 'NS=0', 'GT:GQ:DP'])
                print('InsertionCandidate: ' + inserted_bases + '\t' + str(out_aln[0][2]))
                if len(inserted_bases) > 200 and out_aln[0][2] > 850:
                    return BT.poly(contig, position, reference_bases, inserted_bases, self.ref_o)

class db_node:
    def __init__(self, sequence):
        self.sequence = sequence
        self.num = 1
        self.f_edges = defaultdict(int)
        self.r_edges = defaultdict(int)
                
