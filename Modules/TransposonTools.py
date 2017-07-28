import pysam, os, pyfaidx, sys
from collections import defaultdict
from subprocess import call
from Modules.deBrujin import DBGraph
import time

COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}

def reverse_complement(sequence):
    return ''.join(COMPLEMENT[b.upper()] for b in sequence[::-1])

def create_VCFheader(ref_obj):
    out_string = '##fileformat=VCFv4.2\n'
    out_string += '##fileDate=' + time.strftime('%Y%m%d') + '\n'
    out_string += '##source=BamTools.py\n'
    out_string += '##reference=' + ref_obj.filename.decode('utf8') + '\n'
    for i in range(0, ref_obj.nreferences):
        name = ref_obj.references[i]
        length = ref_obj.lengths[i]
        out_string += '##contig=<ID=' + name + ',length=' +str(length) + ',species="Metraclima zebra">\n'
    out_string += '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n'
    out_string += '##INFO=<ID=FST,Number=1,Type=Float,Description="Fst">\n'
    out_string += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    out_string += '##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference Reads">\n'
    out_string += '##FORMAT=<ID=AR,Number=1,Type=Float,Description="Alt Reads">\n'
    out_string += '##FORMAT=<ID=BR,Number=1,Type=Float,Description="Ambiguous Reads">\n'
    out_string += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\n'

    return out_string

def coverage_calculator(bam_file, contig, start, stop):
    f_reads = 0
    r_reads = 0
    reads = bam_file.fetch(contig, start, stop)
    for read in reads:
        if 'S' not in read.cigarstring and 'H' not in read.cigarstring:
            if read.is_reverse:
                r_reads += 1
            else:
                f_reads += 1
    return (f_reads, r_reads)
def coverage_calculator_cl(bam_file, contig, start, stop):
    n_reads = 0
    reads = bam_file.fetch(contig, start, stop)
    for read in reads:
        n_reads += 1
    return n_reads

def create_candidate_sand_transposons():
    global rock_clipped, sand_clipped, rock_cf, sand_cf
    with open(candidate_sand_transposons, 'w') as f:
        for i in range(0,rock_cf.nreferences):
            contig = rock_cf.references[i]
            length = rock_cf.lengths[i]
            current_hit = False
            for j in range(0, length-bin_size, overlap):
                r_reads = coverage_calculator(rock_cf, contig, j, j+bin_size)
                s_reads = coverage_calculator(sand_cf, contig, j, j+bin_size)
                hit_bin = False
                if s_reads[0] > 30 and r_reads[0] < 5 and s_reads[0] > 5 * s_reads[1]:
                    if j + 4*bin_size < length:
                        for k in range(1,4):
                            r_reads2 = coverage_calculator(rock_cf, contig, j + k*bin_size, j + (k+1)*bin_size)
                            s_reads2 = coverage_calculator(sand_cf, contig, j + k*bin_size, j + (k+1)*bin_size)
                            if s_reads2[1] > 30 and r_reads2[1] < 5 and s_reads2[1] > 5* s_reads2[0]:
                                hit_bin = True
                                #print(contig + '\t' + str(j) + '\t' + str(s_reads) + '\t' + str(r_reads))
                                #print(contig + '\t' + str(j+k*bin_size) + '\t' + str(s_reads2) + '\t' + str(r_reads2))
                if current_hit == True and hit_bin == False:
                    mid_pos = int(j + bin_size/2)
                    r_reads_l = coverage_calculator(rock_cf, contig, mid_pos - 250, mid_pos)
                    s_reads_l = coverage_calculator(sand_cf, contig, mid_pos - 250, mid_pos)
                    r_reads_r = coverage_calculator(rock_cf, contig, mid_pos, mid_pos + 250)
                    s_reads_r = coverage_calculator(sand_cf, contig, mid_pos, mid_pos + 250)
                    r_clipped = coverage_calculator_cl(rock_clipped, contig, mid_pos - 50, mid_pos+50)
                    s_clipped = coverage_calculator_cl(sand_clipped, contig, mid_pos - 50, mid_pos+50)
                
                    print(contig + '\t' + str(mid_pos) + '\tR_CFl_CFr_Cl: ' + str(r_reads_l) + '_' + str(r_reads_r) + '_' + str(r_clipped) + '\tS_CFl_CFr_Cl: ' + str(s_reads_l) + '_' + str(s_reads_r) + '_' + str(s_clipped), file = f)
            
                    current_hit = False
                if hit_bin == True:
                    current_hit = True

def read_candidate_file(infile):
    out_data = []
    with open(infile) as f:
        for line in f:
            contig = line.split()[0]
            pos = int(line.split()[1])
            out_data.append((contig,pos))
    return out_data

def collect_reads(transposon, q_bam_file, or_bam_file):
    all_reads = []
    
    reads = q_bam_file.fetch(transposon[0], transposon[1] - 250, transposon[1])
    for read in reads:
        if not read.is_reverse:
            all_reads.append(read.seq)
            read_name = read.query_name
            mate_contig = read.next_reference_name
            mate_pos = read.mpos
            o_reads = or_bam_file.fetch(mate_contig, max(0,mate_pos - 1), mate_pos + 1)
            count = 0
            for o_read in o_reads:
                if o_read.query_name == read_name:
                    count += 1
                    if o_read.is_reverse:
                        all_reads.append(o_read.seq)
                    else:
                        all_reads.append(reverse_complement(o_read.seq))
                
    reads = q_bam_file.fetch(transposon[0], transposon[1], transposon[1]+250)
    for read in reads:
        if read.is_reverse:
            all_reads.append(read.seq)
            read_name = read.query_name
            mate_contig = read.next_reference_name
            mate_pos = read.mpos
            o_reads = or_bam_file.fetch(mate_contig, max(0,mate_pos - 1), mate_pos + 1)
            count = 0
            for o_read in o_reads:
                if o_read.query_name == read_name:
                    count += 1
                    if not o_read.is_reverse:
                        all_reads.append(o_read.seq)
                    else:
                        all_reads.append(reverse_complement(o_read.seq))

    return(all_reads)
