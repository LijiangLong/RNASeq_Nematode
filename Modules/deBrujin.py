from collections import defaultdict
from graphviz import Digraph
from subprocess import call
import os, pysam, operator, sys
from skbio.alignment import StripedSmithWaterman
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
import pandas
import seaborn as sns

COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}

def reverse_complement(sequence):
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])


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
            
    def indentify_polymorphism(self):
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

            ref_seq = str(self.ref_o.get_seq(contig, pos1-250, pos2+250)).replace('N','')
            #ref_SSW = StripedSmithWaterman(ref_seq)
            self.out_seqs = out_seqs
            self.ref_seq = ref_seq
#            sys.exit()
            if len(out_seqs) == 1:
                
                c = out_seqs[0]
                out_aln = pairwise2.align.globalms(c, ref_seq, 2, -1, -2, -.1)
                #a1 = ref_SSW(a)
                #if abs(a1.query_begin - 250) < abs(a1.query_end - 250):
                #    mid = a1.target_begin - 50
                #else:
                #    mid = a1.target_end_optimal + 50
                #a, b = out_seqs[0][:mid], out_seqs[0][mid:]
                #print('Split')
            else:
                a = out_seqs[0]
                b = out_seqs[1]
                a1 = pairwise2.align.globalms(a+'NNNNN' + b, ref_seq, 2, -1, -2, -.1)
                a2 = pairwise2.align.globalms(b+'NNNNN' + a, ref_seq, 2, -1, -2, -.1)
                if a1[0][2] > a2[0][2]:
                    out_aln = a1
                else:
                    out_aln = a2
            #print(out_aln[0][0])
            #print(out_aln[0][1])
            #print(out_aln[0][2])
            #print(out_aln[0][3])
            #print(out_aln[0][4])

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
                reference_bases = str(self.ref_o.get_seq(contig, position-1, position -1 + len(out_aln[0][1][i1:i2].replace('-',''))))
                vcf_line = '\t'.join([contig, str(position), '.', reference_bases, inserted_bases, '30', 'PASS', 'NS=0', 'GT:GQ:DP'])
                if len(inserted_bases) > 200 and out_aln[0][2] > 850:
                    return(vcf_line, inserted_bases)

class db_node:
    def __init__(self, sequence):
        self.sequence = sequence
        self.num = 1
        self.f_edges = defaultdict(int)
        self.r_edges = defaultdict(int)

class blast_data:
    def __init__(self, blast_file):
        self.data = {}
        with open(blast_file) as f:
            for line in f:
                q_name = line.split('\t')[0]
                h_name = line.split('\t')[1]
                bit_score = float(line.split('\t')[11].rstrip())
                try:
                    self.data[q_name][h_name]+=bit_score
                except KeyError:
                    self.data[q_name] = defaultdict(int)
                    self.data[q_name][h_name]+=bit_score
        good_data = set()
        for key in self.data:
            if len(self.data[key]) - list(self.data[key].values()).count(0) > 0:
                good_data.add(key)
        queries = sorted(good_data)
        with open(blast_file + '.csv', 'w') as f:
            print('Name\t' + '\t'.join(queries), file = f)
            for key in queries:
                print(key, end = '', file = f)
                for key2 in queries:
                    #print('\t' + str(self.data[key][key2]), end = '', file = f)
                    if self.data[key][key2] > 20:
                        print('\t1', end = '', file = f)
                    else:
                        print('\t0', end = '', file = f)
                                                
                print(file = f)
        self.dt = pandas.read_csv(blast_file + '.csv', sep = '\t', header = 0, index_col = 0)
        self.clusterplot = sns.clustermap(self.dt)
        sns.plt.show()
