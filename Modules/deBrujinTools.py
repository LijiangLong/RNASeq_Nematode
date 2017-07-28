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
                reference_bases = str(self.ref_o.get_seq(contig, position-1, position -1 + len(out_aln[0][1][i1:i2].replace('-',''))))
                vcf_line = '\t'.join([contig, str(position), '.', reference_bases, inserted_bases, '30', 'PASS', 'NS=0', 'GT:GQ:DP'])
                if len(inserted_bases) > 200 and out_aln[0][2] > 850:
                    return poly(contig, position, reference_bases, inserted_bases, self.ref_o)

class db_node:
    def __init__(self, sequence):
        self.sequence = sequence
        self.num = 1
        self.f_edges = defaultdict(int)
        self.r_edges = defaultdict(int)
                
