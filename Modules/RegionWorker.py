from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM
from subprocess import call, Popen, PIPE
from multiprocessing import cpu_count
import Modules.SeqOrganizer as SO
import Modules.GtfWorker as GTF
import sys, os, pysam
import pdb
        

            
class RegionMaker():
    def __init__(self, genome_version, species, chrom,coordinates):
        self.genome_version = genome_version
        self.species = species
        self.chrom = chrom
        self.start = int(coordinates.split('-')[0])
        self.end = int(coordinates.split('-')[1])
        self.ref_files = DI.ret_genome_versions(self.species, self.genome_version)
        if int(genome_version.split('WS')[1]) < 253:
            print('current program only work with WS253 or later', file = sys.stderr)
            sys.exit()
        
    def gene_in_region(self):
#            pdb.set_trace()
            
            gtf_file = self.ref_files.annotation_file
            gtf_object = GTF.Gtf_file(gtf_file)
            gtf_lines = gtf_object.lines()
            gene_list = []
            for line in gtf_lines:
                gene_id = line['gene_id']
                chrom = line['seqname']
                start = int(line['start'])
                end = int(line['end'])
                feature = line['feature']
                if feature == 'gene' and chrom == self.chrom and start > self.start and end < self.end:
                    gene_list.append(gene_id)
            out_file = FM.ret_genes_in_region(self.species, self.genome_version,self.chrom,self.start, self.end)
            with open(out_file,'w') as output:
                for gene in gene_list:
                    output.write(gene+'\n')
            
            
            
    

        