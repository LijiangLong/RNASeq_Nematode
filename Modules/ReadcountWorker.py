from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM
from subprocess import call, Popen, PIPE
from multiprocessing import cpu_count
import Modules.SeqOrganizer as SO
import sys, os, pysam
import pdb
        


class ReadcountMaker():
    def __init__(self, genome_version, species, strand,order,samples):
        self.genome_version = genome_version
        self.species = species
        self.strand = strand
        self.order = order
        self.samples = samples
        self.ref_files = DI.ret_genome_versions(self.species, self.genome_version)
        if SI.valid_samples(self.samples) is False:
            print('Problem with one of the samples. Exiting...', file = sys.stderr)
            sys.exit()
        if int(genome_version.split('WS')[1]) < 253:
            print('current program only work with WS253 or later', file = sys.stderr)
            sys.exit()
        
    def count(self):
        for sample in self.samples:
            annotation_file = self.ref_files.annotation_file
            bam_file = FM.ret_asb(sample, self.species, self.genome_version, 'all')
            out_file = FM.ret_srcf(sample, self.species, self.genome_version)
            command = 'htseq-count -f bam -s '+self.strand+' '+bam_file+' '+annotation_file
            command += ' > ' + out_file
            print(command)
            os.system(command)
    

        