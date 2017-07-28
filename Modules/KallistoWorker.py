from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM
from multiprocessing import cpu_count
import Modules.SeqOrganizer as SO
import sys, os, pysam
import pdb
        

class KallistoMaker():
    def __init__(self, genome_version, species, samples):
        self.genome_version = genome_version
        self.species = species
        self.samples = samples
        self.ref_files = DI.ret_genome_versions(self.species, self.genome_version)
        self.idx = self.ref_files.ref_file.split('genomic.fa.gz')[0]+'mRNA_transcripts.idx'
        if SI.valid_samples(self.samples) is False:
            print('Problem with one of the samples. Exiting...', file = sys.stderr)
            sys.exit()
        if int(genome_version.split('WS')[1]) < 253:
            print('current program only work with WS253 or later', file = sys.stderr)
            sys.exit()
        

    def count(self):
#        pdb.set_trace()
#        local_bam = SO.l_bam_directory
        for sample in self.samples:
            idx_file = self.idx
            for fq in SI.ret_fastq_objs(sample):
                mean_frag = str(fq.mean_frag)
                sd_frag = str(fq.sd_frag)
                
                output_dir = '/Volumes/Lijiang_data/packages/Elegans_RNA_seq_LL/Gene_read_count/kallisto/c_elegans/WS255/'+sample
                command = 'mkdir ' + output_dir
                os.system(command)
                if fq.paired_flag:
                    if fq.twofq:
                        print(' '.join(['tophat', '-p', '4', '-G', annotation_file, '-o', local_bam,ref_file, fq.fqfile, fq.fqfile2]), file = sys.stderr)
                        call(['tophat','-p', str(cpu_count()), '-G', annotation_file, '-o', local_bam, ref_file, fq.fqfile, fq.fqfile2], stdout = open(tfile1, 'w'))
                    else:
                        print(' '.join(['tophat', '-p', '4', '-G', annotation_file, '-o', local_bam, ref_file, fq.fqfile, fq.fqfile2]), file = sys.stderr)
                        print(' '.join(['bwa', 'mem', '-t', str(cpu_count()), '-p', '-R', fq.RG_info, '-M', ref_file, fq.fqfile]))
                        call(['bwa', 'mem', '-t', str(cpu_count()), '-p', '-R', fq.RG_info, '-M', ref_file, fq.fqfile], stdout = open(tfile1, 'w'))
                else:
                    command = ' '.join(['kallisto', 'quant','-i',idx_file, '-o',output_dir,'-b 100','--single','-l',mean_frag,'-s',sd_frag, fq.fqfile])
                    print(command)
                    os.system(command)

