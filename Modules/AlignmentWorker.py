from Modules.SeqOrganizer import databaseInfo as DI
from Modules.SeqOrganizer import sampleInfo as SI
from Modules.SeqOrganizer import file_maker as FM
from subprocess import call, Popen, PIPE
from multiprocessing import cpu_count
import Modules.SeqOrganizer as SO
import sys, os, pysam
import pdb
        

class AlignRNAMaker():
    def __init__(self, genome_version, species, samples, split):
        self.genome_version = genome_version
        self.species = species
        self.samples = samples
        self.split = split
#         if len(samples) == 1:
#             if SI.valid_group(samples[0]):
#                 self.samples = SI.valid_group(samples[0])
#                 self.group = samples[0]
#             else:
#                 self.samples = samples
#                 self.group = False
        self.ref_files = DI.ret_genome_versions(self.species, self.genome_version)
        if SI.valid_samples(self.samples) is False:
            print('Problem with one of the samples. Exiting...', file = sys.stderr)
            sys.exit()
        if int(genome_version.split('WS')[1]) < 253:
            print('current program only work with WS253 or later', file = sys.stderr)
            sys.exit()
        

    def align(self):
#        pdb.set_trace()
#        local_bam = SO.l_bam_directory
        if self.split:
            return
        for sample in self.samples:
            bamfiles = []
            ref_file = self.ref_files.ref_file.split('.gz')[0]
            annotation_file = self.ref_files.annotation_file
            hisat2_index_file = ref_file.split('genomic.fa')[0]+'hisat2_index'
            splice_sites_file = ref_file.split('genomic.fa')[0]+'.ss'
            exon_file = ref_file.split('genomic.fa')[0]+'.exon'
            for fq in SI.ret_fastq_objs(sample):
                print(fq.fqfile)
                tfile1 = FM.ret_temp_file()
                if fq.paired_flag:
                    if fq.twofq:
                        out = ' '.join(['hisat2', '-q','--max-intronlen','10000', '-p',str(cpu_count()-1),'-x',hisat2_index_file,'-1', fq.fqfile,'-2',fq.fqfile2])
                        print(out)
                        call(['hisat2', '-q','--max-intronlen','10000', '-p',str(cpu_count()-1),'-x',hisat2_index_file,'-1', fq.fqfile,'-2',fq.fqfile2], stdout = open(tfile1, 'w'))
                        
                else:
                    out = ' '.join(['hisat2', '-q','--max-intronlen','10000', '-p',str(cpu_count()-1),'-x',hisat2_index_file,'-U', fq.fqfile])
                    print(out)
                    call(['hisat2', '-q','--max-intronlen','10000', '-p',str(cpu_count()-1),'-x',hisat2_index_file,'-U', fq.fqfile], stdout = open(tfile1, 'w'))
                tfile2 = FM.ret_temp_file('.bam')
                tfile3 = FM.ret_temp_file('.bam')
                p1 = Popen(['samtools', 'view', '-bh', '-@', str(cpu_count()), tfile1], stdout=PIPE)
                p2 = Popen(['samtools', 'sort','-o', tfile2, '-@', str(cpu_count()), '-'], stdin = p1.stdout)
                p2.communicate()
                call(['samtools','rmdup', tfile2, tfile3])

                #call(['java', '-jar','/usr/local/share/java/picard.jar','SortSam', 'I=' + tfile1, 'O=' + tfile2, 'SORT_ORDER=coordinate'])
                call(['rm', '-f', tfile1, tfile2, 'temp.txt'])
                bamfiles.append(tfile3)
            
            out_file = FM.ret_asb(sample, self.species, self.genome_version, 'all')
            if len(bamfiles) > 1:
                call(['samtools', 'merge', '-f', out_file] + bamfiles)
            else:
                call(['mv', bamfiles[0], out_file])
            call(['samtools', 'index', out_file])
            call(['rm', '-f'] + bamfiles)    

    def split_bamfiles(self):
        for sample in self.samples:
            align_file = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.genome_version, 'all'))
   
            unmapped = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.genome_version, 'unmapped'), mode = 'wb', template = align_file) 
            discordant = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.genome_version, 'discordant'), mode = 'wb', template = align_file)
            inversion = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.genome_version, 'inversion'), mode = 'wb', template = align_file)
            duplication = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.genome_version, 'duplication'), mode = 'wb', template = align_file)
            clipped = pysam.AlignmentFile(FM.ret_asb(sample, self.species, self.genome_version, 'clipped'), mode = 'wb', template = align_file)
    
            for read in align_file.fetch(until_eof=True):
                # Partially mapped reads
                if read.is_unmapped and read.mate_is_unmapped:
                    unmapped.write(read)
                elif read.is_unmapped or read.mate_is_unmapped:
                    discordant.write(read)                
                # Chromosome fusion
                elif read.reference_id!=read.next_reference_id:
                    discordant.write(read)
                # Inversion
                elif read.is_reverse == read.mate_is_reverse:
                    inversion.write(read)
                # Duplication
                elif ((read.pos < read.mpos and read.is_reverse) or (read.pos > read.mpos and read.mate_is_reverse)) and abs(read.isize) > 102:
                    duplication.write(read)
                elif 'S' in read.cigarstring or 'H' in read.cigarstring:
                    clipped.write(read)

            unmapped.close()
            discordant.close()
            inversion.close()
            duplication.close()
            clipped.close()
            align_file.close()

            call(['samtools','index', FM.ret_asb(sample, self.species, self.genome_version, 'unmapped')])                            
            call(['samtools','index', FM.ret_asb(sample, self.species, self.genome_version, 'discordant')])                            
            call(['samtools','index', FM.ret_asb(sample, self.species, self.genome_version, 'inversion')])                            
            call(['samtools','index', FM.ret_asb(sample, self.species, self.genome_version, 'duplication')])                            
            call(['samtools','index', FM.ret_asb(sample, self.species, self.genome_version, 'clipped')])        

    def merge_groups(self):
        if self.group is False:
            return
        [fg,sg] = self.group.split('_')
        trues = []
        falses = []
        bf_True, bf_False = FM.ret_asb(self.group, self.genome_version, 'discordant')
        for sample, grouping in SI.groupings[self.group].items:
            if grouping == True:
                trues.append(FM.ret_asb(sample, self.genome_version, 'discordant'))
            elif grouping == False:
                falses.append(FM.ret_asb(sample, self.genome_version, 'discordant'))
        call(['samtools', 'merge', bf_True, trues])
        call(['samtools', 'merge', bf_False, falses])
