from Modules.SeqOrganizer import databaseInfo as DI
import Modules.SeqOrganizer as SO
from ftplib import FTP
from subprocess import call
import os, sys, xlrd
import pdb


class DBMaker():
    def __init__(self, species, genome_versions, remove_local_flag, remove_remote_flag):
        self.species = species
        self.Genome_version = genome_versions
        self.rlf = remove_local_flag
        self.rrf = remove_remote_flag
        self._ret_wormbase_species()
        if species not in self.available_wormbase_species:
            print(species + ' not found in WormBase ftp site.', file = sys.stderr)
            print('Available species are ' + ', '.join(sorted(self.available_wormbase_species)), file = sys.stderr)
            sys.exit()
        self._ret_wormbase_genome_versions()
        for genome_version in genome_versions:
            if genome_version not in self.available_genome_versions:
                print(genome_version + ' not found in WormBase ftp site.', file = sys.stderr)
                print('Available versions are ' + ', '.join(sorted(self.available_genome_versions)), file = sys.stderr)
                sys.exit()
            if DI.valid_genome_version(species,genome_version):
                print(species + ' ' + genome_version + ' already in SampleDatabase. Remove files first if you would like to redownload this version', file = sys.stderr)
                continue
            if not self.rlf:
                self._download_files(genome_version)
        
    def _ret_wormbase_species(self):
        ftp = FTP('ftp.wormbase.org')
        ftp.login()
        ftp.cwd('pub/wormbase/species')
        self.available_wormbase_species = ftp.nlst()
        ftp.quit()
        
    def _ret_wormbase_genome_versions(self):
        self.available_genome_versions = set()
        ftp = FTP('ftp.wormbase.org')
        ftp.login()
        ftp.cwd('pub/wormbase/species/' + self.species + '/gff')
        t_files = ftp.nlst()
        for t_file in t_files:
            try:
                self.available_genome_versions.add('WS' + t_file.split('.WS')[1].split('.')[0])
            except IndexError:
                continue
        ftp.quit()

    def _download_files(self, gv):
        local_db = SO.l_db_directory
        if not os.path.exists(local_db):
            os.makedirs(local_db)
        if not os.path.exists(local_db + self.species):
            os.makedirs(local_db + self.species)
        if not os.path.exists(local_db + self.species + '/' + gv):
            os.makedirs(local_db + self.species + '/' + gv)
            
        #Download ref_file
        #Find ref file
#         ftp = FTP('ftp.wormbase.org')
#         ftp.login()
#         ftp.cwd('pub/wormbase/species/' + self.species + '/sequence/genomic/')
#         t_files = ftp.nlst()
#         for t_file in t_files:
#             if gv + '.genomic.fa.gz' in t_file:
#                 ref_file = t_file
#                 break
        project='PRJNA13758'
        wormbase_prefix='ftp://ftp.wormbase.org/pub/wormbase/releases/'+gv+'/species/'+self.species+'/'+project+'/'
        file_prefix = self.species+'.'+project+'.'+gv
        ref_file = file_prefix+'.genomic.fa.gz'
        gtf_file = file_prefix+'.canonical_geneset.gtf.gz'
        ss_file = file_prefix+'.ss'
        exon_file = file_prefix+'.exon'
        l_ref = local_db + self.species + '/' + gv + '/' + ref_file
        l_gtf = local_db + self.species + '/' + gv + '/' + gtf_file
        l_ss  = local_db + self.species + '/' + gv + '/' + ss_file
        l_exon = local_db + self.species + '/' + gv + '/' + exon_file
        l_hisat2_index = local_db + self.species + '/' + gv + '/' + file_prefix + '.hisat2_index'
#        curl ${prefix}/c_elegans.${project}.${reference}.genomic.fa.gz 
        wormbase_ref = wormbase_prefix+'/'+ref_file
        wormbase_gtf = wormbase_prefix+'/'+gtf_file
        
        # Download reference fasta file
        print('Downloading ' + ref_file + ' from WormBase', file = sys.stderr)
        command = ' '.join(['curl',wormbase_ref,'>',l_ref])
        os.system(command)
       
        # Download from ftp site
        print('Downloading ' + gtf_file + ' from WormBase', file = sys.stderr)
        command = ' '.join(['curl',wormbase_gtf,'>',l_gtf])
        os.system(command)

        # Update reference database
        command = ' '.join(['gzcat',l_gtf,'|','hisat2_extract_splice_sites.py','-','>',l_ss])
        os.system(command)
        command = ' '.join(['gzcat',l_gtf,'|','hisat2_extract_exons.py','-','>',l_exon])
        os.system(command)
        command = ' '.join(['gunzip', local_db + self.species + '/' + gv + '/' + ref_file])
        os.system(command)
        command = ' '.join(['hisat2-build','-p','6' ,'--ss', l_ss, '--exon', l_exon, l_ref.split('.gz')[0], l_hisat2_index])
        os.system(command)
        DI.add_reference(self.species, gv, ref_file, gtf_file)
       
