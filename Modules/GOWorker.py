import os



class GoMaker():
    def __init__(self, genome_version, species, study_group, population_group):
        self.genome_version = genome_version
        self.species = species
        self.study_group = study_group
        self.population_group = population_group
        if int(genome_version.split('WS')[1]) < 253:
            print('current program only work with WS253 or later', file = sys.stderr)
            sys.exit()
        
    def GO_analysis(self):
        print('not yet done')
'''
        command = 'find_enrichment.py'
        command += ' --indent'
        obo_file = '/Volumes/Lijiang_data/packages/Elegans_RNA_seq_LL/Databases/c_elegans/WS255/gene_ontology.WS255.obo'
        command += ' --obo ' + obo_file
        output_file = '~/Desktop/go_analysis'
        command += ' --outfile '+ output_file
        study_group = '/Volumes/Lijiang_data/RNAseq_analysis/6-5-2017/GO_analysis/study_group.txt'
        population_group = '/Volumes/Lijiang_data/RNAseq_analysis/6-5-2017/GO_analysis/population_group.txt'
        association_file = '/Volumes/Lijiang_data/packages/Elegans_RNA_seq_LL/Databases/c_elegans/WS255/gene_association.WS255.wb.c_elegans'
        command += ' '+study_group + ' ' + population_group + ' ' + association_file
        print(command)
        os.system(command)
'''


