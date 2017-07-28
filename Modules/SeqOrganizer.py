import os, sys, random, inspect, xlrd, pandas
from collections import defaultdict
import pdb
# Where databases are stored
if os.getenv('NGS_ELEGANS_LD') is None:
    l_db_directory = os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/Databases/'
else:
    l_db_directory =  os.getenv('NGS_ELEGANS_LD')

# Where read files are stored
if os.getenv('NGS_ELEGANS_LR') is None:
    l_read_directory = os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/Reads/'
else:
    l_read_directory =  os.getenv('NGS_ELEGANS_LR')

# Where bam files are stored
if os.getenv('NGS_ELEGANS_LB') is None:
    l_bam_directory = os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/Output/'
else:
    l_bam_directory =  os.getenv('NGS_ELEGANS_LB')

# Where polymorphism files are stored
if os.getenv('NGS_ELEGANS_LP') is None:
    l_poly_directory = os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/Polymorphisms/'
else:
    l_poly_directory =  os.getenv('NGS_ELEGANS_LP')
    
# Where read count files are stored
if os.getenv('NGS_ELEGANS_LC') is None:
    l_count_directory = os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/Gene_read_count/'
else:
    l_count_directory =  os.getenv('NGS_ELEGANS_LC')    

if os.getenv('NGS_ELEGANS_RG') is None:
    l_rg_directory = os.path.abspath(inspect.getfile(inspect.currentframe())).split('/Modules')[0] + '/region_genes/'
else:
    l_rg_directory =  os.getenv('NGS_ELEGANS_RG')    



class RefObj():
    def __init__(self, species, genome_version, ref_file, annotation_file):
        self.species = species
        self.genome_version = genome_version
        self.ref_file = l_db_directory + '/' + species + '/' + genome_version + '/' + ref_file
        self.annotation_file = l_db_directory + '/' + species + '/' + genome_version + '/' + annotation_file
        self.unzipped_ref_file = self.ref_file.split('.gz')[0]



class RefDatabase():
    #initiate reference database
    def __init__(self):
        self.ref_database =  l_db_directory + '/AvailableDatabases.xlsx'
        self.ws = 'Databases'
        self.databases = {}
        self.read_database_file()
    #get all available databases
    def read_database_file(self):
        wb = xlrd.open_workbook(self.ref_database)
        sheet = wb.sheet_by_name(self.ws)
        for r_idx in range(1,sheet.nrows):
            ref_species = str(sheet.cell(r_idx,0).value)
            genome_version = str(sheet.cell(r_idx,1).value)
            ref_file = str(sheet.cell(r_idx,2).value)
            annotation_file = str(sheet.cell(r_idx,3).value)
            self.databases[(ref_species,genome_version)] = RefObj(ref_species, genome_version, ref_file, annotation_file)

    def ret_genome_versions(self, species, gv):
        if (species,gv) not in self.databases:
            print('Cant find ' + species + ' ' + gv, file = sys.stderr)
            sys.exit()
        else:
            return self.databases[(species, gv)]

    def valid_genome_version(self, species, gv):
        result = True
        if (species,gv) not in self.databases:
                #print(gv + ' not found', file = sys.stderr)
                result = False
        return result
    
    def add_reference(self, species, gv, ref_file, annotation_file):
        dt = pandas.read_excel(self.ref_database)
        print(dt)
        dt.loc[len(dt)] = [species, gv, ref_file, annotation_file]
        dt.to_excel(self.ref_database, sheet_name = 'Databases', index = False)
    
databaseInfo = RefDatabase()
class SampleObj():
    
    def __init__(self, sample, RG_info, paired_flag, fqfile, fqfile2,mean_frag,sd_frag,twofq = False):
        self.sample = sample
        self.RG_info = RG_info
        self.paired_flag = paired_flag
        self.fqfile = fqfile
        self.fqfile2 = fqfile2
        self.mean_frag = mean_frag
        self.sd_frag = sd_frag
        self.twofq = twofq

        
												
class SampleDatabase():
    def __init__(self):
        self.sample_database = l_read_directory + '/AvailableSamples.xlsx'
        self.ws = 'Samples'
        self.samples = defaultdict(list)
        self.groupings = {}
        self.read_sample_file()
		

    def read_sample_file(self):
        wb = xlrd.open_workbook(self.sample_database)
        sheet = wb.sheet_by_name(self.ws)
#         index_to_grouping = {}
#         for i in range(5, sheet.ncols):
#             self.groupings[str(sheet.cell(0,i).value)] = {}
#             index_to_grouping[i] = str(sheet.cell(0,i).value)
        for r_idx in range(1, sheet.nrows):
            sample = str(sheet.cell(r_idx,0).value)
            RG_info = str(sheet.cell(r_idx,1).value)
            paired_flag = int(sheet.cell(r_idx,2).value)
            fqfile = l_read_directory + str(sheet.cell(r_idx,3).value)
            fqfile2 = l_read_directory + str(sheet.cell(r_idx,4).value)
            mean_frag = sheet.cell(r_idx,7).value
            sd_frag = sheet.cell(r_idx,8).value
            if fqfile2 == l_read_directory:
                twofq = False
            elif fqfile2.split('.')[-1] not in ['fq','fastq', 'gz']:
                print('Error? 2nd fastq file does not seem correct: ' + fqfile2, file = sys.stderr)
            else:
                twofq = True
#             for i in range(5, sheet.ncols):
#                 t_group = str(sheet.cell(r_idx,i).value)
#                 if t_group.lower() in ['true','1']:
#                     self.groupings[index_to_grouping[i]][sample] = True
#                 elif t_group.lower() in ['false','0']:
#                     self.groupings[index_to_grouping[i]][sample] = False
#                 elif t_group.lower() == 'none':
#                     self.groupings[index_to_grouping[i]][sample] = None
#                 else:
#                     print(t_group + ' not an allowed description', file = sys.stderr)
#                     sys.exit()
            self.samples[sample].append(SampleObj(sample, RG_info, paired_flag, fqfile, fqfile2, mean_frag,sd_frag,twofq))


    def valid_group(self, group):
        if group in self.groupings:
            out_samples = []
            for sample, truth in self.groupings[group].items():
                if truth is not None:
                    out_samples.append(sample)
            return sorted(out_samples)
        else:
            return False

    def valid_samples(self, samples):
        result = True
        for sample in samples:
            if sample not in self.samples:
                print(sample + ' not found', file = sys.stderr)
                print('Available options are ' + str(self.samples), file = sys.stderr)
                result = False
        return result

    def ret_fastq_objs(self, sample):
        if sample in self.samples:
            return self.samples[sample]
        
sampleInfo = SampleDatabase()

class FileMaker():
    def __init__(self):
        pass

    def genome_error(self, species, genome_version):
        print(species + '\t' + genome_version + ' does not exits in reference database. Options are:', file = sys.stderr)
        print(', '.join(list(databaseInfo.databases.keys())), file = sys.stderr)
        #sys.exit()

    def sample_error(self, sample):
        print(sample + ' does not exist in sample database. Available samples are:', file = sys.stderr)
        print(', '.join(set(list(sampleInfo.samples.keys()))), file = sys.stderr)
        print('Available groups are:', file = sys.stderr)
        print(', '.join(set(list(sampleInfo.groupings.keys()))), file = sys.stderr)
        #sys.exit()

    def ret_genes_in_region(self,species, genome_version,chrom,start, end):
        file_path = l_rg_directory+ species + '/' + genome_version + '/' + chrom + '_' + str(start) +'-' + str(end)
        if not os.path.isdir(l_rg_directory):
            os.makedirs(l_rg_directory)
        if not os.path.isdir(l_rg_directory + species):
            os.makedirs(l_rg_directory + species)
        if not os.path.isdir(l_rg_directory + species + '/' + genome_version):
            os.makedirs(l_rg_directory + species + '/' + genome_version)
        return file_path
    
    def ret_srcf(self, sample, species, genome_version):
        return self.ret_strain_read_count_file(sample, species, genome_version)
    
    def ret_strain_read_count_file(self, sample, species, genome_version):
        file_path = l_count_directory + species + '/' + genome_version + '/' + sample + '/' + sample
        if not os.path.isdir(l_count_directory):
            os.makedirs(l_count_directory)
        if not os.path.isdir(l_count_directory + species):
            os.makedirs(l_count_directory + species)
        if not os.path.isdir(l_count_directory + species + '/' + genome_version):
            os.makedirs(l_count_directory + species + '/' + genome_version)
        if not os.path.isdir(l_count_directory + species + '/' + genome_version + '/' + sample):
            os.makedirs(l_count_directory + species + '/' + genome_version + '/' + sample)
        
        return file_path
        
    def ret_asb(self, sample, species, genome_version, bam_type):
        return self.ret_aligned_strain_bamfile(sample, species, genome_version, bam_type)

    def ret_aligned_strain_bamfile(self, sample, species, genome_version, bam_type):

        base_bam_dir = l_bam_directory + species + '/' + genome_version + '/' + sample + '/' + sample
        
        type_name = {}       
        type_name['all'] = base_bam_dir + '.merged.bam'
        type_name['discordant'] = base_bam_dir + '.merged.discordant.bam'
        type_name['clipped'] = base_bam_dir + '.merged.clipped.bam'
        type_name['duplication'] = base_bam_dir + '.merged.duplication.bam'
        type_name['inversion'] = base_bam_dir + '.merged.inversion.bam'
        type_name['unmapped'] = base_bam_dir + '.merged.unmapped.bam'
        type_name['region'] = base_bam_dir + '.region.bam'
     
        if databaseInfo.valid_genome_version(species, genome_version) is False:
            self.genome_error(species, genome_version)
        if sampleInfo.valid_samples([sample]) is False and sampleInfo.valid_group(sample) is False:
            self.sample_error(sample)
        if bam_type not in type_name:
            print(bam_type + ' not a legitimate option. Valid options are:')
            print(type_name.keys())
            sys.exit()
  
        if not os.path.isdir(l_bam_directory):
            os.makedirs(l_bam_directory)
        if not os.path.isdir(l_bam_directory + species):
            os.makedirs(l_bam_directory + species)
        if not os.path.isdir(l_bam_directory + species + '/' + genome_version):
            os.makedirs(l_bam_directory + species + '/' + genome_version)
        if not os.path.isdir(l_bam_directory + species + '/' + genome_version + '/' + sample):
            os.makedirs(l_bam_directory + species + '/' + genome_version + '/' + sample)

        if sampleInfo.valid_group(sample):
            bf1 = l_bam_directory + species + '/' + genome_version + '/' + sample + '/' + sample.split('_')[0] + '.merged.' + bam_type + '.bam'
            bf2 = l_bam_directory + species + '/' + genome_version + '/' + sample + '/' + sample.split('_')[1] + '.merged.' + bam_type + '.bam'
            return[bf1,bf2]

        return type_name[bam_type]
   
    def ret_temp_file(self, extension = '.sam'):  
        if not os.path.isdir(l_bam_directory):
            os.makedirs(l_bam_directory)
        if not os.path.isdir(l_bam_directory + 'Temp'):
            os.makedirs(l_bam_directory + 'Temp')
      
        return l_bam_directory + 'Temp/TempFile' + str(random.randint(1, 999999999)) + extension
 
    def ret_vcf_file(self, prefix, species, gv, discovery_samples, excluded_samples, genotype_samples):
        if type(discovery_samples) is str:
            discovery_samples = [discovery_samples]
        if type(excluded_samples) is str:
            excluded_samples = [excluded_samples]
        if type(genotype_samples) is str:
            genotype_samples = [genotype_samples]

        if not os.path.isdir(l_poly_directory):
            os.makedirs(l_poly_directory)
        if not os.path.isdir(l_poly_directory + species):
            os.makedirs(l_poly_directory+species)
        if not os.path.isdir(l_poly_directory + species + '/' + gv):
            os.makedirs(l_poly_directory+species + '/' + gv)
        discovery = ';'.join(sorted(discovery_samples))
        genotyped = ';'.join(sorted(genotype_samples))
        excluded = ';'.join(sorted(excluded_samples))
            
        return l_poly_directory + species + '/' + gv + '/' + prefix + '.Sp:' + species + '.GV:' +gv + '.DS:'+discovery  + '.ES:' + excluded + '.GS:'+genotyped + '.vcf'
            
file_maker = FileMaker()
