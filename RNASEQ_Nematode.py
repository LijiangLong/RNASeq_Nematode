import argparse
#import Modules.SeqOrganizer as SO
import Modules.DatabaseWorker as DBW
import Modules.AlignmentWorker as AW
import Modules.ReadcountWorker as RCW
import Modules.RegionWorker as RW
import Modules.GOWorker as GW
import Modules.KallistoWorker as KW
import pdb

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='Available Commands', dest='command')

db_parser = subparsers.add_parser('db', help='Download or remove databases from Wormbase')
db_parser.add_argument('Species', type = str, help = 'Species you would like to download. Currently only one species can be downloaded at a time')
db_parser.add_argument('Genome_versions', nargs = '+', type = str, help = 'Version of genome(s) you want to download from Wormbase or remove from your computer')
#db_parser.add_argument('-rl', '--removelocal', help = 'Remove database files from your local computer', action="store_true")
#db_parser.add_argument('-rr', '--removeremote', help = 'Remove database files from your remote file server', action="store_true")


alignRNA_parser = subparsers.add_parser('alignRNA', help='Align RNA sequencing reads to a given database')
alignRNA_parser.add_argument('Genome_version', type = str, help = 'Version of the genome you want to use. Assumes appropriate files are installed already with the db command')
alignRNA_parser.add_argument('Species', type = str, help = 'Species you want to align to. Assumes appropriate files are installed already with the db command')
alignRNA_parser.add_argument('Samples', type = str, nargs = '+',  help = 'List of strains you would like to be analyzed. all for all strains.')
alignRNA_parser.add_argument('-s', '--split',  help = 'Just split already existing output bamfile', action="store_true")

readcount_parser = subparsers.add_parser('read_count', help='count reads using mapped bam file')
readcount_parser.add_argument('Genome_version', type = str, help = 'Version of the genome you want to use. Assumes appropriate files are installed already with the db command')
readcount_parser.add_argument('Species', type = str, help = 'Species you want to align to. Assumes appropriate files are installed already with the db command')
readcount_parser.add_argument('Strand', type = str, help = '<yes/no/reverse> whether the data is from a strand-specific assay')
readcount_parser.add_argument('Order', type = str, help = 'For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the samtools sort function of samtools to sort it. Use this option, with name or pos for <order> to indicate how the input data has been sorted. The default is name.')
readcount_parser.add_argument('Samples', type = str, nargs = '+',  help = 'List of samples you would like to be analyzed.')

region_parser = subparsers.add_parser('region', help='analyse a region on one chromosome')
region_parser.add_argument('Genome_version', type = str, help = 'Version of the genome you want to use. Assumes appropriate files are installed already with the db command')
region_parser.add_argument('Species', type = str, help = 'Species you want to align to. Assumes appropriate files are installed already with the db command')
region_parser.add_argument('Chrom', type = str, help = 'chromosome of the region. For example II for chromosomeII')
region_parser.add_argument('Coordinates', type = str, help = 'coodinates of the region. For example 1-100')

Goanalysis_parser = subparsers.add_parser('GO', help='GO analysis of a given list of genes')
Goanalysis_parser.add_argument('Genome_version', type = str, help = 'Version of the genome you want to use. Assumes appropriate files are installed already with the db command')
Goanalysis_parser.add_argument('Species', type = str, help = 'Species you want to analysis')
Goanalysis_parser.add_argument('study_group', type = str, help = 'gene names in a study')
Goanalysis_parser.add_argument('population_group', type = str, help = 'gene names in population (or other study if --compare is specified)')

Kallisto_parser = subparsers.add_parser('kallisto', help='kallisto analysis of isoform expression')
Kallisto_parser.add_argument('Genome_version', type = str, help = 'Version of the genome you want to use. Assumes appropriate files are installed already with the db command')
Kallisto_parser.add_argument('Species', type = str, help = 'Species you want to analysis')
Kallisto_parser.add_argument('Samples', type = str, nargs = '+',  help = 'List of strains you would like to be analyzed. all for all strains.')


args = parser.parse_args()



if args.command is None:
	parser.print_help()

if args.command == 'db':
    db_obj = DBW.DBMaker(args.Species, args.Genome_versions, False, False)
    #db_obj.execute()

if args.command == 'alignRNA':
    aln_obj = AW.AlignRNAMaker(args.Genome_version, args.Species, args.Samples, args.split)
    aln_obj.align()
#    aln_obj.split_bamfiles()

if args.command == 'read_count':
    read_count_obj = RCW.ReadcountMaker(args.Genome_version, args.Species, args.Strand,args.Order,args.Samples)
    read_count_obj.count()
    
if args.command == 'region':
    region_obj = RW.RegionMaker(args.Genome_version, args.Species,args.Chrom, args.Coordinates)
    region_obj.gene_in_region()
    
if args.command == 'GO':
    region_obj = GW.GoMaker(args.Genome_version, args.Species,args.study_group, args.population_group)
    region_obj.GO_analysis()
    
if args.command == 'kallisto':
    kallisto_obj = KW.KallistoMaker(args.Genome_version, args.Species, args.Samples)
    kallisto_obj.count()
    


