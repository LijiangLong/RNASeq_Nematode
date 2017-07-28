# RNASeq_Nematode

six modules are currently there with the following functions:

db: download and pre-load the databases(genomic sequence, idx etc)
alignRNA: align RNAseq reads to genomic sequence using hisat2
read_count: count the number of reads of each gene using HTseq
region: subset genes in a given region
GO: GO analysis(not yet done)
kallisto: count and normalize number of reads of each isoform or gene using kallisto

To run this program, first download the necessary dependencies and then edit the excel in the Reads folder,
run it like:
python3 RNASeq_Nematody.py <module> [args]
