import os
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import multiprocessing
import time
import argparse

rscript = '/'.join(os.path.abspath(__file__).split('/')[:-1]) + '/Helper_Scripts/constructNJDendro.r'
VALID_BASES = set(['A', 'C', 'G', 'T'])
all_4mers = ["AAAA","AAAC","AAAG","AAAT","AACA","AACC","AACG","AACT","AAGA","AAGC","AAGG","AAGT","AATA",
			 "AATC","AATG","AATT","ACAA","ACAC","ACAG","ACAT","ACCA","ACCC","ACCG","ACCT","ACGA","ACGC",
			 "ACGG","ACGT","ACTA","ACTC","ACTG","AGAA","AGAC","AGAG","AGAT","AGCA","AGCC","AGCG","AGCT",
			 "AGGA","AGGC","AGGG","AGTA","AGTC","AGTG","ATAA","ATAC","ATAG","ATAT","ATCA","ATCC","ATCG",
			 "ATGA","ATGC","ATGG","ATTA","ATTC","ATTG","CAAA","CAAC","CAAG","CACA","CACC","CACG","CAGA",
			 "CAGC","CAGG","CATA","CATC","CATG","CCAA","CCAC","CCAG","CCCA","CCCC","CCCG","CCGA","CCGC",
			 "CCGG","CCTA","CCTC","CGAA","CGAC","CGAG","CGCA","CGCC","CGCG","CGGA","CGGC","CGTA","CGTC",
			 "CTAA","CTAC","CTAG","CTCA","CTCC","CTGA","CTGC","CTTA","CTTC","GAAA","GAAC","GACA","GACC",
			 "GAGA","GAGC","GATA","GATC","GCAA","GCAC","GCCA","GCCC","GCGA","GCGC","GCTA","GGAA","GGAC",
			 "GGCA","GGCC","GGGA","GGTA","GTAA","GTAC","GTCA","GTGA","GTTA","TAAA","TACA","TAGA","TATA",
			 "TCAA","TCCA","TCGA","TGAA","TGCA","TTAA"]

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: MetaDownSampler.py
	Author: Rauf Salamzade
    UW Madison MDTP / Rotating in Anantharaman Lab

	This program creates a neighbor-joining dendrogram from TNFs.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--fasta', type=str, help="Multi-FASTA of viral genomes.", required=True)
	parser.add_argument('-o', '--outdir', type=str, help='Path to output directory.', required=True)
	parser.add_argument('-t', '--threads', type=int, help='Number of cores.', required=True)
	args = parser.parse_args()
	return args

def calc_TNF_for_sample(input):
	name, seq, outf = input
	outf_handle = open(outf, 'w')
	genome_kmer_counts = defaultdict(int)
	start = 0
	stop = 4
	while stop <= len(seq):
		kmer = ''
		if stop == len(seq):
			kmer = seq[start:]
		else:
			kmer = seq[start:stop]
		if len(set(kmer).difference(VALID_BASES)) > 0: continue
		rckmer = str(Seq(kmer).reverse_complement())
		kmer_id = sorted([kmer, rckmer])[0]
		genome_kmer_counts[kmer_id] += 1
		start += 1
		stop += 1

	total_kmers = sum(genome_kmer_counts.values())
	printlist = [name]
	for k in all_4mers:
		printlist.append(str(float(genome_kmer_counts[k]) / total_kmers))
	outf_handle.write('\t'.join(printlist) + '\n')

"""
PARSE OPTIONS/INPUT
"""

myargs = create_parser()

input_fasta = os.path.abspath(myargs.fasta)
outdir = os.path.abspath(myargs.outdir) + '/'
threads = myargs.threads

tmpdir = outdir + 'tmp/'

if os.path.isdir(outdir):
	sys.stderr.write("Output directory already exists. Overwriting in five seconds ...\n")
	time.sleep(5)
else:
	os.system('mkdir %s' % outdir)
	os.system('mkdir %s' % tmpdir)

viral_seqs = []
with open(input_fasta) as off:
	for i, rec in enumerate(SeqIO.parse(off, 'fasta')):
		seq = str(rec.seq).upper()
		outf = tmpdir + 'VG_' + str(i) + '.txt'
		viral_seqs.append([rec.id, seq, outf])

sys.stderr.write("Starting to calculate TNFs...\n")
p = multiprocessing.Pool(threads)
p.map(calc_TNF_for_sample, viral_seqs)

tnf_file = outdir + 'NNFs.txt'
tnf_handle = open(tnf_file, 'w')
tnf_handle.write('sample\t' + '\t'.join(all_4mers) + '\n')
tnf_handle.close()
os.system("find %s -name '*.txt' -exec cat '{}' ';' >> %s" % (tmpdir, tnf_file))
sys.stderr.write("Starting to calculate distance matrix ...\n")

sys.stderr.write("Starting neighbor-joining tree construction.\n")
midpoint_tree_file = outdir  + 'TNF_NJ_Phylogeny.midpoint.tre'
os.system('Rscript %s %s %s' % (rscript, tnf_file, midpoint_tree_file))