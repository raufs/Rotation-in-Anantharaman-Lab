import os
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import multiprocessing
import time
import argparse
import itertools
import subprocess

rscript = '/'.join(os.path.abspath(__file__).split('/')[:-1]) + '/Helper_Scripts/constructNJDendroFromDist.r'
VALID_BASES = set(['A', 'C', 'G', 'T'])
all_9mers = set([])

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: MetaDownSampler.py
	Author: Rauf Salamzade
    UW Madison MDTP / Rotating in Anantharaman Lab

	This program creates a neighbor-joining dendrogram from NNFs (Novem-Nucleotide-Frequencies).
	Inspired by Zhang et al. 2017: Viral Phylogenomics Using an Alignment-Free Method: A Three-Step Approach to 
	Determine Optimal Length of k-mer.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--fasta', type=str, help="Multi-FASTA of viral genomes.", required=True)
	parser.add_argument('-o', '--outdir', type=str, help='Path to output directory.', required=True)
	parser.add_argument('-t', '--threads', type=int, help='Number of cores.', required=True)
	args = parser.parse_args()
	return args

def calc_NNF_for_sample(input):
	name, tmp_fasta, outf = input
	outf_handle = open(outf, 'w')
	genome_kmer_counts = defaultdict(int)

	cmd1 = 'jellyfish count -o %s -m %d -s 4000M -t %d -C %s' % (outf + '.jf', 9, 1, tmp_fasta)
	cmd2 = 'jellyfish dump -L %d -c -t %s -o %s' % (1, outf + '.jf', outf + '.kmers.tab')

	subprocess.call(cmd1, shell=True)
	assert (os.path.isfile(outf + '.jf'))
	subprocess.call('rm -f ' + tmp_fasta, shell=True)
	subprocess.call(cmd2, shell=True)
	assert (os.path.isfile(outf + '.kmers.tab'))

	genome_kmer_counts = defaultdict(int)
	with open(outf + '.kmers.tab') as oskf:
		for line in oskf:
			line = line.strip()
			ls = line.split('\t')
			k = ls[0]
			rck = str(Seq(k).reverse_complement())
			alpha = sorted([k, rck])[0]
			genome_kmer_counts[alpha] += int(ls[1])

	total_kmers = sum(genome_kmer_counts.values())
	printlist = [name]
	for k in all_9mers:
		printlist.append(str(float(genome_kmer_counts[k]) / total_kmers))
	outf_handle.write('\t'.join(printlist) + '\n')

"""
PARSE OPTIONS/INPUT
"""

myargs = create_parser()

input_fasta = os.path.abspath(myargs.fasta)
outdir = os.path.abspath(myargs.outdir) + '/'
threads = myargs.threads
os.environ["MKL_NUM_THREADS"] = str(threads)
os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
os.environ["OMP_NUM_THREADS"] = str(threads)
from scipy.spatial import distance


for k in itertools.product('ACGT', repeat=9):
	k = ''.join(list(k))
	rck = str(Seq(k).reverse_complement())
	alpha = sorted([k, rck])[0]
	if not alpha in all_9mers:
		all_9mers.add(alpha)

all_9mers = list(sorted(all_9mers))

tmpdir = outdir + 'tmp/'
if os.path.isdir(outdir):
	sys.stderr.write("Output directory already exists. Overwriting in five seconds ...\n")
	time.sleep(5)
else:
	os.system('mkdir %s' % outdir)
	os.system('mkdir %s' % tmpdir)

viral_seqs = []
sample_results = set([])
with open(input_fasta) as off:
	for i, rec in enumerate(SeqIO.parse(off, 'fasta')):
		seq = str(rec.seq).upper()
		outf_fasta = tmpdir + 'VG_' + str(i) + '.fasta'
		outf_fasta_handle= open(outf_fasta, 'w')
		outf_fasta_handle.write('>' + rec.id + '\n' + seq + '\n')
		outf_fasta_handle.close()
		outf = tmpdir + 'VG_' + str(i) + '.txt'
		sample_results.add(outf)
		viral_seqs.append([rec.id, outf_fasta, outf])

sys.stderr.write("Starting to calculate NNFs...\n")
#p = multiprocessing.Pool(threads)
#p.map(calc_NNF_for_sample, viral_seqs)

samples = set([])
distances = defaultdict(lambda: defaultdict(lambda: "NA"))
for i, s1_f in enumerate(sample_results):
	if i % 100 == 0:
		sys.stderr.write("Have processed distances for %d of %d samples.\n" % (i, len(sample_results)))
	s1_name = None
	s1_freqs = []
	with open(s1_f) as os1f:
		for line in os1f:
			line = line.strip()
			ls = line.split('\t')
			s1_name = ls[0]
			s1_freqs = [float(x) for x in ls[1:]]
	samples.add(s1_name)
	for j, s2_f in enumerate(sample_results):
		if i <= j:
			s2_name = None
			s2_freqs = []
			with open(s2_f) as os2f:
				for line in os2f:
					line = line.strip()
					ls = line.split('\t')
					s2_name = ls[0]
					s2_freqs = [float(x) for x in ls[1:]]
			pair_dist = distance.euclidean(s1_freqs, s2_freqs)
			distances[s1_name][s2_name] = pair_dist

dist_file = outdir + 'sample_distances.txt'
dist_handle = open(dist_file, 'w')

for s1 in sorted(samples):
	printlist = [s1]
	for s2 in sorted(samples):
		dist = distances[s1][s2]
		printlist.append(str(dist))
	dist_handle.write('\t'.join(printlist) + '\n')
dist_handle.close()

sys.stderr.write("Starting neighbor-joining tree construction.\n")
midpoint_tree_file = outdir  + 'NNF_NJ_Dendrogram.midpoint.tre'
os.system('Rscript %s %s %s' % (rscript, dist_file, midpoint_tree_file))