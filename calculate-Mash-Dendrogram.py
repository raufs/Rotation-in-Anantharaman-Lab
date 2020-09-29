import os
import sys
import time
import argparse
from collections import defaultdict
rscript = '/'.join(os.path.abspath(__file__).split('/')[:-1]) + '/Helper_Scripts/constructNJDendroFromDist.r'

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: MetaDownSampler.py
	Author: Rauf Salamzade
    UW Madison MDTP / Rotating in Anantharaman Lab

	This program creates a neighbor-joining dendrogram from MASH distances (Novem-Nucleotide-Frequencies).
	Because it was designed to look at relationships between viral strains/species, we used a k-mer size of 9, 
	which was found to be optimal for looking at viral relationships by 
	Zhang et al. 2017: Viral Phylogenomics Using an Alignment-Free Method: A Three-Step Approach to 
	Determine Optimal Length of k-mer.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--fasta', type=str, help="Multi-FASTA of viral genomes.", required=True)
	parser.add_argument('-o', '--outdir', type=str, help='Path to output directory.', required=True)
	parser.add_argument('-k', '--ksize', type=int, help='K-mer size.', default = 9, required=False)
	parser.add_argument('-s', '--sketch', type=int, help='Sketch size for MASH.', default=10000, required=False)
	parser.add_argument('-t', '--threads', type=int, help='Number of cores.', default=1, required=False)
	args = parser.parse_args()
	return args

"""
PARSE OPTIONS/INPUT
"""

myargs = create_parser()

input_fasta = os.path.abspath(myargs.fasta)
outdir = os.path.abspath(myargs.outdir) + '/'
threads = myargs.threads
ksize = myargs.ksize
sketch = myargs.sketch

if os.path.isdir(outdir):
	sys.stderr.write("Output directory already exists. Overwriting in five seconds ...\n")
	time.sleep(5)
else:
	os.system('mkdir %s' % outdir)

sys.stderr.write('Running MASH to Compute Distances.\n')
mash_sketch = outdir + 'genomes.msh'
mash_distance_file = outdir + 'mash_distances.txt'
os.system('mash sketch -i -s %d -k %d -o %s/genomes %s' % (sketch, ksize, outdir, input_fasta))
os.system('mash dist -p %s %s %s > %s' % (threads, mash_sketch, mash_sketch, mash_distance_file))

sys.stderr.write("Writing Distance Matrix.\n")
distances = defaultdict(lambda: defaultdict(lambda: "NA"))
samples = set([])
with open(mash_distance_file) as omf:
	for line in omf:
		line = line.strip()
		ls = line.split('\t')
		samp1, samp2, dist = ls[:3]
		distances[samp1][samp2] = dist
		samples.add(samp1)
		samples.add(samp2)

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