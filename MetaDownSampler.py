#!/usr/bin/env python

### Program: MetaDownSampler.py
### Author: Rauf Salamzade
### UW Madison MDTP / Rotating in Anantharaman Lab

"""
Special thanks to: Kris Kief and Patricia Tran
"""

import os
import sys
from time import sleep
from Bio import SeqIO
import logging
import time
import argparse
import pysam
from collections import defaultdict
import random
import gzip
import math

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: MetaDownSampler.py
	Author: Rauf Salamzade
    UW Madison MDTP / Rotating in Anantharaman Lab
    This program simulates downsampling reads in a metagenome sample to realistially fragment genomic assemblies.
    
    It takes as input metagenomic assembly in FASTA format, reflexive alignments in sorted and indexed BAM and optionally a 
    set of scaffolds of interest. It will then downsample the number of reads to a user-specified fold and return a
    fragmented assembly and a subset of reads. 
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-m', '--mag', type=str, help="FASTA (single entry) for reference scaffold upon which to call windows on.", required=True)
	parser.add_argument('-a', '--alignment', type=str, help='Reflexive alignment files in BAM format.', required=True)
	parser.add_argument('-o', '--outdir', type=str, help='path to output directory.', required=True)
	parser.add_argument('-b', '--bed', type=str, help='BED file with four columns listing the scaffold, start pos, stop pos, and id. Coordinate ranges on same scaffold should never overlap!', required=False, default=None)
	parser.add_argument('-f', '--downsample_folds', type=float, nargs='+', help='Provide list of downsampling folds. Values should be between 0 (no reads) and 1 (all reads).', required=False, default=[1.0, 0.25, 0.5, 0.75])
	parser.add_argument('-l', '--min_depth', type=int, help='Minimum depth required for scaffolds to be considered covered and prevent fragmenting scaffolds. Should be greater than 1.', required=False, default=1)
	parser.add_argument('-s', '--min_scaff_length', type=int, help='Minimum scaffold length to retain after fragmentation has been performed.', required=False, default=3000)
	args = parser.parse_args()
	return args

def main():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE OPTIONS/INPUT
	"""

	myargs = create_parser()

	mag_fasta = os.path.abspath(myargs.mag)
	alignment_file = os.path.abspath(myargs.alignment)
	bed_file = myargs.bed
	downsample_folds_list = myargs.downsample_folds
	min_depth = myargs.min_depth
	min_scaff_length = myargs.min_scaff_length
	outdir = os.path.abspath(myargs.outdir) + '/'

	downsampling_folds = []
	select_scaffolds = set([])
	select_scaffolds_coords = defaultdict(list)
	select_scaffolds_ids = defaultdict(list)
	try:
		assert(os.path.isfile(mag_fasta) and os.path.isfile(alignment_file))
		all_scaffolds = set([])
		try:
			with open(mag_fasta) as omf:
				for rec in SeqIO.parse(omf, 'fasta'):
					all_scaffolds.add(rec.id)
		except:
			raise RuntimeError("Error reading in Metagenomic Assembly as a FASTA file.")

		try:
			assert(min_depth > 0)
		except:
			raise RuntimeError("Error minimum depth requested for considered a scaffold position covered and likely to be assembled is below 1.")

		try:
			for fold in downsample_folds_list:
				assert(float(fold) > 0.0 and float(fold) <= 1.0)
				downsampling_folds.append(float(fold))
		except:
			raise RuntimeError("Error with downsample folds requested.")

		if bed_file != None:
			with open(bed_file) as obf:
				for line in obf:
					scaff, start, stop, id = line.strip().split('\t')
					try:
						assert (scaff in all_scaffolds)
						select_scaffolds.add(scaff)
						select_scaffolds_coords[scaff].append([int(start), int(stop)])
						select_scaffolds_ids[scaff].append(id)
					except:
						sys.stderr.write("Did not find scaffold %s in metagenomic assembly, ignoring it.\n" % scaff)
		else:
			with open(mag_fasta) as omf:
				for rec in SeqIO.parse(omf, 'fasta'):
					scaff = rec.id
					select_scaffolds.add(scaff)
					start = 1
					stop = len(str(rec.seq))
					select_scaffolds_coords[scaff].append([start, stop])
					select_scaffolds_ids[scaff].append(scaff)
	except:
		raise RuntimeError("Errors with reading in input files, perhaps paths are broken.")

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory already exists. Overwriting in five seconds ...\n")
		#time.sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	# Start program

	# Step 1: Read Alignment File and Get Total Paired-End Reads Per Select Scaffold

	align_handle = pysam.AlignmentFile(alignment_file, 'rb')

	pe_reads_per_scaffold = defaultdict(int)
	for scaff in all_scaffolds:
		total_paired_end = 0
		for read1, read2 in read_pair_generator(align_handle, scaff):
			pe_reads_per_scaffold[scaff] += 1

	total_paired_end_reads = sum(x for x in pe_reads_per_scaffold.values())

	# Step 2: Downsample Paired-End Reads from Total at Requested Folds and

	for fold in downsampling_folds:

		fold_dir = outdir + 'Downsampling_' + str(fold*100.0) + 'X/'
		if not os.path.isdir(fold_dir): os.system('mkdir %s' % fold_dir)

		r1_fastq_file = fold_dir + 'downsampled_R1.fastq'
		r2_fastq_file = fold_dir + 'downsampled_R2.fastq'
		assembly_file = fold_dir + 'downsampled_MAG.fasta'
		outf_r1 = open(r1_fastq_file, 'w')
		outf_r2 = open(r2_fastq_file, 'w')
		outf_assembly = open(assembly_file, 'w')

		downsampled_paired_end_reads = total_paired_end_reads*fold

		for scaff in select_scaffolds:
			scaff_range = set([])
			for coords in select_scaffolds_coords[scaff]:
				scaff_range = scaff_range.union(set(range(coords[0], coords[1] + 1)))
			scaffold_coverage = defaultdict(int)
			scaff_prop_total_pe_reads = pe_reads_per_scaffold[scaff]/float(total_paired_end_reads)
			scaff_downsampled_pe_reads = math.floor(scaff_prop_total_pe_reads*downsampled_paired_end_reads)
			range_of_pe_indices = list(range(0, pe_reads_per_scaffold[scaff]))
			downsampled_pe_read_selection = random.sample(range_of_pe_indices, scaff_downsampled_pe_reads)
			pe_index = 0
			for read1, read2 in read_pair_generator(align_handle, scaff):
				if not pe_index in downsampled_pe_read_selection:
					pe_index += 1
					continue
				query_name = read1.query_name
				r1_sequence = read1.query_sequence
				r1_qual = read1.qual
				r2_sequence = read2.query_sequence
				r2_qual = read2.qual

				location1 = read1.reference_start + 1
				readlen1 = read1.reference_length
				location2 = read1.reference_start + 1
				readlen2 = read1.reference_length
				insertion_location_left_r1 = location1
				insertion_location_right_r1 = location1 + readlen1
				insertion_location_left_r2 = location2
				insertion_location_right_r2 = location2 + readlen2
				read1_pos_range = set(range(insertion_location_left_r1, insertion_location_right_r1 + 1))
				read2_pos_range = set(range(insertion_location_left_r2, insertion_location_right_r2 + 1))

				if len(scaff_range.intersection(read1_pos_range)) > 0 or len(
						scaff_range.intersection(read2_pos_range)) > 0:
					outf_r1.write('@%s 1:N:0:AAAAAA\n%s\n+\n%s\n' % (query_name, r1_sequence, ''.join([str(x) for x in r1_qual])))
					outf_r2.write('@%s 2:N:0:AAAAAA\n%s\n+\n%s\n' % (query_name, r2_sequence, ''.join([str(x) for x in r2_qual])))

					for pos in range(insertion_location_left_r1, insertion_location_right_r1+1):
						scaffold_coverage[pos] += 1
					for pos in range(insertion_location_left_r2, insertion_location_right_r2+1):
						scaffold_coverage[pos] += 1

				pe_index += 1


			scaff_seq = get_scaff_sequence(mag_fasta, scaff)

			split_index = 1
			tmp = []
			for p in range(1, len(scaff_seq)+1):
				if scaffold_coverage[p] >= min_depth and p in scaff_range:
					tmp.append(p)
				else:
					if len(tmp) >= min_scaff_length:
						min_pos = min(tmp)
						max_pos = max(tmp)
						subseq = ""
						if max_pos == len(scaff_seq): subseq = scaff_seq[min_pos-1:]
						else: subseq = scaff_seq[min_pos-1:max_pos]
						scaff_id = select_scaffols_ids[scaff][0]
						if len(select_scaffolds_ids[scaff]) > 1:
							for ci, c in enumerate(select_scaffolds_coords[scaff]):
								cc = set(range(c[0], c[1]))
								if len(cc.intersection(set(tmp))) > 0:
									scaff_id = select_scaffolds_ids[scaff][ci]
						outf_assembly.write('>' + scaff_id + '_' + str(split_index) + '\n' + subseq + '\n')
						split_index += 1
					tmp = []

			if len(tmp) >= min_scaff_length:
				min_pos = min(tmp)
				max_pos = max(tmp)
				subseq = ""
				if max_pos == len(scaff_seq): subseq = scaff_seq[min_pos-1:]
				else: subseq = scaff_seq[min_pos-1:max_pos]
				scaff_id = select_scaffols_ids[scaff][0]
				if len(select_scaffolds_ids[scaff]) > 1:
					for ci, c in enumerate(select_scaffolds_coords[scaff]):
						cc = set(range(c[0], c[1]))
						if len(cc.intersection(set(tmp))) > 0:
							scaff_id = select_scaffolds_ids[scaff][ci]
				outf_assembly.write('>' + scaff_id + '_' + str(split_index) + '\n' + subseq + '\n')
				split_index += 1
			tmp = []
		outf_assembly.close()
		outf_r1.close()
		outf_r2.close()
		os.system('gzip %s %s' % (r1_fastq_file, r2_fastq_file))


def get_scaff_sequence(mag_fasta, scaffold):
	with open(mag_fasta) as omf:
		for rec in SeqIO.parse(omf, 'fasta'):
			if rec.id == scaffold:
				scaff_seq = str(rec.seq)
				return scaff_seq
	return None

def read_pair_generator(align_file, scaffold):
	"""
	Generate read pairs in a BAM file or within a region string.
	Reads are added to read_dict until a pair is found.
	"""
	read_dict = defaultdict(lambda: [None, None])
	for read in align_file.fetch(region=scaffold):
		if not read.is_proper_pair or read.is_secondary or read.is_supplementary: continue
		qname = read.query_name
		if qname not in read_dict:
			if read.is_read1: read_dict[qname][0] = read
			else: read_dict[qname][1] = read
		else:
			if read.is_read1: yield read, read_dict[qname][1]
			else: yield read_dict[qname][0], read
			del read_dict[qname]

main()
