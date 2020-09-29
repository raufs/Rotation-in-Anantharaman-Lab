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
	Program: FilterForReadsOfSpecificScaffolds.py
	Author: Rauf Salamzade
    UW Madison MDTP / Rotating in Anantharaman Lab
    This program gathers scaffolds (or specific coords on them) of interest in MAGs and reads aligning to them.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-m', '--mag', type=str, help="FASTA (single entry) for reference scaffold upon which to call windows on.", required=True)
	parser.add_argument('-a', '--alignment', type=str, help='Reflexive alignment files in BAM format.', required=True)
	parser.add_argument('-o', '--outdir', type=str, help='path to output directory.', required=True)
	parser.add_argument('-b', '--bed', type=str, help='BED file with three columns listing the scaffold, start pos, stop pos, and id.', required=False, default=None)
	parser.add_argument('-n', '--name', type=str, help='Sample name.', required=False, default="sample")
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
	sample_name = myargs.name
	outdir = os.path.abspath(myargs.outdir) + '/'

	select_scaffolds = set([])
	select_scaffolds_coords = defaultdict(list)
	select_scaffolds_ids = defaultdict(list)
	try:
		assert(os.path.isfile(mag_fasta) and os.path.isfile(alignment_file))
		all_scaffolds = set([])
		try:
			with open(mag_fasta) as omf:
				for rec in SeqIO.parse(omf, 'fasta'):
					seqlen = len(str(rec.seq))
					all_scaffolds.add(rec.id)
		except:
			raise RuntimeError("Error reading in Metagenomic Assembly as a FASTA file.")

		if bed_file != None:
			with open(bed_file) as obf:
				for line in obf:
					scaff, start, stop, id = line.strip().split('\t')
					try:
						assert(scaff in all_scaffolds)
						select_scaffolds.add(scaff)
						select_scaffolds_coords[scaff].append([int(start), int(stop)])
						select_scaffolds_ids[scaff].append(id)
					except:
						sys.stderr.write("Did not find scaffold %s in metagenomic assembly, ignoring it.\n" % scaff)
		else:
			select_scaffolds = all_scaffolds
	except:
		raise RuntimeError("Errors with reading in input files, perhaps paths are broken.")

	if os.path.isdir(outdir):
		sys.stderr.write("Output directory already exists. Overwriting in five seconds ...\n")
		#time.sleep(5)
	else:
		os.system('mkdir %s' % outdir)

	# Start program

	align_handle = pysam.AlignmentFile(alignment_file, 'rb')

	r1_fastq_file = outdir + sample_name + '_R1.fastq'
	r2_fastq_file = outdir + sample_name + '_R2.fastq'
	assembly_file = outdir + sample_name + '.fasta'
	outf_r1 = open(r1_fastq_file, 'w')
	outf_r2 = open(r2_fastq_file, 'w')
	outf_assembly = open(assembly_file, 'w')

	for scaff in select_scaffolds:
		scaff_range = set([])
		for coords in select_scaffolds_coords[scaff]:
			scaff_range = scaff_range.union(set(range(coords[0] , coords[1]+1)))
		for read1, read2 in read_pair_generator(align_handle, scaff):
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
			read1_pos_range = set(range(insertion_location_left_r1, insertion_location_right_r1+1))
			read2_pos_range = set(range(insertion_location_left_r2, insertion_location_right_r2+1))

			if len(scaff_range.intersection(read1_pos_range)) > 0 or len(scaff_range.intersection(read2_pos_range)) > 0:
				outf_r1.write('@%s 1:N:0:AAAAAA\n%s\n+\n%s\n' % (query_name, r1_sequence, ''.join([str(x) for x in r1_qual] )))
				outf_r2.write('@%s 2:N:0:AAAAAA\n%s\n+\n%s\n' % (query_name, r2_sequence, ''.join([str(x) for x in r2_qual] )))

	outf_r1.close()
	outf_r2.close()
	os.system('gzip %s %s' % (r1_fastq_file, r2_fastq_file))

	scaff_index = defaultdict(lambda: 1)
	with open(mag_fasta) as omf:
		for rec in SeqIO.parse(omf, 'fasta'):
			if rec.id in select_scaffolds:
				scaff_seq = str(rec.seq)
				for i, coord in enumerate(select_scaffolds_coords[rec.id]):
					start, stop = coord
					id = select_scaffolds_ids[rec.id][i]
					scaff_len = len(str(rec.seq))

					if start == 1 and stop == scaff_len:
						outf_assembly.write('>' + id + '\n' + scaff_seq + '\n')
					else:
						if stop == scaff_len:
							outf_assembly.write('>' + id + '\n' + scaff_seq[start-1:] + '\n')
						else:
							outf_assembly.write('>' + id + '\n' + scaff_seq[start-1:stop] + '\n')
	outf_assembly.close()

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
