import sys
from Bio import SeqIO
from collections import defaultdict

genome_file = sys.argv[1]
te = sys.argv[2]
input_sam = sys.argv[3]
output_sam = open(sys.argv[4], 'w')

genome_id = defaultdict(int)

for rec in SeqIO.parse(genome_file, 'fasta'):
	genome_id[rec.id] = 1


infor_rid = defaultdict(int)

for line in open(input_sam + '.group'):
	lst = line.strip().split('\t')
	rid = lst[0]
	chrom = lst[1].split(',')
	rnext = lst[2].split(',')
	gboo = 0
	tboo = 0
	for i in range(len(chrom)):
		if (chrom[i] in genome_id) or (rnext[i] in genome_id):
			gboo = 1 
		if (chrom[i] == te) or (rnext[i] == te):
			tboo = 1
	if gboo == tboo == 1:
		infor_rid[rid] = 1

for line in open(input_sam):
	if line.startswith('@'):
		output_sam.write(line)
	else:
		lst = line.strip().split('\t')
		if lst[0] in infor_rid:
			output_sam.write(line)
