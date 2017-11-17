import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

genome_file = sys.argv[1]
te_file = sys.argv[2]
blast_file = sys.argv[3]
ref_file = sys.argv[4]
tsv_file = sys.argv[5]

genome = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))

te = [x for x in SeqIO.parse(te_file, 'fasta')]

out_list = open(tsv_file, 'w')

for line in open(blast_file):
	lst = line.strip().split('\t')
	te_id, chrom, start, end = lst[0], lst[1], int(lst[8]), int(lst[9])
	start, end = sorted([start, end])
	out_list.write('%s\t%s\t%s\t%s\n' % (chrom, start, end, te_id))

	te_len = end - start + 1
	chrom_id = genome[chrom].id
	chrom_seq = str(genome[chrom].seq)
	chrom_seq = chrom_seq.replace(chrom_seq[(start-1):end], 'N' * te_len)
	genome[chrom] = SeqRecord(Seq(chrom_seq), id = chrom_id, description =  chrom_id)

out_list.close()

# put masked genome seq and TE seq together
fastout = list(genome.values())
fastout.extend(te)
SeqIO.write(fastout, ref_file, 'fasta')

