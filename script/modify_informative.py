import sys
from Bio.Seq import Seq
from collections import defaultdict

rseq = defaultdict(dict)
rTE = defaultdict(dict)
rTEid = []

input_sam = sys.argv[1]
output_sam = open(sys.argv[2], 'w')
output_fq = open(sys.argv[3], 'w')
TE = sys.argv[4]

with open(input_sam) as infile:
	for line in infile:
		if not line.startswith('@'):
			lst = line.split('\t')
			rid = lst[0]
			flag = lst[1]
			chrom = lst[2]
			cig = lst[5]
			seq = lst[9]
			r = flag[-1]

			if flag.find('r') != -1:
				seq = str(Seq(seq).reverse_complement())
			
			if cig.find('H') == -1:
				rseq[rid][r] = seq

			if chrom == TE:
				if rid not in rTE:
					rTEid.append(rid)
				rTE[rid][r] = 1


with open(input_sam) as infile:
	for line in infile:
		if line.startswith('@'):
			output_sam.write(line)
		else:
			lst = line.split('\t')
			rid = lst[0]
			flag = lst[1]
			r = flag[-1]
				
			cig = lst[5]
			seq = rseq[rid][r]

			if flag.find('r') != -1:
				seq = str(Seq(seq).reverse_complement())

			if cig.find('H') != -1:
				cig = cig.replace('H','S')

			lst[5] = cig
			lst[9] = seq

			output_sam.write('\t'.join(lst))

output_sam.close()

for rid in rTEid:
	r1 = rseq[rid]['1']
	r2 = rseq[rid]['2']
	q1 = 'J' * len(r1)
	q2 = 'J' * len(r2)

	if '1' in rTE[rid]:
		output_fq.write('@%s:1\n%s\n+\n%s\n' % (rid, r1, q1))
	if '2' in rTE[rid]:
		output_fq.write('@%s:2\n%s\n+\n%s\n' % (rid, r2, q2))


