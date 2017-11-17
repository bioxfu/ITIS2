import sys
from Bio.Seq import Seq
from collections import defaultdict

input_sam = sys.argv[1]
output_R1 = sys.argv[2]
output_R2 = sys.argv[3]

r1 = open(output_R1, 'w')
r2 = open(output_R2, 'w')

reads = defaultdict(dict)
rids = []

for line in open(input_sam):
	lst = line.strip().split('\t')
	rid, flag, seq, qual = lst[0], lst[1], lst[9], lst[10]

	# skip if flag contain uU or s
	if flag.find('uU') == -1 and flag.find('s') == -1:
		r = flag[-1]

		# if mapped on the minus strand
		if flag.find('r') != -1:
			seq = str(Seq(seq).reverse_complement())
		
		if rid not in reads:
			rids.append(rid)
		
		reads[rid][r] = '>%s\n%s\n' % (rid, seq) 

for rid in rids:
	r1.write(reads[rid]['1'])
	r2.write(reads[rid]['2'])

