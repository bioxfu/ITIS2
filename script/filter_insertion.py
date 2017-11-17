import sys
from collections import defaultdict


def parse_sr(p):
	paras = {'t': 0, 'CS': 0, 'CE': 0, 'cs': 0, 'ce': 0, 'TS': 0, 'TE': 0}
	for t in p.split('/'):
		if t == '':
			continue
		ts = t.split('=')
		k, v = ts[0], int(ts[1])
		paras[k] = v
	return(paras)


# postion list of te homo in ref genome
lst = sys.argv[1]
ins_file = sys.argv[2]

sr = sys.argv[3]
# /t=3/TS=1/TE=1/ , the minimum requried:
# t:total reads supporting insertion  /3/
# CS:clipped reads cover TE start site /0/
# CE:clipped reads cover TE end site  /0/
# cs:cross reads cover TE start  /0/
# ce:cross reads cover TE end    /0/
# TS:total reads cover TE start  /1/
# TE:total reads cover TE end    /1/
paras = parse_sr(sr)

# the treshhold of NB tag 
nb = int(sys.argv[4])

# the minimum required average mapping quality
MQ = int(sys.argv[5])

# the reads depth range
d = sys.argv[6]
DP_L, DP_H = [int(x) for x in d.split(',')]


filt_bed = open(sys.argv[7], 'w')


homos = defaultdict(lambda : defaultdict(dict))
with open(lst) as file:
	for line in file:
		lst = line.strip().split('\t')
		chrom, s, e, te = lst[0], int(lst[1]), int(lst[2]), lst[3]
		s = s - nb
		e = e + nb
		for i in range(s, e+1):
			homos[te][chrom][i] = 1

# print homos['mping']['Chr1:0-2000000']

with open(ins_file) as file:
	for line in file:
		boo = 1
		lst = line.strip().split('\t')
		chrom, s, e, t, rest = lst[0], int(lst[1]), int(lst[2]), lst[3], '\t'.join(lst[4:])
		tags = t.split(';')

		### parse the rest key and values ###
		other = {}
		for r in tags:
			tmp = r.split('=')
			k, v = tmp[0], tmp[1]
			other[k] = v
		#print other
		#print paras
		# filter support reads number
		if 'SR' in other:
			tmp = [int(x) for x in other['SR'].split(',')]
			cf, tot, r1, r2, r3, r4 = tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5]		
			if tot < paras['t'] or r1 < paras['CS'] or r2 < paras['CE'] or r3 < paras['cs'] or r4 < paras['ce'] or (r1+r3) < paras['TS'] or (r2+r4) < paras['TE']:
				boo = 0
		
		# filter eveage mapping valeu
		if 'MQ' in other and int(other['MQ']) < MQ:
			boo = 0

		# filter bg depth
		if 'DP' in other:
			dp = int(other['DP'])
			if dp < DP_L or dp > DP_H:
				boo = 0

		# mark ins near known site
		near = 0
		for i in range(s, e+1):
			if i in homos[other['NM']][chrom]:
				near = 1
		if near == 1:
			t = t+";NB=Y"
		else:
			t = t+";NB=N"

		if boo == 1:
			filt_bed.write('%s\t%s\t%s\t%s\t%s\n' % (chrom, s, e, t, rest))


