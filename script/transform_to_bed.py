import sys, re, subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, Counter

def deter_ord (aa, bb):
	if aa * bb == 0:     # if juction of one side is determined
		ar = abs(aa) if aa != 0 else abs(bb)
		br = abs(bb) if bb != 0 else abs(aa)
		return([ar, br])
	else:
		ar, br = sorted(abs(aa), abs(bb))
		return ([ar, br])

def mode(num):
	count_dict = dict(Counter(num).most_common())
	max_count = Counter(num).most_common(1)[0][1]
	modes = []
	for x in num:
		if count_dict[x] == max_count: 
			modes.append(x)
	mode = sorted(modes)[0]
	return(mode)

def median(num):
	num = sorted([float(x) for x in num])
	l = len(num)
	me = num[int(l/2)]
	return(me)

def rcd_it(it, chrom, pos, ty, d):
	clus.append(it)
	dirs[d] += 1
	rcder['chrom'] = chrom
	rcder['pos'] = pos
	rcder['ty'] = ty

def compare_two_str(x, y):
	diff_num = 0
	for i in range(len(x)):
		if x[i] != y[i]:
			diff_num += 1
	return diff_num

def check_te(que, sub):
	que = que.upper()
	sub = sub.upper()
	q_l = len(que)
	s_l = len(sub)
	
	for i in range(s_l-q_l):
		tgt = sub[i:(i+q_l)]
		diffcount = compare_two_str(que, tgt)
		if (float(diffcount) / q_l) <= 0.1:
			return(1)
	return(0)

def collect_infor(clu, dirs):
	# the aim of this subroutine is to collect the information of candidate sites
	# use a cluster of support reads to determine if it is a ture TE insertion

	## use ratio to determine the direction of insertion
	if 'S' in dirs and 'R' in dirs:
		if float(dirs['S']) / (dirs['R'] + dirs['S']) > 0.8:
			direc = '+'
		elif float(dirs['R']) / (dirs['R'] + dirs['S']) > 0.8:
			direc = '-'
		else:
			direc = '.'
	else:
		direc = list(dirs.keys())
		direc = '+' if direc[0] == 'S' else '-'
	# direction is now in variable $dir
	
	###		determine the insertion site
	sc = len(clu)	# the number of support reads
	sit_s = []   # exact start site
	sit_e = []   # exact end site
	rou_s = []   # uncertain start site
	rou_e = []	 # uncertain end site
	te_s = []   # the start site of junction at te
	te_e = []   # the end site of junction at te
	in_dict = defaultdict(str)   # used to store read id	
	map_q = 0
	#print clu
	for x in clu:
		lst = x.strip().split('\t')
		rid, d, chrom, pos, ty, mq = lst[0], lst[1], lst[2], int(lst[3]), lst[4], int(lst[5]) 
		m1 = re.search('(GS:|TS:)(-?\d+)', ty)
		m2 = re.search('(GE:|TE:)(-?\d+)', ty)
		if m1 is not None:
			#print m1.group(2)
			sit_s.append(pos)
			te_s.append(int(m1.group(2)))
		elif m2 is not None:
			#print m2.group(2)
			sit_e.append(pos)
			te_e.append(int(m2.group(2)))
		elif re.search('CS', ty) is not None or re.search('ts', ty):
			rou_s.append(pos)
		elif re.search('CE', ty) is not None or re.search('te', ty):
			rou_e.append(pos)
		in_dict[rid] += ty
		map_q += mq
	map_q = int(map_q / sc)
	num_fg = len(in_dict.keys())   # the total num of read pairs suppor insertion
	#print in_dict

	### collet exact insert pos
	
	# site in genome  cover TE start
	if len(sit_s) > 0:
		ss = mode(sit_s)
	elif len(rou_s) > 0:
		ss = median(rou_s)
	else:
		ss = 0

	#  site in genome cover TE end
	if len(sit_e) > 0:
		ee = mode(sit_e)
	elif len(rou_e) > 0:
		ee = median(rou_e)
	else:
		ee = 0

	# start site at TE
	if len(te_s) > 0:
		t_s = mode(te_s)
	else:
		t_s = 'NA'

	# end site at TE
	if len(te_e) > 0:
		t_e = mode(te_e)
	else:
		t_e = 'NA'

	#print [sit_s, sit_e, rou_s, rou_e, te_s, te_e]
	#print [ss, ee, t_s, t_e]

	##################
	clp_s, clp_e, crs_s, crs_e = [len(x) for x in [sit_s, sit_e, rou_s, rou_e]]
	SR = '%s,%s,%s,%s,%s,%s' % (num_fg, sc, clp_s, clp_e, crs_s, crs_e) #  in the order of 'Reads support Start and End'. Fragment suported Start and End
	#print SR
	#print in_dict
	#print '============'

	if direc == '+' and len(sit_s) > 0 and len(sit_e) > 0:
		tsd_l = ss - ee
	elif direc == '-' and len(sit_s) > 0 and len(sit_e) > 0:
		tsd_l = ee - ss
	else:
		tsd_l = -999

	return([tsd_l,[chrom,ss,ee,in_dict,sc,clp_s,clp_e,crs_s,crs_e,"SR=%s;MQ=%s;NM=%s"%(SR,map_q,te),t_s,t_e,direc]])

	# SR : counts of total and every type of supporting reads
	# MQ : the average mapping quality
	# NM : the name of te 
	# TS : the start site of te 
	# TE : the end site of te


def estimate_homo(args):        # check each read pair around the candidate insert sites
	chrom, s_r, e_r, in_dict, direc, t_s, t_e  = args
	# s_r 1 based start site of insertion
	# e_r 1 based end site of insertion
	if t_s == "NA":
		t_s = 0
	if t_e == "NA":
		t_e = len(te_dict[te]) 

	sam_s = s_r - lib_l
	sam_e = e_r + lib_l

	if sam_s < 0:
		sam_s = 0 	
	
	###iterate each read pair
	reads = {}
	cmd = "samtools view -h -X %s %s:%s-%s" % (bam, chrom, int(sam_s), int(sam_e))
	cmd_output = subprocess.check_output(cmd, shell=True, universal_newlines=True)
	#print(cmd_output)
	#print(len(cmd_output.strip().split('\n')))
	for line in cmd_output.strip().split('\n'):
		line = str(line)
		if line.startswith('@'):
			continue
		lst = line.strip().split('\t')
		#print(lst)
		rid, tag, chrom, pos, mq, cig, nchr, npos, tlen, seq  = lst[0], lst[1], lst[2], lst[3], lst[4], lst[5], lst[6], lst[7], lst[8], lst[9] 
		pos = int(pos)
		npos = int(npos)
		tlen = int(tlen)

		if re.search('u', tag) is not None:
			continue

		reads[rid] = 9

		#  1 :  support insertion
		#  2 :  support excision
		#  3 :  flank reads 
		#  4 :  discarded reads
		m1 = re.search('(\d+)S$', cig)
		m2 = re.search('^(\d+)S', cig)

		if rid in in_dict:
			reads[rid] = 1
		elif tlen < 0 and abs(tlen) < 2*lib_l and m1 is not None and int(m1.group(1)) >= 20:
			#### aligned reads with soft clipped ends to test if it support insertion ###
			l = int(m1.group(1))
			que = seq[-l:]
			if direc == "+":
				sub = te_dict[te][int(t_s-1):int(t_s-1+l+5)]
			else:
				sub = str(te_dict[te][int(t_e-l-5):int(t_e-l-5+l+5)].seq.reverse_complement())
			
			if check_te(que, sub) == 1:
				reads[rid] = 1
			
		elif tlen > 0 and abs(tlen) < 2*lib_l and m2 is not None and int(m2.group(1)) >= 20:   
			# at least 20 bp soft clipped to check if it is from TE
			####aligned reads with sfor clipped ends to test if it support insertion###
			l = int(m2.group(1))
			que = seq[:l]
			if direc == "+":
				sub = te_dict[te][int(t_e-l-5):int(t_e-l-5+l+5)]
			else:
				sub = str(te_dict[te][int(t_s-1):int(t_s-1+l+5)].seq.reverse_complement())

			if check_te(que, sub) == 1:
				reads[rid] = 1
		else:
			if nchr == "=" and tlen != 0 and abs(tlen) < 2*lib_l:		
				#  determine the range of pair
				
				if tlen > 0:
					ranges = [pos, pos+tlen-1]
				elif tlen < 0:
					ranges = [npos, npos-tlen-1]

				# check if the range overlap with TE insertion site
				if s_r >= ranges[0] + 20 and ranges[1] >= e_r + 20  and reads[rid] >= 2:
					reads[rid] = 2
					out_sam.write(line+'\n')
	sup = 0
	nsup = 0

	for k, v in reads.items():
		if v == 1:
			sup += 1
		elif v == 2:
			nsup += 1

	return([sup, nsup])


#bam = 'sam/all_reads_aln_ref_and_mping.sort.bam'
#te_file = 'TE/mping.fa'
#te = 'mping'
#ins_file = 'table/mping.ins.loc.sorted.lst'
#lib_l = 300
#window = 150

bam = sys.argv[1]
te_file = sys.argv[2]
te = sys.argv[3]
ins_file = sys.argv[4]
lib_l = int(sys.argv[5])
window = int(sys.argv[6])

out_bed = open(sys.argv[7], 'w')
out_sam = open(sys.argv[8], 'w')

cmd = "samtools view -H %s" % bam
cmd_output = str(subprocess.check_output(cmd, shell=True))
out_sam.write(cmd_output)

te_dict = SeqIO.to_dict(SeqIO.parse(te_file, 'fasta'))

lsts = []
tsds = []

rcder = defaultdict(str)
clus = []
dirs = defaultdict(int)

with open(ins_file) as file:
	for line in file:  # iterate the lst file and get tsd information and candidate insertions 
		#print line
		lst = line.strip().split('\t')
		rid, direc, chrom, pos, ty = lst[0], lst[1], lst[2], int(lst[3]), lst[4] 
		win_s = 0

		if not any(rcder):
			rcd_it(line, chrom, pos, ty, direc)
			continue

		# determine the window step to cluster reads
		#print rcder
		if re.search('C', ty) is not None and re.search('C', rcder['ty']) is not None:
			win_s = window
		elif re.search('G|T', ty) is not None and re.search('G|T', rcder['ty']) is not None:
			win_s = 50
		else:
			win_s = window
		# get a cluter of support reads
		if chrom == rcder['chrom'] and (pos - rcder['pos'] <= win_s):
			rcd_it(line, chrom, pos, ty, direc)
		else:
			#print win_s
			#print ty + '\t' + rcder[ty]
			#print clus
			tsd, inf = collect_infor(clus, dirs)
			if tsd != -999:
				tsds.append(tsd)
			lsts.append(inf)

			rcder = defaultdict(str)
			clus = []
			dirs = defaultdict(int)
			rcd_it(line, chrom, pos, ty, direc)

# EOF ins_file
tsd, inf = collect_infor(clus, dirs)
if tsd != -999:
	tsds.append(tsd)
lsts.append(inf)

#print tsds

######### the final step to print  each candidate#######
#  1, used esimated tsd length to determine the exact insertion site
#
#  2,calculate the bg depth and estimate Genome Type if BAM file is used  #
##########################

# get the most of tsd length
tsd_l = 0
if len(tsds) > 0:
	tsd_l = mode(tsds) 
	tsd_l += 1
	print("Estimate the length of TSD is %s bp" % int(tsd_l))
	tsd_l -= 1
else:
	tsd_l = 1
	print("Can't determine the length TSD, because there are no clipped aligned reads")

#print len(lsts)
#sys.exit()

# check each  candidate insertion
idx = 0
for it in lsts:
	idx += 1
	chrom, ss, ee, in_dict, sc, clp_s, clp_e, crs_s, crs_e, tags, t_s, t_e, direct = it
	if direct == ".":
		print( "Can't determine the direction of insertion at %s\t%s\t%s\n" % (chrom, ss, sc))
		continue     # if a cluster of reads have different diretion, discarded

	sign = 1 if direct == "+" else -1

	### esitimate the exact insertion site  ####
	if clp_s > 0 and clp_e > 0 and sign*(ss-ee) == tsd_l:
		s_p, e_p = ss, ee
	elif clp_s >= 1:
		s_p, e_p = ss, ss-sign*tsd_l
	elif clp_e >= 1: 
		s_p, e_p = ee+sign*tsd_l, ee
	else:
		s_p, e_p = deter_ord(ss, ee)
	
	if sign == 1:
		e_p -= 1
	elif sign == -1:
		s_p -=1
	
	s_p, e_p = sorted([s_p, e_p])

#	if idx == 15:
#		print it
#		print clp_s
#		print clp_e
#		print s_p
#		print e_p
#		print ss
#		print ee
#		print sign
#		print tsd_l

	##  calculate bg depth and esimate homo or hetero
	###### esitmate ratio #####
	sup, nsup, noise = "NA","NA","NA"
	if clp_s > 0 or clp_e > 0:
		sup, nsup = estimate_homo([chrom, s_p, e_p, in_dict, direc, t_s, t_e]) 

	if str(sup) == 'NA' or str(nsup) == 'NA':
		pva = "NA"
	else:
		half = int(sup / 2  + 1)
		cmd = 'Rscript script/genotype_caculator.r %s %s' % (half, nsup)
		pva = float(subprocess.check_output(cmd, shell=True))
		
	if pva == "NA":
		gt = "NA"
	elif pva >= 0.01:
		gt = "Heter"
	elif pva < 0.01:
		gt = "Homo"
	tags += ";GT=%s,%s:%s;PV=%s" % (sup, nsup, gt, pva) 
	######## pick the bg depth ######
		
	pad = int((100- (e_p-s_p))/2)
	s_r, e_r = (s_p - pad, e_p + pad) if pad > 0 else (s_p, e_p)
		
	if s_r < 0:
		s_r = 0 

	cmd = "samtools depth %s -r %s:%s-%s 2>/dev/null" % (bam, chrom, int(s_r), int(e_r))
	cmd_output = subprocess.check_output(cmd, shell=True, universal_newlines=True)
	dep = 0
	m_l = 0

	for line in cmd_output.strip().split('\n'):
		line = str(line)
		lst = line.strip().split('\t')
		rid, pos, d = lst[0], lst[1], int(lst[2])
		dep += d
		m_l += 1

	if m_l > 0 and dep > 0:
		dep = int(dep/m_l)
	else:
		dep = 0
	tags += ";DP=%s" % dep

	out_bed.write("%s\t%s\t%s\t%s;TS=%s;TE=%s\t.\t%s\n" % (chrom, int(s_p), int(e_p), tags, t_s, t_e, direct))

out_bed.close()
out_sam.close()