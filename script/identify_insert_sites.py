import sys, re
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

def read2dict(hit):
	cors = defaultdict(str)
	rid, flag, chrom, pos, mq, cig, nchr, npos, seq = hit[0], hit[1], hit[2], hit[3], hit[4], hit[5], hit[6], hit[7], hit[9]
	rc = -1 if flag.find('r') != -1 else 1
	cs = cigar(cig)
	rid = rid[:-2]
	cors['cig'] = cs
	cors['direc'] = rc
	cors['rid'] =rid
	cors['pos'] = int(pos)
	cors['seq'] = seq
	cors['chrom'] = chrom
	cors['mq'] = mq
	return cors

def cigar(cig):
	if cig.find('*') != -1:
		return 'Z'
	
	read_len = 0
	cs = ''	# simplified cigars
	match_len = 0

	for m in re.findall('(\d+)([MSIH])', cig):
		read_len += int(m[0])
		if m[1] == 'M':
			match_len += int(m[0]) 

	if float(match_len) / read_len > 0.95:
		cs = 'M'
	else:
		m = re.search('^(\d+)[SH]', cig)
		if m is not None:
			d = int(m.group(1))
			if d >= 5:
				cs = cs + 'S:' + str(d)

		m = re.search('(\d+)[SH]$', cig)
		if m is not None:
			d = int(m.group(1))
			if d >= 5:
				cs = cs + 'E:' + str(d)

	if cs == '':
		cs = 'Z'

	return cs

def cross(cors):
	len_tar = len(cors['tar']['seq'])
	len_te = len(cors[TE]['seq'])
	len_inter = int(ins_size / 2)

	if cors[TE]['direc'] == 1 and cors[TE]['pos'] > (TE_len - ins_size - lost):

		#   ---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------
		#                                                               s----->    <--------
		#
		#   ----------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------
		#                 -------->       <------s

		if cors['tar']['direc'] == 1:
			ins_direc = 'R'
			jun = cors['tar']['pos'] + len_tar + len_inter
		else:
			ins_direc = 'S'
			jun = cors['tar']['pos'] - len_inter

		te_jun = cors[TE]['pos'] + len_te
		out_lst.write("%s\t%s\t%s\t%s\tCE:%s\t%s\n" % (cors[TE]['rid'], ins_direc, cors['tar']['chrom'], jun, te_jun, cors['tar']['mq']))
		out_sup_sam.write(rds+'\n')
		return(1)

	elif cors[TE]['direc'] == -1 and cors[TE]['pos'] < (ins_size - len_te + lost)	:

		#   ---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>----------------------------	--
		#                   -------->       s<--------
		#
		#   ----------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------
		#                                                               ------>s	<-----

		len_tar = len(cors['tar']['seq'])
		if cors['tar']['direc'] == 1:
			ins_direc = 'S'
			jun = cors['tar']['pos'] + len_tar + len_inter
		else:
			ins_direc = 'R'
			jun = cors['tar']['pos'] - len_inter

		te_jun = cors[TE]['pos']
		out_lst.write("%s\t%s\t%s\t%s\tCS:%s\t%s\n" % (cors[TE]['rid'], ins_direc, cors['tar']['chrom'], jun, te_jun, cors['tar']['mq']))
		out_sup_sam.write(rds+'\n')
		return(1)

def mat(que, sub):
	que = que.upper()
	sub = sub.upper()
	q_l = len(que)
	s_l = len(sub)

	ref_diff = defaultdict(list)
	for i in range(s_l - q_l):
		tgt = sub[i:(i+q_l)]
		#print que
		#print tgt
		diffcount = compare_two_str(que, tgt)
		ref_diff[diffcount].append(i)

	min_diff = sorted(ref_diff.keys())[0]
	return [min_diff, ref_diff[min_diff]]


def compare_two_str(x, y):
	diff_num = 0
	for i in range(len(x)):
		if x[i] != y[i]:
			diff_num += 1
	return diff_num


def te_start(cors):
	m = re.search('S:(\d+)', cors[TE]['cig'])
	if m is not None and cors[TE]['direc'] == -1:

		# ---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>------------------------------
		#             ------->     <-------
	
		#
		# --------------------------- <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<------------------------------
		#                                                                 ------->   <-----

		l = int(m.group(1))
		que = Seq(cors[TE]['seq'].upper())[:l]
		que_r = que.reverse_complement()
		chrom_t = cors['tar']['chrom']
		pos_t = cors['tar']['pos']
		len_t = len(cors['tar']['seq'])
		jun_te = cors[TE]['pos']
		que = str(que)
		que_r = str(que_r)

		if cors['tar']['direc'] == 1:
			ins_direc = 'S'
			sub = str(genome[chrom_t][(pos_t-1):(pos_t-1+ins_size)].seq)

			diff, juns = mat(que, sub)
			
			if float(diff) / l < 0.05:
				jun = juns[0] + pos_t + l - 1
				mark = 'ts' if len(juns) > 1 else 'TS'
				out_lst.write("%s\t%s\t%s\t%s\t%s:%s\t%s\n" % (cors[TE]['rid'], ins_direc, chrom_t, jun, mark, jun_te, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)
		else:
			ins_direc = 'R'
			sub_start = pos_t - ins_size + len_t
			sub = genome[chrom_t][(sub_start-1):(sub_start-1+ins_size)]

			diff, juns = mat(que_r, sub)

			if float(diff) / l < 0.05:
				jun = juns[0] + sub_start
				mark = 'ts' if len(juns) > 1 else 'TS'
				out_lst.write("%s\t%s\t%s\t%s\t%s:%s\t%s\n" % (cors[TE]['rid'], ins_direc, chrom_t, jun, mark, jun_te, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)


def te_end(cors):
	m = re.search('E:(\d+)', cors[TE]['cig'])
	if m is not None and cors[TE]['direc'] == 1:

		#---------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-------------------------------
		#                                                      -------->     <-------
		#
		#----------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-------------------------------
		#                ------>  <--------      

		l = int(m.group(1))
		end_pos = cors[TE]['pos'] + len(cors[TE]['seq']) - l - 1
		que = Seq(cors[TE]['seq'].upper())[-l:]
		que_r = que.reverse_complement()
		chrom_t = cors['tar']['chrom']
		pos_t = cors['tar']['pos']
		len_t = len(cors['tar']['seq'])
		que = str(que)
		que_r = str(que_r)

		if cors['tar']['direc'] == 1:
			ins_direc = 'R'
			sub = str(genome[chrom_t][(pos_t-1):(pos_t-1+ins_size)].seq)

			diff, juns = mat(que_r, sub)
			
			if float(diff) / l < 0.05:
				jun = juns[0] + pos_t + l - 1
				mark = 'te' if len(juns) > 1 else 'TE'
				out_lst.write("%s\t%s\t%s\t%s\t%s:%s\t%s\n" % (cors[TE]['rid'], ins_direc, chrom_t, jun, mark, end_pos, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)
		else:
			ins_direc = 'S'
			sub_start = pos_t - ins_size + len_t
			sub = str(genome[chrom_t][(sub_start-1):(sub_start-1+ins_size)].seq)

			diff, juns = mat(que, sub)

			if float(diff) / l < 0.05:
				jun = juns[0] + sub_start
				mark = 'te' if len(juns) > 1 else 'TE'
				out_lst.write("%s\t%s\t%s\t%s\t%s:%s\t%s\n" % (cors[TE]['rid'], ins_direc, chrom_t, jun, mark, end_pos, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)

def match(lst):
	# return one hash with key if direction and value is start position of match on subject
	relas = defaultdict(dict)

	for j in range(2):
		for k in range(2,4):
			que = lst[j]
			sub = lst[k]
			diff, locs = mat(que, sub)
			ratio = float(diff) / len(que)
			direc = 1 if j == 0 else -1
			dir_pos = direc * k
			relas[ratio][dir_pos] = locs
			#if lst[4] == 'SRR631734.1559297':
				#print que
				#print sub
				#print ratio		
	# exact the most similarity region and must be small than 0.05
	
	for x in sorted(relas.keys()):
		if x <= 0.05:
			return(relas[x])
		break


def ge_start(cors):
	m = re.search('S:(\d+)', cors['tar']['cig'])
	if m is not None:

		l = int(m.group(1))
		que = Seq(cors['tar']['seq'].upper())[:l]
		ry = Seq(cors['tar']['seq'].upper())[l:]
		que_r = que.reverse_complement()
		ry_r = ry.reverse_complement()
		
		sub_h = "NNNNNNNNNN" + str(genome[TE][:(l + lost + 10)].seq)
		sub_t = str(genome[TE][-(l + lost + 10):].seq) + "NNNNNNNNNN"

		que, que_r, sub_h, sub_t = str(que), str(que_r), str(sub_h), str(sub_t) 
		ma_dict = match([que, que_r, sub_h, sub_t, cors['tar']['rid']])

#		if cors['tar']['rid'] == 'SRR631734.1559297':
#			print que
#			print que_r
#			print sub_h
#			print sub_t

		if ma_dict is not None:
			#print ma_dict
			if 3 in ma_dict and cors[TE]['direc'] == 1 and cors[TE]['pos'] > (TE_len - ins_size - lost - 10):

		#  --------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-------------------------------------------------
		#                                                         -------->        <---------
		#
				adj = ma_dict[3][-1]
				seq_te = sub_t[(adj+l):]
				r = 0
				for i in range(len(seq_te)):
					s = seq_te[i:(i+1)]
					g = ry[i:(i+1)]
					if s == g:
						r = i + 1
					else:
						break

				rr = len(seq_te) - 10 - r
				ins_direc = 'S'
				jun = cors['tar']['pos'] + r
				te_jun = TE_len - rr
				out_lst.write("%s\t%s\t%s\t%s\tGE:%s\t%s\n" % (cors['tar']['rid'], ins_direc, cors['tar']['chrom'], jun, te_jun, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)
			elif -2 in ma_dict and cors[TE]['direc'] == -1 and cors[TE]['pos'] < (ins_size - len(cors[TE]['seq']) + lost + 10):

		#
		#  ---------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-------------------------------------------------
		#                                                         -------->        <---------  
				adj = ma_dict[-2][0]
				seq_te = sub_h[:adj]
				r = 0
				for i in range(len(seq_te)):
					s = seq_te[-(i+1):]
					g = ry_r[-(i+1):]
					if s == g:
						r = i + 1
					else:
						break

				rr = len(seq_te) - 10 - r
				ins_direc = 'R'
				jun = cors['tar']['pos'] + r
				te_jun = rr + 1
#				if cors['tar']['rid'] == 'SRR631734.868315':
#					print adj
#					print seq_te
#					print ry_r
#					print r
#					print rr
#					print jun
#					print te_jun

				out_lst.write("%s\t%s\t%s\t%s\tGS:%s\t%s\n" % (cors['tar']['rid'], ins_direc, cors['tar']['chrom'], jun, te_jun, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)

def ge_end(cors):
	m = re.search('E:(\d+)', cors['tar']['cig'])
	if m is not None:

		l = int(m.group(1))
		tar_l = len(cors['tar']['seq'])

		que = Seq(cors['tar']['seq'].upper())[-l:]
		ry = Seq(cors['tar']['seq'].upper())[:(tar_l-l)]
		que_r = que.reverse_complement()
		ry_r = ry.reverse_complement()
		
		sub_h = "NNNNNNNNNN" + str(genome[TE][:(l + lost + 10)].seq)
		sub_t = str(genome[TE][-(l + lost + 10):].seq) + "NNNNNNNNNN"

		que, que_r, sub_h, sub_t = str(que), str(que_r), str(sub_h), str(sub_t) 
		ma_dict = match([que, que_r, sub_h, sub_t, cors['tar']['rid']])

#		if cors['tar']['rid'] == 'SRR631734.1559297':
#			print que
#			print que_r
#			print sub_h
#			print sub_t
#			print ma_dict

		if ma_dict is not None:

			if -3 in ma_dict and cors[TE]['direc'] == 1 and cors[TE]['pos'] > (TE_len - ins_size - lost - 10):

		#
		#  ---------------------------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<-------------------------------------------------
		#                                     -------->        <----------  
				adj = ma_dict[-3][-1]
				seq_te = sub_t[(adj+l):]
				r = 0
				for i in range(len(seq_te)):
					s = seq_te[i:(i+1)]
					g = ry_r[i:(i+1)]
					if s == g:
						r = i + 1
					else:
						break

				rr = len(seq_te) - 10 - r
				ins_direc = 'R'
				jun = cors['tar']['pos'] + (len(cors['tar']['seq']) - l) - r - 1
				te_jun = TE_len - rr
#				if cors['tar']['rid'] == 'SRR631734.8997453':
#					print adj
#					print seq_te
#					print ry_r
#					print r
#					print rr
#					print jun
#					print te_jun
				out_lst.write("%s\t%s\t%s\t%s\tGE:%s\t%s\n" % (cors['tar']['rid'], ins_direc, cors['tar']['chrom'], jun, te_jun, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)
			elif 2 in ma_dict and cors[TE]['direc'] == -1 and cors[TE]['pos'] < (ins_size - len(cors[TE]['seq']) + lost + 10):

		#  --------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-------------------------------------------------
		#                                    --------->       <---------
		#
				adj = ma_dict[2][0]
				seq_te = sub_h[:adj]
				r = 0
				for i in range(len(seq_te)):
					s = seq_te[-(i+1):]
					g = ry[-(i+1):]
					if s == g:
						r = i + 1
					else:
						break

				rr = len(seq_te) - 10 - r
				ins_direc = 'S'
				jun = cors['tar']['pos'] + (len(cors['tar']['seq']) - l) - r - 1
				te_jun = rr + 1

				out_lst.write("%s\t%s\t%s\t%s\tGS:%s\t%s\n" % (cors['tar']['rid'], ins_direc, cors['tar']['chrom'], jun, te_jun, cors['tar']['mq']))
				out_sup_sam.write(rds+'\n')
				return(1)


#genome_file = 'genome/ref_and_mping.fa'
#sam_file = 'sam/mping.informative.full.sam'
#TE_realign_sam = 'sam/mping.alnte.x.sam'
#TE = 'mping'
#lost = 10
#ins_size = 300

genome_file, sam_file, TE_realign_sam, TE, lost, ins_size, lst_file, sam_out = sys.argv[1:]
lost, ins_size = int(lost), int(ins_size)

genome = SeqIO.to_dict(SeqIO.parse(genome_file, 'fasta'))
TE_len = len(genome[TE].seq)

out_lst = open(lst_file, 'w')
out_sup_sam = open(sam_out, 'w')

## te_aln(TE_realign_sam)
guanxi = defaultdict(list)
te_rcd = {}

with open(TE_realign_sam) as file:
	for line in file:
		lst = line.strip().split('\t')
		rid, flag, pos, cig, seq = lst[0], lst[1], int(lst[3]), lst[5], lst[9]
		cig = cig.replace('H', 'S')

		direc = -1 if flag.find('r') != -1 else 1

		## firstly, save the full infor to $l_seq and $l_as
		if rid not in te_rcd:
			guanxi[rid] = []
			if direc == -1:
				l_seq = str(Seq(seq).reverse_complement())
			else:
				l_seq = seq
			te_rcd[rid] = 1

		if direc == -1:
			seq = str(Seq(l_seq).reverse_complement())
		else:
			seq = l_seq

		lst[5] = cig
		lst[9] = seq

		m1 = re.match('^\d+M$', cig)
		m2 = re.match('^\d+M(\d+)S$', cig)
		m3 = re.match('^\d+S\d+M$', cig)

		if m1 is not None:
			guanxi[rid].append(lst)
		elif m2 is not None:
			ns = int(m2.group(1))
			if (TE_len - (len(seq) - ns) + 1 - lost -10) <= pos:
				guanxi[rid].append(lst)
		elif m3 is not None:
			if pos <= (lost + 10):
				guanxi[rid].append(lst)

## scan_sam($sam_file)
aligns = []

with open(sam_file) as file:
	for line in file:
		if line.startswith('@'):
			out_sup_sam.write(line)
		else:
			lst = line.strip().split('\t')
			rid, flag, chrom = lst[0], lst[1], lst[2]
			if chrom != TE:
				r = '1' if flag.find('1') != -1 else '2'
				lst[0] = rid+':'+r
				r_a_t = '1' if flag.find('2') != -1 else '2'
				r_a_t = rid+':'+r_a_t

				if r_a_t in guanxi:
					for te_aln in guanxi[r_a_t]:
						aligns.append([te_aln, lst])

for hits in aligns:
	te_dict = read2dict(hits[0])
	chrom_dict = read2dict(hits[1])
	rds = '\n'.join(['\t'.join(x) for x in hits])
	cors = {}
	cors['tar'] = chrom_dict
	cors[TE] = te_dict
	tar_cig = cors['tar']['cig']
	te_cig = cors[TE]['cig']
	boo_ns = 0

	if tar_cig.find('M') != -1 and te_cig.find('M') != -1:
		if cross(cors) == 1:
			boo_ns = 1
	if te_cig.find('S') != -1:
		if te_start(cors) == 1:
			boo_ns = 1
	if te_cig.find('E') != -1:
		if te_end(cors) == 1:
			boo_ns = 1
	if tar_cig.find('S') != -1:
		if ge_start(cors) == 1:
			boo_ns = 1
	if tar_cig.find('E') != -1:
		if ge_end(cors) == 1:
			boo_ns = 1

## 
out_lst.close()
out_sup_sam.close()


