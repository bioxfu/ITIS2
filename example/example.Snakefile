configfile: "config.yaml"

rule all:
	input:
		expand('genome/ref_and_{TE}.fa', TE=config['TE']),
		expand('blast/blastn_{TE}_ref_out_TE_position.tsv', TE=config['TE']),
		expand('sam/{sample}_all_reads_aln_ref_and_{TE}.sort.bam', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_reads_aligned_to_{TE}.sam', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_ref_and_{TE}.sam', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_{TE}.informative.sam', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_{TE}.informative.sorted.bam', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_{TE}.informative.full.sam', TE=config['TE'], sample=config['samples']),
		expand('fastq/{sample}_{TE}.informative.full.TE.fastq', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_{TE}.alnte.sam', TE=config['TE'], sample=config['samples']),
		expand('table/{sample}_{TE}.ins.loc.lst.sorted', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_{TE}.support.reads.sam', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_{TE}.support.reads.sorted.bam', TE=config['TE'], sample=config['samples']),
		expand('table/{sample}_{TE}.raw.bed', TE=config['TE'], sample=config['samples']),
		expand('sam/{sample}_{TE}.bg.sam', TE=config['TE'], sample=config['samples']),
		expand('table/{sample}_{TE}.filtered.bed', TE=config['TE'], sample=config['samples'])

rule blastn:
	input:
		TE = 'TE/'+config['TE']+'.fa',
		genome = 'genome/'+config['genome']
	output:
		'blast/blastn_'+config['TE']+'_ref_out'
	shell:
		'blastn -query {input.TE} -subject {input.genome} -outfmt 6 -out {output}'

rule mask_te_homo_in_genome:
	input:
		TE = 'TE/'+config['TE']+'.fa',
		genome = 'genome/'+config['genome'],
		blastn = 'blast/blastn_'+config['TE']+'_ref_out'
	output:
		fa = 'genome/ref_and_{TE}.fa',
		tsv = 'blast/blastn_{TE}_ref_out_TE_position.tsv'
	shell:
		'python script/mask_te_homo_in_genome.py {input.genome} {input.TE} {input.blastn} {output.fa} {output.tsv}'

rule bwa_index:
	input:
		genome = 'genome/ref_and_{TE}.fa',
		TE = 'TE/'+config['TE']+'.fa'
	output:	
		'genome/ref_and_{TE}.fa.bwt',
		'TE/{TE}.fa.bwt'
	shell:
		'bwa index {input.genome}; bwa index {input.TE}' 

rule align_ori_reads:
	input:
		R1 = config['path']+'{sample}_R1.fastq.gz',
		R2 = config['path']+'{sample}_R2.fastq.gz',
		genome = 'genome/ref_and_{TE}.fa',
		genome_idx = 'genome/ref_and_{TE}.fa.bwt'
	output:
		bam = 'sam/{sample}_all_reads_aln_ref_and_{TE}.sort.bam',
		prefix = 'sam/{sample}_all_reads_aln_ref_and_{TE}.sort'
	params:
		cpu = config['cpu']
	shell:
		'touch {output.prefix}; bwa mem -T 20 -t {params.cpu} {input.genome} {input.R1} {input.R2} 2>/dev/null | samtools view -@ 2 -buSh - | samtools sort -@ 2 - {output.prefix}; samtools index {output.bam}; '

rule reads_aligned_to_TE:
	input:
		R1 = config['path']+'{sample}_R1.fastq.gz',
		R2 = config['path']+'{sample}_R2.fastq.gz',
		TE = 'TE/'+config['TE']+'.fa',
		TE_idx = 'TE/'+config['TE']+'.fa.bwt'		
	output:
		sam = 'sam/{sample}_reads_aligned_to_{TE}.sam'
	params:
		cpu = config['cpu']		
	shell:
		'bwa mem -MT 20 -t {params.cpu} {input.TE} {input.R1} {input.R2} 2>/dev/null | samtools view -XS - > {output.sam}'

rule extract_reads_aligned_at_TE:
	input:
		'sam/{sample}_reads_aligned_to_{TE}.sam'
	output:
		R1 = 'fasta/{sample}_{TE}_R1.fa',
		R2 = 'fasta/{sample}_{TE}_R2.fa'
	shell:
		'python script/lean_fq.py {input} {output.R1} {output.R2}'


rule align_reads_associate_with_TE_to_merged_reference:
	input:
		R1 = 'fasta/{sample}_{TE}_R1.fa',
		R2 = 'fasta/{sample}_{TE}_R2.fa',
		genome = 'genome/ref_and_'+config['TE']+'.fa' 
	output:
		sam = 'sam/{sample}_ref_and_{TE}.sam',
		bam = 'sam/{sample}_ref_and_{TE}.sorted.bam'
	params:
		cpu = config['cpu']
	shell:
		'bwa mem -T 20 -v 1 -t {params.cpu} {input.genome} {input.R1} {input.R2} 2>/dev/null | tee {output.sam} | samtools view -@ 2 -buSh - | samtools sort -@ 2 -f - {output.bam}; samtools index {output.bam}'

rule extract_informative_reads:
	input:
		sam = 'sam/{sample}_ref_and_{TE}.sam',
		genome = 'genome/'+config['genome']
	output:
		sam = 'sam/{sample}_{TE}.informative.sam'
	params:
		TE = config['TE']
	shell:
		"grep -v '^@' {input.sam}|cut -f1,3,7|sort|groupBy -g 1 -c 2,3 -o collapse,collapse > {input.sam}.group; python script/extract_informative.py {input.genome} {params.TE} {input.sam} {output.sam}"	

rule sort_informative_reads:
	input:
		'sam/{sample}_{TE}.informative.sam'
	output:
		bam = 'sam/{sample}_{TE}.informative.sorted.bam'
	shell:
		'samtools view -buS {input} | samtools sort -f - {output.bam}; samtools index {output.bam}'	

rule modify_informative_reads:
	input:
		'sam/{sample}_{TE}.informative.sam'
	output:
		x = 'sam/{sample}_{TE}.informative_x.sam',
		full = 'sam/{sample}_{TE}.informative.full.sam',
		fastq = 'fastq/{sample}_{TE}.informative.full.TE.fastq'
	params:
		TE = config['TE']
	shell:
		'samtools view -h -S -X {input} > {output.x}; python script/modify_informative.py {output.x} {output.full} {output.fastq} {params.TE}'


rule realign_informative_reads:
	input:
		TE = 'TE/'+config['TE']+'.fa',
		fastq = 'fastq/{sample}_{TE}.informative.full.TE.fastq'
	output:
		sam = 'sam/{sample}_{TE}.alnte.sam',
		x = 'sam/{sample}_{TE}.alnte.x.sam'
	shell:
		'bwa mem -T 20 {input.TE} {input.fastq} -a 2>/dev/null > {output.sam}; samtools view -S -X {output.sam} > {output.x}'


rule identify_reads_support_insertion:
	input:
		genome = 'genome/ref_and_{TE}.fa',
		full = 'sam/{sample}_{TE}.informative.full.sam',
		x = 'sam/{sample}_{TE}.alnte.x.sam'
	output:
		table = 'table/{sample}_{TE}.ins.loc.lst',
		sort_table = 'table/{sample}_{TE}.ins.loc.lst.sorted',
		sam = 'sam/{sample}_{TE}.support.reads.sam',
		bam = 'sam/{sample}_{TE}.support.reads.sorted.bam'
	params:
		TE = config['TE'],
		lib_l = config['library_insertion_length'],
		lost = config['bases_lost']
	shell:
		'python script/identify_insert_sites.py {input.genome} {input.full} {input.x} {params.TE} {params.lost} {params.lib_l} {output.table} {output.sam}; sort -k 3,3 -k 4,4n {output.table} > {output.sort_table}; samtools view -buS {output.sam} | samtools sort -f - {output.bam};samtools index {output.bam}'

rule generate_bed:
	input:
		bam = 'sam/{sample}_all_reads_aln_ref_and_{TE}.sort.bam',
		TE = 'TE/'+config['TE']+'.fa',
		table = 'table/{sample}_{TE}.ins.loc.lst.sorted'
	output:
		bed = 'table/{sample}_{TE}.raw.bed',
		sam = 'sam/{sample}_{TE}.bg.sam'
	params:
		TE = config['TE'],
		lib_l = config['library_insertion_length'],
		window = config['window_size']
	shell:
		'python script/transform_to_bed.py {input.bam} {input.TE} {params.TE} {input.table} {params.lib_l} {params.window} {output.bed} {output.sam}'

rule filter_bed:
	input:
		blast = 'blast/blastn_{TE}_ref_out_TE_position.tsv',
		table = 'table/{sample}_{TE}.raw.bed'
	output:
		'table/{sample}_{TE}.filtered.bed'
	params:
		filter_str = config['filter_string'],
		NB_cutoff = config['NB_cutoff'],
		MQ_cutoff = config['MQ_cutoff'],
		DP_cutoff = config['DP_cutoff']
	shell:
		'python script/filter_insertion.py {input.blast} {input.table} {params.filter_str} {params.NB_cutoff} {params.MQ_cutoff} {params.DP_cutoff} {output}'

