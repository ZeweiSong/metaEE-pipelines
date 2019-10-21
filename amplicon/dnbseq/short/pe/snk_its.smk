#############################################
# Amplicon sequencing processing pipeline   #
# Used for Pair End reads                   #
# Burst for closed-reference clustering     #
# Winner take all for combining two targets #
#############################################

# This is a Snakemake pipeline for processing amplicon sequencing data based on 
# BGISEQ pair-end method. The insert library is build without controling the sequencing
# direction. So Read1 can start from either the forward primer or reverse primer, 
# randomly and in equal number.

# If you are not familiar with amplicon sequencing, or you have only done amplicon
# sequencing on illumina's platform, I would suggest you to read the following
# instructions carefully (at least once). Most of the illumina based protocol controls
# the direction of amplicons, so Read1 always occurs from the forward primers.
# This is NOT the case for our BGISEQ based protocol, and will result in some differences
# in data processing (and some intersting aftermath).

# Two types of amplicons appears in the PE data:
	# TYPE 1: Read1-forward primer-insert-reverse primer-Read2
	# TYPE 2: Read1-reverse primer-insert-forward primer-Read2
# For PE data, there is much less differnce than that in SE data, but we still distiguish
# these two types since they may different slightly on the reading position.
# For the commands in this pipeline, we refer read1 and read2 as r1 and r2,
# and refer the forward and reverse primers as fw and rv.
# Most of these labels appear as suffixes.
	# These are the two amplicons:
	# r1fw - r2rv
	# r1rv - r2fw

# You will need the sequence of the forward and reverse primers (or the conservative region
# you like to trim off). These primers are different based on the target amplified.
# Change the primer sequences in config.yaml to the one you used (and ask around if you don't know).
# Also check if you have the right Burst database (bacteria or fungi?).
# Check cluster.yaml if you use qsub.

# Get the sample label from its1/samples
# The sample file should be in the format [SampleName].[ReadDirection].fq.gz, i.e. sampleA.r1.fq.gz sampleA.r2.fq.gz
# You can change the file extension in the config file.
# Put all sample files in its1/samples
# A sample should have two corresponding files.

configfile: 'config_its.yaml'

import os
import sys
SAMPLES = []
workpath = config['workpath']
samplepath = workpath + 'samples/'
for file in os.listdir(samplepath): # List all files under the folder. For PE data, one samples correspond to two files.
	string = file.split('.')
	if len(string) != 4: # You have to makesure that "." is not used in naming your samples. All files have the suffix .r1.fq.gz or .r2.fq.gz
		print('{0} contains illegal naming style, please check.'.format(file))
		sys.exit()
	esle:
		SAMPLES.append(string[0])
SAMPLES = tuple(set(SAMPLES))
SAMPLES = {i:i for i in SAMPLES}

# Target of this pipeline: all samples in one profile, either in tsv or biom
rule target:
	input:
		workpath + 'combined.biom',
		workpath + 'combined.taxa.biom',
		workpath + 'combined.taxa.tsv',
		workpath + 'count.txt'
		
# Find the two directions in the PE data		
rule direction:
	input:
		r1 = lambda wildcards: samplepath + SAMPLES[wildcards.sample] + '.r1.fq.gz',
		r2 = lambda wildcards: samplepath + SAMPLES[wildcards.sample] + '.r2.fq.gz'
	output:
		r1fw = workpath + 'split_reads/{sample}.r1fw.fq',
		r2rv = workpath + 'split_reads/{sample}.r2rv.fq',
		r1rv = workpath + 'split_reads/{sample}.r1rv.fq',
		r2fw = workpath + 'split_reads/{sample}.r2fw.fq'
	params:
		fw=config['primers']['fw'],
		rv=config['primers']['rv']
	log:
		workpath + 'logs/split_reads/{sample}.log'
	threads: 2
	run:
		base = {'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','K':'M','M':'K','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N'}
		fwrc = ''.join([base[i] for i in list({params.fw})[0][::-1]])
		rvrc = ''.join([base[i] for i in list({params.rv})[0][::-1]])
		shell('cutadapt {input.r1} {input.r2} -g {params.fw} -a rvrc -G {params.rv} -A fwrc -n 2 --discard-untrimmed -e 0.1 -m 100 --quality-base 64 -j {threads} -o {output.r1fw} -p {output.r2rv} > {log}'.replace('fwrc', fwrc).replace('rvrc', rvrc))
		shell('cutadapt {input.r1} {input.r2} -g {params.rv} -a fwrc -G {params.fw} -A rvrc -n 2 --discard-untrimmed -e 0.1 -m 100 --quality-base 64 -j {threads} -o {output.r1rv} -p {output.r2fw} >> {log}'.replace('rvrc', rvrc).replace('fwrc', fwrc))

# Remove low quality reads		
rule quality:
	input:
		r1fw = workpath + 'split_reads/{sample}.r1fw.fq',
		r2rv = workpath + 'split_reads/{sample}.r2rv.fq',
		r1rv = workpath + 'split_reads/{sample}.r1rv.fq',
		r2fw = workpath + 'split_reads/{sample}.r2fw.fq'
	output:
		r1fw_qc = workpath + 'qc_reads/{sample}.r1fw.qc.fq',
		r2rv_qc = workpath + 'qc_reads/{sample}.r2rv.qc.fq',
		r1rv_qc = workpath + 'qc_reads/{sample}.r1rv.qc.fq',
		r2fw_qc = workpath + 'qc_reads/{sample}.r2fw.qc.fq'
	log:
		workpath + 'logs/qc_reads/{sample}.log'
	params:
		ascii = config['ascii']
	threads: 4
	run:
		shell('vsearch --fastq_filter {input.r1fw} --reverse {input.r2rv} --fastqout {output.r1fw_qc} --fastqout_rev {output.r2rv_qc} --fastq_maxee 1 --fastq_minlen 100 --fastq_ascii {params.ascii} --threads {threads} 2> {log}')
		shell('vsearch --fastq_filter {input.r1rv} --reverse {input.r2fw} --fastqout {output.r1rv_qc} --fastqout_rev {output.r2fw_qc} --fastq_maxee 1 --fastq_minlen 100 --fastq_ascii {params.ascii} --threads {threads} 2>> {log}')
		

# Trim off sequences not in database and convert to FASTA
rule trim_tail:
	input:
		r1fw_qc = workpath + 'qc_reads/{sample}.r1fw.qc.fq',
		r2rv_qc = workpath + 'qc_reads/{sample}.r2rv.qc.fq',
		r1rv_qc = workpath + 'qc_reads/{sample}.r1rv.qc.fq',
		r2fw_qc = workpath + 'qc_reads/{sample}.r2fw.qc.fq'
	output:
		r1fw_qc = workpath + 'trim_tail/{sample}.r1fw.qc.fa',
		r2rv_qc = workpath + 'trim_tail/{sample}.r2rv.qc.fa',
		r1rv_qc = workpath + 'trim_tail/{sample}.r1rv.qc.fa',
		r2fw_qc = workpath + 'trim_tail/{sample}.r2fw.qc.fa'
	params:
		fw_trim = config['trim']['fw'],
		rv_trim = config['trim']['rv']
	run:
		shell('seqtk trimfq -b {params.fw_trim} {input.r1fw_qc} - | seqtk seq -A -L 50 - > {output.r1fw_qc}')
		shell('seqtk trimfq -e {params.rv_trim} {input.r2rv_qc} - | seqtk seq -r -A -L 50 - > {output.r2rv_qc}')
		shell('seqtk trimfq -e {params.rv_trim} {input.r1rv_qc} - | seqtk seq -r -A -L 50 - > {output.r1rv_qc}')
		shell('seqtk trimfq -b {params.fw_trim} {input.r2fw_qc} - | seqtk seq -A -L 50 - > {output.r2fw_qc}')
		
# Align reads to the database
rule align:
	input:
		r1fw_qc = workpath + 'trim_tail/{sample}.r1fw.qc.fa',
		r2rv_qc = workpath + 'trim_tail/{sample}.r2rv.qc.fa',
		r1rv_qc = workpath + 'trim_tail/{sample}.r1rv.qc.fa',
		r2fw_qc = workpath + 'trim_tail/{sample}.r2fw.qc.fa'
	output:
		r1fw_aln = workpath + 'alignments/{sample}.r1fw.b6',
		r2rv_aln = workpath + 'alignments/{sample}.r2rv.b6',
		r1rv_aln = workpath + 'alignments/{sample}.r1rv.b6',
		r2fw_aln = workpath + 'alignments/{sample}.r2fw.b6'
	params:
		acx  = config['database']['acx'],
		edx  = config['database']['edx'],
		mode = config['align']['mode'],
		id   = config['align']['id']
	log:
		workpath + 'logs/alignments/{sample}.log'
	threads: 4
	run:
		shell('burst12 -q {input.r1fw_qc} -a {params.acx} -r {params.edx} -o {output.r1fw_aln} -i {params.id} -m {params.mode} -t {threads} > {log}')
		shell('burst12 -q {input.r2rv_qc} -a {params.acx} -r {params.edx} -o {output.r2rv_aln} -i {params.id} -m {params.mode} -t {threads} >> {log}')
		shell('burst12 -q {input.r1rv_qc} -a {params.acx} -r {params.edx} -o {output.r1rv_aln} -i {params.id} -m {params.mode} -t {threads} >> {log}')
		shell('burst12 -q {input.r2fw_qc} -a {params.acx} -r {params.edx} -o {output.r2fw_aln} -i {params.id} -m {params.mode} -t {threads} >> {log}')

# Keep the paired alignments
rule keep_paired:
	input:
		r1fw_aln = workpath + 'alignments/{sample}.r1fw.b6',
		r2rv_aln = workpath + 'alignments/{sample}.r2rv.b6',
		r1rv_aln = workpath + 'alignments/{sample}.r1rv.b6',
		r2fw_aln = workpath + 'alignments/{sample}.r2fw.b6'
	output:
		fw_aln = workpath + 'align_paired/{sample}.fw.b6',
		rv_aln = workpath + 'align_paired/{sample}.rv.b6'
	log:
		workpath + 'logs/keep_paired/{sample}.log'
	run:
		shell('amplicon_keepPairAln.py -i {input.r1fw_aln},{input.r2rv_aln} -o {output.fw_aln} > {log}')
		shell('amplicon_keepPairAln.py -i {input.r1rv_aln},{input.r2fw_aln} -o {output.rv_aln} >> {log}')

# Count reads retention in the above steps
rule count:
	input:
		direction_fw	=	workpath + 'split_reads/{sample}.r1fw.fq',
		direction_rv	=	workpath + 'split_reads/{sample}.r1rv.fq',
		quality_fw		=	workpath + 'qc_reads/{sample}.r1fw.qc.fq',
		quality_rv		=	workpath + 'qc_reads/{sample}.r1rv.qc.fq',
		align_fw		=	workpath + 'align_paired/{sample}.fw.b6',
		align_rv		=	workpath + 'align_paired/{sample}.rv.b6'
	output: workpath + 'count/{sample}.count'
	params:
		sn ='{sample}'
	run:
		shell('echo {params.sn} > {output}')
		shell('echo -e "forward" >> {output}')
		shell('amplicon_countRetention.py -d {input.direction_fw} -q {input.quality_fw} -a {input.align_fw} >> {output}')
		shell('echo -e "reverse" >> {output}')
		shell('amplicon_countRetention.py -d {input.direction_rv} -q {input.quality_rv} -a {input.align_rv} >> {output}')

# Combine alignments using winner take all method
rule winnerTakeAll:
	input:
		fw_aln = workpath + 'align_paired/{sample}.fw.b6',
		rv_aln = workpath + 'align_paired/{sample}.rv.b6'
	output:
		aln		= workpath + 'align_concat/{sample}.b6',
		tsv		= workpath + 'profiles/{sample}.tsv',
		biom	= workpath + 'profiles/{sample}.biom'
	params:
		sampleName='{sample}'
	log:
		workpath + 'logs/winnerTakeAll/{sample}.log'
	run:
		shell('cat {input.fw_aln} {input.rv_aln} > {output.aln}')
		shell('amplicon_winnerTakeAll.py -i {output.aln} -sn {params.sampleName} -t {output.tsv} -g > {log}')
		shell('biom convert -i {output.tsv} -o {output.biom} --to-json')		

# Concatenate all profiles into one
# Add the taxonomy, and write to biom and tsv file	
rule concat:
	input:
		biom  = expand(workpath + 'profiles/{sample}.biom', sample=SAMPLES),
		count = expand(workpath + 'count/{sample}.count', sample=SAMPLES)
	output:
		combined  = workpath + 'combined.biom',
		biom_taxa = workpath + 'combined.taxa.biom',
		tsv_taxa  = workpath + 'combined.taxa.tsv',
		count 	  = workpath + 'count.txt'
	params:
		taxa=config['database']['tax']
	log:
		workpath + 'logs/concat/concate.log'
	run:
		shell('amplicon_concat.py -i {input.biom} -biom_out {output.combined} > {log}')
		shell('biom add-metadata -i {output.combined} -o {output.biom_taxa} --observation-metadata-fp {params.taxa} --observation-header OTUID,taxonomy --output-as-json --sc-separated taxonomy')
		shell('biom convert -i {output.biom_taxa} -o {output.tsv_taxa} --to-tsv --header-key taxonomy')
		shell('cat {input.count} > {output.count}')