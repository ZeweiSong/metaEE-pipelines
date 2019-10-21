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

configfile: 'config_16s.yaml'

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

# Target of this pipeline: put all samples in one profile, either in tsv or biom
rule target:
	input:
		workpath + 'combined.biom',
		workpath + 'combined.taxa.biom',
		workpath + 'combined.taxa.tsv',
		workpath + 'count.txt'

# We merge R1 and R2 together for a longer asselmby. Usually, variable region in 16S is short enough for effective merging. If it is too long, use the ITS snake instead.
rule pear:
	input:
		r1 = lambda wildcards: samplepath + SAMPLES[wildcards.sample] + '.r1.fq.gz',
		r2 = lambda wildcards: samplepath + SAMPLES[wildcards.sample] + '.r2.fq.gz'
	output:
		merged	= workpath + 'pear/{sample}.assembled.fastq', 
		fw		= workpath + 'pear/{sample}.unassembled.forward.fastq',
		rv		= workpath + 'pear/{sample}.unassembled.reverse.fastq',
		discard	= workpath + 'pear/{sample}.discarded.fastq'
	params:
		sn = workpath + 'pear/{sample}'
	threads: 4
	log: workpath + 'logs/pear/{sample}.log'
	run:
		shell('pear -f {input.r1} -r {input.r2} -o {params.sn} -b 64 -k -j {threads} > {log}')

# Find the two directions in the PE data	
rule direction:
	input:
		merged	= workpath + 'pear/{sample}.assembled.fastq'
	output:
		r1 = workpath + 'split_reads/{sample}.r1.fq',
		r2 = workpath + 'split_reads/{sample}.r2.fq'
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
		shell('cutadapt {input.merged} -g {params.fw} -a rvrc -n 2 --discard-untrimmed -e 0.1 -m 100 --quality-base 64 -j {threads} -o {output.r1} >  {log}'.replace('rvrc', rvrc))
		shell('cutadapt {input.merged} -g {params.rv} -a fwrc -n 2 --discard-untrimmed -e 0.1 -m 100 --quality-base 64 -j {threads} -o {output.r2} >> {log}'.replace('fwrc', fwrc))

# Remove low quality reads, and reverse compliment R2 sequences
rule quality:
	input:
		r1 = workpath + 'split_reads/{sample}.r1.fq',
		r2 = workpath + 'split_reads/{sample}.r2.fq'
	output:
		r1 = workpath + 'qc_reads/{sample}.r1.fa',
		r2_temp = temp(workpath + 'qc_reads/{sample}.r2.temp.fa'),
		r2 = workpath + 'qc_reads/{sample}.r2.fa'
	log:
		workpath + 'logs/qc_reads/{sample}.log'
	params:
		ascii = config['ascii']
	threads: 1
	run:
		shell('vsearch --fastq_filter {input.r1} --fastaout {output.r1} 	 --fastq_maxee 1 --fastq_minlen 100 --fasta_width 0 --fastq_ascii {params.ascii} --threads {threads} 2>  {log}')
		shell('vsearch --fastq_filter {input.r2} --fastaout {output.r2_temp} --fastq_maxee 1 --fastq_minlen 100 --fasta_width 0 --fastq_ascii {params.ascii} --threads {threads} 2>> {log}')
		shell('seqtk seq -r -A {output.r2_temp} > {output.r2}')
	
# Align reads to the database
rule align:
	input:
		r1 = workpath + 'qc_reads/{sample}.r1.fa',
		r2 = workpath + 'qc_reads/{sample}.r2.fa'
	output:
		r1 = workpath + 'alignments/{sample}.r1.b6',
		r2 = workpath + 'alignments/{sample}.r2.b6'
	params:
		acx  = config['database']['acx'],
		edx  = config['database']['edx'],
		mode = config['align']['mode'],
		id   = config['align']['id']
	log:
		workpath + 'logs/alignments/{sample}.log'
	threads: 4
	run:
		shell('burst12 -q {input.r1} -a {params.acx} -r {params.edx} -o {output.r1} -i {params.id} -m {params.mode} -t {threads} >  {log}')
		shell('burst12 -q {input.r2} -a {params.acx} -r {params.edx} -o {output.r2} -i {params.id} -m {params.mode} -t {threads} >> {log}')

# Count reads retention in the above steps
rule count:
	input:
		direction_fw	=	workpath + 'split_reads/{sample}.r1.fq',
		direction_rv	=	workpath + 'split_reads/{sample}.r2.fq',
		quality_fw		=	workpath + 'qc_reads/{sample}.r1.fa',
		quality_rv		=	workpath + 'qc_reads/{sample}.r2.fa',
		align_fw		=	workpath + 'alignments/{sample}.r1.b6',
		align_rv		=	workpath + 'alignments/{sample}.r2.b6'
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
		r1 = workpath + 'alignments/{sample}.r1.b6',
		r2 = workpath + 'alignments/{sample}.r2.b6'
	output:
		aln		= workpath + 'align_concat/{sample}.b6',
		tsv		= workpath + 'profiles/{sample}.tsv',
		biom	= workpath + 'profiles/{sample}.biom'
	params:
		sampleName='{sample}'
	log:
		workpath + 'logs/winnerTakeAll/{sample}.log'
	run:
		shell('cat {input.r1} {input.r2} > {output.aln}')
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