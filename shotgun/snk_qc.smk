# 2019/09/17, 2020/8/18 ZS This is our almost-final pipeline for quality control of sequencing data.
# This pipeline is used generaly on shotgun data where DNA fragments were randomly drawed.
# QCed shotgun data can go down to multilpe analysis pathways. So it is better to
# make a stand alone QC pipeline.
# For amplicon sequencing, QC is often inclued with other pipelines.

# CHECK POINT: You need to first rename your raw data, so the file name looks like:
#				sample-name-XXX_1.fq.gz
#				sample-name-XXX_2.fq.gz
# Remember that "_" is used to identify the correct R1 and R2 file, as well as the sample name.
# You should not use "_" to name your samples.
# Put (or ln -s) all raw data into the path defined in the variable "filepath" below.

configfile: 'config.yaml'

import os
import sys
SAMPLES = {}
filepath = 'data/samples/'
for file in os.listdir(filepath):
		sample_name = file.split('_')[0]
		SAMPLES[sample_name] = SAMPLES.get(sample_name, []) + [filepath + file]
if len(SAMPLES) == 0:
	print('ERROR: no file under {0}.'.format(filepath))
	sys.exit('Pipeline terminated orz.')
else:
	for key, value in SAMPLES.items():
		if len(value) != 2:
			print('{0} does not have the right files.'.format(key))
			sys.exit('Pipeline terminated, please check.')
		else:
			SAMPLES[key] = tuple(sorted(value))

rule target:
	input:
		report='data/count_report.tsv'

# Attemp to merge R1 and R2 reads, reads not assembled is kept.
rule pear:
	input:
		lambda wildcards: SAMPLES[wildcards.sample] # a lambda expression that "sample" is used as wildcards to refer to the value of the dictionary SAMPLES
	output:
		assembled	=	temp('data/pear/{sample}.assembled.fastq'),
		fwd			=	temp('data/pear/{sample}.unassembled.forward.fastq'),
		rev			=	temp('data/pear/{sample}.unassembled.reverse.fastq'),
		discard		=	temp('data/pear/{sample}.discarded.fastq')
	params:
		name = '{sample}',
		phred = config['phred']
	threads: 4
	log: 'log/qc/pear/{sample}.pear.log'
	run:
		shell('pear -f {input[0]} -r {input[1]} -o data/pear/{params.name} -k -b {params.phred} -j {threads} > {log}')

# For unassembled R1 and R2 reads, remove potential adaptors from their tails.
# Convert the two files into one interleaved fastq file.
# Adaptors for DNBSEQ are used here.
rule cutadapt:
	input:
		fwd		= 	'data/pear/{sample}.unassembled.forward.fastq',
		rev		= 	'data/pear/{sample}.unassembled.reverse.fastq'
	output:
		fwd		=	temp('data/cutadapt/{sample}.cutadapt.forward.fq'),
		rev		=	temp('data/cutadapt/{sample}.cutadapt.reverse.fq'),
		inter	=	temp('data/temp/{sample}.interleaved.fq')
	params:
		fw = 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA',
		rv = 'AAGTCGGATCGTAGCCATGTCGTTCTGT',
		phred = config['phred']
	threads: 4
	log: 'log/qc/cutadapt/{sample}.cutadapt.log'
	run:
		
		shell('cutadapt {input.fwd} {input.rev} -a {params.fw} -A {params.rv} -o {output.fwd} -p {output.rev} -m 50 --quality-base {params.phred} -j {threads} > {log}')
		shell('seqtk mergepe {output.fwd} {output.rev} > {output.inter}')

# For assembled reads, discard if maxee > 1
# For unassembled read pair, truncate until maxee <= 1, keep only pair
rule maxee:
	input:
		assembled	='data/pear/{sample}.assembled.fastq',
		pair		='data/temp/{sample}.interleaved.fq'
	output:
		assembled	='data/maxee/{sample}.merged.fa',
		pair		='data/maxee/{sample}.pair.fa'
	params:
		phred = config['phred']
	threads: 4
	log: 'log/qc/maxee/{sample}.maxee.log'
	run:
		shell('vsearch --fastq_filter {input.assembled} --fastq_maxee 1 --fastq_minlen 100 --fastq_maxns 0 --fastaout {output.assembled} --fasta_width 0 --fastq_ascii {params.phred} --threads {threads} 2> {log}')
		shell('vsearch --fastq_filter {input.pair} --fastq_truncee 1 --fastq_minlen 100 --fastq_maxns 0 --fastaout - --fasta_width 0 --fastq_ascii {params.phred} --threads {threads} 2>> {log} | seqtk dropse > {output.pair}')

# Count reads for every step
rule report:
	input:
		raw			='data/samples/{sample}_1.fq.gz',
		merged		='data/pear/{sample}.assembled.fastq',
		merged_qc	='data/maxee/{sample}.merged.fa',
		paired		='data/pear/{sample}.unassembled.forward.fastq',
		paired_qc	='data/maxee/{sample}.pair.fa'
	output:
		counter		='data/count/{sample}.count'
	params:
		name='{sample}'
	run:
		 shell("c1=$(($(gunzip -c {input.raw} | wc -l | awk '{{print $1}}')/4)); \
				c2=$(($(wc -l {input.merged} | awk '{{print $1}}')/4)); \
				c3=$(($(wc -l {input.merged_qc} | awk '{{print $1}}')/2)); \
				c4=$(($(wc -l {input.paired} | awk '{{print $1}}')/4)); \
				c5=$(($(wc -l {input.paired_qc} | awk '{{print $1}}')/4)); \
				echo -e {params.name} '\t' $c1 '\t' $c2 '\t' $c3 '\t' $c4 '\t' $c5 > {output.counter}")

rule concat:
	input:
		counter		=	expand('data/count/{sample}.count', sample=SAMPLES)
	output:
		report='data/count_report.tsv'
	run:
		shell("echo -e sample'\t'raw'\t'merged'\t'merged_QC'\t'paired'\t'paired_QC > {output.report}")
		shell('cat {input.counter} >> {output.report}')