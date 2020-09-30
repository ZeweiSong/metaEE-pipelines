# 2019/09/17 ZS This is our semi-final pipeline for quality control of sequencing data.
# This pipeline is used generaly on shotgun data where DNA fragments were randomly drawed.
# QCed shotgun data can go down to multilpe analysis pathways. So it is better to
# make a stand alone QC pipeline.
# For amplicon sequencing, QC is often inclued with other pipelines.

# CHECK POINT: You need to first rename your raw data, so the file name looks like:
#				sample-name-XXX_1.fq.gz
#				sample-name-XXX_2.fq.gz
# Remember that "_" is used to identify the correct R1 and R2 file, as well as the sample name.
# Put (or ln -s) all raw data into data/sample

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

rule cutadapt:
	input:
		lambda wildcards: SAMPLES[wildcards.sample] # a lambda expression that "sample" is used as wildcards to refer to the value of the dictionary SAMPLES
	output:
		r1 = 'data/cutadapt/{sample}.cutadapt.r1.fq',
		r2 = 'data/cutadapt/{sample}.cutadapt.r2.fq'
	params:
		r1 = config['adaptors']['r1'],
		r2 = config['adaptors']['r2'],
		phred = config['phred']
	threads: 4
	log: 'logs/cutadapt/{sample}.cutadapt.log'
	run:
		shell('cutadapt {input[0]} {input[1]} -a {params.r1} -A {params.r2} -o {output.r1} -p {output.r2} -m 50 --quality-base {params.phred} -j {threads} > {log}')
		
rule truncee:
	input:
		r1 = 'data/cutadapt/{sample}.cutadapt.r1.fq',
		r2 = 'data/cutadapt/{sample}.cutadapt.r2.fq'
	output:
		r1 = 'data/maxee/{sample}.maxee.r1.fa',
		r2 = 'data/maxee/{sample}.maxee.r2.fa'
	params:
		phred = config['phred']
	threads: 4
	log: 'logs/maxee/{sample}.maxee.log'
	run:
		shell('vsearch --fastq_filter {input.r1} --reverse {input.r2} --fastq_ascii {params.phred} --fastq_minlen 50 --fastq_truncee 1 --fasta_width 0 --fastaout {output.r1} --fastaout_rev {output.r2} --threads {threads} 2> {log}')

# Count reads for every step
rule report:
	input:
		raw			='data/samples/{sample}_1.fq.gz',
		cutadapt	='data/cutadapt/{sample}.cutadapt.r1.fq',
		maxee		='data/maxee/{sample}.maxee.r1.fa'
	output:
		count		='data/count/{sample}.count'
	params:
		name='{sample}'
	run:
		 shell("c1=$(($(gunzip -c {input.raw} | wc -l | awk '{{print $1}}')/4));c2=$(($(wc -l {input.cutadapt} | awk '{{print $1}}')/4));c3=$(($(wc -l {input.maxee}| awk '{{print $1}}')/2));echo -e {params.name}'\t'$c1'\t'$c2'\t'$c3 > {output.count}")

rule concat:
	input:
		count=expand('data/count/{sample}.count', sample=SAMPLES)
	output:
		report='data/count_report.tsv'
	run:
		shell("echo -e sample'\t'raw'\t'cutadapt'\t'maxee")
		shell('cat {input.count} > {output.report}')