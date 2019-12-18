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

config = 'config.yaml'

import os
import sys
SAMPLES = {}
filepath = 'data/sample/'
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

rule pear:
	input:
		lambda wildcards: SAMPLES[wildcards.sample] # a lambda expression that "sample" is used as wildcards to refer to the value of the dictionary SAMPLES
	output:
		assembled	='data/pear/{sample}.assembled.fastq',
		forward		='data/pear/{sample}.forward.fastq',
		reverse		='data/pear/{sample}.reverse.fastq',
		discard		='data/pear/{sample}.discard.fastq'
	params:
		name = '{sample}'
		phred = config['phred']
	threads: 4
	log: 'log/pear/{sample}.pear.log'
	run:
		shell('pear -f {input[0]} -r {input[1]} -o {params.name} -k -b {params.phred} -j {threads} > {log}')

# For unassembled files, remove adaptors from their tails.
rule cutadapt:
	input:
		forward		='data/pear/{sample}.forward.fastq',
		reverse		='data/pear/{sample}.reverse.fastq'
	output:
		pair		='data/cutadapt/{sample}.cutadapt.fq'
	params:
		fw = 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA',
		rv = 'AAGTCGGATCGTAGCCATGTCGTTCTGT',
		phred = config['phred']
	threads: 4
	log: 'log/cutadapt/{sample}.cutadapt.log'
	run:
		shell('cutadapt {input.forward} {input.reverse} -a {params.fw} -A {params.rv} --interleaved {output.pair} -m 50 --quality-base {params.phred} -j {threads} > {log}')

# For assembled reads, discard if maxee > 1
# For unassembled read pair, truncate until maxee <= 1, keep only pair
rule maxee:
	input:
		assembled	='data/pear/{sample}.assembled.fastq',
		pair		='data/cutadapt/{sample}.cutadapt.fq'
	output:
		assembled	='data/maxee/{sample}.merged.fa.gz',
		pair		='data/maxee/{sample}.pair.fa.gz'
	params:
		phred = config['phred']
	threads: 4
	log: 'log/maxee/{sample}.maxee.log'
	run:
		shell('vsearch --fastq_filter {input.assembled} --fastq_maxee 1 --fastq_minlen 100 --fastq_maxns 0 --fastaout - --fasta_width 0 --fastq_ascii {params.phred} --threads {threads} 2> {log} | gzip -c > {output.assembled}')
		shell('vsearch --fastq_filter {input.pair} --fastq_truncee 1 --fastq_minlen 100 --fastq_maxns 0 --fastaout - --fasta_width 0 --fastq_ascii {params.phred} --threads {threads} 2>> {log} | seqtk dropse | gzip -c > {output.pair}')

# Count reads for every step
rule report:
	input:
		raw			='data/sample/{sample}.fq.gz',
		merged		='data/pear/{sample}.assembled.fastq',
		merged_qc	='data/maxee/{sample}.merged.fa.gz',
		paired		='data/pear/{sample}.forward.fastq',
		paired_qc	='data/maxee/{sample}.pair.fa.gz'
	output:
		count		='data/count/{sample}.count'
	params:
		name='{sample}'
	run:
		 shell("c1=$(($(gunzip -c {input.raw} | wc -l | awk '{{print $1}}')/4));c2=$(($(wc -l {input.merged} | awk '{{print $1}}')/4));c3=$(($(gunzip -c {input.merged_qc} | wc -l | awk '{{print $1}}')/2));c4=$(($(wc -l {input.paired} | awk '{{print $1}}')/4));c5=$(($(gunzip -c {input.paired_qc} | wc -l | awk '{{print $1}}')/4));echo -e {params.name} '\t' $c1 '\t' $c2 '\t' $c4 '\t' $c3 '\t' $c5 > {output.count}")

rule concat:
	input:
		count=expand('data/count/{sample}.count', sample=SAMPLES)
	output:
		report='data/count_report.tsv'
	run:
		cat {input.count} > {output.report}