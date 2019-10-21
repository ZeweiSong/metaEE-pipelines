#############################################
# Amplicon sequencing processing pipeline   #
# for Pair End reads                        #
# Burst for closed-reference clustering     #
# Winner take all for combining two targets #
#############################################

# This is a Snakemake pipeline for processing amplicon sequencing data based on 
# MiSeq or HiSeq pair-end method. The most popular protocl for illumina's sequencer is
# to use a (very) long primer, which contains the sequence of adaptor to construt the library
# In this case, Read1 will always starting from the forward primer, and vice versa for Read2.
# The pro is that you know what you are expecting in R1 and R2.
# the con is that you have to perfrom PCR with two long primers, and this sometimes can cause
# trouble.

# This protocol goes to the path of "closed-reference", which means that we only align reads
# to known species. We are not about the argue why we prefer closed-ref than de novo here, but
# one short answer is that closed-ref make independent projects accumulable (this is the right word?).

# You will need the sequence of the forward and reverse primers (or the conservative region
# you like to trim off). These primers are different based on the target amplified.
# Change the primer sequences in config.yaml to the one you used (and ask around if you don't know).
# Also check if you have the right Burst database (bacteria or fungi?).
# Check cluster.yaml if you use qsub.

# The sample file should be in the format [SampleName].[ReadDirection].fq.gz, i.e. sampleA.r1.fq.gz sampleA.r2.fq.gz
# Beward that "." should not be used to name your samples as it is a delimiter to identify sample names.
# You can change the file extension in the config file.
# Put all sample files in its1/samples
# A sample should have two corresponding files.

configfile: 'config_illumina_pe.yaml'

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

# Remove low quality reads, and reverse compliment R2 sequences
rule quality:
	input:
		merged	= workpath + 'pear/{sample}.assembled.fastq'
	output:
		merged = workpath + 'qc_reads/{sample}.fa'
	log:
		workpath + 'logs/qc_reads/{sample}.log'
	params:
		ascii = config['ascii']
		trim_fw = config['trim']['fw'],
		trim_rv = config['trim']['rv']
	threads: 1
	run:
		shell('vsearch --fastq_filter {input.merged} --fastaout {output.merged} --fastq_maxee 1 --fastq_minlen 100 --fasta_width 0\
		--fastq_ascii {params.ascii} --fastq_stripleft {params.trim_fw} --fastq_stripright {params.trim_rv}  --threads {threads} 2>  {log}')

# Combine all reads into one
rule combine:
	input:
		samples = expand(workpath + 'qc_reads/{sample}.fa', sample=SAMPLES)
	output:
		merged = workpath + combine/merged.fa'
	log:
		workpath + 'logs/combine/merged.log'	
	threads: 1
	run:
		shell('cat {input.samples} > {output.merged}')
# Dereplicate
rule dereplicate:
	input: 
		merged = workpath + combine/merged.fa'
	output:
		derep = workpath + 'dereplicate/merged.fa'
	log:
		workpath + 'logs/dereplicate/merged.log'
	threads: 2
	run:
		shell('vsearch --derep_fulllength {input.merged} --output {output.derep} --sizeout --fasta_width 0 --threads {threads} 2> {log}')
# Chimera check
rule chimera:
	input: 
		derep = workpath + dereplicate/merged.fa'
	output:
		derep = workpath + 'chimera/merged.fa'
	log:
		workpath + 'logs/chimera/merged.log'
	params:
		db = config['database']['fasta']
	threads: 2
	run:
		shell('vsearch --uchime_ref {input.derep} --db {params.db} --nonchimeras {output.derep} --sizeout --fasta_width 0 --threads {threads} 2> {log}')
# Cluster
	input:
		derep = workpath + 'chimera/merged.fa'
	output:
		cluster = workpath + 'cluster/centroids.fa'
	log:
		workpath + 'logs/cluster/merged.log'
	params:
		id   = config['align']['id']
	threads: 1
	run:
		shell('vsearch --cluster_size {input.derep} --id {params.id} --centroids {output.cluster} --threads {threads} --fasta_width 0')
	
# Align reads to the centroids
rule align:
	input:
		qc_reads = workpath + 'qc_reads/{sample}.fa',
		cluster = workpath + 'cluster/centroids.fa'
	output:
		tsv = workpath + 'otutable.tsv'
	params:
		id   = config['align']['id']
	log:
		workpath + 'logs/align/merged.log'
	threads: 4
	run:
		shell('vsearch --usearch_global {input.qc_reads} --db {input.cluster} --id {params.id} --otutabout {output.tsv} --threads {threads} > {log}')

# Count reads retention in the above steps
rule count:
	input:
		quality_fw	=	workpath + 'qc_reads/{sample}.fa'
		align_fw	=	workpath + 'alignments/{sample}.b6'
	output: workpath + 'count/{sample}.count'
	params:
		sn ='{sample}'
	run:
		shell('echo {params.sn} > {output}')
		shell('amplicon_countRetention.py -d {input.quality_fw} -q {input.quality_fw} -a {input.align_fw} >> {output}')

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