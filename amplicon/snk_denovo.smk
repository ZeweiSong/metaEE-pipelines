#############################################
# Amplicon sequencing processing pipeline   #
# for Pair End reads                        #
# Burst for closed-reference clustering     #
# Winner take all for combining two targets #
#############################################

# This is a Snakemake pipeline for processing amplicon sequencing data based on 
# DNBSeq(TM) pair-end method. 
# Most popular protocl for sequecning amplicon is to use a (very) long primer, 
# which contains the sequence of library adaptor. However, long primer may cause
# problems every long primer inherent, and you need two step PCR for this.
# In our case, we use a short primer with 6 bp barcode at the end. The amplicons
# can be treated as normal insert. In this way, we can use the PCR-free method to
# construct the library. By sacraficing a bit read length, we can achieve a ultra
# high throughput in a single run (~2000 samples).

# This protocol goes to the path of "denovo", which means that we will guess the OTUs first. 
# We will then assign taxonomy to each OTU. 
# We are not about the argue why we prefer closed-ref than de novo here, but
# one short answer is that closed-ref make independent projects accumulable (this is the right word?).

# You will need the sequence of the forward and reverse primers (or the conservative region
# you like to trim off). These primers are different based on the target amplified.
# Change the primer sequences in config.yaml to the one you used (and ask around if you don't know).
# Also check if you have the right Burst database (bacteria or fungi?).
# Check cluster.yaml if you use qsub.

# The sample file should be in the format [SampleName].[ReadDirection].fq.gz, i.e. sampleA.r1.fq.gz sampleA.r2.fq.gz
# Beward that "." should not be used to name your samples as it is a delimiter to identify sample names.
# You can change the file extension in the config file.
# Put all sample files in workpath/samples
# A sample should have two corresponding files.

configfile: 'config.yaml'

import os
import sys
SAMPLES = []
workpath = config['workpath']
project = config['project']
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

fw=config['primers']['fw']
rv=config['primers']['rv']
base = {'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','K':'M','M':'K','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N'}
frc = ''.join([base[i] for i in list(fw)[0][::-1]])
rrc = ''.join([base[i] for i in list(rv)[0][::-1]])

# Target of this pipeline: put all samples in one profile, either in tsv or biom
rule target:
	input:
		biom_taxa = workpath + project + '.taxa.biom',
		tsv_taxa  = workpath + project + '.taxa.tsv'

# We merge R1 and R2 together for a longer asselmby.
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

# Remove the conservative regions. Usually, you just need to remove the two primers. But for some cases, the region has to be longer (like ITS1).
rule cutadapt:
	input:
		merged	= workpath + 'pear/{sample}.assembled.fastq'
	output:
		r1 = workpath + 'cutadapt/{sample}.r1.fq',
		r2 = workpath + 'cutadapt/{sample}.r2.fq'
	params:
		fw=config['primers']['fw'],
		rv=config['primers']['rv']
		fwrc=frc
		rvrc=rrc
		ascii=config['ascii']
	threads: 2
	log: workpath + 'logs/cutadapt/{sample}.log'
	run:
		shell('cutadapt {input.merged} -g {params.fw} -a {params.rvrc} -n 2 --discard-untrimmed -e 0.1 -m 100 --quality-base {params.ascii} -j {threads} -o {output.r1} > {log}')
		shell('cutadapt {input.merged} -g {params.rv} -a {params.fwrc} -n 2 --discard-untrimmed -e 0.1 -m 100 --quality-base {params.ascii} -j {threads} -o {output.r2} >> {log}')

# Remove low quality reads, and reverse compliment R2 sequences
rule quality:
	input:
		r1 = workpath + 'cutadapt/{sample}.r1.fq',
		r2 = workpath + 'cutadapt/{sample}.r2.fq'
	output:
		r1 = workpath + 'quality/{sample}.r1.fa',
		r2 = workpath + 'quality/{sample}.r2.fa',
		combine = temp(workpath + 'relabel/{sample}.fa')
	log:
		workpath + 'logs/quality/{sample}.log'
	params:
		ascii = config['ascii'],
		trim_fw = config['trim']['fw'],
		trim_rv = config['trim']['rv'],
		sn = '{sample}'
	threads: 1
	run:
		shell('vsearch --fastq_filter {input.merged} --fastaout {output.r1} --fastq_maxee 1 --fastq_minlen 100 --fasta_width 0\
		--fastq_ascii {params.ascii} --fastq_stripleft {params.trim_fw} --fastq_stripright {params.trim_rv}  --threads {threads} 2>  {log}')
		shell('vsearch --fastq_filter {input.merged} --fastaout_rev {output.r2} --fastq_maxee 1 --fastq_minlen 100 --fasta_width 0\
		--fastq_ascii {params.ascii} --fastq_stripright {params.trim_fw} --fastq_striplet {params.trim_rv}  --threads {threads} 2>>  {log}')
		shell('cat {output.r1} {output.r2} > - | vsearch --fastq_filter - --fastaout {output.combine} --relabel "sample={params.sn};r1_"')

# Combine all reads into one
rule combine:
	input:
		samples = expand(workpath + 'relabel/{sample}.fa', sample=SAMPLES)
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
		biom = workpath + 'otutable.biom'
	params:
		id   = config['align']['id']
	log:
		workpath + 'logs/align/merged.log'
	threads: 4
	run:
		shell('vsearch --usearch_global {input.qc_reads} --db {input.cluster} --id {params.id} --biomout {output.biom} --threads {threads} > {log}')

# Add the taxonomy, and write to biom and tsv file	
rule taxonomy:
	input:
		biom  = workpath + 'otutable.biom'
	output:
		biom_taxa = workpath + 'combined.taxa.biom',
		tsv_taxa  = workpath + 'combined.taxa.tsv'
	params:
		taxa=config['database']['tax']
	log:
		workpath + 'logs/concat/concate.log'
	run:
		shell('biom add-metadata -i {input.biom} -o {output.biom_taxa} --observation-metadata-fp {params.taxa} --observation-header OTUID,taxonomy --output-as-json --sc-separated taxonomy')
		shell('biom convert -i {output.biom_taxa} -o {output.tsv_taxa} --to-tsv --header-key taxonomy')