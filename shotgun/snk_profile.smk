# Bayesian Inference Of Multi-Alignments on genome SetS
# Snakemake Pipeline

# We recommend subsample reads to 1M per sample for aligning (or lower for testing).
# You can do this using "seqtk sample"
# All files should have extension ".fna", any string before it will be identified
# as sample name.
# Put all subsampled files into the folder defined in the variable "filepath"

# If you use the script "feed_me_the_biomass.sh" to set up this pipeline,
# you should have GTDB_r95 and NCBI fungal genomes softlinked under the 
# folder "database/".
# Remove either gtdb_r95* or ncbi_fungi* if you only want to align to 
# fungal genomes or Prokaryota genomes.

configfile: 'config.yaml'

# Read in sample reads. They should be merged reads.
import os
import sys
SAMPLES = {}
filepath = 'data/maxee_subsamples/'
for file in os.listdir(filepath):
		if file.endswith('.fna'):
			sample_name = file.split('.')[0]
			SAMPLES[sample_name] = filepath + file
if len(SAMPLES) == 0:
	print('ERROR: no file under {0}.'.format(filepath))
	sys.exit('Pipeline terminated orz.')
else:
	pass
DBS = {}
db_path = config['db_path']
for file in os.listdir(db_path):
	if file.endswith('.edx'):
		db_name = file.split('.')[0]
		DBS[db_name] = db_path + file

# What you get is a profile of all samples.
rule target:
	input:
		biom = 'data/profiles.biom',
		tsv = 'data/profiles.tsv'

# Put on sample labels to all sequences.
rule rename:
	input:
		seqs = lambda wildcards: SAMPLES[wildcards.sample]
	output:
		seqs = 'data/rename/{sample}.fna'
	params:
		name = '{sample}'
	log: 'log/profile/rename/{sample}.rename.log'
	run:
		shell("vsearch --fastx_filter {input.seqs} --relabel '{params.name};' --fastaout {output.seqs} --fasta_width 0 2> {log}")

# Concat all sequence files so we only have to align once (per index).
rule concat_sequences:
	input:
		seqs = expand('data/rename/{sample}.fna', sample=SAMPLES)
	output:
		cat = 'data/samples_concat.fna'
	run:
		shell('echo Concat all samples together.')
		shell('cat {input.seqs} > {output.cat}')

# Align to index group
rule align:
	input:
		db = lambda wildcards: DBS[wildcards.db],
		cat = 'data/samples_concat.fna'
	output:
		aln = 'data/alignments/{db}.b6'
	threads: 4
	params:
		mode = config['mode'],
		id = config['id'],
		db_name = db_path + '{db}'
	log: 'log/profile/align/{sample}.align.log'
	run:
		shell('burst15 -q {input.cat} -r {params.db_name}.edx -a {params.db_name}.acx -o {output.aln} -i {params.id} -t {threads} -m {params.mode} -fr > {log}')

# Separate alignments by samples (we tag them at rule rename)
rule split_samples:
	input:
		aln = 'data/alignments/{db}.b6'
	output:
		aln_split = expand('data/aln_split/{sample}.{{db}}.b6', sample=SAMPLES)
	params:
		db_name = '{db}'
	log: 'log/profile/split_samples/{sample}.split_samples.log'
	run:
		shell('touch {output.aln_split}')
		shell('scripts/split_alignments.py -i {input.aln} -n {params.db_name} -o data/aln_split > {log}')

# Solve tied hits for profile
# You can solve at any taxonomic level using the option -l, but species is the most common one.
rule solve_profile:
	input:
		expand('data/aln_split/{{sample}}.{db}.b6', db=DBS)
	output:
		'data/profiles/{sample}.tsv'
	params:
		name = '{sample}',
		taxonomy = config['taxonomy']
	log: 'log/profile/solve_profile/{sample}.solve_profile.log'
	run:
		shell('scripts/solve_profile.py -a {input} -t {params.taxonomy} -o {output} -l species -n {params.name}')

# Convert all tsv profiles into biom and combined'em ALL!
rule convert_biom:
	input:
		tsv = 'data/profiles/{sample}.tsv'
	output:
		biom = 'data/profiles_biom/{sample}.biom'
	run:
		shell('biom convert -i {input.tsv} -o {output.biom} --to-json')

rule concat_profiles:
	input:
		biom = expand('data/profiles_biom/{sample}.biom', sample=SAMPLES)
	output:
		biom = 'data/profiles.biom',
		tsv = 'data/profiles.tsv'
	log: 'log/profile/concat_profiles/{sample}.concat_profiles.log'
	run:
		shell('biom_concat.py -i {input.biom} -biom_out {output.biom} > {log}')
		shell('biom convert -i {output.biom} -o {output.tsv} --to-tsv')