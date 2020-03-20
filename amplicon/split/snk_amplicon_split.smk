# configfile: 'config.yaml'

import fnmatch
import os
import sys
SAMPLES = {}
sampleNames = []
ext = 'fq.gz' # Sometimes people use .fq instead of .fq.gz or other suffix.
for file in os.listdir('raw/'): # List all files under the folder. For PE data, one samples correspond to two files.
	if fnmatch.fnmatch(file, '*'+ext):
		string = file.split('.')[0][:-2]
		sampleNames.append(string)
sampleNames = tuple(set(sampleNames)) # Deduplicate sample names for PE data.
for names in sampleNames:
	SAMPLES[names] = []
for file in os.listdir('raw/'):
	if fnmatch.fnmatch(file, '*'+ext):
		string = file.split('.')[0][:-2]
		SAMPLES[string].append('raw/' + file)
for key in SAMPLES.keys():
	SAMPLES[key].sort()
if len(SAMPLES) == 0:
	print("ERROR: no file found under raw/")
	print("Please check the file extension in the config file.")
	print("Snakemake has terminated.")
	sys.exit("Snake did not make it orz.")

rule target:
	input:
		report='split_report.tsv'

rule split:
	input:
		lambda wildcards: SAMPLES[wildcards.sample]
	output:
		"splitted_files/{sample}.finished"
	params:
		sampleName='{sample}'
	log:
		'log/split/{sample}.log'
	run:
		shell('python split_by_library.py -i {input[0]} {input[1]} -n 10000 -o splitted_files/{params.sampleName} > {log}')
		shell('echo "{params.sampleName} finished." > {output}')

rule summary:
	input:
		report = expand('splitted_files/{sample}.finished', sample=SAMPLES)
	output:
		report='split_report.tsv'
	run:
		shell('cat {input.report} > {output.report}')