configfile: 'config.yaml'

import os
import sys
SAMPLES = []
workpath = config['workpath']
if workpath[-1] != '/':
	workpath = workpath + '/'
project = config['project']
for file in os.listdir(workpath + 'samples/'): # List all files under the folder. For PE data, one samples correspond to two files.
	string = file.split('.')
	if len(string) != 4: # You have to makesure that "." is not used in naming your samples. All files have the suffix .r1.fq.gz or .r2.fq.gz
		print('{0} contains illegal naming style, please check.'.format(file))
		sys.exit()
	else:
		SAMPLES.append(string[0])
SAMPLES = tuple(set(SAMPLES))
SAMPLES = {i:i for i in SAMPLES}

fw=config['primers']['fw']
rv=config['primers']['rv']
base = {'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','K':'M','M':'K','S':'S','W':'W','B':'V','V':'B','D':'H','H':'D','N':'N'}
frc = ''.join([base[i] for i in list(fw)[::-1]])
rrc = ''.join([base[i] for i in list(rv)[::-1]])

rule target:
	input:
		combine_biom			= workpath + 'concat/' + project + '.biom',
		comebin_taxa_biom		= workpath + 'concat/' + project + '.taxa.biom',
		combine_taxa_tsv		= workpath + 'concat/' + project + '.taxa.tsv'

rule cutadapt:
	input:
		r1 = lambda wildcards: workpath + 'samples/' + SAMPLES[wildcards.sample] + '.r1.fq.gz',
		r2 = lambda wildcards: workpath + 'samples/' + SAMPLES[wildcards.sample] + '.r2.fq.gz'
	output:
		r1fw = temp(workpath + 'cutadapt/{sample}.r1fw.fq'),
		r2rv = temp(workpath + 'cutadapt/{sample}.r2rv.fq'),
		r1rv = temp(workpath + 'cutadapt/{sample}.r1rv.fq'),
		r2fw = temp(workpath + 'cutadapt/{sample}.r2fw.fq')
	params:
		fw		= config['primers']['fw'],
		rv		= config['primers']['rv'],
		fwrc	= frc,
		rvrc	= rrc,
		ascii	= config['ascii']
	log:
		workpath + 'logs/cutadapt/{sample}.log'
	threads: 2
	run:
		shell('cutadapt {input.r1} {input.r2} -g {params.fw} -a {params.rvrc} -G {params.rv} -A {params.fwrc} --discard-untrimmed -n 2 -e 0.1 -m 75 --quality-base {params.ascii} -j {threads} -o {output.r1fw} -p {output.r2rv} > {log}')
		shell('cutadapt {input.r1} {input.r2} -g {params.rv} -a {params.fwrc} -G {params.fw} -A {params.rvrc} --discard-untrimmed -n 2 -e 0.1 -m 75 --quality-base {params.ascii} -j {threads} -o {output.r1rv} -p {output.r2fw} >> {log}')

rule quality_control:
	input:
		r1fw 		= workpath + 'cutadapt/{sample}.r1fw.fq',
		r2rv 		= workpath + 'cutadapt/{sample}.r2rv.fq',
		r2fw 		= workpath + 'cutadapt/{sample}.r1rv.fq',
		r1rv 		= workpath + 'cutadapt/{sample}.r2fw.fq'
	output:
		r1fw 		= temp(workpath + 'quality_control/{sample}.r1fw.fq'),
		r2rv 		= temp(workpath + 'quality_control/{sample}.r2rv.fq'),
		r2fw 		= temp(workpath + 'quality_control/{sample}.r2fw.fq'),
		r1rv 		= temp(workpath + 'quality_control/{sample}.r1fv.fq')
	log:
		workpath + 'logs/quality_control/{sample}.log'
	params:
		ascii = config['ascii']
	threads: 4
	run:
		shell('vsearch --fastq_filter {input.r1fw} --reverse {input.r2rv} --fastqout {output.r1fw} --fastqout_rev {output.r2rv} --fastq_maxee 1 --fastq_minlen 50 --fastq_ascii {params.ascii} --threads {threads} 2> {log}')
		shell('vsearch --fastq_filter {input.r1rv} --reverse {input.r2fw} --fastqout {output.r1rv} --fastqout_rev {output.r2fw} --fastq_maxee 1 --fastq_minlen 50 --fastq_ascii {params.ascii} --threads {threads} 2>> {log}')

# Trim off sequences not in database and convert to FASTA
rule trim_tail:
	input:
		r1fw 		= workpath + 'quality_control/{sample}.r1fw.fq',
		r2rv 		= workpath + 'quality_control/{sample}.r2rv.fq',
		r2fw 		= workpath + 'quality_control/{sample}.r2fw.fq',
		r1rv 		= workpath + 'quality_control/{sample}.r1fv.fq'
	output:
		r1fw 		= temp(workpath + 'trim_tail/{sample}.r1fw.fa'),
		r2rv 		= temp(workpath + 'trim_tail/{sample}.r2rv.fa'),
		r2fw 		= temp(workpath + 'trim_tail/{sample}.r2fw.fa'),
		r1rv 		= temp(workpath + 'trim_tail/{sample}.r1fv.fa')
	params:
		fw_trim = config['trim']['fw'],
		rv_trim = config['trim']['rv']
	run:
		shell('seqtk trimfq -b {params.fw_trim} {input.r1fw} - | seqtk seq    -A -L 50 - > {output.r1fw}')
		shell('seqtk trimfq -e {params.rv_trim} {input.r2rv} - | seqtk seq -r -A -L 50 - > {output.r2rv}')
		shell('seqtk trimfq -e {params.rv_trim} {input.r1rv} - | seqtk seq -r -A -L 50 - > {output.r1rv}')
		shell('seqtk trimfq -b {params.fw_trim} {input.r2fw} - | seqtk seq    -A -L 50 - > {output.r2fw}')

rule merge:
	input:
		r1fw 	= workpath + 'trim_tail/{sample}.r1fw.fa',
		r2rv 	= workpath + 'trim_tail/{sample}.r2rv.fa',
		r2fw 	= workpath + 'trim_tail/{sample}.r2fw.fa',
		r1rv 	= workpath + 'trim_tail/{sample}.r1fv.fa'
	output:
		fw 		= workpath + 'merge/{sample}.fw.fa',
		rv		= workpath + 'merge/{sample}.rv.fa'
	run:
		shell('cat {input.r1fw} {input.r2fw} > {output.fw}')
		shell('cat {input.r2rv} {input.r1rv} > {output.rv}')
		
rule alignment:
	input:
		fw 		= workpath + 'merge/{sample}.fw.fa',
		rv		= workpath + 'merge/{sample}.rv.fa'
	output:
		fw 		= temp(workpath + 'alignment/{sample}.fw.b6'),
		rv		= temp(workpath + 'alignment/{sample}.rv.b6'),
		fwrv	= workpath + 'alignment/{sample}.fwrv.b6'
	threads: 4
	params:
		acx  = config['database']['acx'],
		edx  = config['database']['edx'],
		mode = config['align']['mode'],
		id   = config['align']['id']
	run:
		shell('burst12 -q {input.fw} -a {params.acx} -r {params.edx} -o {output.fw} -i {params.id} -m {params.mode} -t {threads} > {log}')
		shell('burst12 -q {input.rv} -a {params.acx} -r {params.edx} -o {output.rv} -i {params.id} -m {params.mode} -t {threads} > {log}')
		shell('amplicon_keepPairAln.py -i {output.fw},{output.rv} -o {output.fwrv} > {log}')
		
rule profiles:
	input:
		fwrv	= workpath + 'alignment/{sample}.fwrv.b6'
	output:
		fwrv	= workpath + 'profiles/{sample}.fwrv.tsv',
		fwrv_biom	= workpath + 'profiles/{sample}.fwrv.biom'
	params:
		sampleName='{sample}'
	run:
		shell('amplicon_winnerTakeAll.py -i {input.fwrv} -sn {params.sampleName} -t {output.fwrv} -g > {log}')
		shell('biom convert -i {output.fwrv} -o {output.fwrv_biom} --to-json')

rule count:
	input:
		r1_preQC 	= workpath + 'cutadapt/{sample}.r1fw.fq',
		r2_preQC 	= workpath + 'cutadapt/{sample}.r2fw.fq',
		r1_aftQC	= workpath + 'quality_control/{sample}.r1fw.fq',
		r2_aftQC	= workpath + 'quality_control/{sample}.r2fw.fq',
		alignment	= workpath + 'alignment/{sample}.fwrv.b6'
	output:
		count = workpath + 'count/{sample}.count'
	run:
		shell("c1=$(($(wc -l {input.r1_preQC} | awk '{{print $1}}')/4));c2=$(($(wc -l {input.r2_preQC} | awk '{{print $1}}')/4));c3=$(($(wc -l {input.r1_aftQC} | awk '{{print $1}}')/2));c4=$(($(wc -l {input.r2_aftQC} | awk '{{print $1}}')/4));c5=$(($(awk '{{print $1}}' {input.alignment} | wc -l  | awk '{{print $1}}')/4));echo -e {params.name} '\t' $c1 '\t' $c2 '\t' $c3 '\t' $c4 '\t' $c5 > {output.count}")
		
rule concat:
	input:
		fwrv	= expand(workpath + 'profiles/{sample}.fwrv.biom', sample=SAMPLES),
		count 	= expand(workpath + 'count/{sample}.count', sample=SAMPLES)
	output:
		combine_biom			= workpath + 'concat/' + project + '.biom',
		comebin_taxa_biom		= workpath + 'concat/' + project + '.taxa.biom',
		combine_taxa_tsv		= workpath + 'concat/' + project + '.taxa.tsv',
		count 						= workpath + 'concat/' + project + '.summary'
	params:
		taxa=config['database']['tax']
	log:
		workpath + 'logs/concat/concate.log'
	run:
		shell('amplicon_concat.py -i {input.fwrc} -biom_out {output.combine_fwrc_biom} > {log}')
		shell('biom add-metadata -i {output.combine_fwrc_biom} -o {comebin_fwrc_taxa_biom} --observation-metadata-fp {params.taxa} --observation-header OTUID,taxonomy --output-as-json --sc-separated taxonomy')
		shell('biom convert -i {comebin_fwrc_taxa_biom} -o {combine_fwrc_taxa_tsv} --to-tsv --header-key taxonomy')
		shell('cat {input.count} > {output.count}')