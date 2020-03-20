#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

De novo split amplicon sequencing data.

This script will search the 21 eeBarcodes in the input FQ file, and output at most 21 sets of files
if they have reads >= threshold.

Coders who love to comment their code are unlikely to have bad luck.

Currently design for PE reads.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%%
from __future__ import print_function
from __future__ import division
import argparse
from metaSeq import io as seqIO
from itertools import combinations
from itertools import permutations
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', nargs='+', help='Input FASTQ files for splitting, use space to seperate twin files.')
parser.add_argument('-n', '--threshold', default=10000, type=int, help='Minimum reads number to keep.')
parser.add_argument('-o', '--prefix', type=str, help='Prefix for output file, include the output path.')
args = parser.parse_args()
inputFile = args.input
threshold = args.threshold
if args.prefix:
    outputPrefix = args.prefix
else:
    outputPrefix = inputFile[0].split('/')
    index = outputPrefix.index('_1')
    outputPrefix = outputPrefix[0:index]
for item in inputFile:
    print('Reading in {0} ...'.format(item))
print('The minimum reads number sets to {0}.'.format(threshold))
print('The output files will have prefix as {0}.'.format(outputPrefix))

#%% Settle the eeBarcode and their mutated deritives
# This section will create a look up dictionary for splitting
eeBarcodesByNumber = {1: 'AACTGC', 2: 'ACGGTT', 3: 'AGCATG', 4: 'AGGACT', 5: 'CAGGAT', 6: 'CAGTTC', 7: 'CCATTG',\
              8: 'CGATAC', 9: 'CTAGAG', 10: 'CTCAGT', 11: 'GAATCC', 12: 'GACTTG', 13: 'GATGAC', 14: 'GCGTTA',\
              15: 'GTACGA', 16: 'GTCCAT', 17: 'TCCTAG', 18: 'TCGATG', 19: 'TGGCAT', 20: 'TTGAGC', 21: 'AACGAG'}
eeBarcodes = {}
for key, value in eeBarcodesByNumber.items():
    eeBarcodes[value] = key

print('There are {0} metaEE barcodes.'.format(len(eeBarcodes)))
print('\tAdding one mutation event ...')

# Add in one mutation barcodes
# metaEE barcode has one mutation tolerance.
mutation = {'A':('T','C','G'), 'T':('A','C','G'), 'C':('A','T','G'), 'G':('A','T','C')}
mutationPool = {1:{i:[] for i in range(1,22,1)}, 2:{i:[] for i in range(1,22,1)}, 3:{i:[] for i in range(1,22,1)}}
#% Create all mutation possibles
for string, number in eeBarcodes.items():
    for i in range(1,4,1): # Mutate from 1 to 3 SNPs, but only one mutation is redundancy free
        for loc in combinations(range(6), i):
            origin = tuple([string[j] for j in loc])
            for mut in permutations(['A','T','C','G']*len(loc), len(loc)):
                if sum([x != y for x, y in zip(origin, mut)]) == i:
                    stringMutate = [s for s in string]
                    for index, item in enumerate(loc):
                        stringMutate[item] = mut[index]
                    mutationPool[i][number].append(''.join(stringMutate))
for i in range(1,4,1):
    for j in range(1,22,1):
        mutationPool[i][j] = tuple(set(mutationPool[i][j]))
# For one mutations case, create another pool
# for a stirng with 6 bp, there should be 6 X 3 = 18 one mutation strings
# for a set of 21 barcodes, there should be 21 X 18 = 378 strings
# Plus the original 21 barcodes, the total candidates will be 399
eeBarcodesMut = {} # This is the final index for splitting
for key, value in mutationPool[1].items():
    for string in value:
        eeBarcodesMut[string] = key
for key, value in eeBarcodes.items(): # add the origin 21 barcodes
    eeBarcodesMut[key] = value
print('The final lookin up dictionary has {0} barcodes.'.format(len(eeBarcodesMut)))

#%% Split
# Reads will be read in in pair, and assign to each barcode if eligible
print('Start splitting {0} ...'.format(' '.join(inputFile)))
count = [0, 0]
split = {i+1:[] for i in range(21)}
seqHandle = seqIO.sequence_twin(inputFile[0], inputFile[1])
for r1, r2 in seqHandle:
    barcodeR1 = eeBarcodesMut.get(r1[1][:6], 999)
    barcodeR2 = eeBarcodesMut.get(r2[1][:6], 998)
    if barcodeR1 == barcodeR2:
        split[barcodeR1].append((r1, r2))
        count[0] += 1
    else:
        count[1] += 1
s = sum([len(i) >= threshold for i in split.values()])
scount = sum([len(i) for i in split.values() if len(i) >= threshold])
print('Finished splitting:')
print('\tTotal reads:\t\t{0}'.format(sum(count)))
print('\tEligible reads:\t\t{0}\t{1:3.2f}%.'.format(count[0], count[0]/sum(count)*100))
print('\tReads passed threshold:\t{0}\t{1:3.2f}%'.format(scount, scount/sum(count)*100))
print('\tFound {0} barcodes with more than {1} reads.'.format(s,  threshold))

# Write to files for barcodes pass the threshold
i = 1
for key, value in split.items():
    if len(value) >= threshold:
        outputR1 = outputPrefix + '_' + str(key) + '.r1.fq.gz'
        outputR2 = outputPrefix + '_' + str(key) + '.r2.fq.gz'
        c = seqIO.write_seqs([r[0] for r in value], outputR1, fastx='q', gz=True)
        print('\tPair{2}: {0} reads wrote to {1}.'.format(c, outputR1, i))
        c = seqIO.write_seqs([r[1] for r in value], outputR2, fastx='q', gz=True)
        print('\tPair{2}: {0} reads wrote to {1}.'.format(c, outputR2, i))
        i += 1
    else:
        pass
print('Finished processing {0}.'.format(' '.join(inputFile)))