#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:36:47 2018

Create softlinks from splitted sequence files to the sample target.

The softlink represents a single target ready for snakemake processing.

You need to create a tab-delimitged file with two columns, and no header.

The first column contains the sample-target name of your project,
the second column contains the prefix of the splitted file.

Coders who love to comment their code are unlikely to have bad luck.

@author: Zewei Song
@email: songzewei@genomics.cn
"""
#%
from __future__ import print_function
from __future__ import division
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_folder', help='Name of the input folder.')
parser.add_argument('-o', '--output_folder', help='Name of the output folder.')
parser.add_argument('-l', '--list', help='Name of the mapping list file.')
args = parser.parse_args()
input_folder = args.input_folder
if input_folder[-1] != '/':
    input_folder += '/'
output_folder = args.output_folder
if output_folder[-1] != '/':
    output_folder += '/'
list_file = args.list

# Get the file list of current target
print('Start creating softlinks from {0} to {1}:'.format(input_folder, output_folder))
file_list = []
with open('list_file', 'r') as f: # this is a two column tab-delimited file
    for line in f:
        line = line.strip('\n').split('\t')
        file_list.append(line)
print('Found {0} candidate targets in {1}.'.format(len(file_list), list_file))

nofile_count = 0
file_count = 0

# Create softlinks for all sample-target
for item in file_list:
    src = ''
    dst = ''
    src1 = input_folder + item[1] + '.r1.fq.gz'
    dst1 = output_folder + item[0] + '.r1.fq.gz'
    src2 = input_folder + item[1] + '.r2.fq.gz'
    dst2 = output_folder + item[0] + '.r2.fq.gz'
    
    if os.path.isfile(src1) and os.path.isfile(src2): # Both raw sequence file exist
        os.symlink(src1, dst1)
        os.symlink(src2, dst2)
        file_count += 1
    else:
        nofile_count += 1

print('Finished.')
print('Fount {0} target without source files.'.format(nofile_count))
print('Softlinks created for {0} targets.'.format(file_count))