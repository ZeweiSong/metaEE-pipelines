# BIOMASS: Bayesian Inference Of Multi-Alignments on genome SetS <h1>

What is it?
It is a pipeline for solving microbial community profile from shotgun sequencing data.
A profile is inferred from multiple alignment of the shotgun reads onto given genome sets, in which aligner reports ALL tied best hits (or all hits pass the threshold).
Tied hits are resolved by calculating their origin probabilities using Bayesian inference.
For example, if a query aligns to target A and B, and the resulting probabilities is 0.8 and 0.2, then the read count 1 is divided into 0.8 to A and 0.2 to B.
It is possible that the origin probability of some target be zero, in which case it gets 0 count.

How to use it?
This pipeline contains two sankemake files, snk_qc.smk and snk_profile.smk
snk_qc.smk: Merge the pair-end read raw sequnence data, quality filter into clean data.
snk_profile.smk: Do the heavy duty aligning, calculate profiles for all samples.

Since aligning to all genomes takes a decent time (in days), it is wise to separate QC and profiling into two steps.
After QC, we recommend to draw 1M reads from each sample for alignments.
We only use merged reads for aligning, so draw from .merged.fa files.
