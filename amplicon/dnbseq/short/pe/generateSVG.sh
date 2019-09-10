snakemake -s snakefile_amplicon_pe --dag | tail -n+3 | dot -Tsvg > dag.svg
# We have to remove the print out line
