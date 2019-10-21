# run on SGE cluster
snakemake \
--snakefile snakefile_amplicon_pe \
--configfile config.yaml \
--cluster-config cluster.yaml \
--jobs 20 \
--cluster "qsub -S /bin/bash -cwd -q {cluster.queue} -P {cluster.project} -l vf={cluster.mem},p={cluster.cores} -binding linear:{cluster.cores} -o {cluster.output} -e {cluster.error}" \
--latency-wait 360 \
--until target
