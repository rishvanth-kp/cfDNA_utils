#!/bin/bash

snakemake --snakefile Snakefile \
--printshellcmds \
--keep-going \
--rerun-incomplete \
--cluster-config cluster.yaml \
--cores 32 \
--cluster 'sbatch --partition={cluster.partition} --ntasks={cluster.cores} --mem={cluster.mem} --time={cluster.time} -o {cluster.logout} -e {cluster.logerror}' 
