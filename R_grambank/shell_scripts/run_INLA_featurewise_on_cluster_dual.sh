#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J INLA_dual
#SBATCH --partition=dlcegpu

#set ..rlibs as r package dir
FILE=../rlibs
if [ -f "$FILE" ]; then
    echo "$FILE exists."
    echo "rlibs folder already exists, won't bother remaking it."
else 
mkdir ../rlibs
fi

export R_LIBS_USER=../rlibs

#set up installation of scripts etc for the cluster.
Rscript cluster_set_up.R

#run INLA featurewise
make INLA_featurewise_dual