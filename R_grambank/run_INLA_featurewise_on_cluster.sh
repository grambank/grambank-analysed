#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J Hedvig_Oceanic_ASR_mcct

##1w:code skirgard$ chmod 755 run_INLA_featurewise_on_cluster.sh
##lingn01w:code skirgard$ ./run_INLA_featurewise_on_cluster.sh

#Step 1
Rscript cluster_set_up.R

FILE=rlibs
if [ -f "$FILE" ]; then
    echo "$FILE exists."
    echo "rlibs folder already exists, won't bother remaking it."
else 
mkdir rlibs
fi

export R_LIBS_USER=rlibs

#Rscript requirements.R 

#Step 2 prep data
#make data
#let's not run this because we don't want to have to move over all the git submodules, let's run make data before moving to cluster.

#Step 3 run INLA featurewise
make INLA_featurewise


