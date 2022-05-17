#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J Hedvig_Oceanic_ASR_mcct

##1w:code skirgard$ chmod 755 run_INLA_featurewise_on_cluster.sh
##lingn01w:code skirgard$ ./run_INLA_featurewise_on_cluster.sh

#Step 1
echo first step, installing necessary packages
Rscript requirements.R 

#Step 2 prep data
make data

#Step 3 run INLA featurewise
make INLA_featurewise


