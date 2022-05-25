#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J INLA_all
#SBATCH --partition=dlcegpu

#set ..rlibs as r package dir
export R_LIBS_USER=../rlibs

#set up installation of scripts etc for the cluster.
Rscript cluster_set_up.R

#run INLA featurewise
make INLA_featurewise


