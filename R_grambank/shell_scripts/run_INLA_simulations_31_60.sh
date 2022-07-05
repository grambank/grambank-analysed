#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J inlasim.61.120
#SBATCH --partition=dlcegpu

source ~/shh_custom_profile 

export R_LIBS_USER=../rlibs

echo "#Running INLA."´
Rscript simulation_INLAmodelscript.R 61 120
