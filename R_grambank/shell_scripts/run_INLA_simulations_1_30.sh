#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J inlasim.1.60
#SBATCH --partition=dlcegpu

source ~/shh_custom_profile 

export R_LIBS_USER=../rlibs

echo "Installing and loading pacakges..."
Rscript pkg_installations.R

echo "#.. done with pkgs"

echo "#simulating data"
#Rscript simulate_data.R

echo "#making prec mats."
#Rscript make_precisionmatrices.R

echo "#Running INLA."´
Rscript simulation_INLAmodelscript.R 1 60
