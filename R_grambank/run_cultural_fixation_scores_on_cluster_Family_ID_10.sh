#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J Hedvig_cultural_fixation-scores_Family_10
#SBATCH --partition=dlcegpu

export R_LIBS_USER=../rlibs

Rscript dist_fixation_scores/fixation_scores_muthukrishna_CLI.R "Family_ID" 10