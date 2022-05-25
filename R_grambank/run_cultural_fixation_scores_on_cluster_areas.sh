#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J Hedvig_cultural_fixation-scores_areas
#SBATCH --partition=dlcegpu

export R_LIBS_USER=../rlibs

Rscript dist_fixation_scores/fixation_scores_muthukrishna_CLI.R "Macroarea" 0

Rscript dist_fixation_scores/fixation_scores_muthukrishna_CLI.R "AUTOTYP-area" 0

Rscript dist_fixation_scores/fixation_scores_muthukrishna_CLI.R "AUTOTYP-area" 20

Rscript dist_fixation_scores/fixation_scores_muthukrishna_CLI.R "AUTOTYP-area" 50

Rscript dist_fixation_scores/fixation_scores_muthukrishna_CLI.R "AUTOTYP-area" 100

