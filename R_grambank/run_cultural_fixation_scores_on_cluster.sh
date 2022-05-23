#!/bin/bash
#SBATCH --cpus-per-task 15
#SBATCH --mem 20G
#SBATCH -J Hedvig_cultural_fixation-scores

##1w:code skirgard$ chmod 755 run_INLA_featurewise_on_cluster.sh
##lingn01w:code skirgard$ ./run_INLA_featurewise_on_cluster.sh

#Step 1
FILE=rlibs
if [ -f "$FILE" ]; then
    echo "$FILE exists."
    echo "rlibs folder already exists, won't bother remaking it."
else 
mkdir rlibs
fi

export R_LIBS_USER=rlibs

Rscript cluster_set_up.R
Rscript dist_fixation_scores/fixation_scores_muthukrishna.R



