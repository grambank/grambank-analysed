# Determine this makefile's path.
THIS_MAKEFILE := $(lastword $(MAKEFILE_LIST))

.PHONY: all clean data pca help pca_process
.DEFAULT_GOAL := help

## Make everything.
all:
	 requirements.log data pca get_unusualness predict_unusualness sp_models sensitivity_testing phylo_signal_per_feature endangerment_analysis test

# runs almost all analysis from the release paper, excluding the analysis that takes a long time and/or more computational resources than a regular domestic computer tends to have (the BRMS-analysis for predicting unusualness score is best run on a cluster and the INLA analysis of the spatiophylogenetic modelling takes a long time)
almost_all_fast: requirements.log data pca get_unusualness phylo_signal_per_feature endangerment_analysis

## clean 
clean:
	rm -r output/non_GB_datasets*
	rm -r output/coverage_plots*
	rm -r output/GB_wide*
	rm -r output/PCA*
	rm -r output/unusualness*
#	rm -r output/spatiophylogenetic_modelling*

#step 0
## Check and install dependencies if needed.
requirements.log:
	@echo "# Installing and loading all the necessary R-packages. #"
	Rscript requirements.R | tee output/requirements.log
 
#step 1
## Fetch data from grambank-cldf into required formats for plotting.
submodules:
	git submodule update --init

data: submodules
	@echo "# Preparing the dataset (including imputation and pruning of global tree. #"
	@mkdir output/non_GB_datasets
	@mkdir output/coverage_plots
	@mkdir output/GB_wide
	@$(MAKE) -f $(THIS_MAKEFILE) output/non_GB_datasets/glottolog-cldf_wide_df.tsv
	@$(MAKE) -f $(THIS_MAKEFILE) output/GB_wide/GB_wide_strict.tsv
	@$(MAKE) -f $(THIS_MAKEFILE) output/GB_wide/GB_wide_binarized.tsv
	@$(MAKE) -f $(THIS_MAKEFILE) output/GB_wide/GB_wide_imputed_binarized.tsv
	@$(MAKE) -f $(THIS_MAKEFILE) output/coverage_plots/WALS_GB_coverage_overlay.png
	@$(MAKE) -f $(THIS_MAKEFILE) output/coverage_plots/coverage_macroarea.png
	@$(MAKE) -f $(THIS_MAKEFILE) output/pruning_EDGE_tree.R
	Rscript unusualness/processing/assigning_AUTOTYP_areas.R

output/non_GB_datasets/glottolog-cldf_wide_df.tsv:
	@echo Fetching Glottolog language table. Note: this script requires the git submodule glottolog-cldf.
	Rscript make_glottolog-cldf_table.R

output/GB_wide/GB_wide_strict.tsv:
	@echo Making the GB data wide and merging dialects
	Rscript make_wide.R

output/GB_wide/GB_wide_binarized.tsv:
	@echo Binarising the GB data.
	Rscript make_wide_binarized.R

output/GB_wide/GB_wide_imputed_binarized.tsv:
	@echo Imputing missing values. This will most likely take a few minutes.
	Rscript impute_missing_values.R | tee impute_missing_values.log
	
output/coverage_plots/WALS_GB_coverage_overlay.png:
	@echo Making plot for comparing GB and WALS. Note: this script requires the git submodule wals
	Rscript compare_coverage_WALS.R

output/coverage_plots/coverage_macroarea.png:
	@echo Making plots for coverage per macroarea.
	Rscript coverage_bar_plots.R

output/pruning_EDGE_tree.R:
	mkdir -p output/spatiophylogenetic_modelling
	mkdir -p output/spatiophylogenetic_modelling/processed_data
	RScript spatiophylogenetic_modelling/processing/pruning_EDGE_tree.R

#step 2
## PCA analysis and plots.
pca:
	@echo "## Running PCA analysis and generating plots. ##"
	mkdir -p output/PCA
	@$(MAKE) -f $(THIS_MAKEFILE) output/PCA/PCA_df.tsv
	@$(MAKE) -f $(THIS_MAKEFILE) output/PCA/PCA_contributions_PCA1_PCA2.png
	@$(MAKE) -f $(THIS_MAKEFILE) output/PCA/PCA_RGB_worldmap.png
	@$(MAKE) -f $(THIS_MAKEFILE) output/PCA/PCA_macroarea.png
	@$(MAKE) -f $(THIS_MAKEFILE) theoretical

output/PCA/PCA_df.tsv:
	Rscript PCA.R

output/PCA/PCA_contributions_PCA1_PCA2.png: 
	Rscript PCA_plot_contributions.R

output/PCA/PCA_RGB_worldmap.png: 
	Rscript PCA_plot_worldmaps.R

output/PCA/PCA_macroarea.png: 
	Rscript PCA_plot_scatter_plots.R

## Comparing PCA to theoretical scores
theoretical:
	@echo "# Comparing PCA to theoretical scores. #"
	Rscript bound_morph.R
	Rscript compare_PCA_loadings_to_theory_scores.R

#step 3
## Unusualness language analysis and generate plots.
get_unusualness:
#	@echo "# Calculating unsualness scores and making plots. Note: these scripts requires internet access.#"

#step 4
## Phylogenetic signal per feature.
phylo_signal_per_feature:
	@echo "## Running analysis of phylogenetic signal per feature. ##"
	Rscript phylo_signal_features_d_stat.R

#step 5
## Run endangerment and structural dissimilarity
endangerment_analysis:
	@echo "## Running analysis of endangerment and structural dissimilarity. ##"
	Rscript diversity_endangerment.R

#step 6
## INLA 
INLA: data INLA_featurewise

INLA_featurewise:
	@echo "## Running INLA-analysis featurewise. #"
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 1 113

INLA_simulation_binary_data:
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "sim" 1 160

#step 7
## BRMS-analysis of unsualness scores
predict_unusualness: 
	@echo "## Running BRMS-analysis of unsualness scores (cluster recommended). #"
	Rscript unusualness/analysis/predict_unusualness.R
	

## Test functions
test:
	Rscript -e 'testthat::test_dir("./tests", reporter="summary")'

make_author_list:
	python3 wrangle_authors/wrangle_author_list.py wrangle_authors/author_list.csv > wrangle_authors/author_list.html

## Display help
help:
	@echo "$$(tput bold)Available commands:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) == Darwin && echo '--no-init --raw-control-chars')

