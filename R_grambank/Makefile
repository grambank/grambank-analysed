# Determine this makefile's path.
THIS_MAKEFILE := $(lastword $(MAKEFILE_LIST))

.PHONY: all clean data pca help pca_process
.DEFAULT_GOAL := help

## Make everything.
all:
	 quick_ones medium_length if_medium_are_finished long_ones if_long_ones_are_finished

# runs almost all analysis from the release paper, excluding the analysis that takes a long time and/or more computational resources than a regular domestic computer tends to have (the BRMS-analysis for predicting unusualness score is best run on a cluster and the INLA analysis of the spatiophylogenetic modelling takes a long time)
quick_ones: requirements.log data pca manhattan_dists

medium_length: get_unusualness make_cultural_fixation_scores

if_medium_are_finished: plot_unusualness make_cultural_fixation_scores_plot

long_ones: INLA get_unusualness predict_unusualness INLA_ASR INLA_simulation_binary_data INLA_predict_BRMS functional_richness_get

if_long_ones_are_finished: INLA_extract_posteriors INLA_compare_settings INLA_plots plot_predict_unusualness INLA_ASR_plot INLA_compare_settings INLA_simulation_binary_data_plot functional_richness_plot


## clean 
clean:
	rm -r output/non_GB_datasets*
	rm -r output/coverage_plots*
	rm -r output/GB_wide*
	rm -r output/PCA*
#	rm -r output/unusualness*
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
	@mkdir -p output/non_GB_datasets
	@mkdir -p output/coverage_plots
	@mkdir -p output/GB_wide
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
	Rscript coverage_world_map.R

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
	Rscript PCA/PCA.R

output/PCA/PCA_contributions_PCA1_PCA2.png: 
	Rscript PCA/PCA_plot_contributions.R

output/PCA/PCA_RGB_worldmap.png: 
	Rscript PCA/PCA_plot_worldmaps.R

output/PCA/PCA_macroarea.png: 
	Rscript PCA/PCA_plot_scatter_plots.R

## Comparing PCA to theoretical scores
theoretical:
	@echo "# Comparing PCA to theoretical scores. #"
	Rscript make_theo_score_fusion.R
	Rscript make_theo_scores.R
	Rscript PCA/theo_scores_compare_PCA_loadings.R

## Unusualness language analysis and generate plots.
get_unusualness: 
#	@echo "# Calculating unsualness scores and making plots. Note: these scripts requires internet access.#"
	Rscript unusualness/analysis/get_unusualness_bayesLCA.R
			@$(MAKE) -f $(THIS_MAKEFILE) plot_unusualness

## BRMS-analysis of unsualness scores
predict_unusualness: 
	@echo "## Running BRMS-analysis of unsualness scores (cluster recommended). #"
	Rscript unusualness/analysis/predict_surprisal.R
		@$(MAKE) -f $(THIS_MAKEFILE) plot_predict_unusualness

plot_unusualness:
	Rscript unusualness/analysis/plot_unusualness_splom.R
	Rscript unusualness/analysis/plot_unusualness_worldmap.R
	Rscript unusualness/analysis/plot_unusualness_ridgeplot.R

plot_predict_unusualness: 
	Rscript unusualness/analysis/predict_surprisal_plot.R
	
## Phylogenetic signal per feature.
#phylo_signal_per_feature:
#	@echo "## Running analysis of phylogenetic signal per feature. ##"
#	Rscript phylo_signal_features_d_stat.R

## Run endangerment and structural dissimilarity
#endangerment_analysis:
#	@echo "## Running analysis of endangerment and structural dissimilarity. ##"
#	Rscript diversity_endangerment.R

## INLA 
INLA: data INLA_featurewise INLA_extract_posteriors

INLA_featurewise:
	@echo "## Running INLA-analysis featurewise. #"
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 1 113 precision_matrices_kappa_2_sigma_1.15.RDS "default"

INLA_featurewise_all_spatial_params_all_priors: INLA_featurewise_spatial_params_2_1_15_all_priors INLA_featurewise_spatial_params_2_2_all_priors INLA_featurewise_spatial_params_2_5_3_all_priors

INLA_featurewise_spatial_params_2_1_15_all_priors:
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 1 113 precision_matrices_kappa_2_sigma_1.15.RDS "loop_all_priors" 

INLA_featurewise_spatial_params_2_2_all_priors:
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 1 113 precision_matrices_kappa_2_sigma_2.RDS "loop_all_priors" 

INLA_featurewise_spatial_params_2_5_3_all_priors_first_part:
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 1 37 precision_matrices_kappa_2.5_sigma_3.RDS "loop_all_priors" 

INLA_featurewise_spatial_params_2_5_3_all_priors_second_part:
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 38 75 precision_matrices_kappa_2.5_sigma_3.RDS "loop_all_priors" 

INLA_featurewise_spatial_params_2_5_3_all_priors_third_part:
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 75 100 precision_matrices_kappa_2.5_sigma_3.RDS "loop_all_priors" 

INLA_featurewise_spatial_params_2_5_3_all_priors_fourth_part:
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "real" 100 113 precision_matrices_kappa_2.5_sigma_3.RDS "loop_all_priors" 

INLA_extract_posteriors:
	Rscript spatiophylogenetic_modelling/analysis/functions/extract_INLA_posteriors.R "real" "quiet"

INLA_plots:
	Rscript spatiophylogenetic_modelling/analysis/featurewise/grambank_dualplots.R

INLA_simulation_binary_data:
	@$(MAKE) -f $(THIS_MAKEFILE) INLA_extract_posteriors_simulated
	Rscript spatiophylogenetic_modelling/analysis/INLA_multi_models.R "sim" 1 160 precision_matrices_kappa_2_sigma_1.15.RDS "default"

INLA_extract_posteriors_simulated:
	Rscript spatiophylogenetic_modelling/analysis/functions/extract_INLA_posteriors.R "sim" "quiet"

INLA_simulation_binary_data_plot:
	Rscript spatiophylogenetic_modelling/analysis/simulations/simulated_INLA_viz.R
	
INLA_compare_settings: 
	Rscript spatiophylogenetic_modelling/analysis/featurewise/summarise_modelprior.R
	Rscript spatiophylogenetic_modelling/analysis/featurewise/summarise_spatialparameters.R
	
INLA_ASR:
	Rscript spatiophylogenetic_modelling/analysis/featurewise/INLA_results_ASR_get.R
	
INLA_ASR_plot: 
	Rscript spatiophylogenetic_modelling/analysis/featurewise/INLA_results_ASR_plot.R

INLA_predict_BRMS:
	Rscript spatiophylogenetic_modelling/analysis/featurewise/predict_sp_effects_BRMS.R

#generate cultural fixation scores and do mds plots
make_cultural_fixation_scores:
	Rscript dists/fixation_scores_muthukrishna_run.R
	
make_cultural_fixation_scores_plot:
	Rscript dists/fixation_scores_muthukrishna_network_americas.R
	Rscript dists/fixation_scores_muthukrishna_barplot_heatmap.R

manhattan_dists:
	Rscript dists/manhattan_dists.R
	
functional_richness_get:
	Rscript functional_richness/make_gower_mfd.R

functional_richness_plot:
	Rscript functional_richness/frischness_plot.R

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

