# core functions for `rare languages` analysis.
# extracted here so we can run tests against it.

#script written by Simon Greenhill and Luke Maurits.

get_col_prob <- function(df, col) {
    counts <- table(as.factor(df[[col]]))
    if(length(counts) == 1) {
      if(names(counts) == "0") {
	      counts[["1"]] <- 0
      } else {
	      counts[["0"]] <- 0
      }
    }
    lapply(as.list(counts), function (x) x / sum(counts))
}

get_probabilities <- function(df) {
    sapply(colnames(df), function (x) get_col_prob(df, x), simplify=FALSE)
}

get_family_probs <- function(gb_data, language_metadata) {
	all_families <- language_metadata %>% pull("Family_group") %>% unique
	family_probs <- list()
	for(family in all_families) {
		langs_in_fam <- language_metadata %>% filter(Family_group == family) %>% pull(ID)
		family_probs[[family]] <- gb_data[rownames(gb_data) %in% langs_in_fam,] %>% get_probabilities
	}
	return(family_probs)
}

get_area_probs <- function(gb_data, language_metadata) {
	all_areas <- language_metadata %>% pull("AUTOTYP_area") %>% unique
	area_probs <- list()
	for(area in all_areas) {
#		print(area)
		langs_in_area <- language_metadata %>% filter(AUTOTYP_area == area) %>% pull(ID)
		area_probs[[area]] <- gb_data[rownames(gb_data) %in% langs_in_area,] %>% get_probabilities
	}
	return(area_probs)
}

#probs[['GB111']][['0']]
#probs[['GB111']][['1']]

get_var_score <- function(var, state, probabilities) {
    p <- probabilities[[var]][[as.character(state)]]
    ifelse(is.null(p), 0.0, p)  # return 0.0 if state is not seen in data
}

get_lang_score <- function(lang, gb_data, language_metadata, family_probabilities, area_probabilities, weights, raw=FALSE) {
    # Validate weights
    stopifnot(length(weights) == 2)
    stopifnot(sum(weights) == 1.0)
    # Pick out relevant family and area probabilities
    my_family <- language_metadata %>% filter(ID == lang) %>% pull(Family_group)
    my_area <- language_metadata %>% filter(ID == lang) %>% pull(AUTOTYP_area)
    my_family_probs <- family_probabilities[[my_family]]
    my_area_probs <- area_probabilities[[my_area]]
    # Family probabilities for Singletons (isolates and lgs which are the only ones repped in their family) don't make sense, so just use the areal probs
    if(my_family == "Singleton") {
	    my_family_probs <- my_area_probs
    }
    # Weight the two
    # No doubt this can be done better by someone with real R-fu...
    my_weighted_probs <- list()
    for(var in colnames(gb_data)) {
	    my_weighted_probs[[var]][["0"]] <- weights[1]*my_family_probs[[var]][["0"]] +
		                      weights[2]*my_area_probs[[var]][["0"]]
	    my_weighted_probs[[var]][["1"]] <- 1 - my_weighted_probs[[var]][["0"]]
    }
    # Use weighted probs to compute data log likelihood
    lang_score <- sapply(
        colnames(gb_data), function (var) get_var_score(var, gb_data[lang, var], my_weighted_probs)
    )
    if (raw) { return(lang_score) }
    return(sum(log(lang_score)))
}

## OLD code
# this is simons original code, which is the same as above except for not controlling for families or area
# core functions for `rare languages` analysis.
# extracted here so we can run tests against it.

get_col_prob_simons_original <- function(df, col) {
  counts <- table(as.factor(df[[col]]))
  lapply(as.list(counts), function (x) x / sum(counts))
}

get_probabilities_simons_original <- function(df) {
  sapply(colnames(df), function (x) get_col_prob(df, x), simplify=FALSE)
}

get_var_score_simons_original <- function(var, state, probabilities) {
  p <- probabilities[[var]][[as.character(state)]]
  ifelse(is.null(p), 0.0, p)  # return 0.0 if state is not seen in data
}

get_lang_score_simons_original <- function(lang, df, probabilities, raw=FALSE) {
  lang_score <- sapply(
    colnames(df), function (var) get_var_score(var, df[lang, var], probabilities)
  )
  if (raw) { return(lang_score) }
  return(sum(log(lang_score)))
}
