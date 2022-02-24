#library(tidyverse)
library(readr)  # read_tsv
library(testthat)
library(tibble)  # column_to_rownames

source("./unusualness_get_lg_score.R")

testdata <- read_tsv(file.path("test_unusualness_get_lg_score.tsv")) %>% column_to_rownames("Language_ID")
# Language GB020    GB021    GB022    GB023
# A            1        0        1        1
# B            1        0        0        1
# C            1        1        0        0
#
#... which means (+ is in most common group, . is not)
# A            +        +        .        +    = +++
# B            +        +        +        +    = ++++
# C            +        .        +        .    = ++
#
# state probabilities are:
#
# .  GB020    GB021    GB022    GB023
# 0    0/3      2/3      2/3      1/3
# 1    3/3      1/3      1/3      2/3
#
# ... and therefore,
#       C has the rarest states:
expected_score_c <- c(3/3, 1/3, 2/3, 1/3)
#       A has a medium level of rare states:
expected_score_a <- c(3/3, 2/3, 1/3, 2/3)
#       B has the most common states:
expected_score_b <- c(3/3, 2/3, 2/3, 2/3)

# .. then we log it, and take the sum:
# A = -1.909543
# B = -1.216395
# C = -2.602690

# ... and B > A > C, or C is the most unusual, while B is least.


context("check test data")
test_that("check test data", {
    expect_equal(nrow(testdata), 3)
    expect_equal(ncol(testdata), 4)
})


context("get_col_prob")
test_that("test get_col_prob", {
    p.gb020 <- get_col_prob(testdata, "GB020")
    p.gb021 <- get_col_prob(testdata, "GB021")
    p.gb022 <- get_col_prob(testdata, "GB022")
    p.gb023 <- get_col_prob(testdata, "GB023")

    expect_equal(p.gb020, list(`1` = 3/3))
    expect_equal(p.gb021, list(`0` = 2/3, `1` = 1/3))
    expect_equal(p.gb022, list(`0` = 2/3, `1` = 1/3))
    expect_equal(p.gb023, list(`0` = 1/3, `1` = 2/3))
})


context("get_probabilities")
test_that("test get_probabilities", {
    probs <- get_probabilities(testdata)
    
    expect_equal(probs[['GB020']][['0']], NULL)
    expect_equal(probs[['GB020']][['1']], 3/3)

    expect_equal(probs[['GB021']][['0']], 2/3)
    expect_equal(probs[['GB021']][['1']], 1/3)

    expect_equal(probs[['GB022']][['0']], 2/3)
    expect_equal(probs[['GB022']][['1']], 1/3)

    expect_equal(probs[['GB023']][['0']], 1/3)
    expect_equal(probs[['GB023']][['1']], 2/3)
})


context("get_var_score")
test_that("test get_var_score", {
    probs <- get_probabilities(testdata)

    expect_equal(get_var_score("GB020", 0, probs), 0)
    expect_equal(get_var_score("GB020", 1, probs), 3/3)
    
    expect_equal(get_var_score("GB021", 0, probs), 2/3)
    expect_equal(get_var_score("GB021", 1, probs), 1/3)

    expect_equal(get_var_score("GB022", 0, probs), 2/3)
    expect_equal(get_var_score("GB022", 1, probs), 1/3)

    expect_equal(get_var_score("GB023", 0, probs), 1/3)
    expect_equal(get_var_score("GB023", 1, probs), 2/3)
    
    # characters rather than integers work correctly
    expect_equal(get_var_score("GB022", '0', probs), 2/3)
    expect_equal(get_var_score("GB022", '1', probs), 1/3)

})



context("get_lang_score")
test_that("test get_lang_score", {
    probs <- get_probabilities(testdata)
    
    a <- get_lang_score('A', testdata, probs, raw=TRUE)
    expect_equal(as.vector(a), expected_score_a)

    b <- get_lang_score('B', testdata, probs, raw=TRUE)
    expect_equal(as.vector(b), expected_score_b)

    c <- get_lang_score('C', testdata, probs, raw=TRUE)
    expect_equal(as.vector(c), expected_score_c)

    # summary
    a <- get_lang_score('A', testdata, probs)
    expect_equal(a, sum(log(expected_score_a)))
    b <- get_lang_score('B', testdata, probs)
    expect_equal(b, sum(log(expected_score_b)))
    c <- get_lang_score('C', testdata, probs)
    expect_equal(c, sum(log(expected_score_c)))
    
    expect_gt(b, a)
    expect_gt(b, c)
    expect_gt(a, c)
})

