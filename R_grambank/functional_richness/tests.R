library(testthat)

# dput(gb[
#     gb$Family_ID %in% c('nucl1708', 'koia1260', 'sepi1257', 'east2433'),
#     c("Language_ID", 'Family_ID', 'GB020', "GB021", "GB022", "GB023", "GB028", 'AUTOTYP_area')
# ])

gb_sample <- structure(list(
    Language_ID = c("abua1245", "auuu1241", "koit1244", "kona1242", "kwom1262",
                    "mehe1243", "moun1252", "nucl1630", "oloo1241",
                    "omie1241", "urim1252", "yapu1240"),
    Family_ID = c("nucl1708", "nucl1708", "koia1260", "east2433", "sepi1257",
                  "sepi1257", "koia1260", "koia1260", "nucl1708", "koia1260",
                  "nucl1708", "nucl1708"),
    GB020 = c(0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L),
    GB021 = c(1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L),
    GB022 = c(1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L),
    GB023 = c(0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L),
    GB028 = c(0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 1L, 0L, 1L, 1L),
    AUTOTYP_area = c("N Coast New Guinea", "N Coast New Guinea",
                     "S New Guinea", "Interior New Guinea", "N Coast New Guinea",
                     "N Coast New Guinea", "S New Guinea", "S New Guinea", "N Coast New Guinea",
                     "S New Guinea", "N Coast New Guinea", "N Coast New Guinea"
    )), row.names = c(3L, 86L, 607L, 618L, 662L, 826L, 869L, 983L, 1011L, 1013L, 1369L, 1467L),
    class = "data.frame")

#
# Language_ID Family_ID       AUTOTYP_area
# 3       abua1245  nucl1708 N Coast New Guinea
# 86      auuu1241  nucl1708 N Coast New Guinea
# 1011    oloo1241  nucl1708 N Coast New Guinea
# 1369    urim1252  nucl1708 N Coast New Guinea
# 1467    yapu1240  nucl1708 N Coast New Guinea
#
# 662    kwom1262  sepi1257 N Coast New Guinea
# 826    mehe1243  sepi1257 N Coast New Guinea
#
# 607     koit1244  koia1260 S New Guinea
# 869     moun1252  koia1260 S New Guinea
# 983     nucl1630  koia1260 S New Guinea
# 1013    omie1241  koia1260 S New Guinea
#
# 618    kona1242  east2433 Interior New Guinea



source("functional_richness/fun_def_make_group_matrix.R")

test_that("make_group_matrix - full with area", {
    mtx <- make_group_matrix(gb_sample, 'AUTOTYP_area', remove=c(), threshold=1)
    expect_equal(rownames(mtx), c("Interior New Guinea", "N Coast New Guinea", "S New Guinea"))

    expect_named(
        which(mtx['Interior New Guinea', ] == 1),
        gb_sample[gb_sample$AUTOTYP_area == 'Interior New Guinea', 'Language_ID']
    )
    expect_named(
        which(mtx['N Coast New Guinea', ] == 1),
        gb_sample[gb_sample$AUTOTYP_area == 'N Coast New Guinea', 'Language_ID']
    )
    expect_named(
        which(mtx['S New Guinea', ] == 1),
        gb_sample[gb_sample$AUTOTYP_area == 'S New Guinea', 'Language_ID']
    )
})



test_that("make_group_matrix - full with family", {
    mtx <- make_group_matrix(gb_sample, 'Family_ID', remove=c(), threshold=1)
    expect_equal(rownames(mtx), c("east2433", "koia1260", "nucl1708", "sepi1257"))

    expect_named(
        which(mtx['nucl1708', ] == 1),
        gb_sample[gb_sample$Family_ID == 'nucl1708', 'Language_ID']
    )
    expect_named(
        which(mtx['koia1260', ] == 1),
        gb_sample[gb_sample$Family_ID == 'koia1260', 'Language_ID']
    )
    expect_named(
        which(mtx['east2433', ] == 1),
        gb_sample[gb_sample$Family_ID == 'east2433', 'Language_ID']
    )
    expect_named(
        which(mtx['sepi1257', ] == 1),
        gb_sample[gb_sample$Family_ID == 'sepi1257', 'Language_ID']
    )
})



test_that("make_group_matrix - subsample with area", {
    remove <- c(
        "kona1242", # ING
        "koit1244", "moun1252", # SNG
        "kwom1262", "mehe1243", "oloo1241" # NCNG
    )
    expected <- list(
        #'Interior New Guinea' =  # Nothing -- have to do this differently
        'N Coast New Guinea' = c('abua1245', 'auuu1241', 'urim1252', 'yapu1240'),
        'S New Guinea' = c('nucl1630', 'omie1241')
    )

    mtx <- make_group_matrix(gb_sample, 'AUTOTYP_area', remove=remove, threshold=1)
    expect_equal(rownames(mtx), c("Interior New Guinea", "N Coast New Guinea", "S New Guinea"))

    # check names
    expect_named(
        which(mtx['N Coast New Guinea', ] == 1), expected[['N Coast New Guinea']]
    )
    expect_named(
        which(mtx['S New Guinea', ] == 1), expected[['S New Guinea']]
    )
    # empty things are tested differently.
    expect_equal(length(which(mtx['Interior New Guinea', ] == 1)), 0)
})




test_that("make_group_matrix - subsample with family", {
    remove <- c(
        "kona1242", # east2433
        "koit1244", "moun1252", # koia1260
        "kwom1262", # sepi1257
        "oloo1241" # nucl1708
    )

    expected <- list(
        'nucl1708' = c('abua1245', 'auuu1241', 'urim1252', 'yapu1240'),
        'sepi1257' = c('mehe1243'),
        'koia1260' = c('nucl1630', 'omie1241'),
        'east2433' = c()
    )
    mtx <- make_group_matrix(gb_sample, 'Family_ID', remove=remove, threshold=1)

    expect_equal(rownames(mtx), c("east2433", "koia1260", "nucl1708", "sepi1257"))

    # check names
    expect_named(which(mtx['nucl1708', ] == 1), expected[['nucl1708']])
    expect_named(which(mtx['koia1260', ] == 1), expected[['koia1260']])
    expect_named(which(mtx['sepi1257', ] == 1), expected[['sepi1257']])
    # empty things are tested differently.
    expect_equal(length(which(mtx['east2433', ] == 1)), 0)

})
