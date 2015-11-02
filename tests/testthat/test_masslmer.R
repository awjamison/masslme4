### tests for masslmer function ###

# use Arabidopsis data from lme4 package
utils::data("Arabidopsis", package="lme4")
data <- Arabidopsis

# define formula; 'amd' will be used as the index variable
formula <- "total.fruits ~ nutrient + (1 | reg)"


test_that("masslmer.multiprocessing=FALSE", {
  res <- masslmer(formula=formula,
                  data=data,
                  indexVars="amd",
                  multiprocessing=FALSE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})


test_that("masslmer.multiprocessing=TRUE", {
  res <- masslmer(formula=formula,
                  data=data,
                  indexVars="amd",
                  multiprocessing=TRUE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})


add.ids <- function(data, indexVars) {
  ## this is a taken from a block in masslmer
  # make sure this is a vector and not a list
  indexVars <- as.character(indexVars)
  # create index table of unique index addedData
  if(length(indexVars) > 1) {
    indexTable <- unique(data[, indexVars])
  } else {
    indexTable <- data.frame(unique(data[, indexVars]))
    # avoid data.frame's non-standard evaluation for names
    names(indexTable) <- indexVars
  }
  # use row number as ID
  indexTable[, "mergeVar"] <- 1:nrow(indexTable)
  # pre-allocate column for index var. "0" means row wasn't found in indexTable
  data[, "mergeVar"] <- 0
  # for i in indexTable, for j in data:
  # This seems ridiculously complicated. Why is there not a hash table object?
  invisible(
  sapply(1:nrow(indexTable), function(i) {
    sapply(1:nrow(data), function(j) {
      if(all(indexTable[i, indexVars] == data[j, indexVars])) {
        data[j, "mergeVar"] <<- indexTable[i, "mergeVar"]
      }
    })
  }))
  # check to make sure all rows have an ID
  if(0 %in% data["mergeVar"]) {
    stop("Found rows in the data that aren't in the index table.")
  }
  return(data)
}

## Split 'data' into two data frames: one with just the DV and ID and another
##   with all the other variables and ID. They will be remerged using ID.
otherVars <- setdiff(names(data), "total.fruits")
data <- add.ids(data, otherVars)
data2 <- data[c("total.fruits", "mergeVar")]
otherVars <- c(otherVars, "mergeVar")
addedData <- data[otherVars]


test_that("masslmer.multiprocessing=FALSE.withaddedData", {
  res <- masslmer(formula=formula,
                  data=data2,
                  indexVars="amd",
                  addedData=addedData,
                  multiprocessing=FALSE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})


test_that("masslmer.multiprocessing=TRUE.withaddedData", {
  res <- masslmer(formula=formula,
                  data=data2,
                  indexVars="amd",
                  addedData=addedData,
                  multiprocessing=TRUE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})


test_that("masslmer.moreDataThanVariables", {
  # add an extra row for an item that doesn't have an entry in addedData
  data.extra = rbind(data2, list(0, 0))
  # dropUnusedData=FALSE should produce an error
  expect_error(masslmer(formula=formula,
                        data=data.extra,
                        indexVars="amd",
                        addedData=addedData,
                        multiprocessing=FALSE,
                        dropUnusedData=FALSE, "Some 'data' entries"))
  # dropUnusedData=TRUE should avoid the error
  expect_is(masslmer(formula=formula,
                     data=data.extra,
                     indexVars="amd",
                     addedData=addedData,
                     multiprocessing=FALSE,
                     dropUnusedData=TRUE), "data.frame")
})


test_that("masslmer.multiprocessing=FALSE.withStarters", {
  res <- masslmer(formula=formula,
                  data=data,
                  indexVars="amd",
                  nStarters=1,
                  multiprocessing=FALSE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})


test_that("masslmer.multiprocessing=TRUE.withStarters", {
  res <- masslmer(formula=formula,
                  data=data,
                  indexVars="amd",
                  nStarters=1,
                  multiprocessing=TRUE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})


test_that("masslmer.multiprocessing=FALSE.argumentPassing", {
  # same as defaults in lmer2 (except onlyCoefs) but passed explicitly
  res <- masslmer(formula=formula,
                  data=data,
                  indexVars="amd",
                  multiprocessing=FALSE,
                  start=NULL,
                  adaptiveStart=TRUE,
                  lmerControls=NULL,
                  warningsToNA=FALSE,
                  errorsToNA=FALSE,
                  triggerConvWarning=FALSE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})


test_that("masslmer.multiprocessing=TRUE.argumentPassing", {
  # same as defaults in lmer2 (except onlyCoefs) but passed explicitly
  res <- masslmer(formula=formula,
                  data=data,
                  indexVars="amd",
                  multiprocessing=TRUE,
                  start=NULL,
                  adaptiveStart=TRUE,
                  lmerControls=NULL,
                  warningsToNA=FALSE,
                  errorsToNA=FALSE,
                  triggerConvWarning=FALSE)
  expect_is(res, "data.frame")
  # columns for: indexVars, IDvar, and coef, se, and t for fixed effects
  expect_equal(ncol(res), 8)
  # nrows == combined number of levels for all indexVars
  expect_equal(nrow(res), 2)
  expect_equal(names(res), c("amd", "IDvar",
                             "coef_(Intercept)", "coef_nutrient",
                             "se_(Intercept)", "se_nutrient",
                             "t_(Intercept)", "t_nutrient"))
  expect_true(all(is.numeric(res[, 'coef_(Intercept)'])))
  expect_false(any(is.na(res[, 'coef_(Intercept)'])))
  expect_true(all(is.numeric(res[, 'coef_nutrient'])))
  expect_false(any(is.na(res[, 'coef_nutrient'])))
})
