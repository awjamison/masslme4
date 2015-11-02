### tests for lmer2 function ###

# use sleep study data from lme4 package
utils::data("sleepstudy", package="lme4")
data <- sleepstudy

# define formula
formula <- "Reaction ~ Days + (1 | Subject)"

# define the optimization routine
controls <- list(lmerControl(optimizer='nloptwrap',
                             optCtrl=list(maxfun=1e9,
                             algorithm="NLOPT_LN_BOBYQA")),

                 lmerControl(optimizer="bobyqa",
                             optCtrl=list(maxfun=1e9)),

                 lmerControl(optimizer="Nelder_Mead",
                             optCtrl=list(maxfun=1e9)))

# expected output to two decimals
test_res <- matrix(c(251.41,10.47, 9.75,0.80, 25.79,13.02), nrow=2, ncol=3)


test_that("lme4.version", {
  # must be at least 1.1.7, which bundles the nloptr package
  expect_true(packageVersion("lme4") >= "1.1.7")
})


test_that("lmer2.formula", {
  res <- lmer2(formula=formula,
               data=data,
               lmerControls=controls,
               onlyCoefs=TRUE)
  expect_is(res, "matrix")
  expect_equal(ncol(res), 3)
  expect_equal(nrow(res), 2)
  expect_true(all(round(res,2) == test_res))
})


test_that("lmer2.onlyCoefs=FALSE", {
  res <- lmer2(formula=formula,
               data=data,
               lmerControls=controls,
               onlyCoefs=FALSE)
  expect_is(res, "list")
  expect_equal(length(res), 3)
  # 1: coef matrix
  expect_is(res[[1]], "matrix")
  expect_equal(ncol(res[[1]]), 3)
  expect_equal(nrow(res[[1]]), 2)
  expect_true(all(round(res[[1]],2) == test_res))
  # 2: lmerMod
  expect_is(res[[2]], "lmerMod")
  # 3: warnings
  expect_is(res[[3]], "list")
  expect_equal(length(res[[3]]), 0)
})


test_that("lmer2.start", {
  res <- lmer2(formula=formula,
               data=data,
               lmerControls=controls,
               onlyCoefs=TRUE,
               start=1)
  expect_is(res, "matrix")
  expect_equal(ncol(res), 3)
  expect_equal(nrow(res), 2)
  expect_true(all(round(res,2) == test_res))
})


test_that("lmer2.errorsToNA=TRUE", {
  # create an error from wrong input value (data as string)
  res <- lmer2(formula=formula,
               data="",
               lmerControls=controls,
               onlyCoefs=FALSE,
               errorsToNA=TRUE)
  expect_is(res, "list")
  expect_equal(length(res), 3)
  # 1: coef matrix
  expect_is(res[[1]], "matrix")
  expect_equal(ncol(res[[1]]), 3)
  expect_equal(nrow(res[[1]]), 2)
  expect_true(all(is.na(res[[1]])))
  expect_equal(colnames(res[[1]]), c("Estimate", "Std. Error", "t value"))
  expect_equal(rownames(res[[1]]), c("(Intercept)", "Days"))
  # 2: lmerMod
  expect_is(res[[2]], "NULL")
  # 3: warnings
  expect_is(res[[3]], "NULL")

  # test if it correctly handles formula with no intercept
  res <- lmer2(formula=paste(formula, "+ 0"),
               data="",
               lmerControls=controls,
               onlyCoefs=FALSE,
               errorsToNA=TRUE)
  expect_equal(nrow(res[[1]]), 1)
  expect_equal(rownames(res[[1]]), "Days")
})


# add a collinear variable
data[, "collinearVar"] <- data[, "Days"]


test_that("lmer2.warning.rankDeficient", {
  newFormula <- paste(formula, "+ collinearVar")

  # test that message is sent
  expect_warning(lmer2(formula=newFormula,
                       data=data,
                       lmerControls=controls,
                       onlyCoefs=FALSE))

  res <- lmer2(formula=newFormula,
               data=data,
               lmerControls=controls,
               onlyCoefs=FALSE)

  # 1: coef matrix
  # lmer will drop the collinear columns from the model, and these need to
  #   be replaced with NA in the coefficients
  expect_is(res[[1]], "matrix")
  expect_equal(ncol(res[[1]]), 3)
  expect_equal(nrow(res[[1]]), 3)
  expect_true(all(is.na(res[[1]][3,])))
  # 2: warnings
  expect_is(res[[3]], "list")
  expect_is(res[[3]][["model1"]][[1]], "warning")
})


test_that("lmer2.warning.scales", {
  # put fixed effects on very different scales to trigger a warning
  data[, "tooBig"] <- data[, "Days"] * 1e7
  newFormula <- "Reaction ~ tooBig + (1 | Subject)"

  expect_warning(lmer2(formula=newFormula,
                       data=data,
                       lmerControls=controls,
                       onlyCoefs=FALSE))

  res <- lmer2(formula=newFormula,
               data=data,
               lmerControls=controls,
               onlyCoefs=FALSE)

  # 1: coef matrix
  expect_is(res[[1]], "matrix")
  expect_equal(ncol(res[[1]]), 3)
  expect_equal(nrow(res[[1]]), 2)
  # 3: warnings
  expect_is(res[[3]], "list")
  expect_is(res[[3]][["model1"]][[1]], "warning")
})


test_that("lmer2.warningsToNA=TRUE", {
  # test using rank deficiency warning
  newFormula <- paste(formula, "+ collinearVar")

  res <- lmer2(formula=newFormula,
               data=data,
               lmerControls=controls,
               onlyCoefs=TRUE,
               warningsToNA=TRUE)

  expect_is(res, "matrix")
  expect_equal(ncol(res), 3)
  expect_equal(nrow(res), 3)
  expect_true(all(is.na(res)))
})


test_that("lmer2.iterate", {
  # trigger a warning so that the optimizer will keep cycling
  res <- lmer2(formula=formula,
               data=data,
               lmerControls=controls,
               onlyCoefs=FALSE,
               triggerConvWarning=TRUE)
  expect_is(res[[3]], "list")
  expect_equal(length(res[[3]]), 3)
  expect_equal(length(res[[3]][["model1"]]), 1)
  expect_is(res[[3]][["model1"]][[1]], "warning")
  expect_equal(length(res[[3]][["model3"]]), 1)
  expect_is(res[[3]][["model3"]][[1]], "warning")
})


test_that("lmer2.iterate.adaptiveStart=FALSE", {
  res <- lmer2(formula=formula,
               data=data,
               lmerControls=controls,
               onlyCoefs=FALSE,
               adaptiveStart=FALSE,
               triggerConvWarning=TRUE)
  expect_is(res[[3]], "list")
  expect_equal(length(res[[3]]), 3)
  expect_equal(length(res[[3]][["model1"]]), 1)
  expect_is(res[[3]][["model1"]][[1]], "warning")
  expect_equal(length(res[[3]][["model3"]]), 1)
  expect_is(res[[3]][["model3"]][[1]], "warning")
})


test_that("lmer2.argumentPassing", {
  # pass the defaults used in lmer explicitly (except control)
  res <- lmer2(formula=formula,
               data=data,
               lmerControls=controls,
               onlyCoefs=TRUE,
               REML=TRUE,
               verbose=0L,
               subset=NULL,
               weights=NULL,
               na.action=na.omit,
               offset=NULL,
               contrasts=NULL,
               devFunOnly=FALSE,)
  expect_is(res, "matrix")
  expect_equal(ncol(res), 3)
  expect_equal(nrow(res), 2)
})


# to execute from R console use test_that::test_dir, test_that::test_file