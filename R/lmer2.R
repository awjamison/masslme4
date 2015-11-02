# lme4 >1.1.8 includes the nloptr optimizers
# lmerControl: rpackages.ianhowson.com/cran/lme4/man/lmerControl.html

# print warnings as they occur
# max out the number of warnings that can be printed (still limited)
# increase the cap on number of prints to near-infinite
options(warn=1, warning.length=8170, max.print=1e9)

#' Wrapper for lme4's \code{lmer} Function
#'
#' Fit linear mixed-effects models using \code{\link[lme4]{lmer}} with access
#' to additional features that support automated refitting and condition
#' handling.
#'
#' @import lme4
#'
#' @param formula same as \code{\link[lme4]{lmer}} \code{formula}.
#' @param data same as \code{\link[lme4]{lmer}} \code{data}.
#' @param onlyCoefs logical: whether to only return model coefficients.
#' @param start same as \code{\link[lme4]{lmer}} \code{start}.
#' @param adaptiveStart logical: whether to use ending theta values of the
#' previous iteration as \code{start} values for the current iteration.
#' @param lmerControls \code{\link[base]{list}} of
#' \code{\link[lme4]{lmerControl}} functions to try when fitting the model.
#' @param warningsToNA logical: whether to return a matrix of NAs instead of
#' coefficients if a warning is raised when fitting the model.
#' @param errorsToNA logical: whether to return a matrix of NAs (with the same
#' shape as a coefficient matrix) if an error is raised when fitting the model.
#' The error is converted to a warning instead of unwinding the stack.
#' @param triggerConvWarning logical: raises a convergence warning.
#' @param ... other potential arguments passed to \code{\link[lme4]{lmer}}.
#'
#' @return If \code{onlyCoef==TRUE} then the output will be a matrix of betas,
#' standard errors, and t-values. Otherwise, the ouput will be a list that
#' contains this matrix as well as the standard ouput of
#' \code{\link[lme4]{lmer}} (an object of class \code{lmerMod}) and a list of
#' any messages, warnings, and errors generated when fitting the model.
#'
#' @examples
#' library(lme4)
#' data(sleepstudy)
#'
#' ## try these optimizers in order
#' controls <- list(lmerControl(optimizer='nloptwrap',
#'                              optCtrl=list(maxfun=1e9,
#'                              algorithm="NLOPT_LN_BOBYQA")),
#'
#'                  lmerControl(optimizer="bobyqa",
#'                              optCtrl=list(maxfun=1e9)),
#'
#'                  lmerControl(optimizer="Nelder_Mead",
#'                              optCtrl=list(maxfun=1e9)))
#'
#' ## fit the model
#' coefs <- lmer2(Reaction ~ Days + (Days | Subject), sleepstudy, TRUE,
#'                lmerControls=controls)
#'
#' @export
lmer2 = function(formula, data, onlyCoefs=FALSE, start=NULL, adaptiveStart=TRUE,
                 lmerControls=NULL, warningsToNA=FALSE, errorsToNA=FALSE,
                 triggerConvWarning=FALSE, ...) {


  # display a message about the number of times the model has been refit
  refitMessage = function(iteration) {
    message(paste("Model finished with convergence warnings. Refitting the",
                  "model with a different optimizer.",
                  paste("(", iteration, ")", sep="")))
  }

  if (is.null(lmerControls)) {
    # default optimizer if no optimizer routine is supplied
    lmerControls <- list(lmerControl(optimizer='nloptwrap',
                                     optCtrl=list(maxfun=1e9,
                                     algorithm="NLOPT_LN_BOBYQA")))
  }

  # prefixes of possible convergence warnings from lme4 as printed by
  #   checkConv or checkHess (checkConv calls checkHess). Right now all of them
  #   will invoke an attempt to refit the model, but this may not be useful.
  #   see github.com/lme4/lme4/blob/master/R/checkConv.R
  convWarnings <- list(checkConv=c("Model failed to converge",
                                   "unable to evaluate scaled gradient"),
                       checkHess=c("Model failed to converge",
                                   "Problem with Hessian check",
                                   "Hessian is numerically singular",
                                   "Model is nearly unidentifiable"))

  make.emptyModel = function(formula) {
    # returns coefficient matrix of NAs, NULL model, and NULL conditions
    # Note: formula is coerced earlier to ensure the error does not originate
    #   from the formula. If it somehow does, this will crash.
    formulaTerms <- terms(formula, simplify=TRUE)
    # extract variables from formula
    vars <- attr(formulaTerms, "term.labels")
    # remove random effects
    vars <- sapply(vars, function(x) if (! grepl(("|"), x, fixed=TRUE)) x)
    vars <- vars[! sapply(vars, is.null)]
    if (attr(formulaTerms, "intercept") == 1) {
      # formula contains an intercept
      vars <- append("(Intercept)", vars)
    }
    # create NA matrix analog to coef(summary(model))
    dummyCoefs <- matrix(NA_real_, length(vars), 3)
    rownames(dummyCoefs) <- vars
    colnames(dummyCoefs) <- c("Estimate", "Std. Error", "t value")
    if (onlyCoefs) {
      return(dummyCoefs)
    } else {
      return(list(coefficients=dummyCoefs, model=NULL, conditions=NULL))
    }
  }

  # store messages, errors, and warnings produced during fitting
  conditions <- list()

  # coerce the input to a formula if not already
  formula <- as.formula(formula)

  # on the first iteration model, is a formula
  model <- formula

  # fitting function
  fit.model <- function(object) {
    if (triggerConvWarning) {
      # used to test the warning handler
      warning("Model failed to converge")
    }
    if(is(object, "formula")) {
      # fit the model for the first time
      lmer(object, data, start=start, control=control, ...)
    } else if(is(object, "lmerMod")) {
      # update the existing model with a new control parameter
      if(adaptiveStart) {
        # start the optimizer with the theta values of the previous iteration
        update(object, start=start, control=control)
      } else {
        refit(object, control=control)
      }
    } else {
      stop("fit.model must take an object of class formula or lmerMod.")
    }
  }

  # count number of iterations
  n <- 1
  # try fitting the model while true
  fit <- TRUE
  while(fit) {
    fit <- FALSE

    # will be updated with conditions triggered during this iteration
    activeConds <- NULL
    activeWarnings <- FALSE
    activeConvWarnings <- FALSE
    activeErrors <- FALSE

    # the lmerControl to use on this iteration
    control <- lmerControls[[n]]

    # promote rank deficiency message to a warning
    control[['checkControl']][['check.rankX']] <- "warn+drop.cols"

    # operate on messages and warnings without unwinding the stack but
    #   allow for the possibility of catching errors
    model <- withCallingHandlers({
          if (errorsToNA) {
            # consistent output structure but no meaningful values
            tryCatch(fit.model(model),
                     error = function(e) {
                       activeErrors <<- TRUE
                       # display the error demoted to a warning
                       warning(as.character(e))
                     })
          } else {
            fit.model(model)
          }
      }, message = function(m) {
          activeConds <<- c(activeConds, list(m))
      }, warning = function(w) {
          activeConds <<- c(activeConds, list(w))
          activeWarnings <<- TRUE
          # scan the warning message for convergence warnings
          wMsg <- w["message"]
          if (any(sapply(unlist(convWarnings),
                         function(x) grepl(x, wMsg, fixed=TRUE)))) {
            activeConvWarnings <<- TRUE
            # got convergence warning: decide whether to refit the model
            if (n < length(lmerControls)) {
              # try refitting the model on the next iteration
              fit <<- TRUE
              # (else) exit the loop; all lmerControl objects exhausted
            }
          }
      })

    # indicate that assignment to model triggered warnings
    if (activeConvWarnings) {
      if (n < length(lmerControls)) {
        refitMessage(n)
      } else {
        warning(paste("All attempts at model convergence failed.",
                      "Parameter estimates may be unreliable."))
      }
    }

    # return if errors are being handled and there was an error
    if (activeErrors) {
      return(make.emptyModel(formula))
    }

    # update with the conditions with those from this iteration
    conditions[[paste("model", n, sep="")]] <- activeConds

    # use theta of previous model as starting value when refitting
    # (this will be ignored if adaptiveStart == FALSE)
    start <- getME(model, "theta")

    n <- n + 1
  }

  # get model coefficients
  coefs <- coef(summary(model))

  if (warningsToNA) {
    if (activeWarnings) {
      # return NAs (probably from convergence warnings)
      coefs[, ] <- NA_real_
    }
  }

  # if variables were dropped from the model (e.g. not full rank) then create
  #   NA entries in coefs for those variables.
  # TODO: maintain the order of the variables in the formula
  # TODO: account for R expanding factors
  startTerms <- attr(terms(model, simplify=TRUE, fixed.only=TRUE),
                     'term.labels')
  diffTerms <- setdiff(startTerms, rownames(coefs))
  if (length(diffTerms) > 0) {
    diffM <- matrix(NA_real_, length(diffTerms), 3)
    rownames(diffM) <- diffTerms
    coefs <- rbind(coefs, diffM)
  }

  if (onlyCoefs) {
    # return coefficient, std. error, and t-value
    return(coefs)
  } else {
    # return coefficient, std. error, and t-value; model object; and conditions
    return(list(coefficients=coefs, model=model, conditions=conditions))
  }
}
