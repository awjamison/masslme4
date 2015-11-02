# it's highly recommended to use the Revolution Analytics version of R for
#   speed (has a better BLAS than native R) and fault tolerance

# print warnings as they occur
# max out the number of warnings that can be printed (still limited)
# increase the cap on number of prints to near-infinite
options(warn=1, warning.length=8170, max.print=1e9)

# JIT compile the code
# TODO: figure out how to suppress the CMD check notes. Or are these cat
#   printed from the compiler package? setCompilerOptions doesn't do anything.
compiler::enableJIT(3)

#' Fit lme4 \code{lmer} Models in Parallel
#'
#' The intended use is to support mass-univarate analysis (frequently used in
#' neuroimaging research) of a large dataset, although other applications
#' are possible. A model is iteratively fit to logically related subsets of
#' the data (e.g. voxels, sensors) and the combined statistics of each dataset
#' are reported. Subsets are determined by index variables ('indexVars'):
#' factors whose levels will be used to split the data.
#'
#' @importFrom foreach %dopar%
#'
#' @param formula same as \code{\link[lme4]{lmer}} \code{formula}.
#' @param data same as \code{\link[lme4]{lmer}} \code{data}.
#' @param indexVars character: fixed effects factors that partition the data
#'   into logically related subsets. Unique combinations of the levels will be
#'   used to split the data, and the subsets will be processed in parallel.
#' @param addedData NULL | data frame: variables to combine with the data.
#'   \code{data} and \code{addedData} must have at least one variable in common,
#'   and will be merged using all shared variables. Merging is done on the data
#'   subsets (not \code{data} itself) at the time they are sent to worker
#'   processes. In repeated measures designs, this may reduce the memory
#'   footprint from having to store all variables in a single data frame
#'   \code{data}.
#' @param dropUnusedData logical: whether to drop rows in \code{data} that
#'   don't have an entry in \code{addedData}.
#' @param saveTo character: path to save the results.
#' @param start same as \code{\link[lme4]{lmer}} \code{start}.
#' @param nStarters integer: number of models to fit in order to get starting
#'   theta values.
#' @param multiprocessing logical: whether to parallelize fitting models.
#' @param maxProcesses integer: maximum number of worker processes to create.
#'   Defaults to the number of CPUs on the host machine.
#' @param logging NULL | <file_path> | 'console': whether to generate logs,
#'   and if so whether to send the output to the console or save it in a file.
#' @param ... other potential arguments passed to \code{\link[masslmer]{lmer2}}.
#'
#' @return data frame containing statistics for each of the index variables.
#'   Currently, only the beta value, standard error, and t-value are reported,
#'   but a future version will allow customizing the output.
#'
#' @examples
#' ## create a model:
#' # people as random intercepts, conditions as fixed effects
#' # fit model separately for different sensors and times
#' # sensor s_k1 ~ N(k1, 1) for k1 = 1,2,3
#' # time t_k2 ~ N(k2, 1) for k2 = 1,2
#' # condition x1 = s*t, x2 = 2*x1
#' # Y_ij = B0 + B1*x_ij + u_j + e_ij for (B0 + u_j) = -2,-1,,0,1,2
#'
#' ## create some fake data
#' # five people, two conditions, two sensors, three times
#' data <- sapply(1:5, function(p) {
#'   p + sapply(1:2, function(x) {
#'   x * sapply(1:2, function(s) {
#'   rnorm(50, s) * sapply(1:3, function(t) {
#'   rnorm(50, t)})})})})
#'
#' data <- data.frame(y = sapply(data, c),
#'                    p = rep(c("A","B","C","D","E"), each = 600),
#'                    x = rep(c("A","B"), each = 300),
#'                    s = rep(c(1,2), each = 150),
#'                    t = rep(c(1,2,3), each = 50))
#'
#' ## fit models in parallel with two worker processes
#' # p is a random effect, x is a fixed effect, and s and t are index variables
#' model <- "y ~ x + (1 | p)"
#' results <- masslmer(model, data, c("s", "t"), maxProcesses = 2)
#'
#' @export
masslmer <- function(formula, data, indexVars, addedData=NULL,
                     dropUnusedData=FALSE, saveTo=NULL, start=NULL,
                     nStarters=NULL, multiprocessing=TRUE, maxProcesses=NULL,
                     logging=NULL, ...) {

  ## make sure "data" is a data frame
  data <- as.data.frame(data)
  # make sure this is a vector and not a list
  indexVars <- as.character(indexVars)

  ## check for errors when merging 'data' and 'addedData'
  if(! is.null(addedData)) {
    # make sure "addedData" is a data frame
    addedData <- as.data.frame(addedData)
    # find common item identifiers
    ids <- intersect(names(data), names(addedData))
    # check that data and addedData have at least one column to merge by
    if(length(ids) < 1) {
      stop(paste("'Data' and 'addedData' must have at least one common column",
                 "to merge by."))
    }
    if(! dropUnusedData) {
      # check that each data entry has a corresponding addedData entry
      dCols <- unique(data[, ids])
      vCols <- unique(addedData[, ids])
      if(length(setdiff(dCols, vCols)) > 0) {
        # merge causes data entries to be dropped
        stop(paste("Some 'data' entries do not have a corresponding",
                   "'addedData' entry. If you want to allow these data entries",
                   "to be dropped, set dropUnusedData = TRUE."))
      }
    }
    if(! all(indexVars %in% names(data))) {
      # merge 'indexVars' now because they will be used for splitting data
      if(! all(indexVars %in% c(names(data), names(addedData)))) {
        stop("indexVars could not be found in either 'data' or 'addedData'.")
      }
      # add indexVars not already in 'data'
      indexVarsToAdd <- setdiff(indexVars, names(data))
      data <- merge(data, addedData[c(indexVarsToAdd, ids)])
    }
  }

  ## split the data into a list of data frames that have unique combinations of
  ##   the index addedData
  # create index table of unique index addedData
  if(length(indexVars) > 1) {
    indexTable <- unique(data[, indexVars])
  } else {
    indexTable <- data.frame(unique(data[, indexVars]))
    # avoid data.frame's non-standard evaluation for names
    names(indexTable) <- indexVars
  }

  # use row number as ID
  indexTable[, "IDvar"] <- 1:nrow(indexTable)
  # pre-allocate column for index var. "0" means row wasn't found in indexTable
  data[, "IDvar"] <- 0
  # for i in indexTable, for j in data:
  # This seems ridiculously complicated. Why is there not a hash table object?
  invisible(
  sapply(1:nrow(indexTable), function(i) {
    sapply(1:nrow(data), function(j) {
      if(all(indexTable[i, indexVars] == data[j, indexVars])) {
        data[j, "IDvar"] <<- indexTable[i, "IDvar"]
      }
    })
  }))
  # check to make sure all rows have an ID
  if(0 %in% data["IDvar"]) {
    stop("Internal error: not all data has an ID from 'indexTable'.")
  }
  # split combined data into a list of data to iterate over
  data <- split(data, data["IDvar"])

  ## fit a quick test model to catch errors early and get column names
  test.model <- function() {
    # get data for all subjects, trials at (time, src) <- (1, 1)
    dat <- data[[1]]
    if(! is.null(addedData)) {
      dat <- merge(dat, addedData)
    }
    suppressWarnings(lmer2(formula, data=dat, onlyCoefs=TRUE,
                           lmerControls=list(lmerControl(
                           optCtrl=list(maxfun=1)))))
  }
  # names of the addedData, starting with "(Intercept)"
  varNames <- rownames(test.model())
  # add prefixes "coef", "se", "t"
  # results in (coef_1, ..., coef_n, se_1, ..., se_n, t_1, ..., t_n)
  varNames <- paste(rep(c("coef", "se", "t"), each=length(varNames)),
                   varNames, sep="_")

  ## Create a wrapper for the main model fitting function that will merge
  ##   'data' and 'addedData' each time it is called. The many small merges
  ##   are intended for when one doesn't want a huge data frame in memory.
  if(! is.null(addedData)) {
    # combine data and addedData before fitting
    fit.model <- function(formula, dat, onlyCoefs=TRUE, ...) {
      dat <- merge(dat, addedData)
      lmer2(formula, dat, onlyCoefs=onlyCoefs, ...)
    }
  } else {
    # necessary addedData are already in the data frame
    fit.model <- function(formula, dat, onlyCoefs=TRUE, ...) {
      lmer2(formula, dat, onlyCoefs=onlyCoefs, ...)
    }
  }

  ## fit some initial models that will be used to get starting values of theta
  # note: this will overwrite 'start' if both values are supplied
  # TODO: get these values from inside the foreach loop
  if(! is.null(nStarters) && nStarters != 0) {
    thetas <- sapply(sample(data, nStarters),
                    function(x) getME(fit.model(formula, x, start=start,
                                      onlyCoefs=FALSE)[["model"]], "theta"))
    # get mean of thetas
    if(is.null(dim(thetas))) {
      start <- mean(thetas)
    } else {
      start <- apply(thetas, 1, mean)
    }
  }

  ## setup the cluster
  if(multiprocessing) {
    # how many worker processes to spawn?
    nCores <- parallel::detectCores()
    if(is.null(maxProcesses)) {
      nWorkers <- nCores
    } else if(maxProcesses > nCores) {
      nWorkers <- nCores
    } else {
      nWorkers <- maxProcesses
    }
    # create socket cluster a la snow (cross-platform)
    # decide how logs will be delivered
    if(is.null(logging)) {
      # no messages
      cl <- parallel::makePSOCKcluster(nWorkers)
    } else if(logging == "console") {
      # messages will go to console of the master (use Rterm on Windows)
      cl <- parallel::makePSOCKcluster(nWorkers, outfile="")
    } else if(class(logging) == "character") {
      # messages will be written to a file
      # remove the log file if it already exists
      if(file.exists(logging)) {
        file.remove(logging)
        }
        cl <- parallel::makePSOCKcluster(nWorkers, outfile=logging)
    } else {
        stop("Unexpected argument for 'logging'.")
    }

    ## fit models in parallel
    # register the cluster with foreach
    doParallel::registerDoParallel(cl)
    # time of cluster initialization
    clusterStartTime <- Sys.time()
    # split the data and fit lmer models in parallel
    results <- foreach::foreach(dat=data, i=iterators::icount(),
                                .packages=c('lme4'), .export="lmer2") %dopar% {
      # TODO: set this initially instead of each iteration
      options(warn=1, warning.length=8170, max.print=1e9)
      motif <- "===================================="
      # print separator 77 characters long + 1 character for each digit of
      # i beyond 1 (can expand to 80 characters with i < 1e4)
      message(paste(motif, i, paste(motif, "\r", sep="")))
      # time since cluster initialization
      message(paste("Runtime:", (Sys.time() - clusterStartTime)[[1]]))
      # process ID
      message(paste("PID:", Sys.getpid()))
      # fit the model
      fit.model(formula, dat, start=start, ...)
      # foreach::.foreachGlobals may be a way to set global envir of workers
    }
    parallel::stopCluster(cl)
    # reshape from a list to a matrix
    results <- sapply(results, c)
  } else {
    # don't use multiprocessing
    results <- sapply(data, function(x) fit.model(formula, x, start=start, ...))
  }

  # transpose; columns after transpose:
  # (est_1, ..., est_n, se_1, ..., se_n, t_1, ..., t_n)
  results <- t(results)
  # add names
  colnames(results) <- varNames
  # combine results with their ID
  # currently this is coerced to a data frame before writing
  results <- cbind(indexTable, results)
  # cache results as R object
  if(! is.null(saveTo)) {
    saveRDS(results, paste(saveTo, "Rds", sep="."))
  }

  return(results)
}
