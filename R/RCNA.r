#' @include correct_gc_bias.r
#' @include estimate_nkr.r
#' @include RCNA.r
#' @include estimate_feature_score.r
#' @include classes.r
#'
#' Run every step of the RCNA workflow
#'
#' @description
#' `run_RCNA` will execute \link{correct_gc_bias}, \link{estimate_nkr}, and \link{estimate_feature_score} in that specific order. For more information, see each of those functions' individual documentation.
#'
#' @title run_RCNA: Perform RCNA copy number detection workflow
#' @param obj A \linkS4class{RCNA_object} class object, or specify each argument (see \link{run_RCNA.default}).
#' @param ... Other arguments
#' @return A \linkS4class{RCNA_object} class object that was used during the workflow, with \link{RCNA_analysis} objects in the `commands` slot that describes the run parameters and results of each step in the workflow.
#' @seealso \linkS4class{RCNA_object}, \linkS4class{RCNA_analysis}, \link{correct_gc_bias}, \link{run_RCNA}, \link{estimate_feature_score}
#' @rdname run_RCNA.default
#' @export run_RCNA

run_RCNA <- function(obj, ...) UseMethod("run_RCNA")

#'
#' Run every step of the RCNA workflow
#'
#' `run_RCNA` will execute \link{correct_gc_bias}, \link{estimate_nkr}, and \link{estimate_feature_score} in that specific order. For more information, see each of those functions' individual documentation, or \link{create_RCNA_object}.
#' @param sample.names Character vector containing names of subjects
#' @param ano.file Character single file path detailing a feature-wise annotation file
#' @param ncpu Numeric value specifying number of cores to use for analysis. Multiple cores will lead to parallel execution.
#' @param out.dir Character vector containing the name of each subject's output directory
#' @param win.size Numeric value detailing the size of the sliding window used to calculate and detect correct GC-content correction.
#' @param gc.step Numeric value detailing the size of each GC-content bin. If providing pre-calculated GC factor file this must match the bins in that file.
#' @param file.raw.coverage Character vector containing the filename of the raw coverage files for GC-content correction. Must be used in combination with `estimate_gc` set to TRUE.
#' @param file.corrected.coverage Character vector containing the filename of the corrected coverage files.
#' @param file.gc.factor Character vector containing the filename of GC factor files. Used if and only if `estimate_gc` is set to FALSE.
#' @param estimate_gc Logical that determines if GC bias factor is calculated. If set to TRUE then GC factor files will be generated for each sample. If set to FALSE then GC factor files must be supplied via `file.gc.factor`.
#' @param nkr Numeric between 0 and 1 which specifies the coverage quantile that should be considered a "normal" karyotype range for each position. Lowering this value may increase sensitivity but also Type I error.
#' @param file.nkr.coverage Character vector containing the filename of the input coverage file for NKR estimation. Defaults to `file coverage` if not specified.
#' @param x.norm Logical vector with length equal to the length of `sample.names`, denoting whether each subject has to be X-normalized. Subjects with an XX karyotype should be set to TRUE to avoid double-counting the coverage on the X chromosome. Set to FALSE if chrX coverage is already normalized.
#' @param score.cutoff Numeric between 0 and 1 which specifies the score filter on the results file. This parameter creates a symmetrical cutoff around 0, filtering all results whose absolute value is less than the specified value. Non-symmetrical cutoffs can be specified using `low.score.cutoff` and `high.score.cutoff`.
#' @param high.score.cutoff Numeric between 0 and 1 which specifies the upper score cutoff. Defaults to `score.cutoff` if not specified.
#' @param low.score.cutoff Numeric between 0 and 1 which specifies the lower score cutoff. Defaults to `score.cutoff` if not specified.
#' @param gcParams Data Frame storing all run parameters for the correct_gc_bias function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.
#' @param nkrParams Data Frame storing all run parameters for the estimate_nkr function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.
#' @param scoreParams Data Frame storing all run parameters for the estimate_feature_score function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.
#' @param commands RCNA_analysis object storing commands and parameters from previous function runs on this object. For more information, see \code{\linkS4class{RCNA_analysis}}.
#' @param verbose Show more messages and warnings. Useful for debugging.
#' @param obj An `RCNA_object` type created by \link{create_RCNA_object}. If specified this takes precedent over all over args.
#' @param ... Additional arguments (unused)
#' @examples
#' ## Run RCNA workflow on example object
#' # See ?example_obj for more information on example
#' example_obj@ano.file <- system.file("examples" ,"annotations-example.csv", package = "RCNA")
#' raw.cov <- system.file("examples", "coverage",
#'                        paste0(example_obj@sample.names, ".txt.gz"), package = "RCNA")
#' example_obj@gcParams$file.raw.coverage <- raw.cov
#' \donttest{example_obj}
#' # Run RCNA workflow
#' result_obj <- run_RCNA(example_obj)
#' \dontshow{system("rm -rf output")}
#' @return  A \linkS4class{RCNA_object} class object that was used during the workflow, with \link{RCNA_analysis} objects in the `commands` slot that describes the run parameters and results of each step in the workflow. For more details on outputs, see \link{estimate_nkr}, \link{correct_gc_bias}, and \link{estimate_feature_score}.
#' @rdname run_RCNA.default
#' @import utils
#' @importFrom methods isClass
#' @export

run_RCNA.default <- function(obj = NULL, sample.names, ano.file, # REQUIRED PARAMETERS
                               out.dir = tempdir(), # Optional arguments
                               gcParams = NULL, win.size = 75, gc.step = 0.01, file.raw.coverage = NULL, file.corrected.coverage = NULL, file.gc.factor = NULL, estimate_gc = TRUE,  # GC ESTIMATION PARAMETERS
                               nkrParams, file.nkr.coverage = NULL, ncpu = 1, nkr = 0.9, x.norm = NULL, # NKR ESTIMATION PARAMETERS
                               scoreParams,  score.cutoff  = 0.5, low.score.cutoff = NULL, high.score.cutoff = NULL, # SCORE ESTIMATION PARAMETERS
                               commands = c(), verbose = FALSE, ...){

  if(missing(obj) | !isClass("RCNA_object", obj)){
    new_RCNA_obj = create_RCNA_object(sample.names = sample.names, ano.file = ano.file, # REQUIRED PARAMETERS
                                      out.dir = out.dir, # Optional arguments
                                      gcParams = gcParams, win.size = win.size, gc.step = gc.step, file.raw.coverage = file.raw.coverage, file.corrected.coverage = file.corrected.coverage, file.gc.factor = file.gc.factor, # GC ESTIMATION PARAMETERS
                                      nkrParams = nkrParams, file.nkr.coverage = file.nkr.coverage, ncpu = ncpu, nkr = nkr, x.norm = x.norm, # NKR ESTIMATION PARAMETERS
                                      scoreParams = scoreParams,  score.cutoff  = score.cutoff, low.score.cutoff = low.score.cutoff, high.score.cutoff = high.score.cutoff, # SCORE ESTIMATION PARAMETERS
                                      commands = commands, verbose = FALSE)
  } else {
    new_RCNA_obj = obj
  }

  run_RCNA.RCNA_object(new_RCNA_obj, estimate_gc = estimate_gc, verbose = verbose)
}

#'
#' Run every step of the RCNA workflow
#'
#' `run_RCNA` will execute \link{correct_gc_bias}, \link{estimate_nkr}, and \link{estimate_feature_score} in that specific order. For more information, see each of those functions' individual documentation.
#' @param obj An `RCNA_object` type created by \link{create_RCNA_object}.
#' @param estimate_gc A logical which determines if GC estimation should be performed. For more information, see \link{correct_gc_bias}.
#' @param verbose If set to TRUE will display more detailed error messages.
#' @param ... Additional arguments (unused).
#' @return  A \linkS4class{RCNA_object} class object that was used during the workflow, with \link{RCNA_analysis} objects in the `commands` slot that describes the run parameters and results of each step in the workflow. For more details on outputs, see \link{estimate_nkr}, \link{correct_gc_bias}, and \link{estimate_feature_score}.
#' @rdname run_RCNA.default
#' @import utils
#' @export

run_RCNA.RCNA_object <- function(obj, estimate_gc = TRUE, verbose = FALSE, ...){

  if(verbose) message("Beginning RCNA analysis...")
  start.time = Sys.time()
  if(verbose) message("[1/3] Beginning GC bias correction")
  obj@commands = append(obj@commands, correct_gc_bias(obj, estimate_gc = estimate_gc, verbose = verbose))
  if(verbose) message("[1/3] GC bias correction completed!")
  correct.gc.time = Sys.time() - start.time
  if(verbose) message("[1/3] Total GC step execution time: ", round(correct.gc.time, 3), " ",units(correct.gc.time))

  start.time.nkr = Sys.time()
  if(verbose) message("[2/3] Beginning NKR Estimation")
  obj@commands = append(obj@commands, estimate_nkr(obj, verbose = verbose))
  if(verbose) message("[2/3] NKR Estimation completed!")
  estimate.nkr.time = Sys.time() - start.time.nkr
  if(verbose) message("[2/3] Total NKR step execution time: ", round(estimate.nkr.time, 3), " ",units(estimate.nkr.time))

  start.time.score = Sys.time()
  if(verbose) message("[3/3] Beginning Feature Score Estimation")
  obj@commands = append(obj@commands, estimate_feature_score(obj, verbose = verbose))
  if(verbose) message("[3/3] Feature Score Estimation completed!")
  estimate.score.time = Sys.time() - start.time.score
  if(verbose) message("[3/3] Total SCORE step execution time: ", round(estimate.score.time,3), " ",units(estimate.score.time))

  if(verbose) message("RCNA analysis completed!")
  final.time = Sys.time() - start.time
  if(verbose) message("Total runtime: ", round(final.time, 3), " ",units(final.time))
  # Return object with commands slot
  return(obj)
}
