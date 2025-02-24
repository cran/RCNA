% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RCNA.r
\name{run_RCNA}
\alias{run_RCNA}
\alias{run_RCNA.default}
\alias{run_RCNA.RCNA_object}
\title{run_RCNA: Perform RCNA copy number detection workflow}
\usage{
run_RCNA(obj, ...)

\method{run_RCNA}{default}(
  obj = NULL,
  sample.names,
  ano.file,
  out.dir = tempdir(),
  gcParams = NULL,
  win.size = 75,
  gc.step = 0.01,
  file.raw.coverage = NULL,
  file.corrected.coverage = NULL,
  file.gc.factor = NULL,
  estimate_gc = TRUE,
  nkrParams,
  file.nkr.coverage = NULL,
  ncpu = 1,
  nkr = 0.9,
  x.norm = NULL,
  scoreParams,
  score.cutoff = 0.5,
  low.score.cutoff = NULL,
  high.score.cutoff = NULL,
  commands = c(),
  verbose = FALSE,
  ...
)

\method{run_RCNA}{RCNA_object}(obj, estimate_gc = TRUE, verbose = FALSE, ...)
}
\arguments{
\item{obj}{An `RCNA_object` type created by \link{create_RCNA_object}.}

\item{...}{Additional arguments (unused).}

\item{sample.names}{Character vector containing names of subjects}

\item{ano.file}{Character single file path detailing a feature-wise annotation file}

\item{out.dir}{Character vector containing the name of each subject's output directory}

\item{gcParams}{Data Frame storing all run parameters for the correct_gc_bias function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.}

\item{win.size}{Numeric value detailing the size of the sliding window used to calculate and detect correct GC-content correction.}

\item{gc.step}{Numeric value detailing the size of each GC-content bin. If providing pre-calculated GC factor file this must match the bins in that file.}

\item{file.raw.coverage}{Character vector containing the filename of the raw coverage files for GC-content correction. Must be used in combination with `estimate_gc` set to TRUE.}

\item{file.corrected.coverage}{Character vector containing the filename of the corrected coverage files.}

\item{file.gc.factor}{Character vector containing the filename of GC factor files. Used if and only if `estimate_gc` is set to FALSE.}

\item{estimate_gc}{A logical which determines if GC estimation should be performed. For more information, see \link{correct_gc_bias}.}

\item{nkrParams}{Data Frame storing all run parameters for the estimate_nkr function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.}

\item{file.nkr.coverage}{Character vector containing the filename of the input coverage file for NKR estimation. Defaults to `file coverage` if not specified.}

\item{ncpu}{Numeric value specifying number of cores to use for analysis. Multiple cores will lead to parallel execution.}

\item{nkr}{Numeric between 0 and 1 which specifies the coverage quantile that should be considered a "normal" karyotype range for each position. Lowering this value may increase sensitivity but also Type I error.}

\item{x.norm}{Logical vector with length equal to the length of `sample.names`, denoting whether each subject has to be X-normalized. Subjects with an XX karyotype should be set to TRUE to avoid double-counting the coverage on the X chromosome. Set to FALSE if chrX coverage is already normalized.}

\item{scoreParams}{Data Frame storing all run parameters for the estimate_feature_score function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.}

\item{score.cutoff}{Numeric between 0 and 1 which specifies the score filter on the results file. This parameter creates a symmetrical cutoff around 0, filtering all results whose absolute value is less than the specified value. Non-symmetrical cutoffs can be specified using `low.score.cutoff` and `high.score.cutoff`.}

\item{low.score.cutoff}{Numeric between 0 and 1 which specifies the lower score cutoff. Defaults to `score.cutoff` if not specified.}

\item{high.score.cutoff}{Numeric between 0 and 1 which specifies the upper score cutoff. Defaults to `score.cutoff` if not specified.}

\item{commands}{RCNA_analysis object storing commands and parameters from previous function runs on this object. For more information, see \code{\linkS4class{RCNA_analysis}}.}

\item{verbose}{If set to TRUE will display more detailed error messages.}
}
\value{
A \linkS4class{RCNA_object} class object that was used during the workflow, with \link{RCNA_analysis} objects in the `commands` slot that describes the run parameters and results of each step in the workflow.

A \linkS4class{RCNA_object} class object that was used during the workflow, with \link{RCNA_analysis} objects in the `commands` slot that describes the run parameters and results of each step in the workflow. For more details on outputs, see \link{estimate_nkr}, \link{correct_gc_bias}, and \link{estimate_feature_score}.

A \linkS4class{RCNA_object} class object that was used during the workflow, with \link{RCNA_analysis} objects in the `commands` slot that describes the run parameters and results of each step in the workflow. For more details on outputs, see \link{estimate_nkr}, \link{correct_gc_bias}, and \link{estimate_feature_score}.
}
\description{
`run_RCNA` will execute \link{correct_gc_bias}, \link{estimate_nkr}, and \link{estimate_feature_score} in that specific order. For more information, see each of those functions' individual documentation.

`run_RCNA` will execute \link{correct_gc_bias}, \link{estimate_nkr}, and \link{estimate_feature_score} in that specific order. For more information, see each of those functions' individual documentation, or \link{create_RCNA_object}.

`run_RCNA` will execute \link{correct_gc_bias}, \link{estimate_nkr}, and \link{estimate_feature_score} in that specific order. For more information, see each of those functions' individual documentation.
}
\examples{
## Run RCNA workflow on example object
# See ?example_obj for more information on example
example_obj@ano.file <- system.file("examples" ,"annotations-example.csv", package = "RCNA")
raw.cov <- system.file("examples", "coverage",
                       paste0(example_obj@sample.names, ".txt.gz"), package = "RCNA")
example_obj@gcParams$file.raw.coverage <- raw.cov
\donttest{example_obj}
# Run RCNA workflow
result_obj <- run_RCNA(example_obj)
\dontshow{system("rm -rf output")}
}
\seealso{
\linkS4class{RCNA_object}, \linkS4class{RCNA_analysis}, \link{correct_gc_bias}, \link{run_RCNA}, \link{estimate_feature_score}
}
