#---------------------------------------#
# RCNA_analysis
#---------------------------------------#

#' @title RCNA_analysis constructor
#' @description An S4 class used to track parameters from a specific RCNA function execution.
#'
#' @slot call A character vector detailing the function call that (`correct_gc_bias`, `estimate_nkr`, `estimate_feature_score`) was performed to produce this S4 object.
#' @slot params A list corresponding to the parameters that were submitted with the associated function call
#' @slot res.files A list containing the names of the flat files that were created by the documented function call.
#' @export RCNA_analysis

# S4 object doesn't need exported constructor since it will be created entirely by other functions in the package
RCNA_analysis = setClass("RCNA_analysis", representation(call = "character", params = "list", res.files = "list"),
		prototype(call = "", params = list(data.frame(matrix(nrow = 0, ncol = 0))), res.files = list()))

#---------------------------------------#
# RCNA.object
#---------------------------------------#

#' @title RCNA.object definition
#'
#' @description An S4 class used to specify parameters for an analysis run
#'
#' @slot sample.names Required. Character vector containing names of subjects
#' @slot ano.file Required. Character single file path detailing a feature-wise annotation file
#' @slot out.dir Required. Character vector containing the name of each subject?s output directory
#' @slot gcParams Data Frame storing all run parameters for the correct_gc_bias function
#' @slot nkrParams Data Frame storing all run parameters for the estimate_nkr function
#' @slot scoreParams Data Frame storing all run parameters for the estimate_feature_score function
#' @slot commands RCNA_analysis object storing commands and parameters from previous function runs on this object
#' @seealso \code{\link{run_RCNA}}
#' \code{\linkS4class{RCNA_analysis}},
#' @export RCNA_object
RCNA_object = setClass("RCNA_object",
		representation(sample.names = "character", ano.file = "character", out.dir = "character", ncpu = "numeric", # GLOBAL ARGUMENTS
				win.size = "numeric", "gc.step" = "numeric", estimate_gc = "logical", # GC correction non-vectorized parameters
				nkr = "numeric", norm.cov.matrix = "character", # NKR estimation non-vectorized parameters
				low.score.cutoff = "numeric", high.score.cutoff = "numeric", # Score estimation non-vectorized parameters
                        gcParams = "data.frame", nkrParams = "data.frame", scoreParams = "data.frame",
						commands = "list"),
                         prototype(out.dir = tempdir(), ncpu = 1, commands = NULL))

#' @title RCNA_object constructor
#'
#' @description An S4 class used to specify parameters for an analysis run
#'
#' @param sample.names Character vector containing names of subjects
#' @param ano.file Character single file path detailing a feature-wise annotation file
#' @param ncpu Numeric value specifying number of cores to use for analysis. Multiple cores will lead to parallel execution.
#' @param out.dir Character vector containing the name of each subject's output directory
#' @param file.coverage Character vector containing the path to the input coverage files for NKR and CNA score estimation.
#' @param win.size Numeric value detailing the size of the sliding window used to calculate and detect correct GC-content correction.
#' @param gc.step Numeric value detailing the size of each GC-content bin. If providing pre-calculated GC factor file this must match the bins in that file.
#' @param file.raw.coverage Character vector containing the filename of the raw coverage files for GC-content correction. Must be used in combination with `estimate_gc` set to TRUE.
#' @param file.corrected.coverage Character vector containing the filename of the corrected coverage files.
#' @param file.gc.factor Character vector containing the filename of GC factor files. Used if and only if `estimate_gc` is set to FALSE.
#' @param estimate_gc Logical that determines if GC bias factor is calculated. If set to TRUE then GC factor files will be generated for each sample. If set to FALSE then GC factor files must be supplied via `file.gc.factor`.
#' @param nkr Numeric between 0 and 1 which specifies the coverage quantile that should be considered a "normal" karyotype range for each position. Lowering this value may increase sensitivity but also Type I error.
#' @param file.nkr.coverage Character vector containing the filename of the input coverage file for NKR estimation. Defaults to `file coverage` if not specified.
#' @param x.norm Logical vector with length equal to the length of `sample.names`, denoting whether each subject has to be X-normalized. Subjects with an XX karyotype should be set to TRUE to avoid double-counting the coverage on the X chromosome. Set to FALSE if chrX coverage is already normalized.
#' @param norm.cov.matrix Character containing the directory or file name of the normalized coverage matrix. Generated by \code{\link{estimate_nkr}} if file doesn't exist.
#' @param file.score.coverage Character vector containing the input coverage file for the scoring function. Defaults to `file.coverage` if not specified.
#' @param score.cutoff Numeric between 0 and 1 which specifies the score filter on the results file. This parameter creates a symmetrical cutoff around 0, filtering all results whose absolute value is less than the specified value. Non-symmetrical cutoffs can be specified using `low.score.cutoff` and `high.score.cutoff`.
#' @param high.score.cutoff Numeric between 0 and 1 which specifies the upper score cutoff. Defaults to `score.cutoff` if not specified.
#' @param low.score.cutoff Numeric between 0 and 1 which specifies the lower score cutoff. Defaults to `score.cutoff` if not specified.
#' @param gcParams Data Frame storing all run parameters for the correct_gc_bias function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.
#' @param nkrParams Data Frame storing all run parameters for the estimate_nkr function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.
#' @param scoreParams Data Frame storing all run parameters for the estimate_feature_score function. Can be specified by a file path to a CSV file, `data.frame`, or (if not specified) will be generated by other arguments.
#' @param commands RCNA_analysis object storing commands and parameters from previous function runs on this object. For more information, see \code{\linkS4class{RCNA_analysis}}.
#' @param verbose Show more messages and warnings. Useful for debugging.
#' @examples
#' # Create an example object - see \link{example_obj} for more information.
#' samples <- c("ex-sample-1", "ex-sample-2", "ex-sample-3")
#' ex.obj <- create_RCNA_object(sample.names = samples,
#'                     ano.file = system.file("examples" ,"annotations-example.csv", package = "RCNA"),
#'                     out.dir = "output",
#'                     file.raw.coverage = system.file("examples", "coverage",
#'                        paste0(samples, ".txt.gz"), package = "RCNA"),
#'                     norm.cov.matrix = file.path("output", "norm-cov-matrix.csv.gz"),
#'                     nkr = 0.9,
#'                     x.norm = "FALSE",
#'                     low.score.cutoff = -0.35,
#'                     high.score.cutoff = 0.35,
#'                     ncpu = 1)
#' class(ex.obj)
#' \dontshow{system("rm -rf output")}
#' @return A \linkS4class{RCNA_object} class object with the specified parameters.
#' @importFrom methods new
#' @seealso \code{\linkS4class{RCNA_analysis}}, \code{\link{run_RCNA}}
#' @export create_RCNA_object

create_RCNA_object = function(sample.names, ano.file, ncpu = 1, out.dir = tempdir(), # GLOBAL PARAMETERS
		file.coverage = NULL, # Include an argument to specify coverage files for both NKR and score estimation
		gcParams = NULL, win.size = 75, gc.step = 0.01, file.raw.coverage = NULL, file.corrected.coverage = NULL, file.gc.factor = NULL, estimate_gc = TRUE,  # GC ESTIMATION PARAMETERS
		nkrParams = NULL, file.nkr.coverage = NULL, nkr = 0.9, x.norm = NULL, norm.cov.matrix = NULL, # NKR ESTIMATION PARAMETERS
		scoreParams = NULL, file.score.coverage = NULL, score.cutoff  = 0.5, low.score.cutoff = NULL, high.score.cutoff = NULL, # SCORE ESTIMATION PARAMETERS
		commands = list(), verbose = FALSE){ # Utility arguments

	# VALIDATE GLOBAL PARAMETERS

	# Validate sample names
	if(!is.character(sample.names)){
		tryCatch({ # Attempt to coerce sample.names to character
					sample.names = as.character(sample.names)
				},
				error = function(cond){
					message("sample.names could not be coerced to character")
				},
				warning = function(cond){
					message("Coercing sample.names caused a warning:")
					message(conditionMessage(cond))
				})
	}
	if(verbose) message("Creating RCNA object with ", length(sample.names), " samples")

	# Throw warning if no output directory specified
	if(missing(out.dir) & verbose){
		warning("No out.dir specified - defaulting to ", tempdir())
	}

  if(!dir.exists(out.dir)) {
    if(verbose) message("Specified output directory (", out.dir, ") not found - creating directory now.")
	  dir.create(out.dir)
	}

	# Check if annotation file is valid
	if(!check_anofile(ano.file)){
		stop("Please specify a valid annotation file")
	}

	# Detect # of cores - if less than specified then throw warning
	if(detectCores() < as.numeric(ncpu)) warning(ncpu, " cores specified but only ", detectCores(), " cores detected.")

	################################################
	# Construct GC estimation/correction parameters
	################################################

	# Validate estimate_gc
	if(!missing(estimate_gc)){
		if(!is.logical(estimate_gc)){
			tryCatch({ # Attempt to coerce estimate_gc to logical
					estimate_gc = as.logical(estimate_gc)
					},
					error = function(cond){
						message("estimate_gc could not be coerced to logical")
					},
					warning = function(cond){
						message("Coercing estimate_gc caused a warning:")
						message(conditionMessage(cond))
			})
		}

		if(length(estimate_gc > 1)){
			stop("Parameter `estimate_gc` has length > 1")
		}
	}

	# Validate win.size
	if(!missing(win.size)){
		if(!is.numeric(win.size)){
			tryCatch({ # Attempt to coerce win.size to numeric
						win.size = as.numeric(win.size)
					},
					error = function(cond){
						message("win.size could not be coerced to numeric")
					},
					warning = function(cond){
						message("Coercing win.size caused a warning:")
						message(conditionMessage(cond))
					})
		}

		if(length(win.size > 1)){
			stop("Parameter `win.size` has length > 1")
		}

		if(win.size < 1){
			stop("Invalid window size - please choose a window size >0. See `?correct_gc_bias()` for more information.")
		}
	}

	# Validate gc.step
	if(!missing(gc.step)){
		if(!is.numeric(gc.step)){
			tryCatch({ # Attempt to coerce gc.step to numeric
						gc.step = as.numeric(gc.step)
					},
					error = function(cond){
						message("gc.step could not be coerced to numeric")
					},
					warning = function(cond){
						message("Coercing gc.step caused a warning:")
						message(conditionMessage(cond))
					})
		}
		if(length(gc.step > 1)){
			stop("Parameter `gc.step` has length > 1")
		}
		if(gc.step <= 0 | gc.step > 1){
			stop("Invalid GC bin size - please choose a valid GC bin size. See `?correct_gc_bias()` for more information.")
		}
	}

	# If gcParams is provided, check for required columns
	if(!missing(gcParams)){

		# Read data
		if(typeof(gcParams)=="character"){ # if config file is a file path, check file path
			if(file.exists(gcParams)){
				gcParams.df = read.csv(gcParams, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
			} else {
				stop("Invalid `gcParams` file path.  See `?correct_gc_bias()` for more information.")
			}
		} else if(typeof(gcParams)=="data.frame"){
			gcParams.df = gcParams
		} else { # Try to coerce gcParams input to data.frame
			tryCatch({
			      if(verbose) ("gcParams parameter ambiguous - trying to coerce to data.frame")
						gcParams.df = as.data.frame(gcParams)
					},
					error = function(cond){
						message("gcParams could not be coerced to data.frame - please provide either a path to the parameter file or a data.frame object in the current workspace.")
					},
					warning = function(cond){
						message("Coercing gcParams caused a warning:")
						message(conditionMessage(cond))
					},
					finally = {
						if(verbose) message("gcParams parameter was coerced to data.frame")
					})
		}

		gcParams = gcParams.df # Assign results of above to gcParams

	} else if(!missing(file.corrected.coverage) | !missing(file.raw.coverage)){

		# Assign gcParams fields to list if specified via argument
		gcParams.list = list()
		if(!missing(file.corrected.coverage)) gcParams.list[["file.corrected.coverage"]] = file.corrected.coverage
		if(!missing(file.gc.factor)) gcParams.list[["file.gc.factor"]] = file.gc.factor
		if(!missing(file.raw.coverage)) gcParams.list[["file.raw.coverage"]] = file.raw.coverage
		if(!missing(sample.names)) gcParams.list[["sample.names"]] = sample.names

		# Attempt to construct a data.frame from list object, pass back error if it fails (e.g. if lengths differ)
		tryCatch({
					if(verbose) message("gcParams missing - trying to construct from specified arguments")
					gcParams = data.frame(gcParams.list, stringsAsFactors = FALSE)
				},
				error = function(cond){
					message("Tried to construct gcParams data.frame from specified conditions. See the following error:")
					message(conditionMessage(cond))
				},
				warning = function(cond){
					message("Constructing gcParams caused a warning:")
					message(conditionMessage(cond))
				},
				finally = {
					if(verbose) message("gcParams data.frame was constructed successfully")
		})
	}

	# VALIDATE GC PARAMS DATA FRAME
	if(!missing(gcParams) | !missing(file.corrected.coverage) | !missing(file.raw.coverage)){

		# Check if raw coverage file exists
		if(!("file.raw.coverage" %in% colnames(gcParams))){
			if(!missing(file.coverage)){ # If raw coverage file is not specified but file.coverage is then use file.coverage
				warning("No raw coverage file specified - using `file.coverage` (", file.coverage, ")")
				gcParams$file.raw.coverage = file.coverage
			} else {
				stop("No gcParams column labeled `file.raw.coverage`. Please provide the raw coverage file. See `?correct_gc_bias()` for more information.")
			}
		}

		# Check if corrected coverage file names exist
		if(!("file.corrected.coverage" %in% colnames(gcParams))){
			# If corrected coverage file isn't provided, construct corrected coverage file names from raw coverage files
			gcParams$file.corrected.coverage = sapply(gcParams$file.raw.coverage, function(x){
				fn = unlist(strsplit(basename(x), split = "\\."))
				fn[1] = paste0(fn[1], ".corrected")
				return(file.path(out.dir, "gc", paste0(fn, collapse = ".")))
			})
			# Re-assign for later use
			file.corrected.coverage = gcParams$file.corrected.coverage
		}

		# Check if file.gc.factor is missing and estimate_gc = FALSE
		if(!("file.gc.factor") %in% colnames(gcParams) & !estimate_gc){
			stop("`estimate_gc` was set to FALSE but no gcParams column labeled `file.gc.factor` used to specify the GC factor file. Please provide the GC factor file or set estimate_gc to TRUE. See `?correct_gc_bias()` for more information.")
		} else if(estimate_gc & "file.gc.factor" %in% colnames(gcParams)){ # Check if GC factor files exist? Do we want to include this?

		}

		# Check if sample.names exist in data.frame, if not assign from mandatory arg
		if(!("sample.names" %in% colnames(gcParams))){
			gcParams$sample.names = sample.names # Assign based on required argument
		}
	}

	#################################
	# Construct NKR estimation config
	#################################

	# If file.nkr.coverage or file.score.coverage are not specified, then try to find a default
	if(missing(file.nkr.coverage)){
		if(!missing(file.coverage)){
			file.nkr.coverage = file.coverage
		} else if(!missing(file.corrected.coverage)) {
			if(verbose) message("Setting normal karyotype range estimation input to corrected coverage file")
			file.nkr.coverage = file.corrected.coverage
		} else if(!missing(gcParams)){
			if(verbose) message("Setting normal karyotype range estimation input to corrected coverage file")
			gcParams$file.corrected.coverage
		}
	}

	# Validate x.norm
	if(!missing(x.norm)){
		if(!is.logical(x.norm)){
			tryCatch({ # Attempt to coerce x.norm to logical
						x.norm = as.logical(x.norm)
					},
					error = function(cond){
						message("x.norm could not be coerced to logical")
					},
					warning = function(cond){
						message("Coercing x.norm caused a warning:")
						message(conditionMessage(cond))
					})
		}
		if(length(x.norm) > 1 & length(x.norm)!=length(sample.names)){
			stop("Parameter `x.norm` has length > 1 and does not match the length of sample.names.")
		}
	}

	# Validate NKR parameter
	if(!missing(nkr)){
		if(!is.numeric(nkr)){
			tryCatch({ # Attempt to coerce nkr to numeric
						nkr = as.numeric(nkr)
					},
					error = function(cond){
						message("nkr could not be coerced to numeric")
					},
					warning = function(cond){
						message("Coercing nkr caused a warning:")
						message(conditionMessage(cond))
					})
		}
		if(length(nkr) > 1){
			stop("Parameter `nkr` has length > 1.")
		}
		if(nkr <= 0 | nkr >= 1){
			stop("Invalid value for `nkr` - must be strictly greater than 0 and less than 1. See `?estimate_nkr()` for more information.")
		}
	}

	# Validate norm.cov.matrix - file existence evaluated at function runtime
	if(!missing(norm.cov.matrix)){
		if(!is.character(norm.cov.matrix)){
			tryCatch({ # Attempt to coerce norm.cov.matrix to character
						norm.cov.matrix = as.character(norm.cov.matrix)
					},
					error = function(cond){
						message("norm.cov.matrix could not be coerced to character")
					},
					warning = function(cond){
						message("Coercing norm.cov.matrix caused a warning:")
						message(conditionMessage(cond))
					})
		}
		if(length(norm.cov.matrix) > 1){
			stop("Parameter `norm.cov.matrix` has length > 1.")
		}
	}

	# If nkrParams is provided, check for required columns
	if(!missing(nkrParams)){

		# Read data
		if(typeof(nkrParams)=="character"){ # if config file is a file path, check file path
			if(file.exists(nkrParams)){
				nkrParams.df = read.csv(nkrParams, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
			} else {
				stop("Invalid `nkrParams` file path.  See `?estimate_nkr()` for more information.")
			}
		} else if(typeof(nkrParams)=="data.frame"){
			nkrParams.df = nkrParams
		} else { # Try to coerce nkrParams input to data.frame
			tryCatch({
			      if(verbose) message("nkrParams parameter ambiguous - trying to coerce to data.frame")
						nkrParams.df = as.data.frame(nkrParams)
					},
					error = function(cond){
						message("nkrParams could not be coerced to data.frame - please provide either a path to the parameter file or a data.frame object in the current workspace.")
					},
					warning = function(cond){
						message("Coercing nkrParams caused a warning:")
						message(conditionMessage(cond))
					},
					finally = {
					  if(verbose) message("nkrParams parameter was coerced to data.frame")
					})
		}

		nkrParams = nkrParams.df # Assign results of above to nkrParams

	} else if(!missing(file.nkr.coverage) | !missing(file.coverage) | !missing(file.corrected.coverage) | !missing(gcParams)){ # as long as there's something to use for file.nkr.coverage

		# Assign nkrParams fields to list if specified via argument
		nkrParams.list = list()
		nkrParams.list[["file.nkr.coverage"]] = file.nkr.coverage
		if(!missing(x.norm)) nkrParams.list[["x.norm"]] = x.norm
		nkrParams.list[["sample.names"]] = sample.names

		# Attempt to construct a data.frame from list object, pass back error if it fails (e.g. if lengths differ)
		tryCatch({
					if(verbose) message("nkrParams missing - trying to construct from specified arguments")
					nkrParams = data.frame(nkrParams.list, stringsAsFactors = FALSE)
				},
				error = function(cond){
					message("Tried to construct nkrParams data.frame from specified conditions. See the following error:")
					message(conditionMessage(cond))
				},
				warning = function(cond){
					message("Constructing nkrParams caused a warning:")
					message(conditionMessage(cond))
				},
				finally = {
					if(verbose) message("nkrParams data.frame was coerced successfully")
				})
	}

	# VALIDATE NKR PARAMS DATA FRAME
	if(!missing(nkrParams) | !missing(file.nkr.coverage) | !missing(file.coverage) | !missing(gcParams)){

		# Check if raw coverage file exists and estimate_gc = TRUE
		if(!("file.nkr.coverage" %in% colnames(nkrParams))){
			stop("No coverage file specified - please use one of the following arguments: file.coverage, file.raw.coverage, file.corected.coverage, file.nkr.coverage")
		}

		if(!"x.norm" %in% colnames(nkrParams)){
			if(missing(x.norm)){
				stop("No chrX normalization specified or found in nkrParams. Please see `?estimate_nkr()` for more information.")
			} else {
				nkrParams$x.norm = x.norm
			}
		}

		# Check if sample.names exist in data.frame, if not assign from mandatory arg
		if(!("sample.names" %in% colnames(nkrParams))){
			nkrParams$sample.names = sample.names # Assign based on required argument
		}
	}

	####################################
	# Construct score estimation config
	####################################

	# Validates score parameter
	if(!missing(score.cutoff)){
		if(!is.numeric(score.cutoff)){
			tryCatch({ # Attempt to coerce score.cutoff to numeric
						score.cutoff = as.numeric(score.cutoff)
					},
					error = function(cond){
						message("score.cutoff could not be coerced to numeric")
					},
					warning = function(cond){
						message("Coercing score.cutoff caused a warning:")
						message(conditionMessage(cond))
					})
		}
		if(length(score.cutoff) > 1){
			stop("Parameter `score.cutoff` has length > 1.")
		}

		if(score.cutoff < -1 | score.cutoff > 1){
			stop("Invalid value for `score.cutoff` - score must be between -1 and 1. See `?estimate_score()` for more information.")
		}
	}

	# Validates high.score.cutoff parameter
	if(!missing(high.score.cutoff)){
		if(!is.numeric(high.score.cutoff)){
			tryCatch({ # Attempt to coerce high.score.cutoff to numeric
						high.score.cutoff = as.numeric(high.score.cutoff)
					},
					error = function(cond){
						message("high.score.cutoff could not be coerced to numeric")
					},
					warning = function(cond){
						message("Coerhigh.score.cutoffng high.score.cutoff caused a warning:")
						message(conditionMessage(cond))
					})
		}
		if(length(high.score.cutoff) > 1){
			stop("Parameter `high.score.cutoff` has length > 1.")
		}
		if(high.score.cutoff > 1){
			stop("Invalid value for `high.score.cutoff` - score must be strictly between -1 and 1. See `?estimate_score()` for more information.")
		}
	} else {
		high.score.cutoff = abs(score.cutoff)
	}

	# Validates low.score.cutoff parameter
	if(!missing(low.score.cutoff)){
		if(!is.numeric(low.score.cutoff)){
			tryCatch({ # Attempt to coerce low.score.cutoff to numeric
						low.score.cutoff = as.numeric(low.score.cutoff)
					},
					error = function(cond){
						message("low.score.cutoff could not be coerced to numeric")
					},
					warning = function(cond){
						message("Coerlow.score.cutoffng low.score.cutoff caused a warning:")
						message(conditionMessage(cond))
					})
		}
		if(length(low.score.cutoff) > 1){
			stop("Parameter `low.score.cutoff` has length > 1.")
		}
		if(low.score.cutoff < -1){
			stop("Invalid value for `low.score.cutoff` - score must be strictly between -1 and 1. See `?estimate_score()` for more information.")
		}
	} else {
		low.score.cutoff = -1*abs(score.cutoff)
	}

	# If file.score.coverage is not specified, then try to find a default
	if(missing(file.score.coverage)){
		if(!missing(file.coverage)){
			file.score.coverage = file.coverage
		} else if(!missing(file.corrected.coverage)) {
			if(verbose) message("Setting score estimation input to corrected coverage file")
			file.score.coverage = file.corrected.coverage
		} else if(!missing(gcParams)){
			if(verbose) message("Setting score estimation input to corrected coverage file")
			gcParams$file.corrected.coverage
		}
	}

	# If scoreParams is provided, check for required columns
	if(!missing(scoreParams)){

		# Read data
		if(typeof(scoreParams)=="character"){ # if config file is a file path, check file path
			if(file.exists(scoreParams)){
				scoreParams.df = read.csv(scoreParams, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
			} else {
				stop("Invalid `scoreParams` file path.  See `?estimate_score()` for more information.")
			}
		} else if(typeof(scoreParams)=="data.frame"){
			scoreParams.df = scoreParams
		} else { # Try to coerce scoreParams input to data.frame
			tryCatch({
			      if(verbose) message("scoreParams parameter ambiguous - trying to coerce to data.frame")
						scoreParams.df = as.data.frame(scoreParams)
					},
					error = function(cond){
						message("scoreParams could not be coerced to data.frame - please provide either a path to the parameter file or a data.frame object in the current workspace.")
					},
					warning = function(cond){
						message("Coercing scoreParams caused a warning:")
						message(conditionMessage(cond))
					},
					finally = {
					  if(verbose) message("scoreParams parameter was coerced to data.frame")
					})
		}

		scoreParams = scoreParams.df # Assign results of above to scoreParams

	} else if(!missing(file.score.coverage) | !missing(file.coverage) | !missing(file.corrected.coverage) | !missing(gcParams)){ # as long as there's something to use for file.nkr.coverage

		# Assign scoreParams fields to list if specified via argument
		scoreParams.list = list()
		scoreParams.list[["file.score.coverage"]] = file.score.coverage
		scoreParams.list[["sample.names"]] = sample.names

		# Attempt to construct a data.frame from list object, pass back error if it fails (e.g. if lengths differ)
		tryCatch({
					if(verbose) message("scoreParams missing - trying to construct from specified arguments")
					scoreParams = data.frame(scoreParams.list, stringsAsFactors = FALSE)
				},
				error = function(cond){
					message("Tried to construct scoreParams data.frame from specified conditions. See the following error:")
					message(conditionMessage(cond))
				},
				warning = function(cond){
					message("Constructing scoreParams caused a warning:")
					message(conditionMessage(cond))
				},
				finally = {
					if(verbose) message("scoreParams data.frame was coerced successfully")
				})
	}

	# VALIDATE NKR PARAMS DATA FRAME
	if(!missing(scoreParams) | !missing(file.nkr.coverage) | !missing(file.coverage) | !missing(gcParams)){

		# Check if raw coverage file exists and estimate_gc = TRUE
		if(!("file.score.coverage" %in% colnames(scoreParams))){
			stop("No coverage file specified - please use one of the following arguments: file.coverage, file.raw.coverage, file.corected.coverage, file.score.coverage")
		}

		# Check if sample.names exist in data.frame, if not assign from mandatory arg
		if(!("sample.names" %in% colnames(scoreParams))){
			scoreParams$sample.names = sample.names # Assign based on required argument
		}
	}

	new("RCNA_object",
			sample.names = sample.names, ano.file = ano.file, out.dir = out.dir, ncpu = ncpu,
			win.size = win.size, gc.step = gc.step, estimate_gc = estimate_gc,
			nkr = nkr, norm.cov.matrix = norm.cov.matrix,
			low.score.cutoff = low.score.cutoff, high.score.cutoff = high.score.cutoff,
			gcParams = gcParams, nkrParams = nkrParams, scoreParams = scoreParams,
			commands = commands)
}

# Helper function to check annotation file
check_anofile = function(ano.file){

	if(is.character(ano.file)){
		ano = read.csv(ano.file, stringsAsFactors = FALSE)
	} else if(!is.data.frame(ano.file)){
		warning("Please provide `ano.file` as either a file path or a `data.frame`.")
	  return(FALSE)
	}

	if(!all(c("feature", "chromosome", "start", "end") %in% colnames(ano))){
		warning("Provided `ano.file` is invalid as an annotation file - should include the following columns:\n feature   chromosome   start   end")
		return(FALSE)
	} else {
		return(TRUE)
	}
}
