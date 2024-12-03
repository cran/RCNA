#' @include classes.r
#' @include RCNA.r

#' Estimate the CNA score for each feature
#'
#' This function estimates the the CNA score for each feature in the annotation file. It creates two flat file text tables with a row for each feature, which is placed in the output directory under `/score` - one with the score filter applied and one with all score results reported.
#'
#' @title estimate_feature_score: Estimate CNV score for each gene in the annotation file
#' @param obj A \linkS4class{RCNA_object} class object, or a data.frame object that contains all the necessary parameters.
#' @param ... Other arguments
#' @details This function can be run as a stand-alone or as part of \link{run_RCNA}.
#' @return A \linkS4class{RCNA_analysis} class object that describes the input parameters and output files generated by this step of the workflow.
#' @seealso \linkS4class{RCNA_object}, \linkS4class{RCNA_analysis}, \link{run_RCNA}
#' @rdname estimate_feature_score.default
#' @export estimate_feature_score

estimate_feature_score <- function(obj, ...) UseMethod("estimate_feature_score")

#' Estimate the CNA score for each feature
#'
#' This function estimates the the CNA score for each feature in the annotation file. It creates two flat file text tables with a row for each feature, which is placed in the output directory under `/score` - one with the score filter applied and one with all score results reported.
#'
#' @param df Path to the config file, or a `data.frame` object containing the valid parameters. Valid column names are `file.score.coverage` and `sample.names`. Additional columns will be ignored.
#' @param sample.names Character vector of sample names. Alternatively can be specified in `df`.
#' @param ano.file Location of the annotation file. This file must be in CSV format and contain the following information (with column headers as specified): "feature,chromosome,start,end".
#' @param out.dir Output directory for results. A subdirectory for results will be created under this + `/nkr/`.
#' @param ncpus Integer number of CPUs to use. Specifying more than one allows this function to be parallelized by feature.
#' @param file.score.coverage Character vector listing the input coverage files. Must be the same length as `sample.names`. Alternatively can be specified in `df`.
#' @param score.cutoff Numeric between 0 and 1 which specifies the score filter on the results file. This parameter creates a symmetrical cutoff around 0, filtering all results whose absolute value is less than the specified value. Non-symmetrical cutoffs can be specified using `low.score.cutoff` and `high.score.cutoff`.
#' @param high.score.cutoff Numeric between 0 and 1 which specifies the upper score cutoff. Defaults to `score.cutoff` if not specified.
#' @param low.score.cutoff Numeric between 0 and 1 which specifies the lower score cutoff. Defaults to `score.cutoff` if not specified.
#' @param verbose When TRUE increases level of detail for command output
#' @param obj An `RCNA_object` type created by \link{create_RCNA_object}. If specified this takes precedent over all over args.
#' @param ... Additional arguments
#' @examples
#' ## Estimate feature scores on example object
#' # See \link{example_obj} for more information on example
#' example_obj@ano.file <- system.file("examples" ,"annotations-example.csv", package = "RCNA")
#' \donttest{example_obj}
#' # Create output directories
#' dir.create(file.path("output", "score"), recursive = TRUE)
#' # Copy example GC-corrected coverage files
#' cov.corrected <- system.file("examples", "gc", package = "RCNA")
#' file.copy(from = cov.corrected, to = "output", recursive = TRUE)
#' # Copy example NKR results for "feature_a"
#' nkr.res <- system.file("examples", "nkr", package = "RCNA")
#' file.copy(from = nkr.res, to = "output", recursive = TRUE)
#' # Run score estimation for "feature_a" and append results
#' estimate_feature_score_analysisObj <- estimate_feature_score(example_obj)
#' example_obj@commands <- c(example_obj@commands, estimate_feature_score_analysisObj)
#' \dontshow{system("rm -rf output")}
#' @return A \linkS4class{RCNA_analysis} class object that describes the input parameters and output files generated by this step of the workflow.
#' @rdname estimate_feature_score.default
#' @import utils parallel doParallel
#' @details
#'
#' The `df` argument corresponds to the `scoreParams` matrix on \linkS4class{RCNA_object}. Valid column names are `sample.names` and `file.score.coverage`. Additional columns will be ignored.
#' @export

estimate_feature_score.default = function(obj = NULL, df = NULL, sample.names = NULL, ano.file, out.dir = NULL, ncpus = 1,
		file.score.coverage = NULL, score.cutoff  = 0.5, low.score.cutoff = NULL, high.score.cutoff = NULL,
		verbose = FALSE, ...){

  if(is.data.frame(obj) | is.character(obj)){
    df = obj
  } else if(!missing(obj)){
    return(estimate_feature_score.RCNA_object(obj, verbose))
  } else if(missing(df)){
		if(verbose) message("No data.frame detected for estimate_feature_score. Constructing run from parameters instead")
		if(!missing(sample.names) & !missing(ano.file)){
			tryCatch({ # Attempt to create an RCNA object
						new_RCNA_obj = create_RCNA_object(sample.names = sample.names, ano.file = ano.file, ncpu = ncpus, ...)
					},
					error = function(cond){
						message("RCNA Object could not be created due to error:")
					},
					warning = function(cond){
						message("Creating an RCNA object caused a warning:")
						message(conditionMessage(cond))
					})
		} else if(missing(sample.names)){
			stop("Missing both `df` argument and `sample.names` argument - please specify at least one.")
		}
	} else {
		if(missing(sample.names)){
			tryCatch({
						if(verbose) message("Pulling sample.names from data.frame")
						sample.names = df[["sample.names"]]
					},
					error = function(cond){
						message("Error while populating sample names from data.frame:")
					},
					warning = function(cond){
						message("Populating sample names produced a warning:")
						message(conditionMessage(cond))
					})
		}

		tryCatch({ # Attempt to create an RCNA object
					new_RCNA_obj = create_RCNA_object(sample.names = sample.names, ano.file = ano.file, ncpu = ncpus, scoreParams = df, ...)
				},
				error = function(cond){
					message("RCNA Object could not be created due to error:")
				},
				warning = function(cond){
					message("Creating an RCNA object caused a warning:")
					message(conditionMessage(cond))
				})
	}

	analysis.obj = estimate_nkr.RCNA_object(obj = new_RCNA_obj, verbose = verbose)

	new_RCNA_obj@commands = analysis.obj
	return(new_RCNA_obj)
}

#' Estimate the CNA score for each feature
#'
#' This function estimates the the CNA score for each feature in the annotation file. It creates two flat file text tables with a row for each feature, which is placed in the output directory under `/score` - one with the score filter applied and one with all score results reported.
#'
#' @rdname estimate_feature_score.default
#' @param obj A RCNA_object type object - parameters will be pulled from the object instead, specifically from the `scoreParams` slot.
#' @param verbose If set to TRUE will display more detail
#' @param ... Additional arguments (unused)
#' @return A \linkS4class{RCNA_analysis} class object that describes the input parameters and output files generated by this step of the workflow.
#' @importFrom R.utils gzip
#' @method estimate_feature_score RCNA_object
#' @details For more parameter information, see \link{estimate_feature_score.default}.
#' @export

estimate_feature_score.RCNA_object = function(obj, verbose = FALSE, ...){

  # Read data
  config.df = obj@scoreParams
  ano.file = obj@ano.file
  lscore = obj@low.score.cutoff
  uscore = obj@high.score.cutoff
  ncpus = obj@ncpu

  nkr.dir = file.path(obj@out.dir, "nkr")

  feature.score.outdir = paste0(obj@out.dir, "/score")
  if(!dir.exists(feature.score.outdir)){
	  if(verbose) message("Creating feature score output directory at ", feature.score.outdir)
	  dir.create(feature.score.outdir, recursive = TRUE)
  }

  if(!dir.exists(nkr.dir)){
	  stop("CI directory not found - expected at ", nkr.dir)
  }

  config.df = config.df[,c("file.score.coverage", "sample.names")]
  do_estimate_feature_score(in.files = config.df$file.score.coverage, sample.names = config.df$sample.names, ano.file = ano.file, nkr.dir = nkr.dir, outfile = feature.score.outdir, lscore = lscore, uscore = uscore, ncpus = ncpus)
  if(verbose) message("Feature score estimation succeeded!")
  return(new("RCNA_analysis", call = "estimate_feature_score", params = list("scoreParams" = config.df, "low.score.cutoff" = lscore, "high.score.cutoff" = uscore), res.files = list("Output" = list.files(feature.score.outdir))))
}

#' @importFrom R.utils gzip
#' @import foreach parallel

do_estimate_feature_score = function(in.files, nkr.dir, ano.file, outfile, lscore = 0.1, uscore = 0.1, sample.names, ncpus = 1) {

	# Coerce numerical args
	lscore = as.numeric(lscore)
	uscore = as.numeric(uscore)
	sample.names = sample.names
	outdir = outfile
	outfile = file.path(outdir, "RCNA-output.csv")
	ncpus = ncpus
	outfile.pass = file.path(outdir, "RCNA-output-scorepass.csv")

	## import feature annotations
	ano = read.csv(unlist(ano.file),sep=',',header=TRUE,stringsAsFactors=FALSE)

	## generate result table
	tbl = as.data.frame(matrix(ncol=10,nrow=0),stringsAsFactors=FALSE)
	tbl.pass = as.data.frame(matrix(ncol=10,nrow=0),stringsAsFactors=FALSE)
	colnames(tbl) = c("sample.name", "feature", "chr", "chr.start", "chr.end", "score", "prc.above", "prc.below", "med.ratio","n.pos")
  res <- 0
	k <- 0
	## for each gene do
	cl=parallel::makeCluster(ncpus)
	registerDoParallel(cl)
	#for(k in 1:length(unique(ano$feature))){
	tbl = foreach(k=1:length(unique(ano$feature)), .combine=rbind) %dopar%{
		feature = unique(ano$feature)[k]
		chr = paste0('chr',ano$chromosome[ano$feature==feature])
		start = ano$start[ano$feature==feature]
		end = ano$end[ano$feature==feature]

		nkr.file = file.path(nkr.dir, paste0(feature,'.RData'))
		if (file.exists(nkr.file)){
			load(nkr.file)
			lci.pos = res$lci.pos
			uci.pos = res$uci.pos
			cov.ratios = res$raw
			cov.ratios[cov.ratios==-Inf] = NA
		} else {
			stop("No feature file detected for ", feature)
		}

		## determine the score by exon
		score = (colSums(cov.ratios>uci.pos, na.rm=TRUE) - colSums(cov.ratios<lci.pos, na.rm=TRUE)) / nrow(cov.ratios)

		## generate result table
		tbl.feature = data.frame(id = names(score))
		tbl.feature$feature = feature
		tbl.feature$chr = chr
		tbl.feature$chr.start = start
		tbl.feature$chr.end = end
		tbl.feature$score = round(score,4)
		tbl.feature$prc.above = round((colSums(cov.ratios>uci.pos)/nrow(cov.ratios)) * 100,2)
		tbl.feature$prc.below = round((colSums(cov.ratios<lci.pos)/nrow(cov.ratios)) * 100,2)
		tbl.feature$med.ratio = round(apply(cov.ratios,2,median,na.rm=TRUE),4)
		tbl.feature$n.pos = nrow(cov.ratios)

		## return tbl.gene
		tbl.feature
	}
	stopCluster(cl)

	tbl.pass = tbl[tbl$score < lscore | tbl$score > uscore,]
	colnames(tbl) = colnames(tbl.pass) = c("Sample Name", "Feature", "Chromosome", "Feature Start Position", "Feature End Position", "Score", "Percent Above NKR",
			"Percent Below NKR", "Median Normalized log2 Ratio","Feature Positions")

	## write result tables to CSV files
	write.table(tbl,outfile,sep=',',quote=TRUE,row.names=FALSE,na='')
	write.table(tbl.pass,outfile.pass,sep=',',quote=TRUE,row.names=FALSE,na='')
}

