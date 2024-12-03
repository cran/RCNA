#' Example RCNA_object
#'
#' An example RCNA object used to run examples and demonstrate the structure of the custom S4 object provided in this package.
#' This example uses a dummy feature ("feature_a") and three coverage files which were subset to the length of the dummy feature to be concise and quick to run.
#' An annotation file has been included in the `inst/` directory along with the coverage files. This object is compiled in the `\code{create_RCNA_object}` function's documentation.
#' In order to use this example, you should make the following replacements:
#' \code{
#' example_obj@ano.file <- system.file("examples" ,"annotations-example.csv", package = "RCNA")
#' raw.cov <- system.file("examples", "coverage", paste0(samples, ".txt.gz"), package = "RCNA")
#' example_obj@gcParams$file.raw.coverage <- raw.cov
#' }
#'
#' This will set the location of the example annotation file and the example raw coverage files to the flat files included with the package.
#'
#' @format An RCNA object created using create_RCNA_object(). See \link{create_RCNA_object} for more details on the slots of this object.
#'
#'
"example_obj"
