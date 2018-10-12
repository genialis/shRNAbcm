#' Create library of shRNA
#'
#' Take a raw file and convert it to a valid fasta format. This file is necessary
#' for performing alignment.
#'
#' @param input A raw file of three columns (\code{AANAT_1522	AANAT	GTATGGGACTCGGGGATCCCAGGTGTGCC})
#' where first value is shRNA, second value is gene name and third value is sequence.
#' @param output Character. Optional. If a character string of length one is supplied, it is
#' stripped of file extension and appended \code{fasta}. If \code{NULL}, name of input is used,
#' where \code{.fasta} extension is appended and written to working directory.
#' @export
#' @author Roman Lustrik (roman@@genialis.com)

makeLibrary <- function(input, output = NULL) {
  # Import data.
  xy <- read.table(input, sep = "", header = FALSE)

  # Define how a single read for fasta is written.
  out <- sprintf("> %s-%s\n%s", xy[, 1], xy[, 2], xy[, 3])

  # Prepare file name. Use default unless specified in output. Strip
  # of file extension and append .fasta when writing to file.
  mycon <- tools::file_path_sans_ext(basename(input))
  if (!is.null(output)) mycon <- tools::file_path_sans_ext(output)

  # Write result to file.
  writeLines(out, con = paste(mycon, ".fasta", sep = ""))

  invisible(paste(mycon, ".fasta", sep = ""))
}
