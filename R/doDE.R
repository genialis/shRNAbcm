#' Make differential expression analysis
#'
#' Input a set of count matrices, define relationships and perform a differential analysis.
#' Return several outputs, including differentil expression result and some extra  comparisons.
#'
#' @param input Character. A relative or absolute path to a parameters file. See Details.
#' @param sample_list Character. Vector of absolute or relative paths to sample expressions or
#' folder containing sample expressions. See \code{importExpressionData} for more info.
#'
#' @details
#' \code{input} file should be an excel (xlsx) file with tabs sample_key, contrasts,
#' overall_contrasts and classification_parameters. \code{sample_key} should have columns
#' sample, <variables> (depends on the experiment) and replicate. \code{contrasts} should
#' have columns group_1 and group_2, and they should hold which groups to test against each
#' other. Sheet \code{overall_contrasts} should have groups of overall contrasts in columns,
#' where rows represent pairs of comparisons. \code{classification_parameters} sheet should
#' have twothree columns, threshold, value and description. See sample file for details.
#' @return Files are returned: \code{deseq_results.txt}, \code{class_results.txt},
#' \code{beneficial_counts.txt}, \code{lethal_counts.txt}.
#' For usage, check out test file.
#' @export

doDE <- function(input, sample_list = NULL) {
  # Import files for sample key and count matrix.
  ## It should look something like the table below. Columns should be sample,
  ## <variables of the sample> and a replicate.
  ## IMPORTANT: First column should be <sample>, second to last column should replicate
  ## and last column should be group, which is a concatenation of variables using
  ## character underscore (_).

  #      sample treatment    time   line  egf replicate                   group
  # 1 TB4693-01       unt initial  gata3 high         1  unt_initial_gata3_high
  # 2 TB4693-02       unt initial  gata3 high         2  unt_initial_gata3_high
  # 3 TB4693-03       unt initial  gata3 high         3  unt_initial_gata3_high
  # 4 TB4693-04       unt initial  gata3 high         4  unt_initial_gata3_high
  # 5 TB4693-05       unt initial ptpn12 high         1 unt_initial_ptpn12_high
  # 6 TB4693-06       unt initial ptpn12 high         2 unt_initial_ptpn12_high
  sam.key <- importSampleKey(x = input)

  ## Import count data based on sample names from sample key.
  ## Fetch only data for specified samples mentioned in sample key constructed above.
  ## Contrasts should be added to sample names like depicted below. Rownames are
  ## shRNA-gene concatenate.

  #                     TB4693_01_unt_initial_gata3_high TB4693_02_unt_initial_gata3_high
  # SLC38A1_1_6391-NAT2                             6575                             5949
  # SLC38A1_1_6392-NAT2                              703                              494
  # SLC38A1_1_6393-NAT2                             1104                             1090
  # SLC38A1_1_6394-NAT2                             1541                              915
  # SLC38A1_1_6395-NAT2                              488                              499
  # SLC38A1_1_6396-NAT2                              175                              169
  expr.data <- importExpressionsData(sample_key = sam.key, sample_list = sample_list)

  # Prepare shRNA and gene names to be used for construction of DE results.
  shrna <- fetchGeneOrHairpin(x = expr.data$shRNA, entity = "shrna")
  gene.names <- fetchGeneOrHairpin(x = expr.data$shRNA, entity = "gene")

  # Prepare list of contrasts to be calculated.
  # Output is a data.frame of contrasts (for each group).
  #   group_1 group_2
  # 1 unt_unt dox_unt
  # 2 unt_tam dox_unt
  # 3 unt_tam dox_tam
  contrast.list <- importContrastList(x = input)

  # Perform differential expression.
  # The result is a list where each list holds results from DESeq2::results().

  # List of 3
  # $ unt_unt.dox_unt:Formal class 'DESeqResults' [package "DESeq2"] with 7 slots
  # .. ..@ priorInfo      :List of 3
  # .. .. ..$ type   : chr "none"
  # .. .. ..$ package: chr "DESeq2"
  # .. .. ..$ version:Classes 'package_version', 'numeric_version'  hidden list of 1
  # .. .. .. ..$ : int [1:3] 1 18 1
  # .. ..@ rownames       : NULL
  # ...
  de.result <- differentialExpression(contrasts = contrast.list, sample_key = sam.key, data = expr.data)

  # Prepare result to contain shRNA, gene and expression result.
  #        shRNA  gene l2fc_unt_unt.dox_unt lfcSE_unt_unt.dox_unt padj_unt_unt.dox_unt
  # 1 AANAT_2131 AANAT           0.61422268             0.8739565            0.9124698
  # 2 AANAT_2740 AANAT          -0.57193511             0.9758090            0.9263345
  # 3  AANAT_304 AANAT           0.91803090             0.6898804            0.7399504
  # 4 AANAT_3349 AANAT           0.02834944             0.3752039            0.9872124
  # 5 AANAT_3958 AANAT           0.28466804             0.4058473            0.9125632
  de.scrape <- scrapeDEresults(dds = de.result, shrna = shrna, gene = gene.names)

  # Write results to a file.
  # TODO: if pool turns out to be an important issue, this information will have to be appended to a filename
  write.table(de.scrape, file = "deseq_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")

  # Trim results of NAs.
  de.trim <- trimResults(de = de.scrape)

  # Import overall contrasts and parameters needed to classify.
  contrast.overall <- importOverallContrasts(x = input)
  param <- importClassificationParams(x = input, trim = TRUE)

  # Classify results based on provided treshold values.
  de.class <- classifyDE(de = de.trim, param = param, ctr = contrast.list, overall_ctr = contrast.overall)

  # Create a count matrix of classified results per each gene.
  ## lethal counts
  ##  beneficial counts

  de.counted.lethal <- countClassified(count = de.class, ctr = contrast.list, overall_ctr = contrast.overall, classy = "lethal")
  de.counted.beneficial <- countClassified(count = de.class, ctr = contrast.list, overall_ctr = contrast.overall, classy = "beneficial")

  write.table(de.class, file = "class_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(de.counted.beneficial, file = "beneficial_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(de.counted.lethal, file = "lethal_counts.txt", sep = "\t", row.names = FALSE, quote = FALSE)

  # This function is used for its side-effects of saving (intermediate) results to text files.
  return(NULL)
}

