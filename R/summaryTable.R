#' Summary and count table of shRNA by sample
#'
#' Produce a summary report. Report will include data for total number of reads, 5' trimmed reads, 3' trimmed reads,
#' number and percent of mapped reads, number and percent of perfectly mapped reads, number and percent of good
#' mapped reads, number and percent of poorly mapped reads and number of unmapped reads.
#'
#' @param samkey Character. A link to sample key file. See source code for how the data should be structured.
#' @param lib Character. A link to library (fasta) file. See source code for how the data should be structured.
#' @param trimreport Character. A link to global trim report from the pipeline. See source code for how the data
#' should be structured.
#' @param mapping_path Character. Path to where mapping files reside.
#' @param mapping_rx Character. Valid regular expression for collating of mapping files (per sample).
#' @param wt_alignscores Integer. Perfect align counts. Align counts that have align scores equal to \code{wt_alignscores} and
#' read lengths equal or higher to \code{wt_readlengths_lower} will be summed.
#' @param wt_readlengths_lower Integer. Perfect align counts. Align counts that have align scores equal to \code{wt_alignscores}
#' and read lengths equal or higher to \code{wt_readlengths_lower} will be summed.
#' @param high_alignscores_lower,high_alignscores_upper,high_readlengths_lower Integer. Good align counts. Align counts with align
#' scores >= \code{high_alignscores_lower} and < \code{high_alignscores_upper} and read lengths equal or higher than
#' \code{high_readlengths_lower} will be summed.
#' @param low_alignscores_upper,low_readlengths_upper Integer. Poor align counts. Align counts below \code{low_alignscores_upper}
#' or read length below \code{low_readlengths_upper} will be excluded from the tally.
#' @param count_alignscores,count_readlengths Integer. For count matrix, align counts for which align scores equal or greater to
#' \code{count_alignscores} and read lengths equal or greater to \code{count_readlengths} will be used in the final tally.
#' @param output_report Character. A name of file (incl. file extension) into which alignment report will be written to.
#' @param output_count Character. A name of file (incl. file extension) into which count matrix (shRNA for rows, sample names
#' for columns) will be written.
#'
#' @return Two files, \code{global_alignment_report.txt} and \code{count_matrix.txt} are written as a side effect. Function
#' returns a list of length 2 where one can access alignment report (\code{output$alignment_report}) or count matrix
#' (\code{output$count_matrix}).
#'
#' @export
#' @importFrom stats aggregate reshape
#' @importFrom utils read.table write.table unstack

summaryTable <- function(samkey, lib, trimreport, mapping_path = ".", mapping_rx = "^.*_mapped_species\\.txt$",
                         wt_alignscores = 0, wt_readlengths_lower = 26,
                         high_alignscores_upper = 0, high_alignscores_lower = -6, high_readlengths_lower = 26,
                         low_alignscores_upper = -6, low_readlengths_upper = 26,
                         count_alignscores = -6, count_readlengths = 26,
                         output_report = "global_alignment_report.txt",
                         output_count = "count_matrix.txt") {

  # Import data ####
  samples <- read.table(samkey, header = TRUE)
  #        sample dox tam rep
  # 1: TB4445-17 unt unt   2
  # 2: TB4445-18 unt unt   3
  # 3: TB4445-19 unt unt   4
  # 4: TB4445-20 dox unt   1

  lib <- read.table(lib, header = FALSE, stringsAsFactors = FALSE)
  lib$seq <- 1:2
  lib <- unstack(lib)
  lib$X1 <- gsub("^>", "", lib$X1)
  names(lib) <- c("shRNA", "sequence")
  #              shRNA                      sequence
  # 1 AANAT_2131-AANAT GTATGAGGCAGCGAAACTCACTGGCTGCC
  # 2 AANAT_2740-AANAT GTATGCCACAGCAGGATGGGGCCCCTGCC
  # 3  AANAT_304-AANAT GTATAGAAGGGTACCAGCGCGTCCTTGCC
  # 4 AANAT_3349-AANAT GTATGAAGCTGAACCTCTCATAGAATGCC

  trim.data <- read.table(trimreport, header = FALSE, stringsAsFactors = FALSE, sep = NULL)
  #          V1
  # 1   sample1
  # 2 2,345,004
  # 3 2,319,009
  # 4 1,718,789

  trim.data$seq <- 1:4
  trim.data <- unstack(trim.data)

  # In case only one sample is handled, unstack will not behave as expected and
  # t() is needed.
  if (all(nrow(trim.data) == 4 & ncol(trim.data) == 1)) {
    trim.data <- data.frame(t(trim.data))
  }

  rownames(trim.data) <- NULL

  trim.data[, 2:4] <- sapply(trim.data[, 2:4], gsub, pattern = ",", replacement = "")
  colnames(trim.data) <- c("sample", "total", "trimmed_5", "trimmed_3")
  trim.data$total <- as.numeric(trim.data$total)
  trim.data$trimmed_5 <- as.numeric(trim.data$trimmed_5)
  trim.data$trimmed_3 <- as.numeric(trim.data$trimmed_3)
  #      sample   total trimmed_5 trimmed_3
  # 1 TB4445-17 2345004   2319009   1718789
  # 2 TB4445-18 2603252   2576740   1927294
  # 3 TB4445-19 3358655   3336617   2670632
  # 4 TB4445-20 3260564   3198344   1671554
  # 5 TB4445-21 2078391   2059003   1527371
  # 6 TB4445-22 2680203   2646307   1779339

  # Import align counts and mapped shRNAs ####
  file.samples <- list.files(path = mapping_path, pattern = mapping_rx, full.names = TRUE)

  xy <- lapply(file.samples, FUN = function(i) {
    ri <- read.table(i, header = FALSE, sep = "", stringsAsFactors = FALSE)
    colnames(ri) <- c("align_counts", "mapped_shRNAs", "sequences", "align_scores")
    ri$read_lengths <- nchar(ri$sequences)
    ri$align_scores <- as.numeric(gsub("^AS:i:", "", ri$align_scores))
    ri
  })
  names(xy) <- samples$sample
  # This result is a list, where each element is mapped shRNAs for a given sample.
  # One sample may look something like this:
  #   align_counts    mapped_shRNAs                     sequences align_scores read_lengths
  # 1           28 AANAT_2131-AANAT GTATGAGGCAGCGAAACTCACTGGCTGCC            0           29
  # 2            3 AANAT_2131-AANAT   GGCAGCGAGTACAAATCGACTGGATAC          -16           27
  # 3            1 AANAT_2131-AANAT GTATGAGGCAGCGAAACCCACTGGCTGCC           -2           29
  # 4            1 AANAT_2131-AANAT GTATGAGGCAGCGAAACGCACTGGCTGCC           -2           29
  # 5           30 AANAT_2740-AANAT      GTATCAGAACATCAGGATGGTGCC          -12           24
  # 6           20 AANAT_2740-AANAT   GTATTTTCCACAGCATGACTCGGTGCC          -14           27

  # Calculate statistics. ####
  trim.data$mapped <- sapply(xy, FUN = function(x) sum(x$align_counts, na.rm = TRUE))
  trim.data["mapped %"] <- signif(with(trim.data, mapped / total), 3)

  trim.data$perfect <- sapply(xy, FUN = function(x) sum(x[x$align_scores == wt_alignscores
                                                          & x$read_lengths >= wt_readlengths_lower, "align_counts"],
                                                        na.rm = TRUE))
  trim.data["perfect %"] <- signif(with(trim.data, perfect / total), 3)

  trim.data$good <- sapply(xy, FUN = function(x) sum(x[x$align_scores < high_alignscores_upper
                                                       & x$align_scores >= high_alignscores_lower
                                                       & x$read_length >= high_readlengths_lower, "align_counts"],
                                                     na.rm = TRUE))
  trim.data["good %"] <- signif(with(trim.data, good / total), 3)

  trim.data$poor <- sapply(xy, FUN = function(x) sum(x[x$align_scores < low_alignscores_upper
                                                       | x$read_lengths < low_readlengths_upper,
                                                       "align_counts"],
                                                     na.rm = TRUE))
  trim.data["poor %"] <- signif(with(trim.data, poor / total), 3)

  trim.data$unmapped <- with(trim.data, total - mapped)
  trim.data["unmapped %"] <- signif(with(trim.data, unmapped / total), 3)

  # sample   total trimmed_5 trimmed_3  mapped mapped % perfect perfect %   good good %   poor poor % unmapped unmapped %
  # 1 TB4445-17 2345004   2319009   1718789 1712494    0.730  970717     0.414 343041 0.1460 398736  0.170   632510      0.270
  # 2 TB4445-18 2603252   2576740   1927294 1897598    0.729 1086360     0.417 378654 0.1450 432584  0.166   705654      0.271
  # 3 TB4445-19 3358655   3336617   2670632 2630073    0.783 1490566     0.444 530136 0.1580 609371  0.181   728582      0.217
  # 4 TB4445-20 3260564   3198344   1671554 1637325    0.502  938649     0.288 321022 0.0985 377654  0.116  1623239      0.498
  # 5 TB4445-21 2078391   2059003   1527371 1503861    0.724  855002     0.411 303352 0.1460 345507  0.166   574530      0.276
  # 6 TB4445-22 2680203   2646307   1779339 1749226    0.653  998534     0.373 353688 0.1320 397004  0.148   930977      0.347

  # Write result of global alignment report. ####
  write.table(trim.data, file = output_report, sep = '\t', row.names = FALSE)

  # Create count matrix ####
  # Create new object which will be one big table of all the data. Add a new column
  # which will designate from which sample individual shRNA comes from. This can now
  # be used for further calculation using R-way logic.
  cnm <- do.call(rbind, xy)
  rownames(cnm) <- NULL

  cnm$sample <- rep(names(xy), times = sapply(xy, nrow))

  # Keep only those shRNAs which satisfy the criteria.
  cnm <- cnm[cnm$align_scores >= count_alignscores & cnm$read_lengths >= count_readlengths, ]

  # Aggreggate align_counts (sum) according to sample and shRNA.
  cnm.agg <- aggregate(align_counts ~ mapped_shRNAs + sample, data = cnm, FUN = sum)
  #      mapped_shRNAs    sample align_counts
  # 1 AANAT_2131-AANAT TB4445-17           30
  # 2 AANAT_2740-AANAT TB4445-17            5
  # 3  AANAT_304-AANAT TB4445-17           47
  # 4 AANAT_3349-AANAT TB4445-17          244
  # 5 AANAT_3958-AANAT TB4445-17          114
  # 6 AANAT_4567-AANAT TB4445-17           12

  # Reshape data to wide format so that rows are shRNAs and columns are samples.
  cnm.agg <- reshape(data = cnm.agg, timevar = "sample", idvar = "mapped_shRNAs", direction = "wide")
  names(cnm.agg) <- gsub("align_counts\\.", "", names(cnm.agg))
  # mapped_shRNAs TB4445-17 TB4445-18 TB4445-19 TB4445-20 ...
  # 1 AANAT_2131-AANAT        30        77         4      ...
  # 2 AANAT_2740-AANAT         5        13        15      ...
  # 3  AANAT_304-AANAT        47        26        73      ...
  # 4 AANAT_3349-AANAT       244       297       181      ...
  # 5 AANAT_3958-AANAT       114        90       237      ...
  # 6 AANAT_4567-AANAT        12        26        39      ...

  # This step is needed so that all shRNAs that appear in output of mapped reads are included.
  cnm.agg <- merge(x = lib[, "shRNA", drop = FALSE], y = cnm.agg,
                   by.x = "shRNA", by.y = "mapped_shRNAs",
                   all.x = TRUE)

  # Note that for some samples, shRNA may be missing, in which case, reshape will add an NA value.
  # This can indicate two things - that there were no reads for this shRNA OR they have been
  # filtered out. In any case, the below line will replace all NAs to keep up with the specs.
  cnm.agg[is.na(cnm.agg)] <- 0

  # Write data to file and return a list of calculated data.
  write.table(cnm.agg, file = output_count, sep = "\t", row.names = FALSE, quote = FALSE)

  return(list(alignment_report = trim.data, count_matrix = cnm.agg))
}
