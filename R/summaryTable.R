#' Summary table
#'
#' Produce a summary report.
#'
#' @param samkey
#' @param lib
#' @param trimreport
#'
#' @retrun Two files, \code{global_alignment_report.txt} and \code{count_matrix.txt}.
#'
#' @export
#' @importFrom data.table fread
#'
#' @author Roman Lustrik (roman@@genialis.com), original idea Nicholas Neill (nicholas.neill@@bcm.edu )

summaryTable <- function(samkey, lib, trimreport, wt_alignscores = 0, wt_readlengths_lower = 26,
                         high_alignscores_upper = 0, high_alignscores_lower = -6,
                         high_readlengths_lower = 26, low_alignscores_upper = -6,
                         low_readlengths_upper = 26, output = "global_alignment_report.txt") {
  # Import data ####
  samples <- fread(samkey)
  #        sample dox tam rep
  # 1: TB4445-17 unt unt   2
  # 2: TB4445-18 unt unt   3
  # 3: TB4445-19 unt unt   4
  # 4: TB4445-20 dox unt   1

  lib <- fread(lib, header = FALSE, stringsAsFactors = FALSE)
  lib$seq <- 1:2
  lib <- unstack(lib)
  lib$X1 <- gsub("^>", "", lib$X1)
  names(lib) <- c("shRNA", "sequence")
  #              shRNA                      sequence
  # 1 AANAT_2131-AANAT GTATGAGGCAGCGAAACTCACTGGCTGCC
  # 2 AANAT_2740-AANAT GTATGCCACAGCAGGATGGGGCCCCTGCC
  # 3  AANAT_304-AANAT GTATAGAAGGGTACCAGCGCGTCCTTGCC
  # 4 AANAT_3349-AANAT GTATGAAGCTGAACCTCTCATAGAATGCC

  trim.data <- fread(trimreport, header = FALSE, stringsAsFactors = FALSE, sep = NULL)
  trim.data$seq <- 1:4
  trim.data <- unstack(trim.data)
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
  file.samples <- paste(samples$sample, "mapped_species.txt", sep = "_")

  xy <- lapply(file.samples, FUN = function(i) {
    ri <- read.table(i, header = FALSE, sep = "", stringsAsFactors = FALSE)
    colnames(ri) <- c("align_counts", "mapped_shRNAs", "sequences", "align_scores")
    ri$read_lengths <- nchar(ri$sequences)
    ri$align_scores <- as.numeric(gsub("^AS:i:", "", ri$align_scores))
    ri
  })
  names(xy) <- samples$sample

  trim.data$mapped <- sapply(xy, FUN = function(x) sum(x$align_counts))
  trim.data["mapped %"] <- signif(with(trim.data, mapped / total), 3)

  trim.data$perfect <- sapply(xy, FUN = function(x) sum(x[x$align_scores == wt_alignscores
                                           & x$read_lengths >= wt_readlengths_lower, "align_counts"]))
  trim.data["perfect %"] <- signif(with(trim.data, perfect / total), 3)

  trim.data$good <- sapply(xy, FUN = function(x) sum(x[x$align_scores < high_alignscores_upper
                                               & x$align_scores >= high_alignscores_lower
                                               & x$read_length >= high_readlengths_lower, "align_counts"]))
  trim.data["good %"] <- signif(with(trim.data, good / total), 3)

  trim.data$poor <- sapply(xy, FUN = function(x) sum(x[x$align_scores < low_alignscores_upper | x$read_lengths < low_readlengths_upper, "align_counts"]))
  trim.data["poor %"] <- signif(with(trim.data, poor / total), 3)

  trim.data$unmapped <- with(trim.data, total - mapped)
  trim.data["unmapped %"] <- signif(with(trim.data, unmapped / total), 3)

  write.table(trim.data, file = output, sep = '\t', row.names = FALSE)

  trim.data
}
#
## generate a table of read species statistics for each hairpin
## $1 is path to sample key
## $2 is path to library file (trimmed fasta)
## $3 is path to global trim report

# args <- c(sample = "sample_key_gata3.txt",
# library = "chromatin_library.fasta",
# trimreport = "global_trim_report.txt"
# )

# samples = as.matrix(read.table(args[1], header = T, sep = '\t', stringsAsFactors = F))
# sample_id = samples[,1]
# library = read.table(args[2], sep = '\t', stringsAsFactors = F)
# library = cbind(library, rep(c(1,2)))
# library = unstack(library)
# shrna = library[,1]
# shrna = gsub('>', '', shrna, fixed = T)

#make alignment report for each sample
# print('Generating alignment report.')
# trim_data = read.table(args[3], stringsAsFactors = F)
# trim_data[[2]] = rep(1:4, times = nrow(samples))
# trim_data = unstack(trim_data)
# total = as.numeric(gsub(',', '', trim_data[[2]], fixed = T))
# trimmed_5 = as.numeric(gsub(',', '', trim_data[[3]], fixed = T))
# trimmed_3 = as.numeric(gsub(',', '', trim_data[[4]], fixed = T))
# align_counts = lapply(sample_id, function(x) as.numeric(read.table(paste(x, '_mapped_species.txt', sep = ''))[[1]]))
# mapped_shrnas = lapply(sample_id, function(x) read.table(paste(x, '_mapped_species.txt', sep = ''))[[2]])
# read_seqs = lapply(sample_id, function(x) as.character(read.table(paste(x, '_mapped_species.txt', sep = ''))[[3]]))
# read_lengths = lapply(read_seqs, nchar)
# align_scores = lapply(sample_id, function(x) read.table(paste(x, '_mapped_species.txt', sep = ''))[[4]])
# align_scores = lapply(align_scores, function(x) as.numeric(gsub('AS:i:', '', x, fixed = T)))
# mapped = sapply(align_counts, sum)
# wt = mapply(function(x,y,z) sum(x[y == 0 & z >= 26]), align_counts, align_scores, read_lengths)
# high_q = mapply(function(x,y,z) sum(x[y < 0 & y >= -6 & z >= 26]), align_counts, align_scores, read_lengths)
# low_q = mapply(function(x,y,z) sum(x[y < -6 | z < 26]), align_counts, align_scores, read_lengths)
# unmapped = total - mapped
# mapped_pc = signif(mapped / total, 3)
# wt_pc = signif(wt / total, 3)
# high_q_pc = signif(high_q / total, 3)
# low_q_pc = signif(low_q / total, 3)
# unmapped_pc = signif(unmapped / total, 3)
# summary = data.frame(sample_id,  total,   trimmed_5  , trimmed_3,   mapped,  mapped_pc,      wt,        wt_pc,   high_q, high_q_pc, low_q, low_q_pc,  unmapped, unmapped_pc)
# colnames(summary) = c('sample', 'total', 'trimmed_5', 'trimmed_3', 'mapped', 'mapped_%', 'perfect', 'perfect_%', 'good', 'good_%', 'poor', 'poor_%', 'unmapped', 'unmapped_%')





# write.table(summary, file = 'original_global_alignment_report.txt', sep = '\t', row.names = F)
#
# #make count matrix
# print('Generating count matrix.')
# good_ind = mapply(function(x,y) which(x >= -6 & y >= 26), align_scores, read_lengths)
# good_counts = mapply(function(x,y) x[y], align_counts, good_ind)
# good_mappings = mapply(function(x,y) x[y], mapped_shrnas, good_ind)
# count_matrix = mapply(function(x,y) sapply(shrna, function(z) sum(x[y == z])), good_counts, good_mappings)
# colnames(count_matrix) = sample_id
#
# # TODO: should this matrix be saved for chromatin, kinase and ubiquitin?
# write.table(count_matrix, file = 'original_count_matrix.txt', sep = '\t')
