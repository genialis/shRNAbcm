context("Testing doDE helper functions")

# Prepare inputs.
inputfile <- system.file("extdata", "template_doDE_inputs.xlsx", package = "shRNAde", mustWork = TRUE)

# Test workflow.
# Simulate three genes, three shRNA each and 12 samples. A simple experimental design is t0 and tend (time end),
# with control and treatments for each group.

set.seed(357)
N <- 12  # number of samples
N.genes <- 30  # number of genes; should be sufficient for tractability of DESeq()
gene <- rep(sprintf("%s%s", "GENE", 1:N.genes), each = 4)  # dictates number of genes
shrna <- sprintf("%s_%s-%s", gene, 1:length(gene), gene)

# Prepare expression data for N samples.
smps <- matrix(rbinom(n = N * length(gene), size = 1000, prob = 0.1), ncol = N)
smps[sample(1:length(smps), size = round(N * length(gene) * 0.1))] <- 0

cts <- data.frame(shRNA = shrna, smps)
colnames(cts) <- c("shRNA", paste("sample", 1:N, sep = "-"))

# Create files and mark them for removal after the test.
dir.create("./tempdata", showWarnings = FALSE)
on.exit({unlink("./tempdata", recursive = TRUE)}, add = TRUE)

for (i in 2:ncol(cts)) {
  samplename <- sprintf("./tempdata/%s_trimmed_trimmed_count_matrix.txt", colnames(cts[i]))
  on.exit({unlink(samplename)}, add = TRUE)  # clean after oneself, but only on exit

  write.table(cts[, c(1, i)], file = samplename, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Remove results that are produced during the workflow.
ben.cnt <- "beneficial_counts.txt"
deseq.res <- "deseq_results.txt"
class.res <- "class_results.txt"
let.cnt <- "lethal_counts.txt"

on.exit(
  {unlink(c(ben.cnt, deseq.res, class.res, let.cnt))},
  add = TRUE)

import.contrast <- shRNAde:::importContrastList(x = inputfile)
import.overall <- shRNAde:::importOverallContrasts(x = inputfile)

out <- doDE(input = inputfile, sample_list = "./tempdata")

beneficial.counts <- read.table(ben.cnt, header = TRUE)
lethal.counts <- read.table(let.cnt, header = TRUE)
class.results <- read.table(class.res, header = TRUE)
deseq.results <- read.table(deseq.res, header = TRUE)

test_that("perform differential expression", {
  # Check DE results.
  expect_equal(nrow(beneficial.counts), N.genes)
  # Number of columns should be equal to number of contrasts + number of overall contrasts.
  # +1 is accounting for gene name.
  expect_equal(ncol(beneficial.counts), sum(nrow(import.contrast), ncol(import.overall)) + 1)
  expect_equal(sum(beneficial.counts[-1]), 12)

  expect_equal(nrow(lethal.counts), N.genes)
  expect_equal(ncol(lethal.counts), sum(nrow(import.contrast), ncol(import.overall)) + 1)
  expect_equal(sum(lethal.counts[-1]), 15)
})
