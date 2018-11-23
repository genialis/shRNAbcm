context("Testing summaryTable")

file.sample <- "sample1_test_spec.txt"
file.samplekey <- "samplekey.txt"
file.lib <- "library.txt"
file.trimreport <- "trimreport.txt"
tar.txt <- "tar.txt"
tcr.txt <- "tcr.txt"

# On exit, the test will clean after itself.
on.exit({
  unlink(c(file.sample, file.samplekey, file.lib, file.trimreport, tar.txt, tcr.txt))
}, add = TRUE)

smp <- read.table(text = "      1 AANAT_2131-AANAT	GTATAAGGCAGCGATGGTGAGCTGCC	AS:i:-14
      1 AANAT_2131-AANAT	GTATGAGGCAGCGAAACTCACTGGCTGCC	AS:i:0
      1 AANAT_2740-AANAT	GTATTACCTTCAACGATGGTGCCCCTGCC	AS:i:-16
      1 AANAT_2740-AANAT	GTATTTATCCACAGCAATGATCGGTGCC	AS:i:-16
      1 AANAT_2740-AANAT	GTATTTCCACAGCAATGACTCGGTGCC	AS:i:-14
      1 AANAT_2740-AANAT	GTTAAACATGATGGGTCCATGCTGCC	AS:i:-18
     60 AANAT_304-AANAT	GTATAGAAGGGTACCAGCGCGTCCTTGCC	AS:i:0
     10 AANAT_304-AANAT	GTATAGAGGGGTACCAGCGCGTCCTTGC	AS:i:-2
      6 AANAT_304-AANAT	GTATAGAGGGGTACCAGCGCGTCCTTGCC	AS:i:-2
      1 AANAT_304-AANAT	GTATAGAAGGGTACCAGCGCGTACTTGCC	AS:i:-2
      1 AANAT_304-AANAT	GTATAGAAGGGTACCAGCGCGTCATTGCC	AS:i:-2
      1 AANAT_304-AANAT	GTATAGAAGGGTACCAGCGCGTCCTTGTC	AS:i:-2
      1 AANAT_304-AANAT	GTATAGAAGGGTACCGGCGCGTCCTTGCC	AS:i:-2
    132 AANAT_3349-AANAT	GTATGAAGCTGAACCTCTCATAGAATGCC	AS:i:0
     24 AANAT_3349-AANAT	GTATGGAATAGAAGTCTCAAAGAGTACC	AS:i:-18
      9 AANAT_3349-AANAT	GTATGGAATAGAAGTCTCAAAGAGTGCC	AS:i:-16
      5 AANAT_3349-AANAT	GTATGAAGCTGAACCTCTCATAGATGCC	AS:i:-2
      3 AANAT_3349-AANAT	GTATGAAGCTGAACCACTCATAGAATGCC	AS:i:-2
      1 AANAT_3349-AANAT	ATATGAAGCTGAACCTCTCATAGAATGCC	AS:i:-2
      1 AANAT_3349-AANAT	GCAGTGCAAGATGAAGGACTCTCATAC	AS:i:-16
      1 AANAT_3349-AANAT	GTATGAAGCGGAACCGCTCATAGAATGCC	AS:i:-4
      1 AANAT_3349-AANAT	GTATGAAGCTGAACCTCTCATGGAATGCC	AS:i:-2
      1 AANAT_3349-AANAT	GTATGAAGCTGAAGCTCTCATAGAATGCC	AS:i:-2",
                 header = FALSE, sep = "\t")

samplekey <- read.table(text = "sample	dox	tam	rep
sample1	unt	unt	2", header = TRUE)

# Write temporary files.
write.table(smp, file = file.sample, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(samplekey, file = file.samplekey, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
cat(file = file.lib, ">AANAT_2131-AANAT
GTATGAGGCAGCGAAACTCACTGGCTGCC
>AANAT_2740-AANAT
GTATGCCACAGCAGGATGGGGCCCCTGCC
>AANAT_304-AANAT
GTATAGAAGGGTACCAGCGCGTCCTTGCC
>AANAT_3349-AANAT
GTATGAAGCTGAACCTCTCATAGAATGCC\n")
cat(file = file.trimreport, "sample1
2,345,004
2,319,009
1,718,789\n")

out <- summaryTable(samkey = file.samplekey,
                    lib = file.lib,
                    mapping_rx = "^.*_test_spec\\.txt$",
                    trimreport = file.trimreport,
                    output_report = tar.txt, # test alignment report
                    output_count = tcr.txt) # tst count matrix

# compare against this
out.gs <- structure(list(alignment_report = structure(list(sample = structure(1L, .Label = "sample1", class = "factor"),
                                                           total = 2345004, trimmed_5 = 2319009, trimmed_3 = 1718789,
                                                           mapped = 264L, `mapped %` = 0.000113, perfect = 193L, `perfect %` = 8.23e-05,
                                                           good = 32L, `good %` = 1.36e-05, poor = 39L, `poor %` = 1.66e-05,
                                                           unmapped = 2344740, `unmapped %` = 1), row.names = c(NA,
                                                                                                                -1L), .Names = c("sample", "total", "trimmed_5", "trimmed_3",
                                                                                                                                 "mapped", "mapped %", "perfect", "perfect %", "good", "good %",
                                                                                                                                 "poor", "poor %", "unmapped", "unmapped %"), class = "data.frame"),
                         count_matrix = structure(list(shRNA = c("AANAT_2131-AANAT",
                                                                 "AANAT_2740-AANAT", "AANAT_304-AANAT", "AANAT_3349-AANAT"
                         ), sample1 = c(1, 0, 80, 144)), .Names = c("shRNA", "sample1"
                         ), row.names = c(NA, -4L), class = "data.frame")), .Names = c("alignment_report",
                                                                                       "count_matrix"))

test_that("check the working of summaryTable", {
  expect_identical(object = out, expected = out.gs)
})
