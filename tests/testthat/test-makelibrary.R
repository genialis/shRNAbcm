context("Testing makeLibrary")

library(shRNAbcm)

# write to these two files (.txt and .fasta)
fl <- "temp_lib.txt"
fl.in <- paste(tools::file_path_sans_ext(fl), ".fasta", sep = "")

td <- data.frame(V1 = c("AANAT_1522", "AANAT_2131", "AANAT_2740"),
                 V2 = c("AANAT", "AANAT", "AANAT"),
                 V3 = c("GTATGGGACTCGGGGATCCCAGGTGTGCC", "GTATGAGGCAGCGAAACTCACTGGCTGCC", "GTATGCCACAGCAGGATGGGGCCCCTGCC"))

write.table(td, file = fl, col.names = FALSE, row.names = FALSE, quote = FALSE)
makeLibrary(input = fl)
unlink(fl)
xy <- readLines(fl.in)
unlink(fl.in)

test_that("correct layout", {
  # test that there are indeed 6 lines read in
  expect_length(xy, 6)
  # test that elements 1, 3 and 5 start with >
  expect_true(all(grepl("^>", xy[seq(from = 1, to = 6, by = 2)])))
  # test that even elements contain a sequence identical to one from td
  expect_true(all(as.character(td$V3) == xy[seq(from = 2, to = 6, by = 2)]))
})