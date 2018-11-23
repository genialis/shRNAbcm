#' Count classified results
#'
#' For each contrast, count how many times neutral, lethal or beneficial species comes up.
#'
#' @param count Data.frame
#' @param ctr Data.frame of contrasts for which counting should be done.
#' @param overall_ctr Data.frame of overall contrasts.
#' @param classy Which classification to count. Possible options are \code{lethal} and \code{beneficial}.
#'
#' @return Matrix of counts where rows are genes and columns are contrasts.
countClassified <- function(count, ctr, overall_ctr, classy) {
  # For each contrast, calculate how many classified values are per gene.
  indiv <- sapply(ctr$contrast, FUN = function(x, count, classy) {
    xy <- count[, colnames(count) %in% c("gene", x)]  # make sure only relevant contrast columns are present
    table(xy, useNA = "ifany")  # for each gene calculates number of beneficial, neutral and lethal
  }, count = count, simplify = FALSE)
  # Subset only specified classy variable (e.g. beneficial).
  indiv <- as.data.frame(sapply(indiv, FUN = function(x, classy) x[, colnames(x) %in% classy], classy = classy))
  indiv$gene <- rownames(indiv)
  rownames(indiv) <- NULL

  # Calculate number of flags per gene for overall contrasts.
  overall <- sapply(names(overall_ctr), FUN = function(x, count) {
    xy <- count[, colnames(count) %in% c("gene", x)]
    table(xy, useNA = "ifany")
  }, count = count, simplify = FALSE)
  overall <- as.data.frame(sapply(overall, FUN = function(x, classy) x[, colnames(x) %in% classy], classy = classy))
  overall$gene <- rownames(overall)
  rownames(overall) <- NULL

  # Construct a matrix of genes onto which values from indiv and overall will be pasted.
  gene.matrix <- data.frame(gene = unique(count$gene))
  out <- Reduce(function(x, y) merge(x, y, by = "gene", sort = FALSE), x = list(gene.matrix, indiv, overall))

  out
}

#' Classify DE result
#'
#' Given DE data for various contrasts (\code{ctr}), each of the contrasts will be classified as lethal, beneficial
#' or neutral from l2fc and p-value based on parameters passed in from \code{param} (\code{l_thr}, \code{p_thr}).
#' Based on a group of (overall) contrasts (\code{overall_ctr}) number of beneficial, lethal or neutral shRNAs in a group of
#' contrasts is counted and classified as lethal, beneficial and neutral using \code{overall_thr} parameter.
#'
#' @param de An object of ("trimmed") results in the form as returned by scrapeDEresults.
#' @param param A data.frame of parameters used during classification.
#' @param ctr A data.frame of contrasts used for constructing DE.
#' @param overall_ctr A data.frame of overall contrasts
#'
#' @return A data.frame with classified results.
#' @importFrom stats na.omit
classifyDE <- function(de, param, ctr, overall_ctr) {

  # For each contrast, classify shRNA values based on param threshold values.
  ctr.classif <- sapply(ctr$contrast, FUN = function(x, de, l_thr, p_thr) {
    # Subset columns for appropriate contrast
    xy <- de[, grepl(x, colnames(de))]
    stopifnot(ncol(xy) == 3)  # make sure that only three columns, l2fc, lfcSE and padj are subsetted

    # Preapare result into which to write classifications.
    out <- rep("neutral", nrow(xy))

    # Subset relevant columns prior to doing any (hand) classification.
    l2fc <- xy[, grepl("l2fc", colnames(xy))]
    pval <- xy[, grepl("padj", colnames(xy))]

    # Classify based on threshold values provided by user.
    out[l2fc < -l_thr & pval < p_thr] <- "lethal"
    out[l2fc > l_thr & pval < p_thr] <- "beneficial"
    # the rest are neutral

    # Set levels explicitly in case one doesn't come up.
    out <- factor(out, levels = c("lethal", "beneficial", "neutral"))

    out
  }, de = de, p_thr = fetchTHvalue(param, th = "p_thr"), l_thr = fetchTHvalue(param, th = "l_thr"), simplify = FALSE)
  # Note that above sapply is set to simplify = FALSE. Coercing to data.frame preserves factor levels.
  ctr.classif <- as.data.frame(ctr.classif)

  # Foc the classified (ctr.classif) values, count how many times they occur by shRNA and classify
  # number of occurrences by a group (overall) of contrasts.
  overall.classif <- sapply(overall_ctr, FUN = function(x, de, count_th) {
    # Prepare objects.
    x <- na.omit(x)  # remove in case there are NAs which would  mess things up
    de <- de[, colnames(de) %in% x, drop = FALSE]

    # Classify based on provided threshold. If some lethal/beneficial modification appears count_th times,
    # consider it as such, otherwise "overall neutral". E.g. if "lethal" appears two times and threshold
    # is three, consider it as "neutral". If it appears three or more times, consider it "lethal".
    out <- rep("neutral", times = nrow(de))

    lethal <- de == "lethal"
    beneficial <- de == "beneficial"

    out[rowSums(lethal) >= count_th] <- "lethal"
    out[rowSums(beneficial) >= count_th] <- "beneficial"
    # the rest are assumed to be neutral

    # Set levels explicitly in case one doesn't come up.
    out <- factor(out, levels = c("lethal", "beneficial", "neutral"))

    out
  }, de = ctr.classif, count_th = fetchTHvalue(param, th  = "overall_thr"), simplify = FALSE)
  # Note that above sapply is set to simplify = FALSE. Coercing to data.frame preserves factor levels.
  overall.classif <- as.data.frame(overall.classif)

  result <- data.frame(shrna = de$shRNA, gene = de$gene, ctr.classif, overall.classif)
  result
}

#' Trim results of NAs.
#'
#' Will remove shRNas which have missing values on at least one contrast.
#'
#' @param de Data.frame. Object from \code{scrapeDEresults}.
#'
#' @return A possible subset of \code{de} where shRNA (rows) for all samples with missing expression have been
#' removed.
trimResults <- function(de) {
  # Account for the fact that columns shRNA and gene (ergo <= 2) should never be missing.
  de[rowSums(is.na(de[, 3:ncol(de)])) <= 2, ]
}

#' Scrape DESeq2 results and format in reusable format
#'
#' A handy function which will extract information from `differentialExpression` result into a more user
#' friendly form.
#'
#' @param dds Data.frame object from \code{differentialExpression}.
#' @param shrna Character. Character of shRNA names.
#' @param gene Character. Character of gene names.
#'
#' @return A data.frame of shRNA, gene names and DESeq2 results.
scrapeDEresults <- function(dds, shrna, gene) {
  # Prepare DE results.
  int.res <- data.frame(do.call(cbind, as.matrix(unlist(dds))))

  # Modify column names. Shorten log2FoldChange to l2fc and append contrast to other columns
  # so that they read e.g. l2fc_<contrast>, lfcSE_<contrast>...
  colnames(int.res) <- gsub("log2FoldChange", replacement = "l2fc", x = colnames(int.res))
  colnames(int.res) <- sapply(strsplit(colnames(int.res), ".", fixed = TRUE), "[", 1)  # split by "." and fetch first element
  colnames(int.res) <- paste(colnames(int.res), rep(names(dds), each = 6), sep = "_")  # append contrast

  # Subset only relevant column names log2fc, lfcSE and padj.
  int.res <- int.res[, grepl("(l2fc)|(lfcSE)|(padj)", colnames(int.res))]

  # Prepare and output final result.
  out <- data.frame(shRNA = shrna, gene = gene, int.res)
  out
}

#' Perform differential expression based on specified design
#'
#' Based on provided contrasts, provide differential expression.
#'
#' @param contrasts A vector of contrasts (length of 2) for group 1 and group 2.
#' @param sample_key Data.frame. A data.frame as returned by \code{importSampleKey}.
#' @param data Count data from \code{importCountData}.
#' @param design Design used to perform differential expression.
#'
#' @return A list of DESeq results as returned by \code{DESeq2::results}.
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results
differentialExpression <- function(contrasts, sample_key, data, design = ~ 0 + group + replicate) {

  out <- apply(contrasts, MARGIN = 1, FUN = function(x, sample_key, data, design) {

    # Based on contrasts, subset data for samples who belong to the contrast (group).
    find.cols <- sprintf("(%s)|(%s)", x[1], x[2])
    data <- data[, grepl(find.cols, colnames(data)), drop = FALSE]

    # Subset sample_key to match (samples) subsetted from data object.
    sample_key <- sample_key[sample_key$group %in% x, ]

    # Perform differential expression analysis
    # TODO: capture errors gracefully
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(data), colData = sample_key[, c("replicate", "group")], design = design)
    dds <- DESeq(object = dds, fitType = "local", minReplicatesForReplace = Inf, betaPrior = FALSE)
    results <- results(dds, contrast = c("group", x[c("group_1", "group_2")]))
  }, sample_key = sample_key, data = data, design = design)

  names(out) <- contrasts$contrast
  out

}

#' Import expression data based on sample key.
#'
#' Imports expression data based on sample name specified in \code{sample_key}. It takes sample
#' name from \code{sample_key} and tries to match the name in \code{sample_list}. For final
#' output, results are merged based on shRNA.
#'
#' @param sample_key A data.frame with columns sample_name, <groups>* and replicate. Number of
#' columns for groups may vary.
#' @param sample_list Character. Vector of absolute or relative paths to sample count matrices. If
#' NULL, subfolder ./data is search for all files ending in .txt.
#'
#' @return A data.frame of shRNA counts for selected samples.
importExpressionsData <- function(sample_key, sample_list = NULL) {
  if (is.null(sample_list)) {
    stop("Please provide a list of sample files or point to a folder where they reside.")
  }

  # The next if checks if provided sample_list is a folder. If it is, it assumes it contains files
  # ending in .txt that need to be scraped.
  if (length(sample_list) == 1 & dir.exists(sample_list) & !is.null(sample_list)) {
    sample_list <- list.files(sample_list, pattern = ".txt$", full.names = TRUE)
  }

  # cm should output a list of count matrices with modified columns that now include group. Final
  # column name is then <sample>_<group>.
  cm <- apply(sample_key, MARGIN = 1, FUN = function(x, samples) {
    # Fetch sample count matrix.
    # This part is needed so that searching for sample-1 doesn't also return sample-10, sample-11...
    pick.samples <- tools::file_path_sans_ext(basename(samples))

    # Having this trimmed_trimmed_mapped_species "snake" behind the sample name may look cumbersome,
    # but it aids in distinguishing samples. This avoids matching sample-1 with sample-10, because
    # it is now sample-1_snake and sample-10_snake.
    sample.regex <- sprintf("%s_trimmed_trimmed_count_matrix", x["sample"])

    sample.to.fetch <- samples[grepl(sample.regex, pick.samples)]

    if (length(sample.to.fetch) != 1) {
      stop("Number of fetched samples is not 1. Check that regex captures correct sample names.")
    }

    new.colname <- paste(x["sample"], x["group"], sep = "_")
    xy <- read.table(sample.to.fetch, header = TRUE, sep = "\t")
    colnames(xy) <- c("shRNA", new.colname)
    xy
  }, samples = sample_list)

  # Merge all count matrices based on shRNA. It is important to merge based on shRNA so that we can
  # relax the assumption that count matrices are all identically structured.
  out <- Reduce(f = function(x, y) merge(x, y, by = "shRNA", sort = FALSE), x = cm)
  out
}

#' Process string to extract shRNA or gene
#'
#' shRNA-gene character is split on "-" and gene or shRNA are returned based on parameter entity.
#'
#' @param x Character. A vector of shRNA-gene names.
#' @param entity Character. Return either shRNA ("shrna") or gene ("gene").
#'
#' @return Depending on which \code{entity} is chosen, shRNA part or gene part are returned.
fetchGeneOrHairpin <- function(x, entity = c("shrna", "gene")) {
  # Split string on - and fetch first element. It is assumed this is shRNA.
  # "shRNA" string holds shRNA-gene combination. Based on which entity user chooses, subset first
  # or second part of the string.
  if (entity == "shrna") get.entity <- 1
  if (entity == "gene") get.entity <- 2

  sp <- strsplit(as.character(x), "-")
  out <- sapply(sp, "[", get.entity)
  out
}

#' Import sample key
#'
#' Take in parameter file (.xlsx) and from tab \code{sample_key} extract information on study design.
#'
#' @param x Character. Absolute or relative path to an .xlsx file. See \code{doDE} description
#' for more information.
#'
#' @return A data.frame with columns sample, replicate and group, used in differential expression.
#' @importFrom readxl read_xlsx
importSampleKey <- function(x) {
  # This function assumes sample, <variables> and replicate columns structure. It creates a
  # new variable called group which concatenates all experiment variables (sample and
  # replicate are assumed not to be part of it).

  # The final data.frame should have the following columns:
  # sample         # this is sample name
  # replicate      # replicate nubmer should always be present
  # group          # this is a variable constructed by concatenating all group variables with _

  xy <- as.data.frame(read_xlsx(x, sheet = "sample_key"))
  # Find colums which are NOT sample and replicate
  exl <- which(colnames(xy) %in% c("sample", "replicate"))
  # Concatenate values from all columns except sample and replicate
  xy$group <- do.call(paste, c(xy[, -exl], sep = "_"))

  # Make sure all variables are factors. This is important when doing DE.
  xy$replicate <- factor(xy$replicate)
  xy$group <- factor(xy$group)

  xy[, c("sample", "replicate", "group")]
}

#' Import contrast list
#'
#' User input should be contrasts (group_1, group_2) which will be used to compare
#' in the DE run.
#'
#' @param x Character. Absolute or relative path to an .xlsx file. See \code{doDE} description
#' for more information.
#'
#' @return A data.frame of contrasts for group_1 and group_2 and a concatenation of
#' group_1.group_2 as column \code{contrast}.
#' @importFrom readxl read_xlsx
importContrastList <- function(x) {
  xy <- as.data.frame(read_xlsx(x, sheet = "contrasts"))
  # Prepare contrasts based on groups.
  xy$contrast <- do.call(paste, c(xy, sep = "."))
  xy
}

#' Import overall contrasts
#'
#' From parameter file's \code{overall_contrasts} fetch all defined overall contrasts.
#'
#' @param x Character. Absolute or relative path to an .xlsx file. See \code{doDE} description
#' for more information.
#'
#' @return A data.frame of overall contrasts. Each column represents a group of contrasts.
#' @importFrom readxl read_xlsx
importOverallContrasts <- function(x) {
  xy <- as.data.frame(read_xlsx(x, sheet = "overall_contrasts"))
  xy
}

#' Import parameters used for classification
#' @param x Character. Absolute or relative path to an .xlsx file. See \code{doDE} description
#' for more information.
#' @param trim Logical. If \code{TRUE} (default), description column will be omitted.
#'
#' @return A data.frame with two columns, one being name of the parameters and second being its value.
#' @importFrom readxl read_xlsx
importClassificationParams <- function(x, trim = TRUE) {
  xy <- as.data.frame(read_xlsx(x, sheet = "classification_parameters"))

  if (trim) {
    return(xy[-which(colnames(xy) %in% c("description"))])
  }
  xy
}

#' Return specified threshold value.
#'
#' A helper function which makes fetching of parameters easy, pretty and extensible.
#'
#' @param x Data.frame. A threshold data.frame as returned by importClassificationParams.
#' @param th Character. Character of length 1. Possible values are \code{l_thr}, \code{p_thr},
#' \code{overall_thr}.
#'
#' @return A value corresponding to \code{th}.
fetchTHvalue <- function(x, th = NULL) {
  if (is.null(th)) stop("You need to specify which threshold value to return.")

  # Attempt to fecth value.
  out <- x$value[x$threshold == th]

  # Fail the query of lookup comes up empty.
  if (length(out) == 0) stop("You have specified a parameters that does not exist in x.")

  out
}