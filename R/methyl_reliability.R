#' Estimate Reliability of Methylation Data from Duplicate Samples
#'
#' The `methyl_reliability()` function estimates probe-level agreement between duplicate samples from the same person. For each probe it computes several intraclass correlation coefficients (ICCs), the Pearson correlation, and summaries of the absolute differences between duplicates, on the beta value scale, the M-value scale, or both.
#'
#' @param methyl_data A numeric matrix containing methylation data, with probe IDs as row names and sample IDs as column names. Each row represents a probe, and each column represents a sample.
#' @param pheno A data frame or matrix of phenotype data with one row per sample.
#' @param per_id The name of the column in `pheno` that identifies people. Samples that share a value are treated as duplicates.
#' @param meth_id The name of the column in `pheno` that identifies samples. Values must match the column names of `methyl_data`.
#' @param methyl_in A character string indicating the scale of `methyl_data`. Must be either `"betas"` or `"mvals"`.
#' @param methyl_out A character string indicating the scale of the estimates to return. Must be `"betas"`, `"mvals"`, or `"both"` (the default).
#' @param cores The number of cores to use. Probes are split into `cores` chunks that are processed in parallel using the base parallel package: forked processes on macOS and Linux, and a socket cluster on Windows. The default of 1 runs without parallel processing.
#'
#' @return A named list with one data frame for each scale requested by `methyl_out`:
#' \describe{
#'   \item{`betas`}{A data frame of reliability estimates on the beta value scale.}
#'   \item{`mvals`}{A data frame of reliability estimates on the M-value scale.}
#' }
#'
#' Each data frame includes the following columns:
#' \itemize{
#'   \item `probe`: The probe ID.
#'   \item `icc_1`: One-way random effects ICC for a single measurement, ICC(1).
#'   \item `icc_c1`: Two-way ICC for consistency of a single measurement, ICC(C,1).
#'   \item `icc_a1`: Two-way ICC for absolute agreement of a single measurement, ICC(A,1).
#'   \item `icc_k`: One-way random effects ICC for the mean of two measurements, ICC(k).
#'   \item `icc_ck`: Two-way ICC for consistency of the mean of two measurements, ICC(C,k).
#'   \item `icc_ak`: Two-way ICC for absolute agreement of the mean of two measurements, ICC(A,k).
#'   \item `pearson`: Pearson correlation between the first and second samples.
#'   \item `diff_mn`: Mean absolute difference between duplicates.
#'   \item `diff_sd`: Standard deviation of the absolute differences.
#'   \item `diff_cv`: Coefficient of variation of the absolute differences.
#'   \item `n_dups`: Number of duplicate pairs with non-missing values for both samples.
#' }
#'
#' @details For each person with more than one sample in `pheno`, the first two samples listed are used as a duplicate pair, and any further samples are ignored. Samples in `pheno` that are not columns of `methyl_data` are dropped with a warning.
#'
#' For each probe, only pairs with both values observed are used. ICCs follow the forms described by Shrout and Fleiss (1979) and McGraw and Wong (1996), with the first and second samples of each pair treated as the two measurements. Metrics are `NA` for probes with fewer than two complete pairs (a single warning reports how many) and where they cannot be computed, such as when a probe has no variance.
#'
#' Data are converted between scales only when `methyl_out` asks for a scale other than `methyl_in`. Before conversion to M-values, beta values of 0 or 1 are moved 1e-6 inside those bounds so that M-values stay finite, and beta values outside 0 to 1 give a warning.
#'
#' @references
#' Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in assessing rater reliability. *Psychological Bulletin*, 86(2), 420-428.
#'
#' McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some intraclass correlation coefficients. *Psychological Methods*, 1(1), 30-46.
#'
#' @seealso [methyl_summarize()] for probe-level summary statistics.
#'
#' @examples
#' data(beta_matrix, package = "MethylSummarizR")
#' data(pheno, package = "MethylSummarizR")
#'
#' # Reliability of beta values on both scales
#' results <- methyl_reliability(beta_matrix, pheno,
#'                               per_id = "person_id", meth_id = "sample_id",
#'                               methyl_in = "betas")
#' head(results$betas)
#'
#' # Reliability on the M-value scale only
#' results <- methyl_reliability(beta_matrix, pheno,
#'                               per_id = "person_id", meth_id = "sample_id",
#'                               methyl_in = "betas", methyl_out = "mvals")
#' head(results$mvals)
#'
#' @export
methyl_reliability <- function(methyl_data, pheno, per_id, meth_id, methyl_in,
                               methyl_out = "both", cores = 1) {
  methyl_in <- check_choice(methyl_in, c("betas", "mvals"), "methyl_in")
  methyl_out <- check_choice(methyl_out, c("betas", "mvals", "both"), "methyl_out")
  cores <- check_cores(cores)
  check_methyl_data(methyl_data)
  if (is.null(colnames(methyl_data))) {
    stop("methyl_data must have sample IDs as column names.")
  }

  # Ensure pheno is either a matrix or data frame with the ID columns
  if (!is.matrix(pheno) && !is.data.frame(pheno)) {
    stop("pheno must be either a matrix or a data frame.")
  }
  pheno <- as.data.frame(pheno)
  check_pheno_column(per_id, pheno, "per_id")
  check_pheno_column(meth_id, pheno, "meth_id")

  # Drop samples that are not in the methylation data
  in_methyl <- as.character(pheno[[meth_id]]) %in% colnames(methyl_data)
  if (!all(in_methyl)) {
    warning(sum(!in_methyl), " sample(s) in pheno were not found in the columns of methyl_data and were dropped.")
  }
  pheno <- pheno[in_methyl & !is.na(pheno[[per_id]]), , drop = FALSE]

  # Keep the first two samples of each person with duplicates
  pheno <- pheno[order(pheno[[per_id]]), , drop = FALSE]
  occurrence <- stats::ave(seq_len(nrow(pheno)), pheno[[per_id]], FUN = seq_along)
  n_samples <- stats::ave(seq_len(nrow(pheno)), pheno[[per_id]], FUN = length)
  first_ids <- as.character(pheno[[meth_id]][occurrence == 1 & n_samples > 1])
  second_ids <- as.character(pheno[[meth_id]][occurrence == 2])

  if (length(first_ids) == 0) {
    stop(paste0("No values of per_id (", per_id, ") are repeated in the phenotype data. Please ensure that replicate/duplicate cases are included in phenotype and methylation data."))
  }

  # Estimate reliability on each requested scale
  n_pairs <- length(first_ids)
  paired <- methyl_data[, c(first_ids, second_ids), drop = FALSE]
  results <- lapply(out_formats(methyl_out), function(format) {
    scaled <- to_scale(paired, methyl_in, format)
    run_chunked(
      reliability_chunk,
      list(scaled[, seq_len(n_pairs), drop = FALSE],
           scaled[, n_pairs + seq_len(n_pairs), drop = FALSE]),
      cores
    )
  })

  n_short <- sum(results[[1]]$n_dups < 2)
  if (n_short > 0) {
    warning(n_short, " probe(s) had fewer than 2 complete duplicate pairs; their metrics are NA.")
  }

  results
}

# Agreement metrics between paired columns by row (base only, see run_chunked)
reliability_chunk <- function(methyl_1, methyl_2) {
  complete <- !is.na(methyl_1) & !is.na(methyl_2)
  methyl_1[!complete] <- NA
  methyl_2[!complete] <- NA
  n <- rowSums(complete)
  k <- 2

  # Mean squares for rows (pairs), columns (duplicates), within rows, and error
  mean_1 <- rowMeans(methyl_1, na.rm = TRUE)
  mean_2 <- rowMeans(methyl_2, na.rm = TRUE)
  grand_mean <- (mean_1 + mean_2) / 2
  pair_mean <- (methyl_1 + methyl_2) / 2
  msr <- k * rowSums((pair_mean - grand_mean)^2, na.rm = TRUE) / (n - 1)
  msw <- rowSums((methyl_1 - methyl_2)^2, na.rm = TRUE) / (2 * n)
  msc <- n * (mean_1 - mean_2)^2 / 2
  sst <- rowSums((methyl_1 - grand_mean)^2, na.rm = TRUE) +
    rowSums((methyl_2 - grand_mean)^2, na.rm = TRUE)
  mse <- (sst - msr * (n - 1) - msc * (k - 1)) / ((n - 1) * (k - 1))

  # Intraclass correlations
  icc_1 <- (msr - msw) / (msr + (k - 1) * msw)
  icc_c1 <- (msr - mse) / (msr + (k - 1) * mse)
  icc_a1 <- (msr - mse) / (msr + (k - 1) * mse + (k / n) * (msc - mse))
  icc_k <- (msr - msw) / msr
  icc_ck <- (msr - mse) / msr
  icc_ak <- (msr - mse) / (msr + (msc - mse) / n)

  # Pearson correlation
  dev_1 <- methyl_1 - mean_1
  dev_2 <- methyl_2 - mean_2
  pearson <- rowSums(dev_1 * dev_2, na.rm = TRUE) /
    sqrt(rowSums(dev_1^2, na.rm = TRUE) * rowSums(dev_2^2, na.rm = TRUE))

  # Absolute differences between duplicates
  abs_diff <- abs(methyl_1 - methyl_2)
  diff_mn <- rowMeans(abs_diff, na.rm = TRUE)
  diff_sd <- sqrt(rowSums((abs_diff - diff_mn)^2, na.rm = TRUE) / (n - 1))
  diff_cv <- diff_sd / diff_mn

  metrics <- cbind(icc_1, icc_c1, icc_a1, icc_k, icc_ck, icc_ak,
                   pearson, diff_mn, diff_sd, diff_cv)
  metrics[n < 2, ] <- NA
  metrics[is.nan(metrics)] <- NA

  data.frame(probe = rownames(methyl_1), metrics, n_dups = as.integer(n), row.names = NULL)
}
