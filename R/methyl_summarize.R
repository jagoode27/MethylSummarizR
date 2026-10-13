#' Summarize Methylation Data
#'
#' The `methyl_summarize()` function provides a comprehensive summary of methylation data on the beta value scale, the M-value scale, or both. The function computes various statistics, including mean, standard deviation, confidence intervals, and more for each probe across samples.
#'
#' @param methyl_data A numeric matrix containing methylation data, with probe IDs as row names. Each row represents a probe, and each column represents a sample.
#' @param methyl_in A character string indicating the scale of `methyl_data`. Must be either `"betas"` or `"mvals"`.
#' @param methyl_out A character string indicating the scale of the summaries to return. Must be `"betas"`, `"mvals"`, or `"both"` (the default).
#' @param cores The number of cores to use. Probes are split into `cores` chunks that are summarized in parallel using the base parallel package: forked processes on macOS and Linux, and a socket cluster on Windows. The default of 1 runs without parallel processing.
#'
#' @return A named list with one data frame for each scale requested by `methyl_out`:
#' \describe{
#'   \item{`betas`}{A data frame summarizing the methylation data on the beta value scale.}
#'   \item{`mvals`}{A data frame summarizing the methylation data on the M-value scale.}
#' }
#'
#' Each data frame includes the following columns:
#' \itemize{
#'   \item `probe`: The probe ID.
#'   \item `mean`: Mean methylation value.
#'   \item `sd`: Standard deviation.
#'   \item `ci95_lb`: Lower bound of the 95% confidence interval.
#'   \item `ci95_ub`: Upper bound of the 95% confidence interval.
#'   \item `coef_var`: Coefficient of variation.
#'   \item `iqr`: Interquartile range.
#'   \item `min`: Minimum value.
#'   \item `p05`: 5th percentile.
#'   \item `p25`: 25th percentile (first quartile).
#'   \item `p50`: Median (50th percentile).
#'   \item `p75`: 75th percentile (third quartile).
#'   \item `p95`: 95th percentile.
#'   \item `max`: Maximum value.
#'   \item `miss_n`: Number of missing values.
#'   \item `nonmiss_n`: Number of non-missing values.
#'   \item `miss_p`: Proportion of missing values.
#'   \item `nonmiss_p`: Proportion of non-missing values.
#' }
#'
#' @details The function checks whether `methyl_data` is a numeric matrix with row names and whether `methyl_in`, `methyl_out`, and `cores` are valid. Data are converted between scales only when `methyl_out` asks for a scale other than `methyl_in`. Before conversion to M-values, beta values of 0 or 1 are moved 1e-6 inside those bounds so that M-values stay finite, and beta values outside 0 to 1 give a warning. Missing values are excluded from each statistic. The 95% confidence interval uses a normal approximation (mean plus or minus 1.96 standard errors) and is `NA` for probes with fewer than two non-missing values. If the number of columns exceeds the number of rows in `methyl_data`, a warning is issued to ensure that the probes are organized as rows and samples as columns.
#'
#' @seealso [methyl_reliability()] for agreement between duplicate samples.
#'
#' @examples
#' # Summarize beta values on both scales
#' data(beta_matrix, package = "MethylSummarizR")
#' results <- methyl_summarize(beta_matrix, methyl_in = "betas")
#' head(results$betas)
#'
#' # Summarize M-values on the beta value scale only
#' data(mval_matrix, package = "MethylSummarizR")
#' results <- methyl_summarize(mval_matrix, methyl_in = "mvals", methyl_out = "betas")
#' head(results$betas)
#'
#' @export
methyl_summarize <- function(methyl_data, methyl_in, methyl_out = "both", cores = 1) {
  methyl_in <- check_choice(methyl_in, c("betas", "mvals"), "methyl_in")
  methyl_out <- check_choice(methyl_out, c("betas", "mvals", "both"), "methyl_out")
  cores <- check_cores(cores)
  check_methyl_data(methyl_data)

  # Summarize each requested scale
  lapply(out_formats(methyl_out), function(format) {
    run_chunked(summarize_chunk, list(to_scale(methyl_data, methyl_in, format)), cores)
  })
}

# Summarize each row of a methylation matrix (base and stats only, see run_chunked)
summarize_chunk <- function(methyl_data) {
  summary <- apply(methyl_data, 1, function(values) {
    n_total <- length(values)
    n_miss <- sum(is.na(values))
    n_nonmiss <- n_total - n_miss

    if (n_nonmiss == 0) {
      return(c(rep(NA_real_, 13), n_miss, n_nonmiss, n_miss / n_total, n_nonmiss / n_total))
    }

    mean_val <- mean(values, na.rm = TRUE)
    sd_val <- stats::sd(values, na.rm = TRUE)
    margin <- if (n_nonmiss > 1) stats::qnorm(0.975) * sd_val / sqrt(n_nonmiss) else NA_real_
    coef_var <- if (mean_val == 0) NA_real_ else sd_val / abs(mean_val)
    iqr_val <- stats::IQR(values, na.rm = TRUE)
    quantiles <- stats::quantile(values, c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00),
                                 na.rm = TRUE, names = FALSE)

    c(mean_val, sd_val, mean_val - margin, mean_val + margin, coef_var, iqr_val,
      quantiles, n_miss, n_nonmiss, n_miss / n_total, n_nonmiss / n_total)
  })

  summary_df <- data.frame(probe = rownames(methyl_data), t(summary), row.names = NULL)
  colnames(summary_df) <- c("probe", "mean", "sd", "ci95_lb", "ci95_ub", "coef_var", "iqr",
                            "min", "p05", "p25", "p50", "p75", "p95", "max",
                            "miss_n", "nonmiss_n", "miss_p", "nonmiss_p")
  summary_df
}
