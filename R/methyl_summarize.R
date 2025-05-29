#' Summarize Methylation Data
#'
#' The `methyl_summarize()` function provides a comprehensive summary of methylation data, either in beta values or M-values format. The function computes various statistics, including mean, standard deviation, confidence intervals, and more for each probe across samples.
#'
#' @param methyl_data A numeric matrix containing methylation data. Each row represents a probe, and each column represents a sample.
#' @param methyl_format A character string indicating the format of the methylation data. Must be either `"betas"` or `"mvals"`.
#'
#' @return A list containing two data frames:
#' \describe{
#'   \item{`betas`}{A data frame summarizing the methylation data in beta values format.}
#'   \item{`mvals`}{A data frame summarizing the methylation data in M-values format.}
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
#' @details The function checks whether the input `methyl_data` is a numeric matrix and whether `methyl_format` is valid. It then summarizes the data, providing extensive statistics for each probe. If the number of columns exceeds the number of rows in `methyl_data`, a warning is issued to ensure that the probes are organized as rows and samples as columns.
#'
#' @examples
#' \dontrun{
#' # example with betas as input
#' data(beta_matrix, package = "MethylSummarizR")
#' results <- methyl_summarize(methyl_data = beta_matrix, methyl_format = "betas")
#' # example with M-values as input
#' data(mval_matrix, package = "MethylSummarizR")
#' results <- methyl_summarize(methyl_data = mval_matrix, methyl_format = "mvals")
#' }
#'
#' @export
methyl_summarize <- function(methyl_data, methyl_format) {
  # Check if methyl_format is valid
  if (!(methyl_format %in% c("betas", "mvals"))) {
    stop("methyl_format must be either 'betas' or 'mvals'.")
  }

  # Ensure methyl_data is a matrix
  if (!is.matrix(methyl_data)) {
    stop("methyl_data must be a matrix.")
  }

  # Ensure methyl_data contains only numeric values
  if (!is.numeric(methyl_data)) {
    stop("methyl_data must contain only numeric values.")
  }

  # Check if matrix orientation is likely to be incorrect
  if (ncol(methyl_data) > nrow(methyl_data)) {
    warning("methyl_data has more columns than rows. Please ensure that sample IDs are in columns and probe IDs are in rows to avoid incorrect results.")
  }

  # Conversion functions
  betas_to_mvals <- function(betas) {
    log2(betas / (1 - betas))
  }

  mvals_to_betas <- function(mvals) {
    2^mvals / (1 + 2^mvals)
  }

  # Define function to calculate confidence interval
  calc_ci95 <- function(data) {
    n <- sum(!is.na(data))  # Calculate n only for non-NA values
    if (n <= 1) return(c(NA, NA))  # CI is not defined for n <= 1
    stderr <- stats::sd(data, na.rm = TRUE) / sqrt(n)
    error_margin <- stats::qnorm(0.975) * stderr
    mean(data, na.rm = TRUE) + c(-error_margin, error_margin)
  }

  # Function to summarize matrix
  summarize_matrix <- function(matrix_data) {
    summary <- apply(matrix_data, 1, function(var) {
      n_miss <- sum(is.na(var))
      n_nonmiss <- sum(!is.na(var))
      total_n <- length(var)

      if (n_nonmiss == 0) {
        return(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, n_miss, n_nonmiss, NA, NA))
      } else {
        mean_val <- mean(var, na.rm = TRUE)
        sd_val <- stats::sd(var, na.rm = TRUE)
        ci95_vals <- calc_ci95(var)
        coef_var <- ifelse(mean_val == 0, NA, sd_val / abs(mean_val))
        iqr_val <- stats::IQR(var, na.rm = TRUE)
        quantiles <- stats::quantile(var, c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1.00), na.rm = TRUE)
        miss_p <- n_miss / total_n
        nonmiss_p <- n_nonmiss / total_n

        c(mean_val, sd_val, ci95_vals[1], ci95_vals[2], coef_var, iqr_val,
          quantiles[1], quantiles[2], quantiles[3], quantiles[4], quantiles[5], quantiles[6], quantiles[7],
          n_miss, n_nonmiss, miss_p, nonmiss_p)
      }
    }, simplify = TRUE)

    # Transpose the summary, convert to dataframe, and add 'probe' as a column
    summary_df <- data.frame(t(summary))
    colnames(summary_df) <- c("mean", "sd", "ci95_lb", "ci95_ub", "coef_var", "iqr",
                              "min", "p05", "p25", "p50", "p75", "p95", "max",
                              "miss_n", "nonmiss_n", "miss_p", "nonmiss_p")
    summary_df$probe <- rownames(matrix_data)

    # Reorder to place 'probe' at the start and remove row names
    summary_df <- summary_df[, c("probe", colnames(summary_df)[-ncol(summary_df)])]
    rownames(summary_df) <- NULL

    return(summary_df)
  }

  # Initialize result list
  results <- list()

  # Perform calculations based on format
  if (methyl_format == "betas") {
    betas_summary <- summarize_matrix(methyl_data)
    mvals_data <- betas_to_mvals(methyl_data)
    mvals_summary <- summarize_matrix(mvals_data)

  } else if (methyl_format == "mvals") {
    mvals_summary <- summarize_matrix(methyl_data)
    betas_data <- mvals_to_betas(methyl_data)
    betas_summary <- summarize_matrix(betas_data)
  }

  # Store results in a list
  results$betas <- betas_summary
  results$mvals <- mvals_summary

  return(results)
}
