#' Convert M-Values to Beta Values
#'
#' Converts methylation M-values to beta values. Beta values represent methylation
#' levels as proportions (0-1), while M-values are log2 ratios that provide better
#' statistical properties for analysis.
#'
#' @param methyl_data A matrix of methylation values with CpG sites as row names and samples as column names. All values must be numeric.
#'
#' @return A matrix of beta values corresponding to the input M-values.
#'
#' @details
#' **Numerical Stability:** Beta values are computed as `1 / (1 + 2^(-M))`, so very large or very small M-values convert to 1 or 0 instead of `NaN`.
#'
#' Missing values stay missing.
#'
#' @seealso [convert_beta_to_m()] for the inverse conversion
#'
#' @examples
#' # Load Methylation M-Values Matrix
#' data(mval_matrix, package = "MethylSummarizR")
#'
#' # Convert M-values to Beta Values
#' beta_values <- convert_m_to_beta(mval_matrix)
#' beta_values[1:5, 1:5]
#'
#' @export
convert_m_to_beta <- function(methyl_data) {
  check_conversion_input(methyl_data)
  1 / (1 + 2^(-methyl_data))
}

#' Convert Beta Values to M-Values
#'
#' Converts methylation beta values to M-values. Beta values represent methylation
#' levels as proportions (0-1), while M-values are log2 ratios that provide better
#' statistical properties for analysis.
#'
#' @param methyl_data A matrix of methylation values with CpG sites as row names and samples as column names. All values must be numeric.
#'
#' @return A matrix of M-values corresponding to the input beta values.
#'
#' @details
#' **Numerical Stability:** The function automatically handles extreme beta values (0 or 1) that would otherwise
#' result in infinite M-values by clamping them to a small epsilon away from the boundaries.
#'
#' **Value Ranges:** Beta values should be between 0 and 1. Values outside this range will generate a warning
#' but the conversion will proceed.
#'
#' Missing values stay missing.
#'
#' @seealso [convert_m_to_beta()] for the inverse conversion
#'
#' @examples
#' # Load Methylation Beta Matrix
#' data(beta_matrix, package = "MethylSummarizR")
#'
#' # Convert Beta Values to M-values
#' m_values <- convert_beta_to_m(beta_matrix)
#' m_values[1:5, 1:5]
#'
#' @export
convert_beta_to_m <- function(methyl_data) {
  check_conversion_input(methyl_data)

  # Check for values outside expected range and warn
  if (any(!is.na(methyl_data) & (methyl_data < 0 | methyl_data > 1))) {
    warning("Beta values outside [0,1] range detected. Results may be unreliable.")
  }

  # Clamp beta values to prevent infinite M-values, then convert
  beta_clamped <- pmax(pmin(methyl_data, 1 - 1e-6), 1e-6)
  m_values <- log2(beta_clamped / (1 - beta_clamped))
  m_values[is.na(methyl_data)] <- NA
  m_values
}

# Validate input to the conversion functions
check_conversion_input <- function(methyl_data) {
  if (!is.matrix(methyl_data)) {
    stop("Input must be a matrix.", call. = FALSE)
  }
  if (!is.numeric(methyl_data)) {
    stop("All values in the methylation matrix must be numeric.", call. = FALSE)
  }
  if (nrow(methyl_data) == 0 || ncol(methyl_data) == 0) {
    stop("Input matrix cannot be empty.", call. = FALSE)
  }
  if (is.null(rownames(methyl_data))) {
    warning("Matrix lacks row names. Consider adding CpG site identifiers as row names.", call. = FALSE)
  }
  invisible(methyl_data)
}
