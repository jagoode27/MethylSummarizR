# Convert M-values to beta values
convert_m_to_beta <- function(methyl_data) {
  check_conversion_input(methyl_data)
  1 / (1 + 2^(-methyl_data))
}

# Convert beta values to M-values
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
