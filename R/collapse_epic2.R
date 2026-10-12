#' Collapse Replicate EPIC v2 Probes
#'
#' The Illumina MethylationEPIC v2 array includes replicate probes for some CpG sites, named with a shared CpG ID followed by an underscore and a suffix (for example, `cg00000029_TC21`). The `collapse_epic2()` function averages replicate probes within each sample so that each CpG ID appears once.
#'
#' @param methyl_data A numeric matrix or data frame of methylation data. Each row represents a probe, and each column represents a sample. Probe IDs come from the row names unless `probe_names` is given.
#' @param probe_names For data frame input, the name of the column that holds probe IDs. Leave as `NULL` (the default) to use row names.
#'
#' @return An object of the same type as `methyl_data` with one row per CpG ID, in the order each ID first appears. For a matrix, or a data frame without `probe_names`, the CpG IDs are the row names. For a data frame with `probe_names`, the CpG IDs are in the first column, named `probe_names`.
#'
#' @details Probe IDs are shortened to the text before the first underscore, and probes that share the shortened ID are averaged. Missing values are ignored, and a CpG ID with no observed values in a sample is `NA`. Probes without a suffix are returned unchanged.
#'
#' @examples
#' epic2 <- matrix(
#'   c(0.10, 0.20, 0.30,
#'     0.12, NA,   0.34,
#'     0.80, 0.70, 0.60),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(c("cg00000029_TC21", "cg00000029_TC22", "cg00000109_BC21"),
#'                   c("s01", "s02", "s03"))
#' )
#' collapse_epic2(epic2)
#'
#' # Probe IDs in a column of a data frame
#' epic2_df <- data.frame(probe = rownames(epic2), epic2, row.names = NULL)
#' collapse_epic2(epic2_df, probe_names = "probe")
#'
#' @export
collapse_epic2 <- function(methyl_data, probe_names = NULL) {
  matrix_input <- is.matrix(methyl_data)
  if (!matrix_input && !is.data.frame(methyl_data)) {
    stop("methyl_data must be a matrix or a data frame.")
  }

  # Get probe IDs from the row names or the probe_names column
  if (is.null(probe_names)) {
    if (is.null(rownames(methyl_data)) ||
        (!matrix_input && .row_names_info(methyl_data) < 0)) {
      stop("Row names are missing in methyl_data.")
    }
    probe_ids <- rownames(methyl_data)
  } else {
    if (matrix_input) {
      stop("probe_names can only be used when methyl_data is a data frame.")
    }
    if (!is.character(probe_names) || length(probe_names) != 1 ||
        !(probe_names %in% colnames(methyl_data))) {
      stop("The specified probe_names column does not exist.")
    }
    probe_ids <- as.character(methyl_data[[probe_names]])
    methyl_data <- methyl_data[, colnames(methyl_data) != probe_names, drop = FALSE]
  }

  # Ensure the remaining columns contain only numeric values
  is_numeric <- if (matrix_input) {
    is.numeric(methyl_data)
  } else {
    all(vapply(methyl_data, is.numeric, logical(1)))
  }
  if (!is_numeric) {
    stop("All columns must contain numeric values (excluding the grouping column if specified).")
  }
  methyl_data <- as.matrix(methyl_data)

  # Average probes that share the text before the first underscore
  cpg_ids <- sub("_.*$", "", probe_ids)
  sums <- rowsum(methyl_data, cpg_ids, reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(1 * !is.na(methyl_data), cpg_ids, reorder = FALSE)
  collapsed <- sums / counts
  collapsed[counts == 0] <- NA

  # Return the same type as the input
  if (matrix_input) {
    return(collapsed)
  }
  result <- as.data.frame(collapsed)
  if (!is.null(probe_names)) {
    result <- cbind(stats::setNames(data.frame(rownames(collapsed)), probe_names), result)
    rownames(result) <- NULL
  }
  result
}
