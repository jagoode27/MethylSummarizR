# Validate a numeric methylation matrix with probes in rows
check_methyl_data <- function(methyl_data) {
  if (!is.matrix(methyl_data)) {
    stop("methyl_data must be a matrix.", call. = FALSE)
  }
  if (!is.numeric(methyl_data)) {
    stop("methyl_data must contain only numeric values.", call. = FALSE)
  }
  if (nrow(methyl_data) == 0 || ncol(methyl_data) == 0) {
    stop("methyl_data cannot be empty.", call. = FALSE)
  }
  if (is.null(rownames(methyl_data))) {
    stop("methyl_data must have probe IDs as row names.", call. = FALSE)
  }
  if (ncol(methyl_data) > nrow(methyl_data)) {
    warning("methyl_data has more columns than rows. Please ensure that sample IDs are in columns and probe IDs are in rows to avoid incorrect results.", call. = FALSE)
  }
  invisible(methyl_data)
}

# Validate a single string against a set of allowed values
check_choice <- function(x, choices, arg) {
  if (!is.character(x) || length(x) != 1 || !(x %in% choices)) {
    quoted <- paste0("'", choices, "'")
    allowed <- if (length(choices) == 2) {
      paste("either", quoted[1], "or", quoted[2])
    } else {
      paste0("one of ", paste(quoted[-length(quoted)], collapse = ", "), ", or ", quoted[length(quoted)])
    }
    stop(arg, " must be ", allowed, ".", call. = FALSE)
  }
  x
}

# Validate the number of cores and return it as an integer
check_cores <- function(cores) {
  if (!is.numeric(cores) || length(cores) != 1 || !is.finite(cores) ||
      cores < 1 || cores != round(cores)) {
    stop("cores must be a single positive whole number.", call. = FALSE)
  }
  as.integer(cores)
}

# Validate a column name in the phenotype data
check_pheno_column <- function(column, pheno, arg) {
  if (!is.character(column) || length(column) != 1 || !(column %in% colnames(pheno))) {
    stop(arg, " must be the name of a column in pheno.", call. = FALSE)
  }
  column
}

# Expand methyl_out into a named vector of output formats
out_formats <- function(methyl_out) {
  formats <- if (methyl_out == "both") c("betas", "mvals") else methyl_out
  stats::setNames(formats, formats)
}

# Convert a methylation matrix to the requested scale
to_scale <- function(methyl_data, from, to) {
  if (from == to) {
    return(methyl_data)
  }
  if (to == "mvals") convert_beta_to_m(methyl_data) else convert_m_to_beta(methyl_data)
}

# Apply fun to row chunks of matrices in parallel and row-bind the results
# (fun must use only base and `pkg::` calls so socket workers can run it)
run_chunked <- function(fun, mats, cores = 1L, fork = .Platform$OS.type != "windows") {
  n_rows <- nrow(mats[[1]])
  n_chunks <- min(cores, n_rows)
  if (n_chunks <= 1) {
    return(do.call(fun, mats))
  }

  chunk_rows <- split(seq_len(n_rows), ceiling(seq_len(n_rows) * n_chunks / n_rows))
  slice_rows <- function(rows) lapply(mats, function(m) m[rows, , drop = FALSE])

  if (fork) {
    results <- suppressWarnings(parallel::mclapply(
      chunk_rows,
      function(rows) do.call(fun, slice_rows(rows)),
      mc.cores = n_chunks
    ))
    failed <- vapply(results, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
    if (any(failed)) {
      first_error <- results[[which(failed)[1]]]
      if (inherits(first_error, "try-error")) {
        stop("A parallel worker failed: ", conditionMessage(attr(first_error, "condition")), call. = FALSE)
      }
      stop("A parallel worker failed to return a result. Try using fewer cores.", call. = FALSE)
    }
  } else {
    worker <- fun
    environment(worker) <- baseenv()
    cl <- parallel::makePSOCKcluster(n_chunks)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # do.call is passed directly so no local environment is sent to the workers
    results <- parallel::parLapply(cl, lapply(chunk_rows, slice_rows), do.call, what = worker)
  }

  combined <- do.call(rbind, unname(results))
  rownames(combined) <- NULL
  combined
}
