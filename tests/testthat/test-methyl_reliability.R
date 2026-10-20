metric_cols <- c("icc_1", "icc_c1", "icc_a1", "icc_k", "icc_ck", "icc_ak",
                 "pearson", "diff_mn", "diff_sd", "diff_cv")

# Per-probe reference using a two-way ANOVA layout of pairs by duplicates
reference_metrics <- function(x, y) {
  complete <- !is.na(x) & !is.na(y)
  ratings <- cbind(x[complete], y[complete])
  n <- nrow(ratings)
  k <- 2
  sst <- stats::var(as.vector(ratings)) * (n * k - 1)
  msr <- stats::var(rowMeans(ratings)) * k
  msw <- mean(apply(ratings, 1, stats::var))
  msc <- stats::var(colMeans(ratings)) * n
  mse <- (sst - msr * (n - 1) - msc * (k - 1)) / ((n - 1) * (k - 1))
  abs_diff <- abs(ratings[, 2] - ratings[, 1])
  c(icc_1 = (msr - msw) / (msr + (k - 1) * msw),
    icc_c1 = (msr - mse) / (msr + (k - 1) * mse),
    icc_a1 = (msr - mse) / (msr + (k - 1) * mse + (k / n) * (msc - mse)),
    icc_k = (msr - msw) / msr,
    icc_ck = (msr - mse) / msr,
    icc_ak = (msr - mse) / (msr + (msc - mse) / n),
    pearson = stats::cor(ratings[, 1], ratings[, 2]),
    diff_mn = mean(abs_diff),
    diff_sd = stats::sd(abs_diff),
    diff_cv = stats::sd(abs_diff) / mean(abs_diff))
}

# Eight probes that each have the same values for three people measured twice
paired_data <- function(first, second) {
  values <- c(first, second)
  matrix(rep(values, each = 8), nrow = 8,
         dimnames = list(paste0("cg", 1:8), paste0("s", seq_along(values))))
}
paired_pheno <- data.frame(sample_id = paste0("s", 1:6),
                           person_id = c("a", "b", "c", "a", "b", "c"))

test_that("methyl_reliability returns the requested scales", {
  both <- methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas")
  expect_type(both, "list")
  expect_named(both, c("betas", "mvals"))
  expect_named(methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas", "mvals"),
               "mvals")
  expect_named(methyl_reliability(mval_matrix, pheno, "person_id", "sample_id", "mvals", "betas"),
               "betas")

  for (result in both) {
    expect_s3_class(result, "data.frame")
    expect_named(result, c("probe", metric_cols, "n_dups"))
    expect_identical(result$probe, rownames(beta_matrix))
    expect_type(result$n_dups, "integer")
    expect_true(all(result$n_dups <= 8))
  }
})

test_that("methyl_reliability matches per-probe reference calculations", {
  samples <- split(pheno$sample_id, pheno$person_id)
  samples <- samples[lengths(samples) > 1]
  first <- vapply(samples, `[`, character(1), 1)
  second <- vapply(samples, `[`, character(1), 2)

  result <- methyl_reliability(mval_matrix, pheno, "person_id", "sample_id", "mvals", "mvals")$mvals
  expected <- t(vapply(rownames(mval_matrix), function(probe) {
    reference_metrics(mval_matrix[probe, first], mval_matrix[probe, second])
  }, numeric(length(metric_cols))))

  expect_equal(as.matrix(result[metric_cols]), expected, ignore_attr = TRUE)
  expect_equal(result$n_dups,
               unname(rowSums(!is.na(mval_matrix[, first]) & !is.na(mval_matrix[, second]))))
})

test_that("methyl_reliability computes known values", {
  methyl_data <- paired_data(c(1, 2, 3), c(2, 3, 4))
  result <- methyl_reliability(methyl_data, paired_pheno, "person_id", "sample_id",
                               "mvals", "mvals")$mvals

  expect_equal(result$icc_1, rep(0.6, 8))
  expect_equal(result$icc_c1, rep(1, 8))
  expect_equal(result$icc_a1, rep(2 / 3, 8))
  expect_equal(result$icc_k, rep(0.75, 8))
  expect_equal(result$icc_ck, rep(1, 8))
  expect_equal(result$icc_ak, rep(0.8, 8))
  expect_equal(result$pearson, rep(1, 8))
  expect_equal(result$diff_mn, rep(1, 8))
  expect_equal(result$diff_sd, rep(0, 8))
  expect_equal(result$diff_cv, rep(0, 8))
  expect_equal(result$n_dups, rep(3L, 8))
})

test_that("methyl_reliability gives perfect agreement for identical duplicates", {
  methyl_data <- paired_data(c(1, 2, 4), c(1, 2, 4))
  result <- methyl_reliability(methyl_data, paired_pheno, "person_id", "sample_id",
                               "mvals", "mvals")$mvals

  for (col in c("icc_1", "icc_c1", "icc_a1", "icc_k", "icc_ck", "icc_ak", "pearson")) {
    expect_equal(result[[col]], rep(1, 8))
  }
  expect_equal(result$diff_mn, rep(0, 8))
  expect_equal(result$diff_sd, rep(0, 8))
  expect_true(all(is.na(result$diff_cv)))
})

test_that("methyl_reliability uses only the first two samples per person", {
  methyl_data <- cbind(paired_data(c(1, 2, 3), c(2, 3, 5)), s7 = 100)
  pheno_extra <- rbind(paired_pheno, data.frame(sample_id = "s7", person_id = "a"))

  expect_equal(
    methyl_reliability(methyl_data, pheno_extra, "person_id", "sample_id", "mvals"),
    methyl_reliability(methyl_data, paired_pheno, "person_id", "sample_id", "mvals")
  )
})

test_that("methyl_reliability drops samples missing from methyl_data", {
  pheno_extra <- rbind(pheno, data.frame(sample_id = "s99", person_id = "p01", dupe = 1L))
  expect_warning(
    result <- methyl_reliability(beta_matrix, pheno_extra, "person_id", "sample_id", "betas"),
    "1 sample(s) in pheno were not found", fixed = TRUE
  )
  expect_equal(result, methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas"))
})

test_that("methyl_reliability accepts factor IDs and matrix pheno", {
  expected <- methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas")

  pheno_factor <- pheno
  pheno_factor$sample_id <- factor(pheno_factor$sample_id, levels = rev(pheno_factor$sample_id))
  pheno_factor$person_id <- factor(pheno_factor$person_id)
  expect_equal(methyl_reliability(beta_matrix, pheno_factor, "person_id", "sample_id", "betas"),
               expected)

  expect_equal(methyl_reliability(beta_matrix, as.matrix(pheno), "person_id", "sample_id", "betas"),
               expected)
})

test_that("methyl_reliability returns NA for probes with too few pairs", {
  methyl_data <- paired_data(c(1, 2, 3), c(2, 3, 5))
  methyl_data["cg1", c("s1", "s2")] <- NA
  expect_warning(
    result <- methyl_reliability(methyl_data, paired_pheno, "person_id", "sample_id",
                                 "mvals", "mvals")$mvals,
    "1 probe(s) had fewer than 2 complete duplicate pairs", fixed = TRUE
  )
  expect_equal(result$n_dups[1], 1L)
  expect_true(all(is.na(result[1, metric_cols])))
  expect_false(anyNA(result[-1, metric_cols]))
})

test_that("methyl_reliability gives the same results from either input scale", {
  expect_equal(
    methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas"),
    methyl_reliability(mval_matrix, pheno, "person_id", "sample_id", "mvals")
  )
})

test_that("methyl_reliability validates its arguments", {
  expect_error(methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "invalid"),
               "methyl_in must be either 'betas' or 'mvals'.", fixed = TRUE)
  expect_error(methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas", "invalid"),
               "methyl_out must be one of 'betas', 'mvals', or 'both'.", fixed = TRUE)
  expect_error(methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas", cores = 1.5),
               "cores must be a single positive whole number.", fixed = TRUE)
  expect_error(methyl_reliability(as.data.frame(beta_matrix), pheno, "person_id", "sample_id", "betas"),
               "methyl_data must be a matrix.", fixed = TRUE)
  expect_error(methyl_reliability(unname(beta_matrix), pheno, "person_id", "sample_id", "betas"),
               "methyl_data must have probe IDs as row names.", fixed = TRUE)

  no_sample_ids <- beta_matrix
  colnames(no_sample_ids) <- NULL
  expect_error(methyl_reliability(no_sample_ids, pheno, "person_id", "sample_id", "betas"),
               "methyl_data must have sample IDs as column names.", fixed = TRUE)

  expect_error(methyl_reliability(beta_matrix, pheno$person_id, "person_id", "sample_id", "betas"),
               "pheno must be either a matrix or a data frame.", fixed = TRUE)
  expect_error(methyl_reliability(beta_matrix, pheno, "missing", "sample_id", "betas"),
               "per_id must be the name of a column in pheno.", fixed = TRUE)
  expect_error(methyl_reliability(beta_matrix, pheno, "person_id", "missing", "betas"),
               "meth_id must be the name of a column in pheno.", fixed = TRUE)

  no_repeats <- pheno[!duplicated(pheno$person_id), ]
  expect_error(methyl_reliability(beta_matrix, no_repeats, "person_id", "sample_id", "betas"),
               "No values of per_id (person_id) are repeated", fixed = TRUE)
})

test_that("methyl_reliability gives identical results in parallel", {
  skip_on_cran()
  expect_identical(
    methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas", cores = 2),
    methyl_reliability(beta_matrix, pheno, "person_id", "sample_id", "betas", cores = 1)
  )
})
