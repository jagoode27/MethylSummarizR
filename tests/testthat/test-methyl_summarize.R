summary_cols <- c("probe", "mean", "sd", "ci95_lb", "ci95_ub", "coef_var", "iqr",
                  "min", "p05", "p25", "p50", "p75", "p95", "max",
                  "miss_n", "nonmiss_n", "miss_p", "nonmiss_p")

test_that("methyl_summarize returns the requested scales", {
  both <- methyl_summarize(mval_matrix, "mvals")
  expect_type(both, "list")
  expect_named(both, c("betas", "mvals"))
  expect_named(methyl_summarize(mval_matrix, "mvals", "betas"), "betas")
  expect_named(methyl_summarize(beta_matrix, "betas", "mvals"), "mvals")

  for (summary in both) {
    expect_s3_class(summary, "data.frame")
    expect_named(summary, summary_cols)
    expect_identical(summary$probe, rownames(mval_matrix))
  }
})

test_that("methyl_summarize computes known values", {
  methyl_data <- matrix(
    c(1,  2,  3,  4,
      NA, 2,  NA, 4,
      NA, NA, NA, NA,
      5,  NA, NA, NA,
      -1, 1,  NA, NA),
    nrow = 5, byrow = TRUE,
    dimnames = list(paste0("cg", 1:5), paste0("s", 1:4))
  )
  result <- methyl_summarize(methyl_data, "mvals", "mvals")$mvals

  full <- result[1, ]
  expect_equal(full$mean, 2.5)
  expect_equal(full$sd, stats::sd(1:4))
  expect_equal(full$ci95_lb, 2.5 - stats::qnorm(0.975) * stats::sd(1:4) / 2)
  expect_equal(full$ci95_ub, 2.5 + stats::qnorm(0.975) * stats::sd(1:4) / 2)
  expect_equal(full$coef_var, stats::sd(1:4) / 2.5)
  expect_equal(full$iqr, 1.5)
  expect_equal(unlist(full[c("min", "p05", "p25", "p50", "p75", "p95", "max")], use.names = FALSE),
               c(1, 1.15, 1.75, 2.5, 3.25, 3.85, 4))
  expect_equal(unlist(full[c("miss_n", "nonmiss_n", "miss_p", "nonmiss_p")], use.names = FALSE),
               c(0, 4, 0, 1))

  expect_equal(result$mean[2], 3)
  expect_equal(unlist(result[2, c("miss_n", "nonmiss_n", "miss_p", "nonmiss_p")], use.names = FALSE),
               c(2, 2, 0.5, 0.5))

  all_missing <- result[3, ]
  expect_true(all(is.na(all_missing[summary_cols[2:14]])))
  expect_equal(unlist(all_missing[c("miss_n", "nonmiss_n", "miss_p", "nonmiss_p")], use.names = FALSE),
               c(4, 0, 1, 0))

  single <- result[4, ]
  expect_equal(single$mean, 5)
  expect_true(all(is.na(single[c("sd", "ci95_lb", "ci95_ub", "coef_var")])))
  expect_equal(single$iqr, 0)

  expect_equal(result$mean[5], 0)
  expect_true(is.na(result$coef_var[5]))
})

test_that("methyl_summarize gives the same results from either input scale", {
  from_betas <- methyl_summarize(beta_matrix, "betas")
  from_mvals <- methyl_summarize(mval_matrix, "mvals")
  expect_equal(from_betas, from_mvals)
})

test_that("methyl_summarize counts missing data", {
  result <- methyl_summarize(beta_matrix, "betas", "betas")$betas
  expect_equal(result$miss_n, unname(rowSums(is.na(beta_matrix))))
  expect_equal(result$miss_n + result$nonmiss_n, rep(ncol(beta_matrix), nrow(beta_matrix)))
})

test_that("methyl_summarize validates its arguments", {
  expect_error(methyl_summarize(beta_matrix, "invalid"),
               "methyl_in must be either 'betas' or 'mvals'.", fixed = TRUE)
  expect_error(methyl_summarize(beta_matrix),
               "methyl_in")
  expect_error(methyl_summarize(beta_matrix, "betas", "invalid"),
               "methyl_out must be one of 'betas', 'mvals', or 'both'.", fixed = TRUE)
  expect_error(methyl_summarize(beta_matrix, "betas", cores = 0),
               "cores must be a single positive whole number.", fixed = TRUE)

  non_numeric_data <- matrix(c("A", "B", "C"), nrow = 1)
  expect_error(methyl_summarize(non_numeric_data, "betas"),
               "methyl_data must contain only numeric values.", fixed = TRUE)

  non_matrix_data <- c(0.1, 0.2, 0.3)
  expect_error(methyl_summarize(non_matrix_data, "betas"),
               "methyl_data must be a matrix.", fixed = TRUE)

  expect_error(methyl_summarize(unname(beta_matrix), "betas"),
               "methyl_data must have probe IDs as row names.", fixed = TRUE)
})

test_that("methyl_summarize issues a warning for likely transposed data", {
  expect_warning(methyl_summarize(t(beta_matrix), "betas"),
                 "more columns than rows")
})

test_that("methyl_summarize gives identical results in parallel", {
  skip_on_cran()
  expect_identical(methyl_summarize(beta_matrix, "betas", cores = 2),
                   methyl_summarize(beta_matrix, "betas", cores = 1))
})
