data(mval_matrix)
data(beta_matrix)

test_that("methyl_summarize handles valid mvals input", {
  result <- methyl_summarize(mval_matrix, "mvals")

  expect_type(result, "list")
  expect_named(result, c("betas", "mvals"))
  expect_s3_class(result$betas, "data.frame")
  expect_s3_class(result$mvals, "data.frame")
  expect_true(all(c("probe", "mean", "sd", "ci95_lb", "ci95_ub", "coef_var", "iqr",
                    "min", "p05", "p25", "p50", "p75", "p95", "max",
                    "miss_n", "nonmiss_n", "miss_p", "nonmiss_p") %in% names(result$mvals)))
})

test_that("methyl_summarize handles valid betas input", {
  result <- methyl_summarize(beta_matrix, "betas")

  expect_type(result, "list")
  expect_named(result, c("betas", "mvals"))
  expect_s3_class(result$betas, "data.frame")
  expect_s3_class(result$mvals, "data.frame")
})

test_that("methyl_summarize correctly handles missing data", {
  result_betas <- methyl_summarize(beta_matrix, "betas")$betas

  expect_true(all(result_betas$miss_n >= 0))
  expect_true(all(result_betas$miss_n <= length(beta_matrix[1, ])))
  expect_true(all(is.na(result_betas$mean) == (result_betas$nonmiss_n == 0)))
})

test_that("methyl_summarize error handling is correct for invalid input", {
  expect_error(methyl_summarize(beta_matrix, "invalid"),
               "methyl_format must be either 'betas' or 'mvals'.")

  non_numeric_data <- matrix(c("A", "B", "C"), nrow = 1)
  expect_error(methyl_summarize(non_numeric_data, "betas"),
               "methyl_data must contain only numeric values.")

  non_matrix_data <- c(0.1, 0.2, 0.3)
  expect_error(methyl_summarize(non_matrix_data, "betas"),
               "methyl_data must be a matrix.")
})

test_that("methyl_summarize issues correct warning for orientation", {
  # Assuming row should be probes, simulate possible orientation warning
  expect_warning(methyl_summarize(t(beta_matrix), "betas"),
                 "more columns than rows")
})
