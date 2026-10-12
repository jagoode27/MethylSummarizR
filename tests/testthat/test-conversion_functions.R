test_that("conversion functions are inverses", {
  expect_equal(convert_beta_to_m(convert_m_to_beta(mval_matrix)), mval_matrix)
  expect_equal(convert_m_to_beta(convert_beta_to_m(beta_matrix)), beta_matrix)
})

test_that("conversion functions compute known values and keep missing values", {
  ids <- list(c("cg1", "cg2"), c("s1", "s2"))
  mvals <- matrix(c(0, 2, -2, NA), nrow = 2, dimnames = ids)
  betas <- matrix(c(0.5, 0.8, 0.2, NA), nrow = 2, dimnames = ids)

  expect_equal(convert_m_to_beta(mvals), betas)
  expect_equal(convert_beta_to_m(betas), mvals)
})

test_that("conversion functions handle extreme values", {
  ids <- list(c("cg1", "cg2"), "s1")
  expect_equal(convert_m_to_beta(matrix(c(2000, -2000), dimnames = ids)),
               matrix(c(1, 0), dimnames = ids))

  m_values <- convert_beta_to_m(matrix(c(1, 0), dimnames = ids))
  expect_true(all(is.finite(m_values)))
  expect_equal(m_values, matrix(c(1, -1) * log2((1 - 1e-6) / 1e-6), dimnames = ids))
})

test_that("conversion functions validate their input", {
  for (convert in list(convert_m_to_beta, convert_beta_to_m)) {
    expect_error(convert(c(0.1, 0.2)), "Input must be a matrix.", fixed = TRUE)
    expect_error(convert(matrix("a")), "All values in the methylation matrix must be numeric.",
                 fixed = TRUE)
    expect_error(convert(matrix(numeric(0), nrow = 0, ncol = 0)), "Input matrix cannot be empty.",
                 fixed = TRUE)
    expect_warning(convert(matrix(0.5)), "Matrix lacks row names.", fixed = TRUE)
  }

  expect_warning(convert_beta_to_m(matrix(1.5, dimnames = list("cg1", "s1"))),
                 "Beta values outside [0,1] range detected.", fixed = TRUE)
})
