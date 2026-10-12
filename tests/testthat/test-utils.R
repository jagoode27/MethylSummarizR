test_that("run_chunked matches sequential results on a socket cluster", {
  skip_on_cran()
  expect_identical(run_chunked(summarize_chunk, list(mval_matrix), cores = 2, fork = FALSE),
                   summarize_chunk(mval_matrix))
})

test_that("run_chunked matches sequential results with forked workers", {
  skip_on_cran()
  skip_on_os("windows")
  expect_identical(run_chunked(summarize_chunk, list(mval_matrix), cores = 2, fork = TRUE),
                   summarize_chunk(mval_matrix))
})

test_that("run_chunked splits several matrices and caps cores at the row count", {
  methyl_1 <- mval_matrix[1:2, 1:8]
  methyl_2 <- mval_matrix[1:2, 9:16]
  expected <- reliability_chunk(methyl_1, methyl_2)
  expect_identical(run_chunked(reliability_chunk, list(methyl_1, methyl_2), cores = 1), expected)

  skip_on_cran()
  expect_identical(run_chunked(reliability_chunk, list(methyl_1, methyl_2), cores = 4, fork = FALSE),
                   expected)
  skip_on_os("windows")
  expect_identical(run_chunked(reliability_chunk, list(methyl_1, methyl_2), cores = 4, fork = TRUE),
                   expected)
})

test_that("run_chunked reports errors from workers", {
  skip_on_cran()
  failing <- function(methyl_data) stop("chunk failed")
  expect_error(run_chunked(failing, list(mval_matrix), cores = 2, fork = FALSE), "chunk failed")
  skip_on_os("windows")
  expect_error(run_chunked(failing, list(mval_matrix), cores = 2, fork = TRUE), "chunk failed")
})

test_that("check_cores accepts only positive whole numbers", {
  expect_identical(check_cores(2), 2L)
  for (bad in list(0, -1, 1.5, NA, Inf, "2", c(1, 2))) {
    expect_error(check_cores(bad), "cores must be a single positive whole number.", fixed = TRUE)
  }
})

test_that("check_choice accepts only listed values", {
  expect_identical(check_choice("both", c("betas", "mvals", "both"), "methyl_out"), "both")
  expect_error(check_choice(c("betas", "mvals"), c("betas", "mvals"), "methyl_in"),
               "methyl_in must be either 'betas' or 'mvals'.", fixed = TRUE)
  expect_error(check_choice(NA, c("betas", "mvals"), "methyl_in"), "methyl_in must be")
})

test_that("check_methyl_data rejects empty input", {
  expect_error(check_methyl_data(matrix(numeric(0), nrow = 0, ncol = 2)),
               "methyl_data cannot be empty.", fixed = TRUE)
})

test_that("out_formats and to_scale map scales", {
  expect_identical(out_formats("both"), c(betas = "betas", mvals = "mvals"))
  expect_identical(out_formats("mvals"), c(mvals = "mvals"))

  expect_identical(to_scale(beta_matrix, "betas", "betas"), beta_matrix)
  expect_equal(to_scale(beta_matrix, "betas", "mvals"), mval_matrix)
  expect_equal(to_scale(mval_matrix, "mvals", "betas"), beta_matrix)
})
