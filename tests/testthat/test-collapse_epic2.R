epic2 <- matrix(
  c(0.80, 0.70, 0.60,
    0.10, 0.20, 0.30,
    0.30, 0.30, 0.30,
    0.20, NA,   0.40,
    NA,   0.50, NA,
    NA,   NA,   0.90),
  nrow = 6, byrow = TRUE,
  dimnames = list(c("cg02_BC21", "cg01_TC21", "cg04", "cg01_TC22", "cg03_TC11", "cg03_TC12"),
                  c("s1", "s2", "s3"))
)

expected <- matrix(
  c(0.80, 0.70, 0.60,
    0.15, 0.20, 0.35,
    0.30, 0.30, 0.30,
    NA,   0.50, 0.90),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("cg02", "cg01", "cg04", "cg03"), c("s1", "s2", "s3"))
)

test_that("collapse_epic2 averages replicate probes in a matrix", {
  expect_equal(collapse_epic2(epic2), expected)
})

test_that("collapse_epic2 uses row names of a data frame", {
  expect_equal(collapse_epic2(as.data.frame(epic2)), as.data.frame(expected))
})

test_that("collapse_epic2 uses a probe ID column of a data frame", {
  epic2_df <- data.frame(s1 = epic2[, "s1"], probe = rownames(epic2),
                         s2 = epic2[, "s2"], s3 = epic2[, "s3"], row.names = NULL)
  result <- collapse_epic2(epic2_df, probe_names = "probe")

  expect_named(result, c("probe", "s1", "s2", "s3"))
  expect_identical(result$probe, rownames(expected))
  expect_identical(rownames(result), as.character(seq_len(nrow(expected))))
  expect_equal(as.matrix(result[c("s1", "s2", "s3")]), expected, ignore_attr = TRUE)
})

test_that("collapse_epic2 validates its arguments", {
  expect_error(collapse_epic2(epic2[, 1]),
               "methyl_data must be a matrix or a data frame.", fixed = TRUE)
  expect_error(collapse_epic2(unname(epic2)),
               "Row names are missing in methyl_data.", fixed = TRUE)
  expect_error(collapse_epic2(data.frame(s1 = 1:2)),
               "Row names are missing in methyl_data.", fixed = TRUE)
  expect_error(collapse_epic2(epic2, probe_names = "probe"),
               "probe_names can only be used when methyl_data is a data frame.", fixed = TRUE)
  expect_error(collapse_epic2(as.data.frame(epic2), probe_names = "probe"),
               "The specified probe_names column does not exist.", fixed = TRUE)

  non_numeric <- data.frame(probe = c("cg01_TC21", "cg01_TC22"), s1 = c("a", "b"))
  expect_error(collapse_epic2(non_numeric, probe_names = "probe"),
               "All columns must contain numeric values", fixed = TRUE)
})
