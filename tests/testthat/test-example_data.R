test_that("pheno flags the second and later samples from each person", {
  expect_identical(pheno$dupe, as.integer(duplicated(pheno$person_id)))
})

test_that("epic2_betas collapses to the CpGs in beta_matrix", {
  collapsed <- collapse_epic2(epic2_betas)
  expect_identical(dimnames(collapsed), dimnames(beta_matrix))

  cpg_ids <- sub("_.*$", "", rownames(epic2_betas))
  single <- rownames(beta_matrix) %in% names(which(table(cpg_ids) == 1))
  expect_identical(collapsed[single, ], beta_matrix[single, ])
})
