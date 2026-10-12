# set up
rm(list = ls())
gc()
devtools::load_all()
set.seed(20261012)

# sample and person IDs: 7 people with two samples, 1 with three, and 3 with one
n_probes <- 30
person_ids <- sprintf("p%02d", 1:11)
pheno <- data.frame(
  sample_id = sprintf("s%02d", 1:20),
  person_id = sample(c(rep(person_ids[1:7], each = 2), rep(person_ids[8], 3), person_ids[9:11]))
)

# m-values from a probe mean, a person effect, and probe-specific sample noise
probe_ids <- sprintf("cg%08d", sort(sample.int(99999999, n_probes)))
probe_mean <- stats::rnorm(n_probes, mean = 0, sd = 2)
person_effect <- matrix(
  stats::rnorm(n_probes * length(person_ids), mean = 0, sd = 0.8),
  nrow = n_probes,
  dimnames = list(NULL, person_ids)
)
noise_sd <- stats::runif(n_probes, min = 0.1, max = 0.8)
noise <- matrix(stats::rnorm(n_probes * nrow(pheno)), nrow = n_probes) * noise_sd
mval_matrix <- (probe_mean + person_effect[, pheno$person_id] + noise) |>
  round(digits = 3) |>
  `dimnames<-`(list(probe_ids, pheno$sample_id))

# add missingness
prop_miss <- 0.05
mval_matrix[sample(length(mval_matrix), round(prop_miss * length(mval_matrix)))] <- NA

# create beta matrix
beta_matrix <- convert_m_to_beta(mval_matrix)

# save data
usethis::use_data(mval_matrix, beta_matrix, pheno, overwrite = TRUE)
