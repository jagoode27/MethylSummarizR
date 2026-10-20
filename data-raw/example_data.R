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

# flag the second and later samples from each person
pheno$dupe <- as.integer(duplicated(pheno$person_id))

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

# epic v2 version: suffixed probe IDs, with replicate probes for 9 CpGs
n_reps <- rep(1L, n_probes)
n_reps[sample(n_probes, 9)] <- c(rep(2L, 8), 3L)
cpg_rows <- rep(seq_len(n_probes), n_reps)
rep_num <- sequence(n_reps)
design <- sample(c("TC1", "TC2", "BC1", "BC2"), n_probes, replace = TRUE)

# replicate probes share the probe mean and person effect, with new noise and missingness
is_rep <- rep_num > 1
rep_noise <- matrix(stats::rnorm(sum(is_rep) * nrow(pheno)), nrow = sum(is_rep)) *
  noise_sd[cpg_rows[is_rep]]
rep_mvals <- (probe_mean[cpg_rows[is_rep]] +
                person_effect[cpg_rows[is_rep], pheno$person_id] + rep_noise) |>
  round(digits = 3)
rep_mvals[sample(length(rep_mvals), round(prop_miss * length(rep_mvals)))] <- NA
epic2_mvals <- mval_matrix[cpg_rows, ]
epic2_mvals[is_rep, ] <- rep_mvals
rownames(epic2_mvals) <- paste0(probe_ids[cpg_rows], "_", design[cpg_rows], rep_num)
epic2_betas <- convert_m_to_beta(epic2_mvals)

# save data
usethis::use_data(mval_matrix, beta_matrix, epic2_betas, pheno, overwrite = TRUE)
