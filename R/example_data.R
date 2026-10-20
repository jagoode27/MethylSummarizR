#' Beta Matrix with Missing Values
#'
#' A 30 x 20 matrix of synthetic methylation beta values with some missing data.
#' (Probes are in rows and samples are in columns.) Some samples are duplicates
#' from the same person; see [pheno].
#'
#' @format A matrix with 30 rows and 20 columns. Row names are probe IDs and
#'   column names are sample IDs. About 5 percent of values are missing.
#' @seealso [mval_matrix] for the same data as M-values.
#' @examples
#' data(beta_matrix)
#' str(beta_matrix)
"beta_matrix"

#' M-Values Matrix with Missing Values
#'
#' A 30 x 20 matrix of synthetic methylation M-values with some missing data.
#' (Probes are in rows and samples are in columns.) Some samples are duplicates
#' from the same person; see [pheno].
#'
#' @format A matrix with 30 rows and 20 columns. Row names are probe IDs and
#'   column names are sample IDs. About 5 percent of values are missing.
#' @seealso [beta_matrix] for the same data as beta values.
#' @examples
#' data(mval_matrix)
#' str(mval_matrix)
"mval_matrix"

#' EPIC v2 Beta Matrix with Replicate Probes
#'
#' A 40 x 20 matrix of synthetic methylation beta values with MethylationEPIC v2
#' probe IDs. It holds the 30 CpGs in [beta_matrix], each with an EPIC v2
#' suffix, and 9 CpGs are measured by two or three replicate probes. The first
#' probe of each CpG has the values in [beta_matrix], so [collapse_epic2()]
#' returns a matrix with the same row and column names as [beta_matrix].
#'
#' @format A matrix with 40 rows and 20 columns. Row names are probe IDs (for
#'   example, `cg29302050_BC21`) and column names are sample IDs matching
#'   [pheno]. About 5 percent of values are missing.
#' @seealso [collapse_epic2()] to average the replicate probes, and
#'   [beta_matrix] for the data with one probe per CpG.
#' @examples
#' data(epic2_betas)
#' dim(epic2_betas)
#' dim(collapse_epic2(epic2_betas))
"epic2_betas"

#' Phenotype Data for the Example Matrices
#'
#' Synthetic phenotype data linking the 20 samples in [beta_matrix],
#' [mval_matrix], and [epic2_betas] to 11 people. Seven people have two
#' samples, one person has three samples, and three people have one sample,
#' giving 8 duplicate pairs for [methyl_reliability()]. Keeping the rows with
#' `dupe` equal to 0 gives one sample per person.
#'
#' @format A data frame with 20 rows and 3 columns:
#' \describe{
#'   \item{sample_id}{Sample ID, matching the column names of the example matrices.}
#'   \item{person_id}{Person ID, shared by samples from the same person.}
#'   \item{dupe}{1 for the second and later samples from a person, in the order listed, and 0 otherwise.}
#' }
#' @examples
#' data(pheno)
#' table(pheno$person_id)
#' table(pheno$dupe)
"pheno"
