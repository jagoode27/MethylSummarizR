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

#' Phenotype Data for the Example Matrices
#'
#' Synthetic phenotype data linking the 20 samples in [beta_matrix] and
#' [mval_matrix] to 11 people. Seven people have two samples, one person has
#' three samples, and three people have one sample, giving 8 duplicate pairs
#' for [methyl_reliability()].
#'
#' @format A data frame with 20 rows and 2 columns:
#' \describe{
#'   \item{sample_id}{Sample ID, matching the column names of the example matrices.}
#'   \item{person_id}{Person ID, shared by samples from the same person.}
#' }
#' @examples
#' data(pheno)
#' table(pheno$person_id)
"pheno"
