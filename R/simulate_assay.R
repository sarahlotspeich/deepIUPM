#' Simulate a dilution assay with partial deep sequencing information
#' @name simulate_assay
#' @param M Total number of wells.
#' @param n Total number of underlying distinct viral lineages (DVL) to detect
#' @param lambda A vector of length \code{n} with the expected number of infected cells per DVL
#' @param q Proportion of p24-positive wells to be deep sequenced
#' @param remove_undetected (Logical) If \code{remove_undetected = TRUE} (the default), DVL that were not detected in any wells are removed
#' @return A matrix with \code{n} rows (one per DVL) and \code{M} columns (one per well)
#' @export
#'

simulate_assay <- function(M, n, lambda, q, remove_undetected = TRUE) {
  # Generate Z = I(Xij >= 1) for all i, j
  Z <- rbinom(n = M * n,
              size = 1,
              prob = (1 - exp(- rep(x = lambda, each = M))))
  data_mat <- matrix(data = Z, nrow = n, ncol = M, byrow = TRUE)

  # Calculate M_P: number of wells with >= 1 DVL (assumed to be p24+)
  MP <- sum(colSums(data_mat) >= 1)

  # Calculate the number of p24 positive wells for VOA-USDA
  if (q == 0) {
    m <- MP
  } else {
    m <- ceiling((1 - q) * MP)
  }

  # Randomly select columns to make missing
  p24_pos <- which(colSums(data_mat) >= 1)
  p24_neg <- which(colSums(data_mat) == 0)
  if (q > 0) {
    make_miss <- p24_pos[1:(MP - m)]
    data_mat[, make_miss] <- NA
  }

  if (remove_undetected) {
    # Check for which DVL detected
    detected <- rowSums(data_mat, na.rm = TRUE) >= 1
    # Number of DVL detected
    n_ <- sum(detected)
    data_mat <- data_mat[detected, ]
    if (n_ == 0) {
      return(NULL)
    } else if (n_ == 1) {
      data_mat <- matrix(data_mat, ncol = M, byrow = TRUE)
    }
  }

  # Reorder to look like the paper
  ## With non-missing data first
  #data_mat <- matrix(data_mat, ncol = M, byrow = TRUE)
  if (q > 0) {
    data_mat <- data_mat[, c(setdiff(p24_pos, make_miss), p24_neg, make_miss)]
  } else {
    data_mat <- data_mat[, c(p24_pos, p24_neg)]
  }
  if (n_ == 1) {
    data_mat <- matrix(data_mat, ncol = M, byrow = TRUE)
  }
  return(data_mat)
}
