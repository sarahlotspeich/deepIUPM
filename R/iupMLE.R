#' Maximum likelihood estimator for infectious units per million (IUPM) to incorporate partial deep sequencing information
#' @name iupMLE
#' @param data Assay data, with rows representing the distinct viral lineages (DVL) and columns representing the wells. Any columns with \code{NA} values will be treated as positive but missing deep sequencing information.
#' @param maxit The maximum number of iterations (passed to \code{optim}). Default is \code{maxit = 1E4}.
#' @param lb Lower-bound on the IUPM (passed to \code{optim}). Default is \code{lb = 1E-6}.
#' @param ub Upper-bound on the IUPM (passed to \code{optim}). Default is \code{ub = Inf}.
#' @return A named list with the following slots:
#' \item{mle}{MLE of the IUPM}
#' \item{bc_mle}{Bias-corrected MLE of the IUPM}
#' \item{cov}{Covariance matrix for the MLE of the IUPM}
#' \item{cov_bc}{Covariance matrix for the bias-corrected MLE of the IUPM}
#' @export
#'
iupMLE <- function(data, maxit = 1E4, lb = 1E-6, ub = Inf) {
# Helpful constants
  M = ncol(assay) # number of wells (total)
  n = nrow(data) # number of DVL detected
  MN = sum(colSums(assay) == 0, na.rm = TRUE) # number of wells (p24-negative)
  MP = M - MN # number of wells (p24-positive) = xM
  m = MP - sum(is.na(colSums(assay))) # number of deep-sequenced p24-positive wells = m1
  q = m / MP # proportion of p24-positive wells deep sequenced
  Y = rowSums(assay, na.rm = TRUE) # number of infected wells per DVL = count

  # Fit MLE
  optimization = optim(par = - log(1 - Y / M),
                       fn = loglik,
                       gr = gloglik,
                       n = n,
                       M = M,
                       xM = MP,
                       m1 = m,
                       count = Y,
                       method = "L-BFGS-B",
                       control = list(maxit = maxit),
                       lower = rep(lb, n),
                       upper = rep(ub, n),
                       hessian = T)
  lambda_hat = optimization$par
  Lambda_hat = sum(lambda_hat) # MLE of the IUPM
  cov = solve(optimization$hessian) #Variance of IUPM

  # Bias correction
  ## Calculate Fisher information
  I = fisher(
    lambda = lambda_hat,
    M = M,
    q = q
  )
  ## Standard error (from inverting Fisher's information)
  V = 1 / I
  SE = sqrt(sum(V))
  ## Derivative of the Fisher information
  dI = deriv_fisher(
    lambda = lambda_hat,
    M = M,
    q = q
  )
  ## Expectation of third Derivative
  expd3 = exp_third_deriv(
    lambda = lambda_hat,
    M = M,
    q = q
  )
  ## Bias correction
  BC = - (2 * dI + expd3) / (2 * I ^ 2)

  ## Bias corrected MLE
  Lambda_hat_bc = sum(lambda_hat - BC)

  return(list("mle" = Lambda_hat,
              "bc_mle" = Lambda_hat_bc,
              "cov" = cov,
              "cov_bc" = matrix(NA, nrow = n, ncol = n)
              )
         )
}
