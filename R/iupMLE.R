#' Maximum likelihood estimator for infectious units per million (IUPM) to incorporate partial deep sequencing information
#' @name iupMLE
#' @param data Assay data, with rows representing the distinct viral lineages (DVL) and columns representing the wells. Any columns with \code{NA} values will be treated as positive but missing deep sequencing information.
#' @param maxit The maximum number of iterations (passed to \code{optim}). Default is \code{maxit = 1E4}.
#' @param lb Lower-bound on the IUPM (passed to \code{optim}). Default is \code{lb = 1E-6}.
#' @param ub Upper-bound on the IUPM (passed to \code{optim}). Default is \code{ub = Inf}.
#' @return A named list with the following slots:
#' \item{mle}{MLE of the IUPM}
#' \item{bc_mle}{Bias-corrected MLE of the IUPM}
#' \item{se}{Standard error for the MLE of the IUPM}
#' \item{se_bc}{Standard error for the bias-corrected MLE of the IUPM}
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
  se = sqrt(diag(cov))

  # Bias corrected MLE
  lambda_hat_bc = BC(lambda = lambda_hat,
                     M = M,
                     q = q)
  Lambda_hat_bc = sum(lambda_hat_bc)

  ## Variance for BC-MLE
  nabla = grad_bc(lambda = lambda_hat,
                  M = M,
                  q = q)

  var_bc = vector()
  for (i in 1:n) {
    nabla_i = matrix(data = nabla[, i],
                     ncol = 1)
    var_bc = append(var_bc,
                    t(nabla_i) %*% cov %*% nabla_i)
  }
  se_bc = sqrt(var_bc)

  return(list("mle" = Lambda_hat,
              "bc_mle" = Lambda_hat_bc,
              "se" = sqrt(sum(se ^ 2)),
              "se_bc" = sqrt(sum(se_bc ^ 2))
              )
         )
}
