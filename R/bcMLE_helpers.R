fisher = function(lambda, M, q) {
  Lambda = sum(lambda)

  # Calculate Fisher's information
  I = M * ( (1 - q) * exp(- Lambda) + q) / (exp( lambda) - 1)
  I = I + M * (1 - q) / (exp(Lambda) - 1)

  # Return
  return(I)
}

deriv_fisher = function(lambda, M, q) {
  Lambda = sum(lambda)

  # Calculate Fisher's information
  dI = - M * (1 - q) * exp(Lambda) / (exp(Lambda) - 1) ** 2
  dI = dI - M * (1 - q) * exp(- Lambda) / (exp(lambda) - 1)
  dI = dI - M * (q + (1 - q) * exp(- Lambda)) * exp(lambda) / (exp(lambda) - 1) ^ 2

  # Return
  return(dI)
}

exp_third_deriv = function(lambda, M, q) {
  Lambda = sum(lambda)

  exp_d3 = M * (q + (1 - q) * exp(- Lambda)) * (exp(lambda) + 1) / (exp(lambda) - 1) ^ 2
  exp_d3 = exp_d3 + M * (1 - q) * (exp(Lambda) + 1) / (exp(Lambda) - 1) ^ 2

  # Return
  return(exp_d3)
}

BC = function(lambda, M, q) {
  # Bias correction
  ## Calculate Fisher information
  I = fisher(
    lambda = lambda,
    M = M,
    q = q
  )

  ## Derivative of the Fisher information
  dI = deriv_fisher(
    lambda = lambda,
    M = M,
    q = q
  )

  ## Expectation of third Derivative
  expd3 = exp_third_deriv(
    lambda = lambda,
    M = M,
    q = q
  )

  ## Bias correction
  BC = - (2 * dI + expd3) / (2 * I ^ 2)

  ## Bias corrected MLE
  lambda = lambda - BC

  #Lambda_hat_bc = sum(lambda_hat - BC)
  return(lambda)
}

fisher_mat = function(lambda, M, q) {
  Lambda = sum(lambda)

  # Calculate off-diagonal elements (a constant)
  off_diag = M * (1 - q) / (exp(Lambda) - 1)

  I = matrix(data = off_diag,
             nrow = length(lambda),
             ncol = length(lambda))

  # Calculate diagonal elements
  diag(I) = fisher(lambda = lambda,
                   M = M,
                   q = q)

  # Return
  return(I)
}

grad_bc = function(lambda, M, q) {
  # BC-MLE (observed)
  lambda_bc = BC(lambda = lambda,
                 M = M,
                 q = q)

  # Fill the gradient to start with BC-MLE (observed)
  g = matrix(data = - lambda_bc,
             nrow = length(lambda),
             ncol = length(lambda),
             byrow = TRUE)

  # Size of the perturbation for derivatives
  hM = M ^ (- 1 / 2)

  # Calculate element-specific gradient terms
  for (i in 1:length(lambda)) {
    # Perturb the ith element of lambda
    lambda_pert = lambda
    lambda_pert[i] = lambda_pert[i] + hM

    # BC-MLE (perturbed)
    lambda_pert_bc = BC(lambda = lambda_pert,
                        M = M,
                        q = q)

    # Add it to the ith elements of the gradient
    g[i, ] = g[i, ] + lambda_pert_bc
  }
  g = 1 / hM * g

  return(g)
}
