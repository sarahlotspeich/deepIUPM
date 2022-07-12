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
