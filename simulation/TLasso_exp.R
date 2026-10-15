################################################################################
# Simulation Experiments for TLasso
# ------------------------------------------------------------------------------
# This script reproduces the results in Table 1 and Figure 4 in the main text
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Setup
source("helper/Graph Difference Operator.R")
source("helper/Beta Structure Generation.R")
source("../function/TLasso.R")
source("../function/other methods/HD_SplitKnockoff.R")
library(genlasso)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions
gamma_estimation_error <- function(gamma_est, gamma_true) {
  if (length(gamma_est) != length(gamma_true)) {
    stop("Dimension mismatch.")
  }
  diff <- gamma_est - gamma_true
  list(
    l2_norm = sqrt(sum(diff^2)),
    linf_norm = max(abs(diff))
  )
}

mean.se = function(x, ns = 3) { paste0(format(round(mean(x), ns), nsmall = ns), ' (', format(round(sd(x)/sqrt(length(x)), ns), nsmall = ns), ')')}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Simulation settings
set.seed(2025)
nIterations = 100
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Table 1-1: piecewise constant over Chain
################################################################################
################################################################################
################################################################################
### Simulation settings
set.seed(2025)
p = 150; n = 300; A = 0.5; rho = 0.5
p1 = 150
m1 = 20
Sigma = toeplitz(rho^(1:p-1))
D <- create.chain(p)
m <- nrow(D)
Gamma_list_ChainConstant = sapply(1:nIterations, function(it) {
  ################################################################################
  ################################################################################
  ################################################################################
  ### Generate the Data
  Beta0 = Chain_PiecewiseConstant(p,p1,m1,A)
  X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)
  signal_ind = which(abs(D%*%Beta0) > 1e-4)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (1)
  TLasso_fit1 <- TLasso(X, y, D, X_N_option = 1)
  gamma_TLasso1 <- TLasso_fit1$gamma
  l2_11 = gamma_estimation_error(gamma_TLasso1, D %*% Beta0)$l2_norm
  linf_11 = gamma_estimation_error(gamma_TLasso1, D %*% Beta0)$linf_norm
  TP_11 = length(intersect(which(abs(gamma_TLasso1) > 1e-4),signal_ind))
  MS_11 = length(which(abs(gamma_TLasso1) > 1e-4))
  print(l2_11)
  print(linf_11)
  print(TP_11)
  print(MS_11)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (2)
  TLasso_fit2 <- TLasso(X, y, D, X_N_option = 2)
  gamma_TLasso2 <- TLasso_fit2$gamma
  l2_12 = gamma_estimation_error(gamma_TLasso2, D %*% Beta0)$l2_norm
  linf_12 = gamma_estimation_error(gamma_TLasso2, D %*% Beta0)$linf_norm
  TP_12 = length(intersect(which(abs(gamma_TLasso2) > 1e-4),signal_ind))
  MS_12 = length(which(abs(gamma_TLasso2) > 1e-4))
  print(l2_12)
  print(linf_12)
  print(TP_12)
  print(MS_12)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (3)
  TLasso_fit3 <- TLasso(X, y, D, X_N_option = 3)
  gamma_TLasso3 <- TLasso_fit3$gamma
  l2_13 = gamma_estimation_error(gamma_TLasso3, D %*% Beta0)$l2_norm
  linf_13 = gamma_estimation_error(gamma_TLasso3, D %*% Beta0)$linf_norm
  TP_13 = length(intersect(which(abs(gamma_TLasso3) > 1e-4),signal_ind))
  MS_13 = length(which(abs(gamma_TLasso3) > 1e-4))
  print(l2_13)
  print(linf_13)
  print(TP_13)
  print(MS_13)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### GenLasso
  lambda_GenLasso = sqrt(n*log(m))
  GenLasso_fit <- genlasso(y, X, D, minlam = lambda_GenLasso)
  beta_GenLasso <- coef(GenLasso_fit, lambda = lambda_GenLasso)$beta
  gamma_GenLasso <- D %*% beta_GenLasso
  l2_2 = gamma_estimation_error(gamma_GenLasso, D %*% Beta0)$l2_norm
  linf_2 = gamma_estimation_error(gamma_GenLasso, D %*% Beta0)$linf_norm
  TP_2 = length(intersect(which(abs(gamma_GenLasso) > 1e-4),signal_ind))
  MS_2 = length(which(abs(gamma_GenLasso) > 1e-4))
  print(l2_2)
  print(linf_2)
  print(TP_2)
  print(MS_2)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### SplitLasso
  SplitLasso_fit <- SplitLasso(X, y, D, nu_values = 10)
  gamma_SplitLasso <- SplitLasso_fit$gammas
  l2_3 = gamma_estimation_error(gamma_SplitLasso, D %*% Beta0)$l2_norm
  linf_3 = gamma_estimation_error(gamma_SplitLasso, D %*% Beta0)$linf_norm
  TP_3 = length(intersect(which(abs(gamma_SplitLasso) > 1e-4),signal_ind))
  MS_3 = length(which(abs(gamma_SplitLasso) > 1e-4))
  print(l2_3)
  print(linf_3)
  print(TP_3)
  print(MS_3)
  
  r = c(l2_11, l2_12, l2_13, l2_2, l2_3,
        linf_11, linf_12, linf_13, linf_2, linf_3,
        TP_11, TP_12, TP_13, TP_2, TP_3,
        MS_11, MS_12, MS_13, MS_2, MS_3)
})
################################################################################
################################################################################
################################################################################
### Results
apply(Gamma_list_ChainConstant, 1, mean.se)
save(Gamma_list_ChainConstant, file = "Gamma_list_ChainConstant.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Table 1-2: piecewise constant over Lattice
################################################################################
################################################################################
################################################################################
### Simulation settings
set.seed(2025)
p = 150; n = 300; A = 0.5; rho = 0.5
s = 15
t = 10
p1 = 150
m1 = 20
Sigma = toeplitz(rho^(1:p-1))
D <- create.lattice(s,t)
m <- nrow(D)

Gamma_list_LatticeConstant = sapply(1:nIterations, function(it) {
  ################################################################################
  ################################################################################
  ################################################################################
  ### Generate the Data
  Beta0 = Lattice_PiecewiseConstant(s, t, p1, m1, A)
  X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)
  signal_ind = which(abs(D%*%Beta0) > 1e-4)
  print(min(abs(D %*% Beta0)[signal_ind]))
  plot(D %*% Beta0)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (1)
  TLasso_fit1 <- TLasso(X, y, D, X_N_option = 1)
  gamma_TLasso1 <- TLasso_fit1$gamma
  l2_11 = gamma_estimation_error(gamma_TLasso1, D %*% Beta0)$l2_norm
  linf_11 = gamma_estimation_error(gamma_TLasso1, D %*% Beta0)$linf_norm
  TP_11 = length(intersect(which(abs(gamma_TLasso1) > 1e-4),signal_ind))
  MS_11 = length(which(abs(gamma_TLasso1) > 1e-4))
  print(l2_11)
  print(linf_11)
  print(TP_11)
  print(MS_11)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (2)
  TLasso_fit2 <- TLasso(X, y, D, X_N_option = 2)
  gamma_TLasso2 <- TLasso_fit2$gamma
  l2_12 = gamma_estimation_error(gamma_TLasso2, D %*% Beta0)$l2_norm
  linf_12 = gamma_estimation_error(gamma_TLasso2, D %*% Beta0)$linf_norm
  TP_12 = length(intersect(which(abs(gamma_TLasso2) > 1e-4),signal_ind))
  MS_12 = length(which(abs(gamma_TLasso2) > 1e-4))
  print(l2_12)
  print(linf_12)
  print(TP_12)
  print(MS_12)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (3)
  TLasso_fit3 <- TLasso(X, y, D, X_N_option = 3)
  gamma_TLasso3 <- TLasso_fit3$gamma
  l2_13 = gamma_estimation_error(gamma_TLasso3, D %*% Beta0)$l2_norm
  linf_13 = gamma_estimation_error(gamma_TLasso3, D %*% Beta0)$linf_norm
  TP_13 = length(intersect(which(abs(gamma_TLasso3) > 1e-4),signal_ind))
  MS_13 = length(which(abs(gamma_TLasso3) > 1e-4))
  print(l2_13)
  print(linf_13)
  print(TP_13)
  print(MS_13)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### GenLasso
  lambda_GenLasso = sqrt(n*log(m))
  GenLasso_fit <- genlasso(y, X, D, minlam = lambda_GenLasso)
  beta_GenLasso <- coef(GenLasso_fit, lambda = lambda_GenLasso)$beta
  gamma_GenLasso <- D %*% beta_GenLasso
  l2_2 = gamma_estimation_error(gamma_GenLasso, D %*% Beta0)$l2_norm
  linf_2 = gamma_estimation_error(gamma_GenLasso, D %*% Beta0)$linf_norm
  TP_2 = length(intersect(which(abs(gamma_GenLasso) > 1e-4),signal_ind))
  MS_2 = length(which(abs(gamma_GenLasso) > 1e-4))
  print(l2_2)
  print(linf_2)
  print(TP_2)
  print(MS_2)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### SplitLasso
  SplitLasso_fit <- SplitLasso(X, y, D, nu_values = 10)
  gamma_SplitLasso <- SplitLasso_fit$gammas
  l2_3 = gamma_estimation_error(gamma_SplitLasso, D %*% Beta0)$l2_norm
  linf_3 = gamma_estimation_error(gamma_SplitLasso, D %*% Beta0)$linf_norm
  TP_3 = length(intersect(which(abs(gamma_SplitLasso) > 1e-4),signal_ind))
  MS_3 = length(which(abs(gamma_SplitLasso) > 1e-4))
  print(l2_3)
  print(linf_3)
  print(TP_3)
  print(MS_3)
  
  r = c(l2_11, l2_12, l2_13, l2_2, l2_3,
        linf_11, linf_12, linf_13, linf_2, linf_3,
        TP_11, TP_12, TP_13, TP_2, TP_3,
        MS_11, MS_12, MS_13, MS_2, MS_3)
})
################################################################################
################################################################################
################################################################################
### Results
apply(Gamma_list_LatticeConstant, 1, mean.se)
save(Gamma_list_LatticeConstant, file = "Gamma_list_LatticeConstant.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Table 1-3: piecewise linear over Lattice
################################################################################
################################################################################
################################################################################
### Simulation settings
set.seed(2025)
p = 150; n = 300; A = 0.5; rho = 0.5
s = 15
t = 10
p1 = 150
m1 = 20
Sigma = toeplitz(rho^(1:p-1))
D <- create.lattice(s,t)
D <- t(D) %*% D
m <- nrow(D)

Gamma_list_LatticeLinear = sapply(1:nIterations, function(it) {
  ################################################################################
  ################################################################################
  ################################################################################
  ### Generate the Data
  Beta0 = Lattice_PiecewiseConstant(s, t, p1, m1/5*4, A)
  X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)
  signal_ind = which(abs(D%*%Beta0) > 1e-4)
  print(min(abs(D %*% Beta0)[signal_ind]))
  plot(D %*% Beta0)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (1)
  TLasso_fit1 <- TLasso(X, y, D, X_N_option = 1)
  gamma_TLasso1 <- TLasso_fit1$gamma
  l2_11 = gamma_estimation_error(gamma_TLasso1, D %*% Beta0)$l2_norm
  linf_11 = gamma_estimation_error(gamma_TLasso1, D %*% Beta0)$linf_norm
  TP_11 = length(intersect(which(abs(gamma_TLasso1) > 1e-4),signal_ind))
  MS_11 = length(which(abs(gamma_TLasso1) > 1e-4))
  print(l2_11)
  print(linf_11)
  print(TP_11)
  print(MS_11)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (2)
  TLasso_fit2 <- TLasso(X, y, D, X_N_option = 2)
  gamma_TLasso2 <- TLasso_fit2$gamma
  l2_12 = gamma_estimation_error(gamma_TLasso2, D %*% Beta0)$l2_norm
  linf_12 = gamma_estimation_error(gamma_TLasso2, D %*% Beta0)$linf_norm
  TP_12 = length(intersect(which(abs(gamma_TLasso2) > 1e-4),signal_ind))
  MS_12 = length(which(abs(gamma_TLasso2) > 1e-4))
  print(l2_12)
  print(linf_12)
  print(TP_12)
  print(MS_12)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso (3)
  TLasso_fit3 <- TLasso(X, y, D, X_N_option = 3)
  gamma_TLasso3 <- TLasso_fit3$gamma
  l2_13 = gamma_estimation_error(gamma_TLasso3, D %*% Beta0)$l2_norm
  linf_13 = gamma_estimation_error(gamma_TLasso3, D %*% Beta0)$linf_norm
  TP_13 = length(intersect(which(abs(gamma_TLasso3) > 1e-4),signal_ind))
  MS_13 = length(which(abs(gamma_TLasso3) > 1e-4))
  print(l2_13)
  print(linf_13)
  print(TP_13)
  print(MS_13)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### GenLasso
  lambda_GenLasso = sqrt(n*log(m))
  GenLasso_fit <- genlasso(y, X, D, minlam = lambda_GenLasso)
  beta_GenLasso <- coef(GenLasso_fit, lambda = lambda_GenLasso)$beta
  gamma_GenLasso <- D %*% beta_GenLasso
  l2_2 = gamma_estimation_error(gamma_GenLasso, D %*% Beta0)$l2_norm
  linf_2 = gamma_estimation_error(gamma_GenLasso, D %*% Beta0)$linf_norm
  TP_2 = length(intersect(which(abs(gamma_GenLasso) > 1e-4),signal_ind))
  MS_2 = length(which(abs(gamma_GenLasso) > 1e-4))
  print(l2_2)
  print(linf_2)
  print(TP_2)
  print(MS_2)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### SplitLasso
  SplitLasso_fit <- SplitLasso(X, y, D, nu_values = 10)
  gamma_SplitLasso <- SplitLasso_fit$gammas
  l2_3 = gamma_estimation_error(gamma_SplitLasso, D %*% Beta0)$l2_norm
  linf_3 = gamma_estimation_error(gamma_SplitLasso, D %*% Beta0)$linf_norm
  TP_3 = length(intersect(which(abs(gamma_SplitLasso) > 1e-4),signal_ind))
  MS_3 = length(which(abs(gamma_SplitLasso) > 1e-4))
  print(l2_3)
  print(linf_3)
  print(TP_3)
  print(MS_3)
  
  r = c(l2_11, l2_12, l2_13, l2_2, l2_3,
        linf_11, linf_12, linf_13, linf_2, linf_3,
        TP_11, TP_12, TP_13, TP_2, TP_3,
        MS_11, MS_12, MS_13, MS_2, MS_3)
})
################################################################################
################################################################################
################################################################################
### Results
apply(Gamma_list_LatticeLinear, 1, mean.se)
save(Gamma_list_LatticeLinear, file = "Gamma_list_LatticeLinear.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure 4: 1D sparse fused Lasso problem
################################################################################
################################################################################
################################################################################
### Simulation settings
set.seed(2025)
p = 150; n = 300; A = 1; rho = 0.2
p1 = 30
m1 = 8
Sigma = toeplitz(rho^(1:p-1))

D_1 <- diag(p)
m_1 <- nrow(D_1)
D_2 <- create.chain(p)
m_2 <- nrow(D_2)
D_3 <- rbind(diag(p), create.chain(p))
m_3 <- nrow(D_3)

Gamma_list_1D_SFL = sapply(1:nIterations, function(it) {
  ################################################################################
  ################################################################################
  ################################################################################
  ### Generate the Data
  Beta0 = Chain_PiecewiseConstant(p,p1,m1,A)
  X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)
  signal_ind_1 = which(abs(D_1%*%Beta0) > 1e-4)
  signal_ind_2 = which(abs(D_2%*%Beta0) > 1e-4)
  signal_ind_3 = which(abs(D_3%*%Beta0) > 1e-4)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso
  TLasso_fit1 <- TLasso(X, y, D_1)
  gamma_TLasso1 <- TLasso_fit1$gamma
  l2_11 = gamma_estimation_error(gamma_TLasso1, D_1 %*% Beta0)$l2_norm
  linf_11 = gamma_estimation_error(gamma_TLasso1, D_1 %*% Beta0)$linf_norm
  TP_11 = length(intersect(which(abs(gamma_TLasso1) > 1e-4),signal_ind_1))
  MS_11 = length(which(abs(gamma_TLasso1) > 1e-4))
  print(l2_11)
  print(linf_11)
  print(TP_11)
  print(MS_11)
  TLasso_fit2 <- TLasso(X, y, D_2)
  gamma_TLasso2 <- TLasso_fit2$gamma
  l2_12 = gamma_estimation_error(gamma_TLasso2, D_2 %*% Beta0)$l2_norm
  linf_12 = gamma_estimation_error(gamma_TLasso2, D_2 %*% Beta0)$linf_norm
  TP_12 = length(intersect(which(abs(gamma_TLasso2) > 1e-4),signal_ind_2))
  MS_12 = length(which(abs(gamma_TLasso2) > 1e-4))
  print(l2_12)
  print(linf_12)
  print(TP_12)
  print(MS_12)
  TLasso_fit3 <- TLasso(X, y, D_3)
  gamma_TLasso3 <- TLasso_fit3$gamma
  l2_13 = gamma_estimation_error(gamma_TLasso3, D_3 %*% Beta0)$l2_norm
  linf_13 = gamma_estimation_error(gamma_TLasso3, D_3 %*% Beta0)$linf_norm
  TP_13 = length(intersect(which(abs(gamma_TLasso3) > 1e-4),signal_ind_3))
  MS_13 = length(which(abs(gamma_TLasso3) > 1e-4))
  print(l2_13)
  print(linf_13)
  print(TP_13)
  print(MS_13)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### GenLasso
  lambda_GenLasso = sqrt(n*log(m_1))
  GenLasso_fit1 <- genlasso(y, X, D_1, minlam = lambda_GenLasso)
  beta_GenLasso1 <- coef(GenLasso_fit1, lambda = lambda_GenLasso)$beta
  gamma_GenLasso1 <- D_1 %*% beta_GenLasso1
  l2_21 = gamma_estimation_error(gamma_GenLasso1, D_1 %*% Beta0)$l2_norm
  linf_21 = gamma_estimation_error(gamma_GenLasso1, D_1 %*% Beta0)$linf_norm
  TP_21 = length(intersect(which(abs(gamma_GenLasso1) > 1e-4),signal_ind_1))
  MS_21 = length(which(abs(gamma_GenLasso1) > 1e-4))
  print(l2_21)
  print(linf_21)
  print(TP_21)
  print(MS_21)
  lambda_GenLasso = sqrt(n*log(m_2))
  GenLasso_fit2 <- genlasso(y, X, D_2, minlam = lambda_GenLasso)
  beta_GenLasso2 <- coef(GenLasso_fit2, lambda = lambda_GenLasso)$beta
  gamma_GenLasso2 <- D_2 %*% beta_GenLasso2
  l2_22 = gamma_estimation_error(gamma_GenLasso2, D_2 %*% Beta0)$l2_norm
  linf_22 = gamma_estimation_error(gamma_GenLasso2, D_2 %*% Beta0)$linf_norm
  TP_22 = length(intersect(which(abs(gamma_GenLasso2) > 1e-4),signal_ind_2))
  MS_22 = length(which(abs(gamma_GenLasso2) > 1e-4))
  print(l2_22)
  print(linf_22)
  print(TP_22)
  print(MS_22)
  lambda_GenLasso = sqrt(n*log(m_3))
  GenLasso_fit3 <- genlasso(y, X, D_3, minlam = lambda_GenLasso)
  beta_GenLasso3 <- coef(GenLasso_fit3, lambda = lambda_GenLasso)$beta
  gamma_GenLasso3 <- D_3 %*% beta_GenLasso3
  l2_23 = gamma_estimation_error(gamma_GenLasso3, D_3 %*% Beta0)$l2_norm
  linf_23 = gamma_estimation_error(gamma_GenLasso3, D_3 %*% Beta0)$linf_norm
  TP_23 = length(intersect(which(abs(gamma_GenLasso3) > 1e-4),signal_ind_3))
  MS_23 = length(which(abs(gamma_GenLasso3) > 1e-4))
  print(l2_23)
  print(linf_23)
  print(TP_23)
  print(MS_23)
  
  ################################################################################
  ################################################################################
  ################################################################################
  ### SplitLasso
  SplitLasso_fit1 <- SplitLasso(X, y, D_1, nu_values = 10)
  gamma_SplitLasso1 <- SplitLasso_fit1$gammas
  l2_31 = gamma_estimation_error(gamma_SplitLasso1, D_1 %*% Beta0)$l2_norm
  linf_31 = gamma_estimation_error(gamma_SplitLasso1, D_1 %*% Beta0)$linf_norm
  TP_31 = length(intersect(which(abs(gamma_SplitLasso1) > 1e-4),signal_ind_1))
  MS_31 = length(which(abs(gamma_SplitLasso1) > 1e-4))
  print(l2_31)
  print(linf_31)
  print(TP_31)
  print(MS_31)
  SplitLasso_fit2 <- SplitLasso(X, y, D_2, nu_values = 10)
  gamma_SplitLasso2 <- SplitLasso_fit2$gammas
  l2_32 = gamma_estimation_error(gamma_SplitLasso2, D_2 %*% Beta0)$l2_norm
  linf_32 = gamma_estimation_error(gamma_SplitLasso2, D_2 %*% Beta0)$linf_norm
  TP_32 = length(intersect(which(abs(gamma_SplitLasso2) > 1e-4),signal_ind_2))
  MS_32 = length(which(abs(gamma_SplitLasso2) > 1e-4))
  print(l2_32)
  print(linf_32)
  print(TP_32)
  print(MS_32)
  SplitLasso_fit3 <- SplitLasso(X, y, D_3, nu_values = 10)
  gamma_SplitLasso3 <- SplitLasso_fit3$gammas
  l2_33 = gamma_estimation_error(gamma_SplitLasso3, D_3 %*% Beta0)$l2_norm
  linf_33 = gamma_estimation_error(gamma_SplitLasso3, D_3 %*% Beta0)$linf_norm
  TP_33 = length(intersect(which(abs(gamma_SplitLasso3) > 1e-4),signal_ind_3))
  MS_33 = length(which(abs(gamma_SplitLasso3) > 1e-4))
  print(l2_33)
  print(linf_33)
  print(TP_33)
  print(MS_33)
  
  r = c(l2_11, l2_12, l2_13, l2_21, l2_22, l2_23, l2_31, l2_32, l2_33,
        linf_11, linf_12, linf_13, linf_21, linf_22, linf_23, linf_31, linf_32, linf_33,
        TP_11, TP_12, TP_13, TP_21, TP_22, TP_23, TP_31, TP_32, TP_33,
        MS_11, MS_12, MS_13, MS_21, MS_22, MS_23, MS_31, MS_32, MS_33)
})
################################################################################
################################################################################
################################################################################
### Results
apply(Gamma_list_1D_SFL, 1, mean.se)
save(Gamma_list_1D_SFL, file = "Gamma_list_1D_SFL.RData")
