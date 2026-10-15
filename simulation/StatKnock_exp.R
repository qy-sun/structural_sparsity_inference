################################################################################
# Simulation Experiments for StatKnock
# ------------------------------------------------------------------------------
# This script reproduces Figures S.8-S.11 in the Supplement
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
source("../function/StatKnock.R")
source("../function/other methods/BY_Structural.R")
source("../function/other methods/DS_Structural.R")
source("../function/other methods/HD_SplitKnockoff.R")
library(GKnockoff)     # package version 0.1.0
library(SplitKnockoff) # package version 1.2
library(genlasso)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions
fdp_power <- function(selected_index, signal_index){
  fdp = 0
  power = 0
  ssr = 0
  if (length(selected_index) != 0){
    num_selected <- length(selected_index)
    tp <- length(intersect(selected_index, signal_index))
    fp <- num_selected - tp
    fdp <- fp / max(num_selected, 1)
    power <- tp / length(signal_index)
    if (power == 1){
      ssr = 1
    }else{
      ssr = 0
    }
  }  
  return(list(fdp = fdp, power = power, ssr = ssr))
}

mean.se = function(x, ns = 1) { 
  paste0(format(round(mean(x)*100, ns), nsmall = ns), ' (', 
         format(round(100*sd(x)/sqrt(length(x)), ns), nsmall = ns), ')')}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Simulation settings
set.seed(2025)
nIterations = 100
q = 0.2
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.8 (a): piecewise constant over chain (n=500, p=120)
set.seed(2025)
p = 120; n = 500
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- c(0.2, 0.5, 0.8)
ChainPiecewiseConstant_result_all_n500 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ################################################################################
      ################################################################################
      ################################################################################
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      y <- X %*% Beta0 + rnorm(n, sd = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-4)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### DS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel <- DS_Fit$DS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### Split Knockoff
      nu_hat <- SplitLasso(normc(X), normc(y), D)$best_nu
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = nu_hat,
        sign = 'disabled'
      )
      SplitKnockoff_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_Sel <- SplitKnockoff_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_result_all_n500[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_result_all_n500[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_result_all_n500, file = "StatKnock_ChainPiecewiseConstant_n500signal12.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.8 (b): piecewise constant over chain (n=250, p=120)
set.seed(2025)
p = 120; n = 250
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- c(0.2, 0.5, 0.8)
ChainPiecewiseConstant_result_all_n250 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ################################################################################
      ################################################################################
      ################################################################################
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      y <- X %*% Beta0 + rnorm(n, sd = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-4)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### DS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel <- DS_Fit$DS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### Split Knockoff
      nu_hat <- SplitLasso(normc(X), normc(y), D)$best_nu
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = nu_hat,
        sign = 'disabled'
      )
      SplitKnockoff_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_Sel <- SplitKnockoff_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_result_all_n250[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_result_all_n250[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_result_all_n250, file = "StatKnock_ChainPiecewiseConstant_n250signal12.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.9: 1D sparse fused Lasso

set.seed(2025)
p = 150; n = 600
p1 = 40
m1 = 12
D <- rbind(diag(p), create.chain(p))
m <- nrow(D)

rho_values <- c(0.5, 0.8)
A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)

SFL_1D_result_all <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))

  SFL_1D_list <- list()
  
  for (A in A_values) {
    cat("Running for A =", A, "\n")
    
    fdp_power_list_SFL_1D = sapply(1:nIterations, function(it) {
      ################################################################################
      ################################################################################
      ################################################################################
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      y <- X %*% Beta0 + rnorm(n, sd = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-4)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### DS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel <- DS_Fit$DS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      ################################################################################
      ################################################################################
      ################################################################################
      ### Split Knockoff
      nu_hat <- SplitLasso(normc(X), normc(y), D)$best_nu
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = nu_hat,
        sign = 'disabled'
      )
      SplitKnockoff_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_Sel <- SplitKnockoff_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    SFL_1D_list[[paste0("A=", A)]] <- fdp_power_list_SFL_1D
  }
  
  SFL_1D_result_all[[paste0("rho=", rho)]] <- SFL_1D_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- SFL_1D_result_all[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(SFL_1D_result_all, file = "StatKnock_SFL_1D.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.10: heavy-tailed t(3) design (n=250, p=120)
set.seed(2025)
df_t <- 3
p = 120; n = 250
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- 0.2
ChainPiecewiseConstant_t3_result_all_n250 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ################################################################################
      ################################################################################
      ################################################################################
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rt(n * p, df = df_t) / sqrt(df_t / (df_t - 2)), n) %*% chol(Sigma)
      y <- X %*% Beta0 + rnorm(n, sd = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-4)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### DS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel <- DS_Fit$DS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### Split Knockoff
      nu_hat <- SplitLasso(normc(X), normc(y), D)$best_nu
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = nu_hat,
        sign = 'disabled'
      )
      SplitKnockoff_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_Sel <- SplitKnockoff_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_t3_result_all_n250[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_t3_result_all_n250[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_t3_result_all_n250, file = "StatKnock_ChainPiecewiseConstant_n250signal12_t3.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.11 (a): heavy-tailed errors (n=250, p=120)
set.seed(2025)
df_t <- 3
p = 120; n = 250
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- 0.2
ChainPiecewiseConstant_heavytailerr_result_all_n250 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ################################################################################
      ################################################################################
      ################################################################################
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      err <- rt(n, df = df_t) / sqrt(df_t / (df_t - 2))
      y <- X %*% Beta0 + err
      signal_ind = which(abs(D %*% Beta0) > 1e-4)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### DS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel <- DS_Fit$DS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### Split Knockoff
      nu_hat <- SplitLasso(normc(X), normc(y), D)$best_nu
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = nu_hat,
        sign = 'disabled'
      )
      SplitKnockoff_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_Sel <- SplitKnockoff_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_heavytailerr_result_all_n250[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_heavytailerr_result_all_n250[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_heavytailerr_result_all_n250, file = "StatKnock_ChainPiecewiseConstant_n250signal12_heavytailerr.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.11 (b): skewed errors (n=250, p=120)
set.seed(2025)
p = 120; n = 250
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- 0.2
ChainPiecewiseConstant_skewerr_result_all_n250 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ################################################################################
      ################################################################################
      ################################################################################
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      err <- rexp(n, rate = 1) - 1
      y <- X %*% Beta0 + err
      signal_ind = which(abs(D %*% Beta0) > 1e-4)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### DS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel <- DS_Fit$DS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### Split Knockoff
      nu_hat <- SplitLasso(normc(X), normc(y), D)$best_nu
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = nu_hat,
        sign = 'disabled'
      )
      SplitKnockoff_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_Sel <- SplitKnockoff_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_skewerr_result_all_n250[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_skewerr_result_all_n250[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_skewerr_result_all_n250, file = "StatKnock_ChainPiecewiseConstant_n250signal12_skewerr.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.11 (c): heteroskedastic errors (n=250, p=120)
set.seed(2025)
p = 120; n = 250
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- 0.2
ChainPiecewiseConstant_heteroerr_result_all_n250 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ################################################################################
      ################################################################################
      ################################################################################
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      mu <- as.vector(X %*% Beta0)
      w <- 1 + abs(mu) / stats::sd(mu)
      w <- w / sqrt(mean(w^2))
      y <- matrix(mu + rnorm(n, sd = 1) * w, ncol = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-4)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### DS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel <- DS_Fit$DS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ################################################################################
      ################################################################################
      ################################################################################
      ### Split Knockoff
      nu_hat <- SplitLasso(normc(X), normc(y), D)$best_nu
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = nu_hat,
        sign = 'disabled'
      )
      SplitKnockoff_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_Sel <- SplitKnockoff_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_heteroerr_result_all_n250[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_heteroerr_result_all_n250[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_heteroerr_result_all_n250, file = "StatKnock_ChainPiecewiseConstant_n250signal12_heteroerr.RData")
