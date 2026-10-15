################################################################################
# Simulation Experiments for SATLasso+StatKnock
# ------------------------------------------------------------------------------
# This script reproduces the results in Figure 5 in the main text
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
source("../function/SATLasso+StatKnock.R")
source("../function/other methods/HD_SplitKnockoff.R")
source("../function/other methods/HD_GKnockoff.R")
source("../function/other methods/BY_Structural.R")
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
## Figure 5 (a): piecewise constant over chain (n=500, m=1000)

set.seed(2025)
p = 1001; n = 500
m1 = 5
D <- create.chain(p)
m <- nrow(D)
A = 1
rho = 0.2
Sigma = toeplitz(rho^(1:p - 1))

p1_values <- c(20, 40, 80, 160)
ChainPiecewiseConstant_list <- list()

for (p1 in p1_values) {
  cat("  Running for p1 =", p1, "\n")
    
  fdp_power_list_ChainPiecewiseConstant_m1000 = sapply(1:nIterations, function(it) {
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
    ### SATLasso+StatKnock
    SATLasso_StatKnock_Fit <- SATLasso_StatKnock(X = X, y = y, D = D, q = q)
    SATLasso_StatKnock_Sel <- SATLasso_StatKnock_Fit$S
    SATLasso_StatKnock_plus_Sel <- SATLasso_StatKnock_Fit$S_plus
    fdp_power1 = fdp_power(selected_index = SATLasso_StatKnock_Sel, signal_index = signal_ind)
    fdp_power2 = fdp_power(selected_index = SATLasso_StatKnock_plus_Sel, signal_index = signal_ind)
      
    ################################################################################
    ################################################################################
    ################################################################################
    ### B-Y
    BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
    fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      
    ################################################################################
    ################################################################################
    ################################################################################
    ### GKnockoff
    GKnockoff_Fit = HGKnockoff(X, y, D = D, q = q)
    GKnockoff_Sel <- GKnockoff_Fit$selected
    fdp_power4 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
    fdp_power4
      
    ################################################################################
    ################################################################################
    ################################################################################
    ### Split Knockoff
    SplitKnockoff_Fit = HD_SplitKnockoff(X, y, D = D, q = q)
    SplitKnockoff_Sel <- SplitKnockoff_Fit$Sel
    fdp_power5 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
    fdp_power5
    
    r = c(fdp_power1$fdp, fdp_power1$power,
          fdp_power2$fdp, fdp_power2$power,
          fdp_power3$fdp, fdp_power3$power,
          fdp_power4$fdp, fdp_power4$power,
          fdp_power5$fdp, fdp_power5$power
    )
    return(r)
  })
    
  ChainPiecewiseConstant_list[[paste0("p1=", p1)]] <- fdp_power_list_ChainPiecewiseConstant_m1000
}

for (p1 in p1_values) {
  result_matrix <- ChainPiecewiseConstant_list[[paste0("p1=", p1)]]
  cat("p1 =", p1, "\n")
  print(apply(result_matrix, 1, mean.se))
}

save(ChainPiecewiseConstant_list, file = "SATLasso_StatKnock_ChainPiecewiseConstant_m1000.RData")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure 5 (b): piecewise constant over chain (n=500, m=5000)

set.seed(2025)
p = 5001; n = 500
m1 = 5
D <- create.chain(p)
m <- nrow(D)
A = 1
rho = 0.2
Sigma = toeplitz(rho^(1:p - 1))

p1_values <- c(20, 40, 80, 160)
ChainPiecewiseConstant_list <- list()

for (p1 in p1_values) {
  cat("  Running for p1 =", p1, "\n")
  
  fdp_power_list_ChainPiecewiseConstant_m5000 = sapply(1:nIterations, function(it) {
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
    ### SATLasso+StatKnock
    SATLasso_StatKnock_Fit <- SATLasso_StatKnock(X = X, y = y, D = D, q = q)
    SATLasso_StatKnock_Sel <- SATLasso_StatKnock_Fit$S
    SATLasso_StatKnock_plus_Sel <- SATLasso_StatKnock_Fit$S_plus
    fdp_power1 = fdp_power(selected_index = SATLasso_StatKnock_Sel, signal_index = signal_ind)
    fdp_power2 = fdp_power(selected_index = SATLasso_StatKnock_plus_Sel, signal_index = signal_ind)
    
    ################################################################################
    ################################################################################
    ################################################################################
    ### B-Y
    BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
    fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
    
    ################################################################################
    ################################################################################
    ################################################################################
    ### GKnockoff
    GKnockoff_Fit = HGKnockoff(X, y, D = D, q = q)
    GKnockoff_Sel <- GKnockoff_Fit$selected
    fdp_power4 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
    fdp_power4
    
    ################################################################################
    ################################################################################
    ################################################################################
    ### Split Knockoff
    SplitKnockoff_Fit = HD_SplitKnockoff(X, y, D = D, q = q)
    SplitKnockoff_Sel <- SplitKnockoff_Fit$Sel
    fdp_power5 = fdp_power(selected_index = SplitKnockoff_Sel, signal_index = signal_ind)
    fdp_power5
    
    r = c(fdp_power1$fdp, fdp_power1$power,
          fdp_power2$fdp, fdp_power2$power,
          fdp_power3$fdp, fdp_power3$power,
          fdp_power4$fdp, fdp_power4$power,
          fdp_power5$fdp, fdp_power5$power
    )
    return(r)
  })
  
  ChainPiecewiseConstant_list[[paste0("p1=", p1)]] <- fdp_power_list_ChainPiecewiseConstant_m5000
}

for (p1 in p1_values) {
  result_matrix <- ChainPiecewiseConstant_list[[paste0("p1=", p1)]]
  cat("p1 =", p1, "\n")
  print(apply(result_matrix, 1, mean.se))
}

save(ChainPiecewiseConstant_list, file = "SATLasso_StatKnock_ChainPiecewiseConstant_m5000.RData")
