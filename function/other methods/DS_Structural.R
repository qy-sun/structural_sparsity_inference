################################################################################
# Data Splitting (DS) for Structural Sparsity
# ------------------------------------------------------------------------------
# This script extends the DS procedure to structural sparsity problems.
#
# Reference:
# Dai, C., Lin, B., Xing, X., et al. (2023).
# False discovery rate control via data splitting.
# Journal of the American Statistical Association, 118(544): 2503-2520.
################################################################################
library(glmnet)
library(genlasso)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## DS_D
DS_D <- function(X, y, D, num_split = 50, q = 0.1){
  n <- dim(X)[1]; p <- dim(X)[2]
  m <- dim(D)[1]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = m)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)

  for(iter in 1:num_split){
    # randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)

    # run GenLasso on the first half of the data
    lambda_GenLasso = sqrt(length(sample_index1)*log(m))
    GenLasso_fit <- genlasso(y[sample_index1], X[sample_index1, ], D, minlam = lambda_GenLasso)
    beta1 <- coef(GenLasso_fit, lambda = lambda_GenLasso)$beta
    gamma1 <- D %*% beta1
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index) != 0){
      # run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      gamma2 <- D %*% beta2

      # calculate the mirror statistics
      M <- sign(gamma1 * gamma2) * (abs(gamma1) + abs(gamma2))

      selected_index = analys(M, abs(M), q)
      DS_selected_index = selected_index

      # number of selected variables
      if(length(selected_index) != 0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
      }
    }else{
      DS_selected_index = NULL
    }
  }

  # multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)

  # rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank) != 0){
    null_feature <- numeric()

    # backtracking
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    MDS_selected_index <- setdiff(feature_rank, null_feature)
  }else{
    MDS_selected_index = NULL
  }
  return(list(DS_feature = DS_selected_index, MDS_feature = MDS_selected_index))
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions
analys <- function(mm, ww, q){
  # mm: mirror statistics
  # ww: absolute value of mirror statistics
  # q:  FDR control level
  cutoff_set <- max(ww)
  for(t in ww){
    ps <- length(mm[mm > t])
    ng <- length(na.omit(mm[mm < -t]))
    rto <- (ng + 1)/max(ps, 1)
    if(rto <= q){
      cutoff_set <- c(cutoff_set, t)
    }
  }
  cutoff <- min(cutoff_set)
  selected_index <- which(mm > cutoff)

  return(selected_index)
}
