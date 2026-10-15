################################################################################
# High-dimensional GKnockoff and FuSIS
# ------------------------------------------------------------------------------
# This script extends two procedures proposed in Liu et al. (2024)
#   - FuSIS():  Fused Sure Independence Screening extended to linear 
#               combinations defined by D.
#   - HGKnockoff(): High-dim Generalized Knockoff combining SIS, FuSIS, 
#                   and GKnockoff for FDR-controlled selection.
#
# Dependencies: GKnockoff (version 0.1.0)
# 
# Reference: 
# Liu, J., Sun, A., & Ke, Y. (2024). A generalized knockoff procedure for FDR control 
# in structural change detection. Journal of Econometrics, 239(2) 105331.
################################################################################
library(GKnockoff)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## (1) FuSIS
FuSIS <- function(X, y, D, stop = dim(D)[1]){
  y <- matrix(y, ncol = 1)
  m <- dim(D)[1]
  
  corr_gamma <- cor(X %*% t(D), y)
  k <- min(stop, m)
  fusis_sel <- order(abs(corr_gamma), decreasing = TRUE)[1:k]
  
  return(fusis_sel)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## (2) HGKnockoff
HGKnockoff <- function(X, y, D, q = 0.1){
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(D)[1]
  if (n >= 2 * m) {
    r <- qr(D)$rank
    if (r == m) {
      message("Low-dimensional setting: running GKnockoff directly.")
      return(GKnockoff(X, y, D, fdr = q))
    }
    if (r < m) {
      message("Low-dimensional setting and not full row rank.")
      return(list(selected = c()))
    }
  }
  
  # Data Splitting
  indices <- sample(1:n)
  n1 <- floor(n/2)
  n2 <- n - n1
  index1 <- indices[1:n1]
  index2 <- indices[(n1 + 1):n]
  
  # 1) Screen beta&gamma: SIS then FuSIS
  corr_beta <- cor(X[index1, ],y[index1])
  k_beta <- min(p, floor(n2/2))
  S_beta <- order(abs(corr_beta), decreasing = TRUE)[1:k_beta]
  D_sub <- D[, S_beta, drop = FALSE]
  keep <- which(rowSums(abs(D_sub)) > 0)
  corr_gamma <- cor(X[index1, S_beta] %*% t(D_sub[keep, , drop = FALSE]), y[index1])
  k_gamma <- min(length(keep), floor(n2/2))
  S_gamma <- keep[order(abs(corr_gamma), decreasing = TRUE)[seq_len(k_gamma)]]
  
  # 2) Selection: GKnockoff
  D_scr = D[S_gamma,S_beta, drop = F]
  
  # Ensure full row rank
  if (qr(D_scr)$rank < nrow(D_scr)) {
    for (g in rev(S_gamma)) {
      pos <- which(S_gamma == g)
      D_scr <- D_scr[-pos, , drop = FALSE]
      S_gamma <- S_gamma[-pos]
      
      if (qr(D_scr)$rank == nrow(D_scr)) {
        break
      }
    }
  }
  
  HGKnockoff_Fit <- tryCatch({
    result <- GKnockoff(X[index2, S_beta], y[index2], D = D_scr, fdr = q)
    result
  }, error = function(e) {
    message("Error in GKnockoff: ", conditionMessage(e))
    list(selected = NULL)
  })
  HGKnockoff_Sel <- S_gamma[HGKnockoff_Fit$selected]
  
  return(list(selected = HGKnockoff_Sel))
}
