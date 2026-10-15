################################################################################
# TLasso: Transformational Lasso
# ------------------------------------------------------------------------------
# This script implements the TLasso method for structural sparsity recovery.
#
# Inputs:
#   X (n x p), y (n), D (m x p);
#   X_N_option in {1,2,3} for auxiliary design generation; proj in {TRUE,FALSE}.
#
# Outputs:
#   list(
#     gamma        = estimated Dbeta,
#     beta         = recovered beta on original p coordinates,
#     beta_N       = coefficients on auxiliary columns (if any),
#     best_lam     = CV-selected lambda,
#     X_star, y_star = transformed design/response used for Lasso
#   )
################################################################################
library(glmnet)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## TLasso
TLasso <- function(X, y, D, X_N_option = 1, proj = TRUE) {
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- nrow(D)
  qrD <- qr(D)
  r <- qrD$rank

  resid_proj <- function(A, U) {
    A - U %*% (MASS::ginv(crossprod(U)) %*% crossprod(U, A))
  }

  if (r == m){
    D_plus = D
    X_plus = X
  } else {
    Q_full <- qr.Q(qrD, complete = TRUE)
    N <- Q_full[, (r + 1):m, drop = FALSE]
    D_plus <- cbind(D, N)

    if (qr(X)$rank < p) stop("rank(X) < p")
    if (n < p + m - r) stop("n < p + m - r")
    tries <- 0
    repeat {
      if (X_N_option == 1) {
        X_N <- matrix(rnorm(n * (m - r)), nrow = n, ncol = m - r)
      } else if (X_N_option == 2) {
        X_N <- matrix(sample(c(-1, 1), n * (m - r), replace = TRUE), nrow = n, ncol = m - r)
      } else if (X_N_option == 3) {
        if ((m - r) > p) stop("Option 3 requires (m - r) <= p")
        X_N <- sapply(sample(ncol(X), m - r, replace = FALSE), function(j) sample(X[, j]))
      } else {
        stop("X_N_option must be 1, 2, or 3.")
      }

      if (proj) {
        if (n >= p) {
          X_N <- resid_proj(X_N, X)
        } else {
          message("Projection skipped: n < p")
        }
      }

      X_plus = cbind(X, X_N)
      if (qr(X_plus)$rank == p + m - r) break
      tries <- tries + 1
      if (tries >= 100) stop("resampling failed after 100 attempts")
    }
  }

  if (r == p){
    X_star = X_plus %*% solve(D_plus)
    y_star = y
  } else {
    E <- pracma::nullspace(D_plus)
    D_tilde <- rbind(D_plus, t(E))
    D_tilde_inverse <- solve(D_tilde)
    D_dagger <- D_tilde_inverse[, 1:m, drop = FALSE]
    D_star <- D_tilde_inverse[, (m + 1):(m + (p - r)), drop = FALSE]
    XDstar <- X_plus %*% D_star
    y_star <- resid_proj(y, XDstar)
    X_star <- resid_proj(X_plus %*% D_dagger, XDstar)
  }

  col_scale <- sqrt(colSums(X_star^2) / n)
  col_scale[col_scale < 1e-10] <- 1
  X_star_scaled <- sweep(X_star, 2, col_scale, "/")
  X_fit <- if (m == 1L) cbind(X_star_scaled, 0) else X_star_scaled
  cv_fit <- cv.glmnet(X_fit, y_star, family = "gaussian", alpha = 1, nfolds = 5, standardize = F, intercept = F)
  best_lambda <- cv_fit$lambda.min
  final_fit <- glmnet(X_fit, y_star, family = "gaussian", alpha = 1,
                      lambda = best_lambda, standardize = F, intercept = F,
                      thresh = 1e-14, maxit = 1e7)
  gamma <- as.vector(final_fit$beta[seq_len(m), , drop = FALSE]) / col_scale

  # recover beta
  if (r == p) {
    beta_plus <- solve(D_plus) %*% gamma
  } else {
    beta_partial <- D_dagger %*% gamma
    residual_y <- y - X_plus %*% beta_partial
    alpha <- MASS::ginv(t(X_plus %*% D_star) %*% (X_plus %*% D_star)) %*% t(X_plus %*% D_star) %*% residual_y
    beta_plus <- beta_partial + D_star %*% alpha
  }
  beta = beta_plus[1:p]
  beta_N = if (r == m) numeric(0) else beta_plus[(p+1):length(beta_plus)]

  return(list(gamma = gamma,
              beta = beta,
              beta_N = beta_N,
              best_lam = best_lambda,
              X_star = X_star,
              y_star = y_star))
}
