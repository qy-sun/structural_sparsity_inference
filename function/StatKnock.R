################################################################################
# StatKnock: Stability Transformational Knockoff
# ------------------------------------------------------------------------------
# This script implements the StatKnock method for selection in structural sparsity problems.
#
# Inputs:
#   X (n x p), y (n), D (m x p);
#   L (number of subsamplings), q (target mFDR/FDR for knockoff/knockoff+),
#   cv/cvfolds/lambda_cr (Lasso tuning), X_N_option/proj (augmentation choices),
#   knockoff_args (e.g., list(fixed_method = "equi")),
#   V_extra (optional nuisance-space basis).
#
# Outputs:
#   list(S       = selected set via knockoff (offset = 0; mFDR),
#        S_plus  = selected set via knockoff+ (offset = 1; FDR),
#        W       = knockoff statistics,
#        T1, T_plus = thresholds,
#        X_star, y_star = transformed design/response)
################################################################################
library(knockoff)
library(glmnet)
library(MASS)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## StatKnock
StatKnock <- function(X, y, D, L = 100, q = 0.1,
                      cv = TRUE, cvfolds = 5, lambda_cr = c("min","1se"),
                      X_N_option = 1, proj = TRUE,
                      knockoff_args = list(), V_extra = NULL){
  lambda_cr <- match.arg(lambda_cr)
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(D)[1]
  qrD <- qr(D)
  r <- qrD$rank

  resid_proj <- function(A, U) {
    A - U %*% (MASS::ginv(crossprod(U)) %*% crossprod(U, A))
  }

  if (r == m){
    D_plus = D
    X_plus = X
  } else {
    message("rank(D) < m, auxiliary noises introduced!")
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

  if (r == p && is.null(V_extra)){
    if (n < 2 * m) {
      stop(sprintf("Insufficient effective samples for fixed-X knockoffs: %d < %d.", n, 2 * m))
    }
    X_star = X_plus %*% solve(D_plus)
    y_star = y
  } else {
    if (r < p) {
      message("rank(D) < p, projection used!")
      E <- pracma::nullspace(D_plus)
      D_tilde <- rbind(D_plus, t(E))
      D_tilde_inverse <- solve(D_tilde)
      D_dagger <- D_tilde_inverse[, 1:m]
      D_star <- D_tilde_inverse[, (m + 1):(m + (p - r))]
      XDstar <- X_plus %*% D_star
    } else {
      D_dagger <- solve(D_plus)
      XDstar <- NULL
    }
    if (!is.null(V_extra)) {
      V_extra <- as.matrix(V_extra)
      if (nrow(V_extra) != n) stop("nrow(V_extra) must equal nrow(X)")
      XDstar <- cbind(XDstar, V_extra)
    }
    qrXD <- qr(XDstar, tol = 1e-10)
    rk <- qrXD$rank
    if (n - rk < 2 * m) {
      stop(sprintf("Insufficient effective samples for fixed-X knockoffs: %d < %d.", n - rk, 2 * m))
    }
    y_star <- qr.qty(qrXD, y)[(rk + 1):n, , drop = FALSE]
    X_star <- qr.qty(
      qrXD,
      X_plus %*% D_dagger
    )[(rk + 1):n, , drop = FALSE]
    if (!is.null(V_extra) && qr(X_star)$rank < m) {
      stop("rank(X_star) < m after nuisance projection")
    }
  }
  res <- generate_knockoff(
    X_star,
    method_args = knockoff_args
  )
  X_org <- res$X_org
  X_tilde <- res$X_tilde
  X_aug <- cbind(X_org, X_tilde)

  if (r != p || !is.null(V_extra)) {
    lift_QW <- function(Z) {
      Z_full <- matrix(0, n, ncol(Z))
      Z_full[(rk + 1):n, ] <- Z
      qr.qy(qrXD, Z_full)
    }
    X_aug <- lift_QW(X_aug)
    X_star <- lift_QW(X_star)
    y_star <- lift_QW(y_star)
  }

  V_pair <- X_aug[, seq_len(m), drop = FALSE] -
            X_aug[, m + seq_len(m), drop = FALSE]
  v_norm2 <- colSums(V_pair^2)
  if (any(!is.finite(v_norm2) | v_norm2 <= 0)) {
    stop("Response reflection requires every s_j > 0.")
  }
  R_set <- which(runif(m) < 0.5)
  if (length(R_set)) {
    dd <- V_pair[, R_set, drop = FALSE]
    y_star <- y_star - 2 * dd %*%
              (crossprod(dd, y_star) / v_norm2[R_set])
  }

  S_hat <- rep(0, 2*m)
  if(cv == TRUE){
    cv.fit = cv.glmnet(x = X_aug, y = y_star, family = "gaussian", intercept = F, alpha = 1, nfolds = cvfolds, standardize = F)
    if (lambda_cr == "min") {
      best_lambda = cv.fit$lambda.min
    } else {
      best_lambda = cv.fit$lambda.1se
    }
  } else {
    best_lambda = sqrt(log(2*m)/n) / sqrt(n)

  }

  for(i in 1:L){
    index <- sample((1:n), floor(n/2), replace = F)

    fit1 = glmnet(x = X_aug[index, ], y = y_star[index], family = "gaussian", intercept = F, alpha = 1, lambda = best_lambda, standardize = F,
                  thresh = 1e-14, maxit = 1e7)
    bhat1 <- as.vector(stats::coef(fit1)[-1])
    S_index1 <- which(bhat1 != 0)

    fit2 = glmnet(x = X_aug[-index, ],y = y_star[-index], family = "gaussian", intercept = F, alpha = 1, lambda = best_lambda, standardize = F,
                  thresh = 1e-14, maxit = 1e7)
    bhat2 <- as.vector(stats::coef(fit2)[-1])
    S_index2 <- which(bhat2 != 0)

    S_index <- intersect(S_index1, S_index2)

    S = rep(0, 2*m)
    S[S_index] = 1

    S_hat = S_hat + S
  }
  S_hat[c(R_set, R_set + m)] <- S_hat[c(R_set + m, R_set)]

  # compute Knockoff statistics
  W <- (S_hat[1:m] - S_hat[(m+1):(2*m)]) / L

  # find the knockoff threshold T
  t = sort(unique(abs(W)[abs(W) > 0]))
  ratio = numeric(length(t))

  for (j in seq_along(t)) {
    ratio[j] = (sum(W <= -t[j]))/max(1, sum(W >= t[j]))
  }
  id = which(ratio <= q)[1]
  if (is.na(id)) {
    T1 = Inf
  } else {
    T1 = t[id]
  }

  # find the knockoff_plus threshold T_plus
  t_plus = sort(unique(abs(W)[abs(W) > 0]))
  ratio_plus = numeric(length(t_plus))

  for (j in seq_along(t_plus)) {
    ratio_plus[j] = (1 + sum(W <= -t_plus[j]))/max(1, sum(W >= t_plus[j]))
  }
  id_plus = which(ratio_plus <= q)[1]
  if (is.na(id_plus)) {
    T_plus = Inf
  } else {
    T_plus = t_plus[id_plus]
  }

  # set of discovered variables
  S <- which(W >= T1)
  S_plus <- which(W >= T_plus)

  return(list(S = S, S_plus = S_plus, W = W, T1 = T1, T_plus = T_plus,
              X_star = X_star, y_star = y_star))
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions
haar_uperp <- function(X, m) {
  G <- matrix(rnorm(nrow(X) * m), nrow(X), m)
  G <- G - X %*% solve(crossprod(X), crossprod(X, G))
  qz <- qr(G)
  qr.Q(qz) %*% diag(sign(diag(qr.R(qz))), m)
}

create.fixed.haar <- function(X_org, method = "equi") {
  if (ncol(X_org) == 1L && nrow(X_org) > 1L && method %in% c("equi", "sdp")) {
    X <- X_org / sqrt(sum(X_org^2))
    return(list(X = X, Xk = haar_uperp(X, 1L)))
  }
  kn <- create.fixed(X_org, method = method)
  X <- kn$X
  m <- ncol(X)
  G <- crossprod(X)
  s <- diag(G - crossprod(X, kn$Xk))
  A <- X %*% (diag(m) - solve(G, diag(s, m)))
  R <- kn$Xk - A
  ev <- eigen(crossprod(R), symmetric = TRUE)
  ev$values[ev$values < 0] <- 0
  Cm <- diag(sqrt(ev$values), m) %*% t(ev$vectors)
  list(X = X, Xk = A + haar_uperp(X, m) %*% Cm)
}

generate_knockoff <- function(X_org, method_args = list()) {
  X_org <- scale(X_org, center = FALSE, scale = sqrt(colSums(X_org^2)))
  fixed_method <- if (!is.null(method_args$fixed_method)) {
    method_args$fixed_method
  } else {
    "equi"
  }
  knock <- create.fixed.haar(X_org, method = fixed_method)

  list(
    X_org = knock$X,
    X_tilde = knock$Xk
  )
}
