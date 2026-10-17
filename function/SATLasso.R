################################################################################
# SATLasso: Sequential Adjusted Transformational Lasso
# ------------------------------------------------------------------------------
# This script implements the SATLasso method for screening in structural sparsity problems.
#
# Reference:
#   Luo, S., & Chen, Z. (2014).
#   Sequential Lasso cum EBIC for feature selection with ultra-high dimensional
#   feature space. Journal of the American Statistical Association, 109(507), 1229-1240.
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
SATLasso <- function(X,y,D, alpha_grid = 10^seq(-1, 1, length.out = 10),
                     stop = NULL, signal_ind = NULL, eta0 = 2.1, gbic = TRUE){
  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(D)
  qrD <- qr(D)
  r <- qrD$rank
  n_eff <- n - p + r

  if (r == m){
    final_res <- STLasso(X, y, D, stop = stop, signal_ind = signal_ind, eta0 = eta0, gbic = gbic)
    final_res$alpha <- NA_real_
    return(final_res)
  }

  y <- matrix(y, ncol = 1)

  if (!gbic) alpha_grid <- alpha_grid[ceiling(length(alpha_grid) / 2)]

  alpha_score <- numeric(length(alpha_grid))
  fits <- vector("list", length(alpha_grid))

  Q_full <- qr.Q(qrD, complete = TRUE)
  Nmat <- Q_full[, (r + 1):m, drop = FALSE]
  D_plus <- cbind(D, Nmat)
  
  if (r == p) {
    D_plus_inv <- solve(D_plus)
  } else {
    E <- pracma::nullspace(D_plus)
    D_tilde <- rbind(D_plus, t(E))
    D_inv <- solve(D_tilde)
    
    D_dagger <- D_inv[, 1:m]
    D_star <- D_inv[, (m + 1):(m + (p - r))]
  }
  
  X_N <- matrix(0, n, (m - r))
  y_aug <- rbind(1/sqrt(2*n) * y, matrix(0, m - r, 1))
  X_aug_obs <- cbind(1/sqrt(2*n) * X, 1/sqrt(2*n) * X_N)

  cb <- NULL
  if (!is.null(signal_ind)) {
    cb <- function(path) {
      length(intersect(path, signal_ind)) == length(signal_ind)
    }
  }

  for (i in seq_along(alpha_grid)) {
    alpha <- alpha_grid[i]
    X_aug <- rbind(X_aug_obs, cbind(matrix(0, m - r, p), sqrt(alpha) * diag(m - r)))

    if (r == p) {
      X_star <- X_aug %*% D_plus_inv
      y_star <- y_aug
    } else {
      XDstar <- X_aug %*% D_star
      Mmat <- diag(n + m - r) -
        XDstar %*% MASS::ginv(t(XDstar) %*% XDstar) %*% t(XDstar)
      y_star <- Mmat %*% y_aug
      X_star <- Mmat %*% X_aug %*% D_dagger
    }

    fits[[i]] <- SLasso(X_star, y_star, stop = stop, callback = cb,
                             n_eff = n_eff, eta0 = eta0, gbic = gbic)
    if (gbic) alpha_score[i] <- refit_gbic(fits[[i]]$em, X, y, D, eta0)
  }

  best_alpha_idx <- which.min(alpha_score)
  best_alpha <- alpha_grid[best_alpha_idx]
  print(paste("best alpha =", best_alpha))

  final_res <- fits[[best_alpha_idx]]
  final_res$alpha <- best_alpha
  if (!gbic) {
    final_res$em <- NA_integer_
    final_res$bic <- NULL
  }
  return(final_res)
}
################################################################################
################################################################################
################################################################################
### Helper functions
## (1) STLasso
STLasso <- function(X,y,D, X_N_option = 1, proj = TRUE, stop = NULL, signal_ind = NULL, eta0 = 2.1, gbic = TRUE){
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- nrow(D)
  qrD <- qr(D)
  r <- qrD$rank
  n_eff <- n - p + r
  
  if (r == m){
    D_plus = D
    X_plus = X
  } else {
    Q_full <- qr.Q(qrD, complete = TRUE)
    N <- Q_full[, (r + 1):m, drop = FALSE]
    D_plus <- cbind(D, N)
    
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
        P_X <- diag(n) - X %*% MASS::ginv(t(X) %*% X) %*% t(X)
        X_N <- P_X %*% X_N
      } else {
        message("Projection skipped: n < p")
      }
    }
    
    X_plus = cbind(X, X_N) 
  }
  
  if (r == p){
    X_star = X_plus %*% solve(D_plus)
    y_star = y
  } else {
    E <- pracma::nullspace(D_plus)
    D_tilde <- rbind(D_plus, t(E))
    D_tilde_inverse <- solve(D_tilde)
    D_dagger <- D_tilde_inverse[, 1:m]
    D_star <- D_tilde_inverse[, (m + 1):(m + (p - r))]
    XDstar <- X_plus %*% D_star
    M_XDstar <- diag(n) - XDstar %*% MASS::ginv(t(XDstar) %*% XDstar) %*% t(XDstar)
    y_star <- M_XDstar %*% y
    X_star <- M_XDstar %*% X_plus %*% D_dagger
  }
  
  cb <- NULL
  if (!is.null(signal_ind)) {
    cb <- function(path) {
      length(intersect(path, signal_ind)) == length(signal_ind)
    }
  }
  
  final_res <- SLasso(X_star, y_star, stop = stop, callback = cb, n_eff = n_eff, eta0 = eta0, gbic = gbic)
  return(final_res)
}
################################################################################
################################################################################
################################################################################
## (2) SATLasso_alpha (fix alpha)
SATLasso_alpha <- function(X,y,D, alpha = 0.5, stop = NULL, signal_ind = NULL, eta0 = 2.1, gbic = TRUE){
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- nrow(D)
  qrD <- qr(D)
  r <- qrD$rank
  n_eff <- n - p + r
  
  if (r == m){
    final_res <- STLasso(X, y, D, stop = stop, signal_ind = signal_ind, eta0 = eta0, gbic = gbic)
    return(final_res)
  }
  
  Q_full <- qr.Q(qrD, complete = TRUE)
  N <- Q_full[, (r + 1):m, drop = FALSE]
  D_plus <- cbind(D, N)
  
  X_N <- matrix(0, n, (m-r))
  
  y_aug <- rbind(1/sqrt(2*n) * y, matrix(0, m-r, 1))
  X_aug <- rbind(
    cbind(1/sqrt(2*n) * X,     1/sqrt(2*n) * X_N),
    cbind(matrix(0, m-r, p),  sqrt(alpha) * diag(m-r))
  )
  
  if (r == p) {
    X_star <- X_aug %*% solve(D_plus)
    y_star <- y_aug
  } else {
    E <- pracma::nullspace(D_plus)
    D_tilde <- rbind(D_plus, t(E))
    D_inv <- solve(D_tilde)
    
    D_dagger <- D_inv[, 1:m]
    D_star <- D_inv[, (m + 1):(m + (p - r))]
    
    XDstar <- X_aug %*% D_star
    M <- diag(n+m-r) - XDstar %*% MASS::ginv(t(XDstar) %*% XDstar) %*% t(XDstar)
    
    y_star <- M %*% y_aug
    X_star <- M %*% X_aug %*% D_dagger
  }
  
  cb <- NULL
  if (!is.null(signal_ind)) {
    cb <- function(path) {
      length(intersect(path, signal_ind)) == length(signal_ind)
    }
  }
  
  final_res <- SLasso(X_star, y_star, stop = stop, callback = cb, n_eff = n_eff, eta0 = eta0, gbic = gbic)
  return(final_res)
}
################################################################################
################################################################################
################################################################################
## (3) SLasso (Luo & Chen, 2014)
SLasso <- function(X,y, stop = NULL, callback = NULL, n_eff = NULL, eta0 = 2.1, gbic = TRUE, tol = 1e-10){
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(n_eff)) n_eff <- n
  K_max <- if (is.null(stop)) p else stop
  if (!is.numeric(X) || length(y) != n || any(!is.finite(X)) || any(!is.finite(y)))
    stop("SLasso requires a finite numeric design and matching response")
  if (length(K_max) != 1L || !is.finite(K_max) || K_max < 0 || K_max != floor(K_max))
    stop("SLasso path limit must be a nonnegative finite integer")
  if (length(tol) != 1L || !is.finite(tol) || tol < 0 || tol >= 1)
    stop("SLasso tol must be in [0, 1)")
  K_max <- min(K_max, p)
  
  ## Initial model
  Xcheck <- X
  Ycheck <- y
  original_norm <- sqrt(colSums(X^2))
  if (any(!is.finite(original_norm))) stop("Non-finite SLasso column norm")
  active <- original_norm > 0
  path <- integer(K_max)
  RSS <- numeric(K_max)

  k <- 0L
  for (i in seq_len(K_max)){
    ## Sequential Select
    sel <- sequential_select(Xcheck, Ycheck, active, original_norm, tol = tol)
    if (is.null(sel)) break
    index.m <- sel$index.m
    active <- sel$active

    ## Update Projector
    u <- Xcheck[, index.m, drop = FALSE] / sel$col_norm
    Xcheck <- Xcheck - u %*% crossprod(u, Xcheck)
    Ycheck <- as.numeric(Ycheck - u %*% crossprod(u, Ycheck))
    next_rss <- sum(Ycheck^2)
    if (!is.finite(next_rss) || any(!is.finite(Xcheck))) stop("Non-finite SLasso residual update")
    active[index.m] <- FALSE
    k <- k+1L
    path[k] <- index.m
    RSS[k] <- next_rss

    ## Early Stop
    if (!is.null(callback) && callback(path[seq_len(k)])) break
  }

  path <- path[seq_len(k)]
  RSS <- RSS[seq_len(k)]
  EBIC <- NULL
  em <- NA_integer_
  if (gbic) {
    if (length(n_eff) != 1L || !is.finite(n_eff) || n_eff <= 1)
      stop("GBIC requires explicit n_eff > 1 (n - p + rank(D) for SATLasso)")
    eta <- if (p > 1L) eta0 -log(n_eff)/log(p) else eta0
    if (length(eta) != 1L || !is.finite(eta)) stop("GBIC requires finite eta")
    sz <- seq_along(path)
    EBIC <- n_eff*log(RSS[sz]/n_eff) + sz*log(n_eff) + eta*lchoose(p,sz)
    em <- if (length(sz)) path[seq_len(sz[which.min(EBIC)])] else integer(0)
  }

  return(list(es = path,
              em = em,
              bic = EBIC,
              rss = RSS))
}

sequential_select <- function(X, y, active, original_norm, tol = 1e-10){
  original_indices <- which(active)
  if (!length(original_indices)) return(NULL)
  X_new <- X[, original_indices, drop = FALSE]
  col_norm <- sqrt(colSums(X_new^2))
  if (any(!is.finite(col_norm))) stop("Non-finite SLasso residual norm")
  valid <- col_norm > tol * original_norm[original_indices] & col_norm > 0
  active[original_indices[!valid]] <- FALSE
  if (!any(valid)) return(NULL)
  original_indices <- original_indices[valid]
  col_norm <- col_norm[valid]
  score <- as.numeric(abs(crossprod(X_new[, valid, drop = FALSE], y)) / col_norm)
  if (any(!is.finite(score))) stop("Non-finite SLasso selection score")
  if (max(score) == 0) return(NULL)
  best <- which.max(score)
  return(list(index.m = original_indices[best],
              col_norm = col_norm[best],
              active = active))
}

################################################################################
################################################################################
################################################################################
## (4) refit_gbic: GBIC of the constrained least squares refit on the original scale
##
## The augmented residual equals the data fit plus alpha*||t(N)%*%gamma||^2, a penalty
## that grows with alpha and shrinks the deviation of gamma from col(D), so augmented
## GBIC scores rank models within a single alpha. The refit uses the n original
## observations and c_S parameters, while the combinatorial penalty uses eta0.
refit_gbic <- function(S, X, y, D, eta0) {
  m <- nrow(D)
  n <- nrow(X)
  Sc <- setdiff(1:m, S)
  Z <- if (length(Sc) == 0) diag(ncol(X)) else pracma::nullspace(D[Sc, , drop = FALSE])
  if (is.null(Z) || ncol(Z) == 0) {
    rss <- sum(y^2)
    c_S <- 0
  } else if (ncol(Z) >= n) {
    return(Inf)
  } else {
    beta_hat <- Z %*% (MASS::ginv(X %*% Z) %*% y)
    rss <- sum((y - X %*% beta_hat)^2)
    c_S <- ncol(Z)
  }
  n*log(rss/n) + c_S*log(n) + eta0*lchoose(m, length(S))
}
