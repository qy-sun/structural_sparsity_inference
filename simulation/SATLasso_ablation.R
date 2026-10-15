################################################################################
# Simulation Ablations for SATLasso
# ------------------------------------------------------------------------------
# This script reproduces Figures S.2 and S.3 in the Supplement
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
library(glmnet)
library(MASS)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions
fdp_power <- function(selected_index, signal_index) {
  fdp <- 0
  power <- 0
  ssr <- 0
  if (length(selected_index) != 0) {
    num_selected <- length(selected_index)
    tp <- length(intersect(selected_index, signal_index))
    fp <- num_selected - tp
    fdp <- fp / max(num_selected, 1)
    power <- tp / length(signal_index)
    ssr <- ifelse(power == 1, 1, 0)
  }
  return(list(fdp = fdp, power = power, ssr = ssr))
}

precompute_D <- function(D) {
  m <- nrow(D)
  p <- ncol(D)
  qrD <- qr(D)
  r <- qrD$rank

  if (r == m) {
    N <- NULL
    D_plus <- D
  } else {
    Q_full <- qr.Q(qrD, complete = TRUE)
    N <- Q_full[, (r + 1):m, drop = FALSE]
    D_plus <- cbind(D, N)
  }

  if (r == p) {
    D_dagger <- solve(D_plus)
    D_star <- NULL
  } else {
    QN <- qr.Q(qr(t(D_plus)), complete = TRUE)
    E <- QN[, (m + 1):(p + m - r), drop = FALSE]
    D_tilde <- rbind(D_plus, t(E))
    D_inv <- solve(D_tilde)
    D_dagger <- D_inv[, 1:m, drop = FALSE]
    D_star <- D_inv[, (m + 1):(m + (p - r)), drop = FALSE]
  }

  D_dagger_1 <- D_dagger[1:p, , drop = FALSE]
  D_dagger_2 <- if (r == m) NULL else D_dagger[(p + 1):(p + m - r), , drop = FALSE]
  D_star_1 <- if (r == p) NULL else D_star[1:p, , drop = FALSE]
  D_star_2 <- if (r == p || r == m) NULL else D_star[(p + 1):(p + m - r), , drop = FALSE]

  return(list(p = p, m = m, r = r,
              D_dagger_1 = D_dagger_1, D_dagger_2 = D_dagger_2,
              D_star_1 = D_star_1, D_star_2 = D_star_2))
}

ATLasso_pre <- function(X, D_info) {
  list(XD1 = X %*% D_info$D_dagger_1,
       XS1 = if (is.null(D_info$D_star_1)) NULL else X %*% D_info$D_star_1)
}

ATLasso_star <- function(X, y, D_info, alpha, pre = NULL) {
  y <- as.numeric(y)
  n <- nrow(X)
  p <- D_info$p
  m <- D_info$m
  r <- D_info$r
  mr <- m - r
  c0 <- 1 / sqrt(2 * n)
  if (is.null(pre)) pre <- ATLasso_pre(X, D_info)

  y_aug <- c(c0 * y, rep(0, mr))

  XD <- rbind(c0 * pre$XD1, sqrt(alpha) * D_info$D_dagger_2)

  if (D_info$r == p) {
    return(list(X_star = XD, y_star = y_aug))
  }

  U <- rbind(c0 * pre$XS1, sqrt(alpha) * D_info$D_star_2)
  UtU_inv <- MASS::ginv(crossprod(U))
  X_star <- XD - U %*% (UtU_inv %*% crossprod(U, XD))
  y_star <- as.numeric(y_aug - U %*% (UtU_inv %*% crossprod(U, y_aug)))

  return(list(X_star = X_star, y_star = y_star))
}

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

  # Initial model
  Xcheck <- X
  Ycheck <- y
  original_norm <- sqrt(colSums(X^2))
  if (any(!is.finite(original_norm))) stop("Non-finite SLasso column norm")
  active <- original_norm > 0
  path <- integer(K_max)
  RSS <- numeric(K_max)

  k <- 0L
  for (i in seq_len(K_max)){
    # Sequential Select
    sel <- sequential_select(Xcheck, Ycheck, active, original_norm, tol = tol)
    if (is.null(sel)) break
    index.m <- sel$index.m
    active <- sel$active

    # Update Projector
    u <- Xcheck[, index.m, drop = FALSE] / sel$col_norm
    Xcheck <- Xcheck - u %*% crossprod(u, Xcheck)
    Ycheck <- as.numeric(Ycheck - u %*% crossprod(u, Ycheck))
    next_rss <- sum(Ycheck^2)
    if (!is.finite(next_rss) || any(!is.finite(Xcheck))) stop("Non-finite SLasso residual update")
    active[index.m] <- FALSE
    k <- k+1L
    path[k] <- index.m
    RSS[k] <- next_rss

    # Early Stop
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

  return(list(es=path,
              em=em,
              bic=EBIC,
              rss=RSS))
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
mms_from_path <- function(fit_path, signal_ind) {
  hit <- match(signal_ind, fit_path$es)
  if (any(is.na(hit))) {
    return(list(mms = NA_integer_, censored = 1L))
  }
  return(list(mms = as.integer(max(hit)), censored = 0L))
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Simulation settings
nIterations = 100

p <- 1000
n <- 200
p1 <- 20
m1 <- 5
A <- 1
rho_x <- 0.2
Sigma <- toeplitz(rho_x^(0:(p - 1)))

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.2: ATLasso sensitivity to alpha and lambda
D <- rbind(diag(p), create.chain(p))
m <- nrow(D)
D_info <- precompute_D(D)
r <- D_info$r

alpha_grid <- 10^seq(-1, 4, length.out = 26)
nlam       <- 40
rho_grid   <- 10^seq(0, -3, length.out = nlam)

Sens_AL <- list()

for (iter in 1:nIterations) {
  set.seed(iter)
  Beta0 <- Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
  X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)

  gamma_star <- as.numeric(D %*% Beta0)
  sig_gamma  <- which(abs(gamma_star) > 1e-10)
  norm_gs    <- sqrt(sum(gamma_star^2))
  pre <- ATLasso_pre(X, D_info)

  for (ia in seq_along(alpha_grid)) {
    alpha <- alpha_grid[ia]
    st <- ATLasso_star(X, y, D_info, alpha = alpha, pre = pre)
    Xs <- st$X_star
    ys <- st$y_star

    fit0    <- glmnet(Xs, ys, family = "gaussian", alpha = 1,
                      standardize = TRUE, intercept = FALSE, nlambda = 3)
    lam_max <- max(fit0$lambda)
    lam_grid <- lam_max * rho_grid

    fit <- glmnet(Xs, ys, family = "gaussian", alpha = 1,
                  lambda = lam_grid, standardize = TRUE, intercept = FALSE,
                  maxit = 1e6)
    G <- as.matrix(fit$beta)
    lam_used <- fit$lambda
    nL <- ncol(G)

    for (jl in seq_len(nL)) {
      gh   <- G[, jl]
      diff <- gh - gamma_star
      sel  <- which(abs(gh) > 1e-10)
      tp   <- length(intersect(sel, sig_gamma))
      Sens_AL[[length(Sens_AL) + 1]] <- data.frame(
        iter = iter,
        alpha = alpha, lambda = lam_used[jl], rho = lam_used[jl] / lam_max,
        rel_err = sqrt(sum(diff^2)) / norm_gs,
        abs_err = sqrt(sum(diff^2)),
        nnz = length(sel),
        TPR = tp / length(sig_gamma),
        FDR_supp = (length(sel) - tp) / max(length(sel), 1),
        n = n, m = m, r = r, n_signal = length(sig_gamma),
        stringsAsFactors = FALSE
      )
    }
  }
}

Sens_AL <- do.call(rbind, Sens_AL)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.3: SATLasso sensitivity to eta
D_1 <- diag(p)
D_2 <- create.chain(p)
D_3 <- rbind(diag(p), create.chain(p))
D_list <- list(D_1 = D_1, D_2 = D_2, D_3 = D_3)

alpha_fixed <- 0.5

K_max <- floor(n / log(n))

use_n_eff <- "n-p+r"
Sens_eta <- list()

for (iter in 1:nIterations) {
  set.seed(iter)
  Beta0 <- Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
  X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)

  for (dd in seq_along(D_list)) {
    Dname <- names(D_list)[dd]
    D <- D_list[[dd]]
    m <- nrow(D)
    signal_ind <- which(abs(D %*% Beta0) > 1e-10)

    D_info <- precompute_D(D)
    r <- D_info$r
    n_star <- n + m - r
    n_eff  <- n - p + r

    alpha_used <- alpha_fixed

    eta_ref     <- 2 - log(n_eff) / log(m)
    eta_default <- 2.1 - log(n_eff) / log(m)
    eta_common  <- c(0, 0.5, 0.75, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0)
    eta_vals    <- sort(unique(c(eta_common, eta_default)))
    eta_tab <- data.frame(eta = eta_vals, stringsAsFactors = FALSE)
    eta_tab$is_bic     <- as.integer(abs(eta_tab$eta - 0)           < 1e-9)
    eta_tab$is_default <- as.integer(abs(eta_tab$eta - eta_default) < 1e-9)
    eta_tab$eta_label  <- ifelse(eta_tab$is_bic == 1,     "BIC (eta=0)",
                          ifelse(eta_tab$is_default == 1, "default (2.1 - log_m n_eff)",
                                 sprintf("eta=%.2f", eta_tab$eta)))

    pre <- ATLasso_pre(X, D_info)
    st  <- ATLasso_star(X, y, D_info, alpha = alpha_used, pre = pre)
    fit_path <- SLasso(st$X_star, st$y_star, stop = K_max, gbic = FALSE)

    mms_res <- mms_from_path(fit_path, signal_ind)
    screened_Kmax <- fit_path$es[seq_len(min(K_max, length(fit_path$es)))]
    ssr_Kmax <- as.integer(all(signal_ind %in% screened_Kmax))

    model_size <- seq_len(min(K_max, length(fit_path$es)))
    RSS <- fit_path$rss[model_size]
    GBIC_base <- n_eff * log(RSS / n_eff) + model_size * log(n_eff)
    GBIC_penalty <- lchoose(m, model_size)
    for (e in 1:nrow(eta_tab)) {
      EBIC <- GBIC_base + eta_tab$eta[e] * GBIC_penalty
      size_hat <- if (length(model_size)) model_size[which.min(EBIC)] else 0L
      sel <- list(selected = fit_path$es[seq_len(size_hat)], size_hat = size_hat)
      fp  <- fdp_power(selected_index = sel$selected, signal_index = signal_ind)
      Sens_eta[[length(Sens_eta) + 1]] <- data.frame(
        iter = iter, D = Dname,
        eta_label = eta_tab$eta_label[e], eta = eta_tab$eta[e],
        is_bic = eta_tab$is_bic[e], is_default_eta = eta_tab$is_default[e],
        size_hat = sel$size_hat, FDP = fp$fdp, PSR = fp$power, SSR = fp$ssr,
        MMS = mms_res$mms, MMS_censored = mms_res$censored, SSR_Kmax = ssr_Kmax,
        n = n, n_star = n_star, n_eff = n_eff, use_n_eff = use_n_eff,
        m = m, r = r, n_signal = length(signal_ind),
        eta_ref = eta_ref, eta_default = eta_default, K_max = K_max,
        alpha = alpha_used,
        stringsAsFactors = FALSE
      )
    }
  }
}

Sens_eta <- do.call(rbind, Sens_eta)

save(Sens_AL, Sens_eta, file = "SATLasso_ablation.RData")
