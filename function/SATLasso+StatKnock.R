################################################################################
# SATLasso + StatKnock
# ------------------------------------------------------------------------------
# This script implements the two-stage SATLasso and StatKnock procedure.
################################################################################
source("../function/StatKnock.R")
source("../function/SATLasso.R")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## SATLasso_StatKnock
SATLasso_StatKnock <- function(X, y, D,
                               L = 100, q = 0.1,
                               cv = TRUE, cvfolds = 5,
                               lambda_cr = c("min", "1se"),
                               knockoff_args = list(),
                               gbic = TRUE,
                               screening_fraction = 0.5) {
  lambda_cr <- match.arg(lambda_cr)
  y <- matrix(y, ncol = 1)
  n <- nrow(X); p <- ncol(X); m <- nrow(D)

  # If low-dimensional, skip screening
  if (n - (p - qr(D)$rank) >= 2 * m) {
    message("Low-dimensional setting: running StatKnock directly.")
    fit0 <- StatKnock(X, y, D, L = L, q = q,
                      cv = cv, cvfolds = cvfolds,
                      lambda_cr = lambda_cr,
                      knockoff_args = knockoff_args)
    return(list(S = fit0$S, S_plus = fit0$S_plus))
  }

  if (length(screening_fraction) != 1L || !is.finite(screening_fraction) ||
      screening_fraction <= 0 || screening_fraction >= 1) {
    stop("screening_fraction must be between 0 and 1")
  }

  ################################################################################
  # 1) Random split into (X1, y1) and (X2, y2)
  ################################################################################
  indices <- sample.int(n)
  n1 <- floor(n * screening_fraction); n2 <- n - n1
  idx1 <- indices[1:n1]
  idx2 <- indices[(n1 + 1):n]
  X1 <- X[idx1, , drop = FALSE]; y1 <- y[idx1, , drop = FALSE]
  X2 <- X[idx2, , drop = FALSE]; y2 <- y[idx2, , drop = FALSE]

  ################################################################################
  # 2) Screening on (X1, y1)
  ################################################################################
  if (gbic) {
    S_gamma <- SATLasso(X1, y1, D, stop = floor(n2 / log(n2)), gbic = TRUE)$em
    S_beta <- which(colSums(abs(D[S_gamma, , drop = FALSE]) != 0) > 0)
  } else {
    path <- SATLasso(X1, y1, D, stop = floor(n2 / log(n2)), gbic = FALSE)$es
    S_gamma <- integer(0); S_beta <- integer(0)
    for (g in path) {
      cand_g <- c(S_gamma, g)
      cand_b <- which(colSums(abs(D[cand_g, , drop = FALSE]) != 0) > 0)
      if (2 * length(cand_g) + length(cand_b) -
          qr(D[cand_g, cand_b, drop = FALSE])$rank > n2) break
      S_gamma <- cand_g; S_beta <- cand_b
    }
  }

  D_scr <- D[S_gamma, S_beta, drop = FALSE]

  ################################################################################
  # 3) Nuisance projection on (X2, y2)
  ################################################################################
  S_beta_c <- setdiff(seq_len(p), S_beta)
  edge <- seq_len(m)
  is_chain <- m == p - 1L &&
    all(D[cbind(edge, edge)] == -1) &&
    all(D[cbind(edge, edge + 1L)] == 1) &&
    all(rowSums(D != 0) == 2L)

  if (is_chain) {
    cuts <- sort(S_gamma)
    starts <- c(1L, cuts + 1L)
    ends <- c(cuts, p)
    Vcols <- vector("list", length(starts))
    for (k in seq_along(starts)) {
      ii <- intersect(starts[k]:ends[k], S_beta_c)
      if (length(ii)) Vcols[[k]] <- rowSums(X2[, ii, drop = FALSE])
    }
    Vcols <- Vcols[!vapply(Vcols, is.null, logical(1))]
    V_extra <- if (length(Vcols)) do.call(cbind, Vcols) else NULL
  } else {
    D_comp <- D[setdiff(seq_len(m), S_gamma), , drop = FALSE]
    basis <- if (nrow(D_comp)) pracma::nullspace(D_comp) else diag(p)
    if (is.null(basis) || !length(S_beta_c)) {
      V_extra <- NULL
    } else {
      basis <- as.matrix(basis)
      V_extra <- X2[, S_beta_c, drop = FALSE] %*%
        basis[S_beta_c, , drop = FALSE]
    }
  }

  dimV <- if (is.null(V_extra)) 0L else qr(V_extra)$rank
  if (dimV == 0L) V_extra <- NULL

  ################################################################################
  # 4) Selection on (X2, y2)
  ################################################################################
  m_scr <- length(S_gamma)
  if (m_scr < 2) {
    Sel <- integer(0)
    Sel_plus <- integer(0)
  } else {
    stat_fit <- tryCatch(
      StatKnock(X2[, S_beta, drop = FALSE],
                y2, D = D_scr, V_extra = V_extra,
                L = L, q = q,
                cv = cv, cvfolds = cvfolds,
                lambda_cr = lambda_cr,
                knockoff_args = knockoff_args),
      error = function(e) {
        if (gbic) stop(e)
        warning("StatKnock failed on the screened problem: ", conditionMessage(e))
        NULL
      })
    if (is.null(stat_fit)) {
      Sel <- integer(0)
      Sel_plus <- integer(0)
    } else {
      Sel <- S_gamma[stat_fit$S]
      Sel_plus <- S_gamma[stat_fit$S_plus]
    }
  }

  return(list(S = Sel, S_plus = Sel_plus))
}
