################################################################################
# High-dimensional SplitKnockoff and SplitLasso
# ------------------------------------------------------------------------------
# This script implements the high-dimensional SplitKnockoff and SplitLasso.
#
# Dependencies: SplitKnockoff (version 1.2)
#
# Reference:
#   Cao, Y., Sun, X., & Yao, Y. (2024).
#   Controlling the false discovery rate in transformational sparsity:
#   Split knockoffs.
#   Journal of the Royal Statistical Society: Series B (Statistical Methodology),
#   86(2), 386-410.
################################################################################
library(glmnet)
library(SplitKnockoff)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## (1) SplitLasso
SplitLasso <- function(X, y, D, nu_values = 10^seq(0, 2, by = 0.4),
                       lambda_values = 10^seq(0, -8, by = -0.4), nfolds = 5) {
  y <- as.vector(y)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- nrow(D)

  penalty <- matrix(1, p + m, 1)
  penalty[1:p, 1] <- 0

  test_size <- floor(n / nfolds)
  rand_rank <- sample(n)
  loss_sl <- array(Inf, dim = c(nfolds, length(nu_values), length(lambda_values)))

  for (k in seq_len(nfolds)) {
    test <- rand_rank[(test_size * (k - 1) + 1):(test_size * k)]
    X_train <- sk.center(X[-test, , drop = FALSE])
    y_train <- y[-test] - mean(y[-test])
    X_test <- sk.center(X[test, , drop = FALSE])
    y_test <- y[test] - mean(y[test])

    for (i in seq_along(nu_values)) {
      creat.result <- sk.create(X_train, y_train, D, nu_values[i])
      fit <- glmnet(cbind(creat.result$A_beta, creat.result$A_gamma), creat.result$tilde_y,
                    alpha = 1, lambda = lambda_values, penalty.factor = penalty,
                    standardize = FALSE, intercept = FALSE)
      coefs <- as.matrix(fit$beta)
      for (j in seq_len(ncol(coefs))) {
        y_sl <- X_test %*% coefs[1:p, j]
        loss_sl[k, i, j] <- sum((y_sl - y_test)^2) / test_size
      }
    }
  }

  mean_loss_sl <- apply(loss_sl, c(2, 3), mean)
  best_index <- which(mean_loss_sl == min(mean_loss_sl), arr.ind = TRUE)[1, ]
  best_nu <- nu_values[best_index[1]]
  best_lambda <- lambda_values[best_index[2]]
  print(paste("best nu =", best_nu))

  creat.result <- sk.create(X, y, D, best_nu)
  screen_step <- glmnet(cbind(creat.result$A_beta, creat.result$A_gamma), creat.result$tilde_y,
                        alpha = 1, lambda = best_lambda, penalty.factor = penalty,
                        standardize = FALSE, intercept = FALSE)
  coefs <- screen_step$beta
  betas <- coefs[1:p, ]
  gammas <- coefs[(p + 1):(p + m), ]

  S <- which(abs(gammas) > 1e-4)

  list(S = S, gammas = gammas, betas = betas, best_nu = best_nu, best_lam = best_lambda)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## (2) SplitLasso_screen
SplitLasso_screen <- function(X, y, D, n2, n_2, nu_values = 10^seq(2, 0, by = -0.4),
                              lambda_values = 10^seq(0, -8, by = -0.4), nfolds = 10,
                              screen_lambda = 10^seq(0, -6, by = -0.1)) {
  y <- as.vector(y)
  n <- dim(X)[1]

  cvfit <- cv.glmnet(X, y, lambda = screen_lambda)
  index <- which(cvfit$lambda == cvfit$lambda.1se)[1]
  beta_hat <- as.vector(cvfit$glmnet.fit$beta[, index])
  beta_supp <- which(abs(beta_hat) >= 1e-2)

  if (!length(beta_supp)) {
    return(list(beta_supp = integer(0), gamma_supp = integer(0),
                gammas = numeric(0), betas = numeric(0),
                best_nu = NA_real_, best_lam = NA_real_))
  }

  X <- X[, beta_supp, drop = FALSE]
  D <- D[, beta_supp, drop = FALSE]
  p <- ncol(D)
  m <- nrow(D)

  penalty <- matrix(1, p + m, 1)
  penalty[1:p, 1] <- 0

  test_size <- floor(n / nfolds)
  rand_rank <- sample(n)
  loss_sl <- array(Inf, dim = c(nfolds, length(nu_values), length(lambda_values)))

  for (k in seq_len(nfolds)) {
    test <- rand_rank[(test_size * (k - 1) + 1):(test_size * k)]
    X_train <- sk.center(X[-test, , drop = FALSE])
    y_train <- y[-test] - mean(y[-test])
    X_test <- sk.center(X[test, , drop = FALSE])
    y_test <- y[test] - mean(y[test])

    for (i in seq_along(nu_values)) {
      creat.result <- sk.create(X_train, y_train, D, nu_values[i])
      fit <- glmnet(cbind(creat.result$A_beta, creat.result$A_gamma), creat.result$tilde_y,
                    alpha = 1, lambda = lambda_values, penalty.factor = penalty,
                    standardize = FALSE, intercept = FALSE)
      coefs <- as.matrix(fit$beta)
      nonzeros <- colSums(coefs != 0)
      for (j in seq_len(ncol(coefs))) {
        y_sl <- X_test %*% coefs[1:p, j]
        loss_sl[k, i, j] <- sum((y_sl - y_test)^2) / test_size
        if (nonzeros[j] > n2) {
          loss_sl[k, i, j] <- Inf
        }
      }
    }
  }

  mean_loss_sl <- apply(loss_sl, c(2, 3), mean)
  sd_loss_sl <- apply(loss_sl, c(2, 3), sk.sd.pop)
  candidate <- which(mean_loss_sl <= min(mean_loss_sl + sd_loss_sl, na.rm = TRUE))
  last <- candidate[length(candidate)]
  best_nu <- nu_values[((last - 1) %% nrow(mean_loss_sl)) + 1]
  best_lambda <- lambda_values[((last - 1) %/% nrow(mean_loss_sl)) + 1]
  print(paste("best nu =", best_nu))

  creat.result <- sk.create(X, y, D, best_nu)
  screen_step <- glmnet(cbind(creat.result$A_beta, creat.result$A_gamma), creat.result$tilde_y,
                        alpha = 1, lambda = best_lambda, penalty.factor = penalty,
                        standardize = FALSE, intercept = FALSE)
  coefs <- screen_step$beta
  betas <- coefs[1:p, ]
  gammas <- coefs[(p + 1):(p + m), ]

  gamma_supp <- which(abs(gammas) > 1e-4)
  gamma_room <- n_2 - length(beta_supp)
  if (gamma_room <= 0L) {
    gamma_supp <- integer(0)
  } else if (length(gamma_supp) > gamma_room) {
    gamma_sort <- sort(abs(gammas), decreasing = TRUE)
    threshold <- gamma_sort[gamma_room + 1L]
    gamma_supp <- which(abs(gammas) > threshold)
  }

  list(beta_supp = beta_supp, gamma_supp = gamma_supp, gammas = gammas, betas = betas,
       best_nu = best_nu, best_lam = best_lambda)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## (3) HD_SplitKnockoff
HD_SplitKnockoff <- function(
    X, y, D, q = 0.1,
    seed = NULL,
    fold = 10,
    eta = 0.1,
    method = "knockoff+",
    normalize = TRUE,
    copy = TRUE,
    lambda_seq = 10^seq(0, -6, by = -0.01),
    sign = "disabled"
) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)

  # Sample splitting
  sample_index1 <- sample(seq_len(n), size = floor(n / 2), replace = FALSE)
  sample_index2 <- setdiff(seq_len(n), sample_index1)

  # Screening
  n2 <- length(sample_index2)
  screen_res <- SplitLasso_screen(
    X[sample_index1, , drop = FALSE],
    y[sample_index1],
    D,
    n2 = n2,
    n_2 = n2,
    nfolds = fold
  )
  S_beta <- screen_res$beta_supp
  S_gamma <- screen_res$gamma_supp
  D_scr <- D[S_gamma, S_beta, drop = FALSE]
  if (length(S_gamma) == 0) {
    message("No gamma selected in the screening step! Return Nothing.")
    return(list(
      Sel = integer(0),
      S_beta = S_beta,
      S_gamma = S_gamma,
      D_scr = D_scr,
      Fit = NULL
    ))
  }

  # Selection
  option <- list(data = NA, dim = length(data), dimnames = NULL)
  option$q <- q
  option$method <- method
  option$eta <- eta
  option$normalize <- ifelse(normalize, "true", "false")
  option$copy <- ifelse(copy, "true", "false")
  option$lambda <- lambda_seq
  option$nu <- screen_res$best_nu
  option$sign <- sign
  option <- option[-1]

  SplitKnockoff_Fit <- tryCatch({
    result <- sk.filter(
      X[sample_index2, S_beta, drop = FALSE],
      D_scr,
      matrix(y[sample_index2], ncol = 1),
      option
    )
    result
  }, error = function(e) {
    warning("SplitKnockoff failed; returning an empty selection: ", conditionMessage(e))
    list(selected = NULL)
  })

  SplitKnockoff_Sel <- S_gamma[SplitKnockoff_Fit$results$data]

  return(list(
    Sel = SplitKnockoff_Sel,
    S_beta = S_beta,
    S_gamma = S_gamma,
    D_scr = D_scr,
    Fit = SplitKnockoff_Fit
  ))
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions

################################################################################
################################################################################
################################################################################
### (1) Center matrix columns
sk.center <- function(A) sweep(A, 2, colMeans(A), "-")

################################################################################
################################################################################
################################################################################
### (2) Population standard deviation
sk.sd.pop <- function(z) sqrt(mean((z - mean(z))^2))

################################################################################
################################################################################
################################################################################
### (3) Construct the augmented Split Lasso design
sk.create <- function(X, y, D, nu) {
  n <- nrow(X)
  m <- nrow(D)

  A_beta <- rbind(X / sqrt(n), D / sqrt(nu))
  A_gamma <- rbind(matrix(0, n, m), -diag(m) / sqrt(nu))

  tilde_y <- c(y / sqrt(n), rep(0, m))

  structure(list(call = match.call(),
                 A_beta = A_beta,
                 A_gamma = A_gamma,
                 tilde_y = tilde_y))
}

################################################################################
################################################################################
################################################################################
### (4) Normalize matrix columns
normc <- function(X) {
  if (is.null(dim(X)))
    X <- matrix(X, ncol = 1)

  n = nrow(X)
  p = ncol(X)
  X <- X - matrix(rep(colMeans(X), n), nrow = n, ncol = p, byrow = TRUE)
  factors <- 1 / sqrt(colSums(X^2))
  factors <- matrix(rep(factors, n), nrow = n, ncol = p, byrow = TRUE)
  Y = X * factors * sqrt(n - 1)
  return(Y)
}

################################################################################
################################################################################
################################################################################
### (5) Patch SplitKnockoff normalization
sk.patch.normc <- function() {
  if (!requireNamespace("SplitKnockoff", quietly = TRUE)) {
    warning("Package SplitKnockoff is unavailable; normc was left unpatched.")
    return(invisible(FALSE))
  }
  ok <- tryCatch({
    assignInNamespace("normc", normc, ns = "SplitKnockoff")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) warning("Could not install the aligned normc into the SplitKnockoff namespace.")
  invisible(ok)
}
sk.patch.normc()
