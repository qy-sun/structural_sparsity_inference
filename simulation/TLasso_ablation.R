################################################################################
# Simulation Ablation for TLasso
# ------------------------------------------------------------------------------
# This script reproduces Table S.2 in the Supplement
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Setup
source("helper/Graph Difference Operator.R", encoding = "UTF-8")
source("helper/Beta Structure Generation.R", encoding = "UTF-8")
source("../function/TLasso.R", encoding = "UTF-8")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions
popcount <- function(x, width) {
  vapply(x, function(value) {
    sum(as.integer(intToBits(value))[seq_len(width)])
  }, integer(1))
}

enumerate_lattice_cut5 <- function(height = 15L, width = 10L, budget = 5L) {
  masks <- 0:(2^width - 1L)
  horizontal <- popcount(bitwAnd(bitwXor(masks, bitwShiftR(masks, 1L)),
                                 2^(width - 1L) - 1L), width)
  transitions <- lapply(masks, function(previous) {
    cost <- popcount(bitwXor(previous, masks), width) + horizontal
    cbind(mask = masks[cost <= budget], cost = cost[cost <= budget])
  })

  dp <- array(0, dim = c(height, length(masks), budget + 1L))
  dp[1L, , 1L] <- 1
  for (remaining in 2:height) {
    for (left in 0:budget) {
      for (previous in masks) {
        edge <- transitions[[previous + 1L]]
        keep <- edge[, "cost"] <= left
        if (any(keep)) {
          dp[remaining, previous + 1L, left + 1L] <-
            sum(dp[remaining - 1L,
                   edge[keep, "mask"] + 1L,
                   left - edge[keep, "cost"] + 1L])
        }
      }
    }
  }

  partitions <- list()
  visit <- function(path, left) {
    remaining <- height - length(path)
    if (remaining == 0L) {
      stopifnot(left == 0L)
      partitions[[length(partitions) + 1L]] <<- unlist(
        lapply(path, function(mask) {
          bitwAnd(bitwShiftR(mask, 0:(width - 1L)), 1L)
        }), use.names = FALSE)
      return(invisible(NULL))
    }

    edge <- transitions[[tail(path, 1L) + 1L]]
    for (j in seq_len(nrow(edge))) {
      next_mask <- edge[j, "mask"]
      cost <- edge[j, "cost"]
      if (cost <= left &&
          dp[remaining, next_mask + 1L, left - cost + 1L] > 0) {
        visit(c(path, next_mask), left - cost)
      }
    }
    invisible(NULL)
  }

  for (first_mask in masks[bitwAnd(masks, 1L) == 0L]) {
    left <- budget - horizontal[first_mask + 1L]
    if (left >= 0L && dp[height, first_mask + 1L, left + 1L] > 0) {
      visit(first_mask, left)
    }
  }
  do.call(rbind, partitions)
}

diagnose_design <- function(Z, signal_index) {
  G <- crossprod(Z) / nrow(Z)
  G_signal <- G[signal_index, signal_index, drop = FALSE]
  G_inverse <- solve(G_signal)
  c(omega = 1 - max(rowSums(abs(G[-signal_index, signal_index, drop = FALSE] %*%
                                  G_inverse))),
    Cmin = min(eigen(G_signal, symmetric = TRUE, only.values = TRUE)$values),
    invnorm = max(rowSums(abs(G_inverse))))
}

fit_at_lambda <- function(Z, y, lambda, analytical = FALSE) {
  lambda_path <- lambda
  if (analytical) {
    lambda_max <- max(max(abs(crossprod(Z, y))) / nrow(Z) * 1.01,
                      lambda * 1.01)
    lambda_path <- exp(seq(log(lambda_max), log(lambda), length.out = 60))
    lambda_path[length(lambda_path)] <- lambda
  }

  fit <- glmnet::glmnet(Z, y, lambda = lambda_path, standardize = FALSE,
                        intercept = FALSE, thresh = 1e-14, maxit = 1e7)
  stopifnot(fit$jerr == 0, length(fit$lambda) == length(lambda_path),
            abs(tail(fit$lambda, 1) - lambda) < 1e-12)

  theta <- as.numeric(fit$beta[, ncol(fit$beta)])
  gradient <- as.numeric(crossprod(Z, y - Z %*% theta)) / nrow(Z)
  error <- pmax(abs(gradient) - lambda, 0)
  active <- which(theta != 0)
  error[active] <- abs(gradient[active] - lambda * sign(theta[active]))
  kkt <- max(error) / max(1, lambda)
  stopifnot(kkt < 1e-4)
  list(theta = theta, kkt = kkt)
}

summarize_results <- function(result) {
  valid <- result[result$defined, ]
  out <- data.frame(
    design = result$design[1], arm = result$arm[1],
    planned = nrow(result), defined = nrow(valid),
    undefined = sum(!result$defined), exact_count = sum(valid$exact),
    exact_rate = mean(valid$exact), subset_rate = mean(valid$subset),
    sign_rate = mean(valid$sign_exact), max_kkt = max(valid$kkt))
  for (name in c("lambda", "L2", "Linf", "TP", "MS")) {
    out[[name]] <- mean(valid[[name]])
    out[[paste0(name, "_se")]] <- sd(valid[[name]]) / sqrt(nrow(valid))
  }
  out
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Table S.2: selection behavior under different penalty choices
n <- 300L
p <- 150L
A <- 4
sigma <- 1
nIterations <- 100L
n_signal <- 5L
designs <- c("chain", "lattice1", "lattice2")

B <- create.lattice(15, 10)
e <- rep(1 / sqrt(p), p)
patterns <- enumerate_lattice_cut5()
stopifnot(identical(dim(patterns), c(300L, 150L)),
          !anyDuplicated(as.data.frame(patterns)),
          all(patterns[, 1] == 0),
          all(rowSums(abs(patterns %*% t(B))) == n_signal))

project <- TLasso
expressions <- as.list(body(TLasso))
cut <- which(vapply(expressions, function(expression) {
  is.call(expression) && identical(expression[[1]], as.name("<-")) &&
    identical(expression[[2]], as.name("cv_fit"))
}, logical(1)))[1]
stopifnot(is.finite(cut))
body(project) <- as.call(c(list(as.name("{")), expressions[2:(cut - 1)],
  list(quote(list(X = X_star, Z = X_star_scaled, scale = col_scale,
                  W = if (r == p) diag(n) else resid_proj(diag(n), XDstar))))))

all_results <- list()
all_diagnostics <- list()

################################################################################
################################################################################
################################################################################
### Run the three designs
for (design in designs) {
  D <- switch(design,
              chain = create.chain(p),
              lattice1 = B,
              lattice2 = crossprod(B))
  m <- nrow(D)
  Sigma <- crossprod(D) + tcrossprod(e)
  Sigma_factor <- chol(Sigma)
  rows <- list()
  diagnostics <- list()

  for (iter in seq_len(nIterations)) {
    design_index <- match(design, designs)
    data_seed <- 2091000 + 100 * design_index + iter
    position_seed <- 2090000 + 100 * design_index + iter

    set.seed(data_seed)
    X <- matrix(rnorm(n * p), n) %*% Sigma_factor
    noise <- rnorm(n, sd = sigma)
    projection <- project(X, noise, D)
    Z <- projection$Z
    scales <- projection$scale
    foldid <- sample(rep(1:5, length.out = n))

    set.seed(position_seed)
    partition_id <- orientation <- NA_integer_
    if (design == "lattice1") {
      partition_id <- sample.int(nrow(patterns), 1)
      orientation <- sample(c(-1, 1), 1)
      Beta0 <- orientation * A * (patterns[partition_id, ] - 0.5)
    } else {
      invisible(capture.output(
        Beta0 <- if (design == "chain") {
          Chain_PiecewiseConstant(p, p, n_signal, A)
        } else {
          Lattice_PiecewiseConstant(15, 10, p, n_signal - 1L, A)
        }))
    }

    truth <- as.numeric(D %*% Beta0)
    signal_ind <- which(abs(truth) > 1e-4)
    stopifnot(length(signal_ind) == n_signal)
    y_original <- as.numeric(X %*% Beta0 + noise)
    y <- as.numeric(projection$W %*% y_original)

    diagnostic <- diagnose_design(Z, signal_ind)
    omega <- unname(diagnostic["omega"])
    defined <- is.finite(omega) && omega > 0
    lambda_analytical <- if (defined) {
      2 * sigma / omega * sqrt(2 * log(m) / n)
    } else {
      NA_real_
    }
    bound <- if (defined) {
      (omega / sqrt(diagnostic["Cmin"]) + diagnostic["invnorm"]) *
        lambda_analytical
    } else {
      NA_real_
    }

    cv <- glmnet::cv.glmnet(Z, y, foldid = foldid, standardize = FALSE,
                            intercept = FALSE)
    for (arm in c("analytical", "CV")) {
      lambda <- if (arm == "analytical") lambda_analytical else cv$lambda.min
      L2 <- Linf <- TP <- MS <- exact <- subset <- sign_exact <- kkt <- NA_real_
      if (is.finite(lambda)) {
        fit <- fit_at_lambda(Z, y, lambda, analytical = arm == "analytical")
        gamma <- fit$theta / scales
        selected <- which(abs(gamma) > 1e-4)
        L2 <- sqrt(sum((gamma - truth)^2))
        Linf <- max(abs(gamma - truth))
        TP <- length(intersect(signal_ind, selected))
        MS <- length(selected)
        exact <- as.integer(setequal(signal_ind, selected))
        subset <- as.integer(all(selected %in% signal_ind))
        sign_exact <- as.integer(exact &&
                                 all(sign(gamma[signal_ind]) == sign(truth[signal_ind])))
        kkt <- fit$kkt
      }
      rows[[length(rows) + 1L]] <- data.frame(
        design = design, rep = iter, n = n, p = p, m = m, A = A,
        s = n_signal, arm = arm, defined = is.finite(lambda), lambda = lambda,
        L2 = L2, Linf = Linf, TP = TP, MS = MS, exact = exact,
        subset = subset, sign_exact = sign_exact, kkt = kkt)
    }

    diagnostics[[iter]] <- data.frame(
      design = design, rep = iter, data_seed = data_seed,
      position_seed = position_seed, partition_id = partition_id,
      orientation = orientation, omega = omega,
      Cmin = unname(diagnostic["Cmin"]),
      defined = defined,
      min_signal = min(abs(truth[signal_ind])),
      max_signal = max(abs(truth[signal_ind])),
      beta_min = min(abs(scales[signal_ind] * truth[signal_ind])),
      bound = unname(bound),
      beta_min_pass = if (defined) {
        min(abs(scales[signal_ind] * truth[signal_ind])) > bound
      } else {
        NA
      },
      support = paste(signal_ind, collapse = ";"),
      negative_vertices = paste(which(Beta0 < 0), collapse = ";"))

  }

  replicates <- do.call(rbind, rows)
  diagnostics <- do.call(rbind, diagnostics)

  all_results[[design]] <- replicates
  all_diagnostics[[design]] <- diagnostics
}

################################################################################
################################################################################
################################################################################
### Combine results
replicates <- do.call(rbind, all_results)
diagnostics <- do.call(rbind, all_diagnostics)
summary <- do.call(rbind, lapply(
  split(replicates, list(replicates$design, replicates$arm), drop = TRUE),
  summarize_results))
summary <- summary[
  order(match(summary$design, designs),
        match(summary$arm, c("analytical", "CV"))), ]
rownames(summary) <- NULL

save(replicates, diagnostics, summary, file = "TLasso_ablation.RData")
