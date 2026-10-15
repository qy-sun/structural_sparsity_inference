################################################################################
# Simulation Ablation for StatKnock
# ------------------------------------------------------------------------------
# This script reproduces Figure S.1 and Table S.1 in the Supplement
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
library(knockoff)
library(glmnet)
library(MASS)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Helper functions
haar_uperp <- function(X, m) {
  G  <- matrix(rnorm(nrow(X) * m), nrow(X), m)
  G  <- G - X %*% solve(crossprod(X), crossprod(X, G))
  qz <- qr(G)
  qr.Q(qz) %*% diag(sign(diag(qr.R(qz))), m)
}

create.fixed.haar <- function(X_org, method = "equi") {
  if (ncol(X_org) == 1L && nrow(X_org) > 1L && method %in% c("equi", "sdp")) {
    X <- X_org / sqrt(sum(X_org^2))
    return(list(X = X, Xk = haar_uperp(X, 1L)))
  }
  kn <- create.fixed(X_org, method = method)
  X  <- kn$X
  m  <- ncol(X)
  G  <- crossprod(X)
  s  <- diag(G - crossprod(X, kn$Xk))
  A  <- X %*% (diag(m) - solve(G, diag(s, m)))
  R  <- kn$Xk - A
  ev <- eigen(crossprod(R), symmetric = TRUE)
  ev$values[ev$values < 0] <- 0
  Cm <- diag(sqrt(ev$values), m) %*% t(ev$vectors)
  list(X = X, Xk = A + haar_uperp(X, m) %*% Cm)
}

generate_knockoff_fixed <- function(X_org) {
  X_org <- scale(X_org, center = FALSE, scale = sqrt(colSums(X_org^2)))
  knock <- create.fixed.haar(X_org)
  list(X_org = knock$X, X_tilde = knock$Xk)
}

knockoff_threshold <- function(W, q, offset) {
  ts <- sort(unique(abs(W)[abs(W) > 0]))
  ratio <- numeric(length(ts))
  for (j in seq_along(ts)) {
    ratio[j] <- (offset + sum(W <= -ts[j])) / max(1, sum(W >= ts[j]))
  }
  id <- which(ratio <= q)[1]
  if (is.na(id)) Inf else ts[id]
}


StatKnock_CombRules <- function(X, y, D, L = 100, q = 0.1,
                               cvfolds = 5, X_N_option = 1, proj = TRUE) {
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(D)[1]
  qrD <- qr(D)
  r <- qrD$rank

  # Construct the transformed design
  if (r == m) {
    D_plus <- D
    X_plus <- X
  } else {
    message("rank(D) < m, auxiliary noises introduced!")
    Q_full <- qr.Q(qrD, complete = TRUE)
    N      <- Q_full[, (r + 1):m, drop = FALSE]
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
        P_X <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
        X_N <- P_X %*% X_N
      } else {
        message("Projection skipped: n < p")
      }
    }
    X_plus <- cbind(X, X_N)
  }

  if (r == p) {
    if (n < 2 * m) {
      stop(sprintf("Insufficient effective samples for fixed-X knockoffs: %d < %d.", n, 2 * m))
    }
    X_star <- X_plus %*% solve(D_plus)
    y_star <- y
  } else {
    message("rank(D) < p, projection used!")
    E                <- pracma::nullspace(D_plus)
    D_tilde          <- rbind(D_plus, t(E))
    D_tilde_inverse  <- solve(D_tilde)
    D_dagger         <- D_tilde_inverse[, 1:m]
    D_star           <- D_tilde_inverse[, (m + 1):(m + (p - r))]
    XDstar           <- X_plus %*% D_star
    sv               <- svd(XDstar, nu = n)
    rk               <- sum(sv$d > 1e-10 * max(sv$d))
    if (n - rk < 2 * m) {
      stop(sprintf("Insufficient effective samples for fixed-X knockoffs: %d < %d.", n - rk, 2 * m))
    }
    Q_W              <- sv$u[, (rk + 1):n, drop = FALSE]
    y_star           <- t(Q_W) %*% y
    X_star           <- t(Q_W) %*% X_plus %*% D_dagger
  }

  # Generate one shared knockoff copy
  res     <- generate_knockoff_fixed(X_star)
  X_org   <- res$X_org
  X_tilde <- res$X_tilde
  X_aug   <- cbind(X_org, X_tilde)

  if (r != p) {
    X_aug  <- Q_W %*% X_aug
    X_star <- Q_W %*% X_star
    y_star <- Q_W %*% y_star
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

  cv.fit <- cv.glmnet(x = X_aug, y = y_star, family = "gaussian",
                      intercept = FALSE, alpha = 1, nfolds = cvfolds,
                      standardize = FALSE)
  best_lambda <- cv.fit$lambda.min

  # Aggregate complementary splits
  n1  <- rep(0, 2 * m)
  n2  <- rep(0, 2 * m)
  n12 <- rep(0, 2 * m)
  nor <- rep(0, 2 * m)

  for (i in 1:L) {
    index <- sample((1:n), floor(n / 2), replace = FALSE)

    fit1     <- glmnet(x = X_aug[index, ], y = y_star[index], family = "gaussian",
                       intercept = FALSE, alpha = 1, lambda = best_lambda,
                       standardize = FALSE, thresh = 1e-14, maxit = 1e7)
    bhat1    <- as.vector(stats::coef(fit1)[-1])
    S_index1 <- which(bhat1 != 0)

    fit2     <- glmnet(x = X_aug[-index, ], y = y_star[-index], family = "gaussian",
                       intercept = FALSE, alpha = 1, lambda = best_lambda,
                       standardize = FALSE, thresh = 1e-14, maxit = 1e7)
    bhat2    <- as.vector(stats::coef(fit2)[-1])
    S_index2 <- which(bhat2 != 0)

    e1 <- rep(0, 2 * m)
    e1[S_index1] <- 1
    e2 <- rep(0, 2 * m)
    e2[S_index2] <- 1

    n1  <- n1  + e1
    n2  <- n2  + e2
    n12 <- n12 + e1 * e2
    nor <- nor + pmax(e1, e2)
  }
  sw <- seq_len(2 * m); sw[c(R_set, R_set + m)] <- c(R_set + m, R_set)
  n1 <- n1[sw]; n2 <- n2[sw]; n12 <- n12[sw]; nor <- nor[sw]

  pi1     <- n1  / L
  Pi_none <- pi1
  Pi_avg  <- (n1 + n2) / (2 * L)
  Pi_and  <- n12 / L
  Pi_or   <- nor / L

  # Compute the aggregation statistics and thresholds
  spd <- function(P) P[1:m] - P[(m + 1):(2 * m)]

  W <- list(
    none = spd(Pi_none),
    avg  = spd(Pi_avg),
    or   = spd(Pi_or),
    and  = spd(Pi_and)
  )


  T_mfdr <- sapply(W, knockoff_threshold, q = q, offset = 0)
  T_fdr  <- sapply(W, knockoff_threshold, q = q, offset = 1)

  S_mfdr <- lapply(names(W), function(k) which(W[[k]] >= T_mfdr[[k]]))
  S_fdr  <- lapply(names(W), function(k) which(W[[k]] >= T_fdr[[k]]))
  names(S_mfdr) <- names(S_fdr) <- names(W)

  list(W = W, T_mfdr = T_mfdr, T_fdr = T_fdr,
       S_mfdr = S_mfdr, S_fdr = S_fdr,
       best_lambda = best_lambda)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Figure S.1 and Table S.1: StatKnock aggregation rules
set.seed(2025)
nIterations = 100
q = 0.2
L = 100

p = 120
n = 500
p1 = 120
m1 = 12
rho = 0.2
D <- create.chain(p)
m <- nrow(D)
Sigma = toeplitz(rho^(1:p - 1))
A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rules <- c("none", "avg", "or", "and")
figure_seed <- 198
figure_A <- 0.5

fdp_power <- function(selected_index, signal_index) {
  fdp = 0
  power = 0
  if (length(selected_index) != 0) {
    num_selected <- length(selected_index)
    tp <- length(intersect(selected_index, signal_index))
    fp <- num_selected - tp
    fdp <- fp / max(num_selected, 1)
    power <- tp / length(signal_index)
  }
  return(list(fdp = fdp, power = power))
}

mean.se = function(x, ns = 3) {
  paste0(format(round(mean(x), ns), nsmall = ns), " (",
         format(round(sd(x) / sqrt(length(x)), ns), nsmall = ns), ")")
}

StatKnock_ablation_result <- list()

for (A in A_values) {
  cat("Running for A =", A, "\n")

  result_matrix <- sapply(1:nIterations, function(it) {
    ################################################################################
    ################################################################################
    ################################################################################
    ### Generate the Data
    Beta0 = Chain_PiecewiseConstant(p, p1, m1, A,
                                    min_interval = 3, edge_buffer = 3)
    X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
    y <- X %*% Beta0 + rnorm(n, sd = 1)
    signal_ind = which(abs(D %*% Beta0) > 1e-4)

    ################################################################################
    ################################################################################
    ################################################################################
    ### StatKnock aggregation rules
    fit <- StatKnock_CombRules(X = X, y = y, D = D, L = L, q = q)

    r <- unlist(setNames(lapply(rules, function(rule) {
      result <- fdp_power(fit$S_fdr[[rule]], signal_ind)
      c(FDR = result$fdp, Power = result$power)
    }), rules))
    return(r)
  })

  StatKnock_ablation_result[[paste0("A=", A)]] <- result_matrix
}

for (A in A_values) {
  cat("\nResults for A =", A, "\n")
  print(apply(StatKnock_ablation_result[[paste0("A=", A)]], 1, mean.se))
}

set.seed(figure_seed)
Beta0 <- Chain_PiecewiseConstant(p, p1, m1, figure_A,
                                 min_interval = 3, edge_buffer = 3)
X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
y <- X %*% Beta0 + rnorm(n, sd = 1)
signal_ind <- which(abs(D %*% Beta0) > 1e-4)
StatKnock_ablation_figure_result <- StatKnock_CombRules(
  X = X, y = y, D = D, L = L, q = q
)
StatKnock_ablation_figure_result$seed <- figure_seed
StatKnock_ablation_figure_result$setting <- c(
  n = n, p = p, m = m, m1 = m1, A = figure_A, rho = rho, q = q, L = L
)
StatKnock_ablation_figure_result$signal_ind <- signal_ind

save(StatKnock_ablation_result, StatKnock_ablation_figure_result,
     file = "StatKnock_ablation.RData")

################################################################################
################################################################################
################################################################################
### Figure S.1
library(ggplot2)
library(gridExtra)

rule_labels <- c(none = "(a) original rule", avg = "(b) average rule",
                 or = "(c) \u2018OR\u2019 rule", and = "(d) \u2018AND\u2019 rule")
fit <- StatKnock_ablation_figure_result

make_panel <- function(rule) {
  W <- fit$W[[rule]]
  T0 <- fit$T_mfdr[[rule]]
  T1 <- fit$T_fdr[[rule]]
  statistics <- data.frame(
    Index = seq_along(W),
    W = W,
    Signal = ifelse(seq_along(W) %in% fit$signal_ind, "Signal", "Null")
  )
  statistics$Signal <- factor(statistics$Signal,
                              levels = c("Null", "Signal"))

  breaks <- c(-1, 0, 1)
  labels <- c("-1.0", "0.0", "1.0")
  if (is.finite(T1) && T1 > 0.08 && T1 < 0.92) {
    breaks <- sort(c(breaks, -T1, T1))
    labels <- c("-1.0", "-T", "0.0", "T", "1.0")
  }

  panel <- ggplot(statistics, aes(Index, W, color = Signal,
                                  shape = Signal)) +
    geom_point(size = 1.5, show.legend = FALSE) +
    scale_color_manual(values = c(Null = "black", Signal = "green3")) +
    scale_shape_manual(values = c(Null = 19, Signal = 15)) +
    scale_y_continuous(breaks = breaks, labels = labels,
                       limits = c(-1.2, 1.2)) +
    labs(x = "Index", y = "W", caption = rule_labels[[rule]]) +
    theme(
      axis.title = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      plot.caption = element_text(size = 11, hjust = 0.5,
                                  family = "serif",
                                  margin = margin(t = 7, b = 2))
    )

  if (is.finite(T0)) {
    panel <- panel +
      geom_hline(yintercept = c(-T0, T0), linetype = "dashed",
                 color = "skyblue", linewidth = 0.5)
  }
  if (is.finite(T1)) {
    panel <- panel +
      geom_hline(yintercept = c(-T1, T1), linetype = "dashed",
                 color = "blue", linewidth = 0.5)
  }
  panel
}

panels <- lapply(rules, make_panel)
grDevices::pdf(file = NULL)
figure <- arrangeGrob(grobs = panels, nrow = 2, ncol = 2)
ggsave("AND_rule_ablation.png", plot = figure,
       width = 8.5, height = 6, dpi = 1000, bg = "white")
invisible(grDevices::dev.off())
