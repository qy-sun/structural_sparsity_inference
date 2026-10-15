################################################################################
# Benjamini-Yekutieli (BY) for Structural Sparsity
# ------------------------------------------------------------------------------
# This script extends the BY procedure to structural sparsity problems.
#
# Reference:
#   Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate
#   in multiple testing under dependency. Annals of Statistics, 29(4), 1165-1188.
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## BY_GLasso
BY_GLasso <- function(X, y, D, q, method = "auto") {
  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(D)

  # Check if n < p and choose method accordingly
  if (method == "auto") {
    if (n < p) {
      method <- "univariate_high"
    } else {
      method <- "OLS"
    }
  }

  # Define methods
  if (method == "univariate_high") {
    # Univariate regression method for high dimensions
    p_values <- univariate_high_estimate(X, y, D, p, m)

  } else if (method == "OLS") {
    # OLS method
    OLS_Fit <- lm(y ~ X - 1)
    coefficients <- summary(OLS_Fit)$coefficients[,"Estimate"]
    noise_variance <- summary(OLS_Fit)$sigma^2
    t_statistics <- numeric(m)
    XtX_inv <- solve(t(X) %*% X)
    for (j in 1:m) {
      t_statistics[j] <- D[j, ] %*% coefficients /
        sqrt(noise_variance * t(D[j, ]) %*% XtX_inv %*% D[j, ])
    }
    p_values <- 2 * pt(-abs(t_statistics), n - p)
  }
  # BY selection
  BY_Sel <- BY_selection(p_values, q)

  return(BY_Sel)
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
### (1) Univariate Regression Method for high dimensions
univariate_high_estimate <- function(X, y, D, p, m) {
  uni_beta <- numeric(p)
  uni_std <- numeric(p)
  for (j in 1:p) {
    x <- X[, j]
    model <- lm(y ~ x)
    uni_beta[j] <- summary(model)$coefficients["x", "Estimate"]
    uni_std[j] <- summary(model)$coefficients["x", "Std. Error"]
  }
  uni_beta_scaled <- uni_beta / uni_std
  uni_D <- D %*% uni_beta_scaled
  uni_D_rescaled <- numeric(m)
  for (j in 1:m) {
    uni_D_rescaled[j] <- uni_D[j] / sqrt(sum(D[j, ]^2))
  }
  p_values <- 2 * pnorm(abs(uni_D_rescaled), lower.tail = FALSE)
  return(p_values)
}
################################################################################
################################################################################
################################################################################
### (2) Function to perform BY adjustment and selection
BY_selection <- function(p_values, q) {
  BY <- (p.adjust(p_values, method = "BY") < q)
  BY_Sel <- which(as.vector(BY) == TRUE)
  return(BY_Sel)
}
