################################################################################
# Real Data Experiment (FF5)
# ------------------------------------------------------------------------------
# This script applies TLasso and StatKnock to the Banks industry portfolio and
# the Fama-French five factors from the Kenneth French Data Library.
#
# The two source archives are public. Run "Real Data Preprocess_FF5.R" once to
# download them and build financial_banks_ff5_model_data.csv.
#
# Workflow:
#   1) Load the daily excess-return and factor panel.
#   2) Split the time axis into quarterly blocks, build the block-diagonal
#      design X, and construct the first-order transformation matrix D.
#   3) Apply TLasso for estimation and StatKnock for FDR control.
#   4) Count transformation selection frequencies across 200 iterations
#      for each method.
#
# Outputs:
#   - frequency_FF5.csv: combined summary of selection frequencies
#     (columns: Method, Index, Count)
#   - TLasso_Fits_FF5.Rdata: TLasso beta, gamma, and beta_N from all iterations
#   - lookup/FF5_row_labels.csv: Index to (block, coefficient, type)
################################################################################

source("../../function/TLasso.R")
source("../../function/StatKnock.R")
library(readr)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Real Data Settings
set.seed(2025)
nIterations = 200
q = 0.2
FREQ = "quarter"
FREE_INTERCEPT = TRUE

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Load Data
data <- read_csv("financial_banks_ff5_model_data.csv")
data <- as.data.frame(data)
data$date <- as.Date(data$date)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Create the Block Design
FF5_FACTORS <- c("Mkt-RF", "SMB", "HML", "RMW", "CMA")
FF5_COEF_NAMES <- c("alpha", "beta_MKT", "beta_SMB", "beta_HML", "beta_RMW", "beta_CMA")

make_block_id <- function(dates, freq = c("year", "half", "quarter", "bimonth", "month", "week")) {
  freq <- match.arg(freq)
  dates <- as.Date(dates)
  yr <- as.integer(format(dates, "%Y"))
  mo <- as.integer(format(dates, "%m"))

  key <- switch(freq,
    year    = sprintf("%04d", yr),
    half    = sprintf("%04d-H%d", yr, (mo - 1) %/% 6 + 1),
    quarter = sprintf("%04d-Q%d", yr, (mo - 1) %/% 3 + 1),
    bimonth = sprintf("%04d-B%d", yr, (mo - 1) %/% 2 + 1),
    month   = sprintf("%04d-%02d", yr, mo),
    # ISO week uses the ISO year at the turn of the year
    week    = format(dates, "%G-W%V")
  )

  lab <- sort(unique(key))
  list(id = match(key, lab), label = lab, B = length(lab), freq = freq)
}

make_ff5_design <- function(dat, freq = "year", response = "Banks_excess") {
  blk <- make_block_id(dat$date, freq)
  B <- blk$B
  K <- length(FF5_FACTORS)
  Kp1 <- K + 1L
  n <- nrow(dat)

  y <- as.numeric(dat[[response]])
  F_mat <- cbind(1, as.matrix(dat[, FF5_FACTORS, drop = FALSE]))
  colnames(F_mat) <- FF5_COEF_NAMES

  # X is block diagonal: row t is nonzero only in the columns of block b(t)
  X <- matrix(0, nrow = n, ncol = B * Kp1)
  col_off <- (blk$id - 1L) * Kp1
  for (k in seq_len(Kp1)) {
    X[cbind(seq_len(n), col_off + k)] <- F_mat[, k]
  }
  colnames(X) <- paste0(rep(blk$label, each = Kp1), ":", rep(FF5_COEF_NAMES, times = B))

  list(X = X, y = y, block = blk, B = B, K = K, p = B * Kp1,
       n = n, n_per_block = as.integer(table(blk$id)))
}

design <- make_ff5_design(data, freq = FREQ)
X <- design$X
y <- design$y

n = dim(X)[1]
p = dim(X)[2]
B = design$B

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Transformation Matrix
create.chain <- function(p) {
  assoc_matrix <- matrix(0, nrow = p - 1, ncol = p)
  for (i in 1:(p - 1)) {
    assoc_matrix[i, i] <- -1
    assoc_matrix[i, i + 1] <- 1
  }
  return(assoc_matrix)
}

make_ff5_D <- function(B, free_intercept = TRUE) {
  K <- length(FF5_FACTORS)
  E <- if (free_intercept) cbind(matrix(0, K, 1), diag(K)) else diag(K + 1)

  G <- create.chain(B)
  kronecker(G, E)
}

D <- make_ff5_D(B, free_intercept = FREE_INTERCEPT)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Create Row Labels
# The boundary is reported at the later block
label_D_rows <- function(block_label, free_intercept = TRUE) {
  nm <- if (free_intercept) FF5_COEF_NAMES[-1] else FF5_COEF_NAMES
  B <- length(block_label)
  data.frame(
    row = seq_len((B - 1) * length(nm)),
    block = rep(block_label[-1], each = length(nm)),
    coefficient = rep(nm, times = B - 1),
    type = "break",
    stringsAsFactors = FALSE
  )
}

row_labels <- label_D_rows(design$block$label, FREE_INTERCEPT)
names(row_labels)[names(row_labels) == "row"] <- "Index"
write.csv(row_labels, "lookup/FF5_row_labels.csv", row.names = FALSE)

m <- nrow(D)
print(m)

stopifnot(n - p + qr(D)$rank >= 2 * m)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Initialize Selection Frequencies
methods <- c("TLasso", "StatKnock", "StatKnock_plus")

for (method in methods) {
  assign(paste0(method, "_frequency"), numeric(m))
}

TLasso_Fits <- vector("list", nIterations)

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Run methods
for (i in 1:nIterations) {
  print(i)

  ################################################################################
  ################################################################################
  ################################################################################
  ### TLasso
  TLasso_Fit <- TLasso(X = X, y = y, D = D)
  TLasso_Fits[[i]] <- TLasso_Fit[c("beta", "gamma", "beta_N")]
  TLasso_Sel <- which(abs(TLasso_Fit$gamma) > 1e-4)

  ################################################################################
  ################################################################################
  ################################################################################
  ### StatKnock
  StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
  StatKnock_Sel <- StatKnock_Fit$S
  StatKnock_plus_Sel <- StatKnock_Fit$S_plus

  for (method in methods) {
    freq_name <- paste0(method, "_frequency")
    sel_name <- paste0(method, "_Sel")

    if (exists(freq_name) && exists(sel_name)) {
      freq <- get(freq_name)
      sel <- get(sel_name)

      if (length(sel)) stopifnot(all(sel >= 1), all(sel <= m))
      freq[sel] <- freq[sel] + 1

      assign(freq_name, freq)
    }
  }
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Summarize Results
summary_list <- list()

for (method in methods) {
  freq_name <- paste0(method, "_frequency")
  if (exists(freq_name)) {
    freq <- get(freq_name)
    nonzero_index <- which(freq > 0)
    nonzero_freq <- freq[nonzero_index]

    if (length(nonzero_index) > 0) {
      sorted_order <- order(nonzero_freq, decreasing = TRUE)
      sorted_index <- nonzero_index[sorted_order]
      sorted_freq <- nonzero_freq[sorted_order]

      df <- data.frame(
        Method = method,
        Index = sorted_index,
        Count = sorted_freq
      )

      summary_list[[method]] <- df
    }
  }
}

final_summary <- do.call(rbind, summary_list)

write.csv(final_summary, file = "frequency_FF5.csv", row.names = FALSE)
save(TLasso_Fits, file = "TLasso_Fits_FF5.Rdata")
