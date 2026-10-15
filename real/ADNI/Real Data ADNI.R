################################################################################
# Real Data Experiment (ADNI)
# ------------------------------------------------------------------------------
# This script applies multiple structural FDR-control methods to brain imaging
# data from the ADNI (Alzheimer's Disease Neuroimaging Initiative) dataset.
# The data include subject-level sMRI features (Level 5 parcellation).
# 
# Note: Access to the ADNI dataset is governed by privacy and data-sharing terms;
# obtain and preprocess the data before running this script.
#
# Workflow:
#   1) Load sMRI feature matrix (X) and response (y).
#   2) Construct transformational matrices D1, D2, D3:
#        - D1: identity (zero-order)
#        - D2: first-order transformation
#        - D3: second-order transformation (trimmed to remove zero rows)
#   3) Apply multiple methods:
#        - SATLasso + StatKnock
#        - BY for structural sparsity
#        - High-dimensional SplitKnockoff
#   4) Count feature/edge selection frequencies across 200 iterations
#        for each method-structure combination.
#
# Outputs:
#   - frequency_ADNI.csv: combined summary of selection frequencies
#     (columns: Method, Group [1-3], Index, Count)
################################################################################

source("../../function/SATLasso+StatKnock.R", chdir = TRUE)
source("../../function/other methods/HD_SplitKnockoff.R")
source("../../function/other methods/BY_Structural.R")
source("../../simulation/helper/Graph Difference Operator.R")
library(SplitKnockoff) # package version 1.2
library(readr)
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Real Data Settings
nIterations = 200
q = 0.2
screening_fraction = 0.4
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Load Data
data <- read_csv("ADNI_sMRI_Level5.csv")
metadata_columns <- c("...1", "X", "id", "mri_scan_date", "cognitive_exam_date", "ADNI_MEM")
feature_columns <- setdiff(colnames(data), metadata_columns)
y <- as.matrix(data$ADNI_MEM); y <- y - mean(y)
X <- as.matrix(data[, feature_columns])
stopifnot(ncol(X) == 279L)
X <- scale(X, center = TRUE,  scale = FALSE)
X <- scale(X, center = FALSE, scale = sqrt(colSums(X^2) / nrow(X)))[,]

n = dim(X)[1]
p = dim(X)[2]
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Create the Graph
data_lookup <- read.csv("lookup/multilevel_lookup_table.csv")
pairs <- list()
# Loop over each possible pair
for (i in 1:(nrow(data_lookup) - 1)) {
  for (j in (i + 1):nrow(data_lookup)) {
    # Check if levels match
    if (data_lookup$Level_3[i] == data_lookup$Level_3[j]) {
      # Store indices of matching pairs
      pairs[[length(pairs) + 1]] <- c(i, j)
    }
  }
}
pair_indices <- do.call(rbind, pairs)
pair_indices <- as.matrix(pair_indices)
stopifnot(nrow(pair_indices) == 1122L)

col1 <- data_lookup[[1]]
col2 <- data_lookup[[2]]
src_endpoint <- paste(col1[pair_indices[,1]], col2[pair_indices[,1]], sep = "_")
dst_endpoint <- paste(col1[pair_indices[,2]], col2[pair_indices[,2]], sep = "_")
edge_df <- data.frame(
  source = src_endpoint,
  target = dst_endpoint,
  stringsAsFactors = FALSE
)
edge_df <- read.csv("lookup/edge_endpoints.csv")
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Different Transformation Matrices

D_1 <- diag(p)

D_2 <- matrix(0, nrow = nrow(pair_indices), ncol = nrow(data_lookup))
for (k in 1:nrow(pair_indices)) {
  D_2[k, pair_indices[k, 1]] <- -1
  D_2[k, pair_indices[k, 2]] <- 1
}

D_3_raw <- t(D_2) %*% D_2

process_matrix <- function(D) {
  keep_rows <- which(!apply(D, 1, function(row) all(row == 0)))
  
  D_trimmed <- D[keep_rows, , drop = FALSE]
  
  list(
    matrix = D_trimmed,
    rows_kept = keep_rows
  )
}

D_3_processed <- process_matrix(D_3_raw)
D_3 <- D_3_processed$matrix
D_3_index <- D_3_processed$rows_kept
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Initialize Selection Frequencies
m_1 <- nrow(D_1)
print(m_1)
m_2 <- nrow(D_2)
print(m_2)
m_3 <- nrow(D_3)
print(m_3)

methods <- c("SATLasso_StatKnock", "SATLasso_StatKnock_plus",
             "BY", "SplitKnockoff", "SplitKnockoff_plus")

freq_len <- function(g) if (g == 3) p else get(paste0("m_", g))
for (i in 1:3) {
  m_val <- freq_len(i)
  for (method in methods) {
    assign(paste0(method, "_frequency_", i), numeric(m_val))
  }
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## Run methods
################################################################################
################################################################################
################################################################################
### BY
BY_Sel_1 = BY_GLasso(X = X, y = y, D = D_1, q = q)
BY_Sel_2 = BY_GLasso(X = X, y = y, D = D_2, q = q)
BY_Sel_3 = D_3_index[BY_GLasso(X = X, y = y, D = D_3, q = q)]

for (i in 1:nIterations) {
  print(i)
  set.seed(i)
  ################################################################################
  ################################################################################
  ################################################################################
  ### SATLasso+StatKnock
  SATLasso_StatKnock_Fit_1 <- SATLasso_StatKnock(X = X, y = y, D = D_1, q = q, gbic = FALSE,
                                                 screening_fraction = screening_fraction)
  SATLasso_StatKnock_Sel_1 <- SATLasso_StatKnock_Fit_1$S
  SATLasso_StatKnock_plus_Sel_1 <- SATLasso_StatKnock_Fit_1$S_plus
  SATLasso_StatKnock_Fit_2 <- SATLasso_StatKnock(X = X, y = y, D = D_2, q = q, gbic = FALSE,
                                                 screening_fraction = screening_fraction)
  SATLasso_StatKnock_Sel_2 <- SATLasso_StatKnock_Fit_2$S
  SATLasso_StatKnock_plus_Sel_2 <- SATLasso_StatKnock_Fit_2$S_plus
  SATLasso_StatKnock_Fit_3 <- SATLasso_StatKnock(X = X, y = y, D = D_3, q = q, gbic = FALSE,
                                                 screening_fraction = screening_fraction)
  SATLasso_StatKnock_Sel_3 <- D_3_index[SATLasso_StatKnock_Fit_3$S]
  SATLasso_StatKnock_plus_Sel_3 <- D_3_index[SATLasso_StatKnock_Fit_3$S_plus]
  ################################################################################
  ################################################################################
  ################################################################################
  ### Split Knockoff
  s1 <- sample.int(1e9, 1); s2 <- sample.int(1e9, 1); s3 <- sample.int(1e9, 1)
  SplitKnockoff_fit <- HD_SplitKnockoff(X, y, D = D_1, q = q, method = "knockoff", seed = s1)
  SplitKnockoff_Sel_1 <- SplitKnockoff_fit$Sel
  SplitKnockoff_plus_Sel_1 <- if (is.null(SplitKnockoff_fit$Fit$W)) SplitKnockoff_fit$Sel else
    SplitKnockoff_fit$S_gamma[SplitKnockoff::sk.select(
      SplitKnockoff_fit$Fit$W, q, list(method = "knockoff+"))$S]
  SplitKnockoff_fit <- HD_SplitKnockoff(X, y, D = D_2, q = q, method = "knockoff", seed = s2)
  SplitKnockoff_Sel_2 <- SplitKnockoff_fit$Sel
  SplitKnockoff_plus_Sel_2 <- if (is.null(SplitKnockoff_fit$Fit$W)) SplitKnockoff_fit$Sel else
    SplitKnockoff_fit$S_gamma[SplitKnockoff::sk.select(
      SplitKnockoff_fit$Fit$W, q, list(method = "knockoff+"))$S]
  SplitKnockoff_fit <- HD_SplitKnockoff(X, y, D = D_3, q = q, method = "knockoff", seed = s3)
  SplitKnockoff_Sel_3 <- D_3_index[SplitKnockoff_fit$Sel]
  SplitKnockoff_plus_Sel_3 <- if (is.null(SplitKnockoff_fit$Fit$W)) SplitKnockoff_fit$Sel else
    SplitKnockoff_fit$S_gamma[SplitKnockoff::sk.select(
      SplitKnockoff_fit$Fit$W, q, list(method = "knockoff+"))$S]
  SplitKnockoff_plus_Sel_3 <- D_3_index[SplitKnockoff_plus_Sel_3]
  rm(SplitKnockoff_fit)
  
  for (method in methods) {
    for (g in 1:3) {
      freq_name <- paste0(method, "_frequency_", g)
      sel_name <- paste0(method, "_Sel_", g)
      
      if (exists(freq_name) && exists(sel_name)) {
        freq <- get(freq_name)
        sel <- get(sel_name)
        
        if (length(sel)) stopifnot(all(sel >= 1), all(sel <= freq_len(g)))
        freq[sel] <- freq[sel] + 1
        
        assign(freq_name, freq)
      }
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
  for (i in 1:3) {
    freq_name <- paste0(method, "_frequency_", i)
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
          Group = i,
          Index = sorted_index,
          Count = sorted_freq
        )
        
        summary_list[[paste0(method, "_", i)]] <- df
      }
    }
  }
}

final_summary <- do.call(rbind, summary_list)

write.csv(final_summary, file = "frequency_ADNI.csv", row.names = FALSE)
