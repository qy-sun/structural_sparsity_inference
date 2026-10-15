################################################################################
# Beta Structure Generation
# ------------------------------------------------------------------------------
# This script synthesizes structured coefficient vectors beta (Beta0) for
# simulation studies under structural sparsity:
#
#   1. Chain graph
#      - Chain_PiecewiseConstant(p, p1, m1, A, min_interval=1, edge_buffer=1)
#          Alternating +/-A/2 blocks; m1 counts internal changepoints plus
#          the boundary to the zero tail when p1 < p. Returns length-p vector.
#      - Chain_PiecewiseLinear(p, p1, m1, A_target, D, min_interval=1, edge_buffer=1)
#          Piecewise linear segments scaled so the smallest magnitude of
#          D %*% beta above 1e-8 equals A_target.
#
#   2. Lattice graph
#      - Lattice_PiecewiseConstant(s, t, p1, m1, A)
#          Creates an sxt matrix with A/2 blocks and sparse -A/2 inserts
#          to induce edges; returns vectorized length-(s*t) Beta0.
#
#   3. Complete graph
#      - CompleteGraph_PiecewiseConstant(p, p1, m1, A)
#          First p1 entries start at A, then m1 positions are randomly
#          adjusted (centered to keep mean near A); rest set to 0.
#
# Conventions:
#   - p = total length; p1 = requested prefix size; m1 = # of changepoints/edits.
#   - For lattice graphs, the initial block has floor(p1/t)*t vertices (complete rows).
#   - A (or A_target) controls magnitude of structural transformations.
#   - Functions may plot Beta0 (and D %*% Beta0) for quick inspection.
################################################################################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## 1. Chain Graph
################################################################################
################################################################################
################################################################################
### (1) Piecewise Constant over Chain Graph
Chain_PiecewiseConstant <- function(p, p1, m1, A, min_interval = 1, edge_buffer = 1) {
  if (p1 < p) {
    m1 = m1 - 1
  }

  generate_changepoints <- function(p1, m1, min_interval, edge_buffer) {
    if (m1 < 0) stop("m1 is negative after accounting for the inactive tail.")
    if (m1 == 0) return(integer(0))
    if (edge_buffer > p1 - edge_buffer) {
      stop("No candidate positions satisfy the edge buffer constraint.")
    }
    candidate_range <- (edge_buffer):(p1 - edge_buffer)
    max_possible <- floor((length(candidate_range) - 1) / min_interval) + 1

    if (m1 > max_possible) {
      stop("m1 is too large to satisfy the minimum interval and edge buffer constraints.")
    }

    for (attempt in seq_len(1000L)) {
      positions <- c()
      available <- candidate_range

      while (length(positions) < m1) {
        if (length(available) == 0) {
          break
        }
        pos <- available[sample.int(length(available), 1)]
        positions <- c(positions, pos)
        available <- available[abs(available - pos) >= min_interval]
      }

      if (length(positions) == m1) {
        return(sort(positions))
      }
    }
    stop(sprintf("Could not generate %s changepoints after 1000 attempts (min_interval=%s, edge_buffer=%s).",
                 m1, min_interval, edge_buffer))
  }

  changepos <- generate_changepoints(p1, m1, min_interval, edge_buffer)
  changepos_ex <- c(0, changepos, p1)
  print(changepos)

  Beta0 <- numeric(p)
  A1 = A/2

  for (i in 1:(length(changepos_ex) - 1)) {
    start <- changepos_ex[i] + 1
    end <- changepos_ex[i + 1]
    Beta0[start:end] <- A1
    A1 <- -A1
  }

  if (p1 < p) {
    Beta0[(p1 + 1):p] <- 0
  }

  plot(Beta0, main = "Beta0 Structure")

  return(Beta0)
}
################################################################################
################################################################################
################################################################################
### (2) Piecewise Linear over Chain Graph
Chain_PiecewiseLinear <- function(p, p1, m1, A_target, D,
                                  min_interval = 1, edge_buffer = 1) {
  if (p1 < p) {
    m1 = m1 - 2
  }

  generate_changepoints <- function(p1, m1, min_interval, edge_buffer) {
    if (m1 < 0) stop("m1 is negative after accounting for the inactive tail.")
    if (m1 == 0) return(integer(0))
    if (edge_buffer > p1 - edge_buffer) {
      stop("No candidate positions satisfy the edge buffer constraint.")
    }
    candidate_range <- edge_buffer:(p1 - edge_buffer)
    max_possible <- floor((length(candidate_range) - 1) / min_interval) + 1
    if (m1 > max_possible) {
      stop("m1 is too large to satisfy the minimum interval and edge buffer constraints.")
    }
    for (attempt in seq_len(1000L)) {
      positions <- c()
      available <- candidate_range
      while (length(positions) < m1) {
        if (length(available) == 0) {
          break
        }
        pos <- available[sample.int(length(available), 1)]
        if (all(abs(positions - pos) >= min_interval)) {
          positions <- c(positions, pos)
        }
        available <- setdiff(available, pos)
      }
      if (length(positions) == m1) {
        return(sort(positions))
      }
    }
    stop(sprintf("Could not generate %s changepoints after 1000 attempts (min_interval=%s, edge_buffer=%s).",
                 m1, min_interval, edge_buffer))
  }

  changepos <- generate_changepoints(p1, m1, min_interval, edge_buffer)
  changepos_ex <- c(1, changepos, p1)
  print(changepos)

  Beta0_raw <- numeric(p)
  current_value <- 1
  for (i in 1:(length(changepos_ex) - 1)) {
    start <- changepos_ex[i]
    end <- changepos_ex[i + 1]
    next_value <- ifelse(current_value == 1, -1, 1)
    Beta0_raw[start:end] <- seq(from = current_value, to = next_value, length.out = end - start + 1)
    current_value <- next_value
  }
  if (p1 < p) {
    Beta0_raw[(p1 + 1):p] <- 0
  }

  D_beta <- as.vector(D %*% Beta0_raw)
  min_val_nonzero <- min(abs(D_beta[abs(D_beta) > 1e-8]))
  scale_factor <- A_target / min_val_nonzero
  Beta0 <- scale_factor * Beta0_raw

  plot(Beta0)
  plot(D %*% Beta0)
  return(Beta0)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## 2. Lattice Graph
Lattice_PiecewiseConstant <- function(s, t, p1, m1, A) {
  # s         : number of rows
  # t         : number of columns
  # p1        : initial block size parameter; fills floor(p1/t) complete rows
  # m1        : edge-count parameter; subtracts t if p1 < s*t, then inserts floor(m1/4) points
  # A         : twice the absolute value of the nonzero coefficients (+/-A/2)

  if (p1 < s*t) {
    m1 = m1 - t
  }
  if (m1 < 0) stop("m1 is negative after accounting for boundary edges.")

  A1 = A/2                # absolute value of the nonzero coefficients
  p1_rows = floor(p1/t)    # number of rows with non-zero entries (A1)
  m0 = floor(m1/4)         # number of -A1 points to insert, sparsely

  for (attempt in seq_len(1000L)) {
    Beta0 <- matrix(0, nrow = s, ncol = t)
    Beta0[1:p1_rows, ] <- A1
    if (m0 == 0) break

    valid_positions <- which(
      Beta0 == A1 &
        rbind(0, Beta0[-s,]) != 0 &
        rbind(Beta0[-1,], 0) != 0 &
        cbind(0, Beta0[,-t]) != 0 &
        cbind(Beta0[,-1], 0) != 0,
      arr.ind = TRUE
    )
    if (nrow(valid_positions) < m0) {
      stop("Too few candidate positions for the requested lattice inserts.")
    }

    chosen <- matrix(numeric(0), ncol = 2)
    for (i in sample.int(nrow(valid_positions))) {
      pos <- valid_positions[i, ]
      row_range <- max(1, pos[1]-2):min(s, pos[1]+2)
      col_range <- max(1, pos[2]-2):min(t, pos[2]+2)

      if (all(Beta0[row_range, col_range] != -A1)) {
        Beta0[pos[1], pos[2]] <- -A1
        chosen <- rbind(chosen, pos)
        if (nrow(chosen) == m0) break
      }
    }

    if (nrow(chosen) == m0) break
  }
  if (m0 > 0 && nrow(chosen) != m0) {
    stop(sprintf("Could not place %s points with the 2-cell exclusion buffer after 1000 attempts (s=%s, t=%s, p1=%s).",
                 m0, s, t, p1))
  }

  # Plot
  image(1:t, 1:s, t(Beta0[s:1, ]), col = c("white", "red", "blue"), axes = FALSE,
        main = "Sparse Grid Beta0")
  text(rep(1:t, each = s), rep(s:1, t), labels = as.vector(Beta0))

  Beta0 = as.vector(t(Beta0))

  return(Beta0)
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
## 3. Complete Graph
CompleteGraph_PiecewiseConstant <- function(p, p1, m1, A) {
  # p    : total length of Beta0
  # p1   : number of initial values set to A (non-zero block)
  # m1   : number of values to randomly adjust among the first p1 entries
  # A    : base value for initial non-zero entries

  # Initialize Beta0
  Beta0 <- rep(0, p)
  Beta0[1:p1] <- A

  # Choose m1 positions to adjust
  changepos <- sort(sample(1:p1, m1))

  # Generate random adjustments and normalize to maintain mean A
  adjustment <- rnorm(m1)
  adjustment <- adjustment - mean(adjustment) + A

  # Apply adjustments
  Beta0[changepos] <- adjustment

  # Plot
  plot(Beta0, type = "o", main = "Beta0 with Adjustments",
       xlab = "Index", ylab = "Value")

  return(Beta0)
}
