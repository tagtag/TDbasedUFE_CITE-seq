library(Matrix)
library(RSpectra)

# -----------------------------
# 0) ユーティリティ
# -----------------------------
split_sparse_by_col <- function(X, block_size = 2000L) {
  p <- ncol(X)
  idx <- split(seq_len(p), ceiling(seq_len(p) / block_size))
  lapply(idx, function(jj) X[, jj, drop = FALSE])
}

block_col_ranges <- function(blocks) {
  lens <- vapply(blocks, ncol, integer(1))
  ends <- cumsum(lens)
  starts <- ends - lens + 1L
  data.frame(start = starts, end = ends, len = lens)
}

# -----------------------------
# 1) Mode-1 unfolding の左特異ベクトル U1
#    X_(1) X_(1)^T = A diag(||b_j||^2) A^T
# -----------------------------
make_mode1_op <- function(A_blocks, B_blocks) {
  stopifnot(length(A_blocks) == length(B_blocks))
  N <- nrow(A_blocks[[1]])

  # w_j = ||b_j||^2 をブロックで保持
  w_blocks <- lapply(B_blocks, function(Bb) colSums(Bb^2))

  f <- function(x, args) {
    y <- numeric(N)
    for (b in seq_along(A_blocks)) {
      Ab <- A_blocks[[b]]
      wb <- w_blocks[[b]]
      t1 <- as.numeric(crossprod(Ab, x))   # len_b
      t1 <- wb * t1
      y <- y + as.numeric(Ab %*% t1)
    }
    y
  }
  list(f = f, n = N)
}

# -----------------------------
# 2) Mode-2 unfolding の左特異ベクトル U2
#    X_(2) X_(2)^T = B diag(||a_j||^2) B^T
# -----------------------------
make_mode2_op <- function(A_blocks, B_blocks) {
  stopifnot(length(A_blocks) == length(B_blocks))
  K <- nrow(B_blocks[[1]])

  w_blocks <- lapply(A_blocks, function(Ab) colSums(Ab^2))

  f <- function(x, args) {
    y <- numeric(K)
    for (b in seq_along(B_blocks)) {
      Bb <- B_blocks[[b]]
      wb <- w_blocks[[b]]
      t1 <- as.numeric(crossprod(Bb, x))
      t1 <- wb * t1
      y <- y + as.numeric(Bb %*% t1)
    }
    y
  }
  list(f = f, n = K)
}

# -----------------------------
# 3) Mode-3 unfolding の左特異ベクトル U3
#    X_(3) X_(3)^T = (A^T A) .* (B^T B)
#    ※巨大なのでブロック演算で y = G3 x を実装
# -----------------------------
make_mode3_op <- function(A_blocks, B_blocks) {
  stopifnot(length(A_blocks) == length(B_blocks))
  info <- block_col_ranges(A_blocks)
  M <- sum(info$len)

  f <- function(x, args) {
    y <- numeric(M)

    for (p in seq_along(A_blocks)) {
      Ap <- A_blocks[[p]]; Bp <- B_blocks[[p]]
      rp <- info[p, ]
      xp <- x[rp$start:rp$end]

      # 対角ブロック
      GAp <- as.matrix(crossprod(Ap, Ap))
      GBp <- as.matrix(crossprod(Bp, Bp))
      Hp <- GAp * GBp
      y[rp$start:rp$end] <- y[rp$start:rp$end] + as.numeric(Hp %*% xp)

      # 非対角ブロック
      if (p < length(A_blocks)) {
        for (q in (p + 1L):length(A_blocks)) {
          Aq <- A_blocks[[q]]; Bq <- B_blocks[[q]]
          rq <- info[q, ]
          xq <- x[rq$start:rq$end]

          GA_pq <- as.matrix(crossprod(Ap, Aq))
          GB_pq <- as.matrix(crossprod(Bp, Bq))
          Hpq <- GA_pq * GB_pq  # len_p x len_q

          y[rp$start:rp$end] <- y[rp$start:rp$end] + as.numeric(Hpq %*% xq)
          y[rq$start:rq$end] <- y[rq$start:rq$end] + as.numeric(t(Hpq) %*% xp)
        }
      }
    }

    y
  }

  list(f = f, n = M)
}

# -----------------------------
# 4) HOSVD（上位r1,r2,r3）
# -----------------------------
hosvd_rank_trunc <- function(A, B, r1, r2, r3, block_size = 2000L, eigs_tol = 1e-6) {
  stopifnot(inherits(A, "sparseMatrix"),
            inherits(B, "sparseMatrix"),
            ncol(A) == ncol(B))

  A_blocks <- split_sparse_by_col(A, block_size)
  B_blocks <- split_sparse_by_col(B, block_size)

  # Mode-1
  op1 <- make_mode1_op(A_blocks, B_blocks)
  e1 <- eigs_sym(op1$f, k = r1, n = op1$n, which = "LM", opts = list(tol = eigs_tol))
  U1 <- e1$vectors
  s1 <- sqrt(pmax(e1$values, 0))

  # Mode-2
  op2 <- make_mode2_op(A_blocks, B_blocks)
  e2 <- eigs_sym(op2$f, k = r2, n = op2$n, which = "LM", opts = list(tol = eigs_tol))
  U2 <- e2$vectors
  s2 <- sqrt(pmax(e2$values, 0))

  # Mode-3
  op3 <- make_mode3_op(A_blocks, B_blocks)
  e3 <- eigs_sym(op3$f, k = r3, n = op3$n, which = "LM", opts = list(tol = eigs_tol))
  U3 <- e3$vectors
  s3 <- sqrt(pmax(e3$values, 0))

  # コア G = X ×1 U1^T ×2 U2^T ×3 U3^T
  # X_{i,k,j} = A_{i,j} B_{k,j}
  # => jごとに rank-1 外積で加算
  G <- array(0, dim = c(r1, r2, r3))
  col_offset <- 0L

  for (b in seq_along(A_blocks)) {
    Ab <- A_blocks[[b]]
    Bb <- B_blocks[[b]]
    mb <- ncol(Ab)

    PA <- as.matrix(crossprod(U1, Ab))  # r1 x mb
    PB <- as.matrix(crossprod(U2, Bb))  # r2 x mb

    for (t in seq_len(mb)) {
      j_global <- col_offset + t
      cvec <- U3[j_global, ]            # length r3
      M12 <- PA[, t, drop = FALSE] %*% t(PB[, t, drop = FALSE])  # r1 x r2
      for (r in seq_len(r3)) {
        G[, , r] <- G[, , r] + M12 * cvec[r]
      }
    }
    col_offset <- col_offset + mb
  }

  list(
    U1 = U1, U2 = U2, U3 = U3,
    s1 = s1, s2 = s2, s3 = s3,
    core = G
  )
}