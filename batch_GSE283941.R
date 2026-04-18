source("HOSVD.R")
x <- readMM("GSE283941_matrix_protein.mtx.gz")
x1 <- readMM("GSE283941_matrix_rna.mtx.gz")
fit <- hosvd_rank_trunc(
  x, x1,
  r1 = 10, r2 = 10, r3 = 10,
  block_size = 1000L,
  eigs_tol = 1e-5
)
save(file="fit_GSE283941",fit)
