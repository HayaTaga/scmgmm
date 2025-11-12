#========================
# 比較手法の実装
#========================

# 便利関数：事前・事後のインデックス
.get_masks <- function(T, T0) {
  pre <- seq_len(T) <= T0
  post <- !pre
  list(pre = pre, post = post)
}

# 便利関数：事後平均のギャップから tau を集約
.aggregate_tau <- function(Y0, Yd, w) {
  as.numeric(Y0 - as.numeric(Yd %*% w))
}

#---------------------------------
# 1) 古典的 SCM（単純二乗誤差最小化, w>=0, sum w=1）
#---------------------------------
fit_scm_classic <- function(data, ridge = 0) {
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  J <- ncol(Yd)
  T0 <- data$T0
  ms <- .get_masks(T, T0)

  Yp <- Y0[ms$pre]
  Xp <- Yd[ms$pre, , drop = FALSE] # デザイン行列（事前のドナー）
  # QP: min_w ||Yp - Xp w||^2 + ridge*||w||^2  s.t. 1'w=1, w>=0
  D <- 2 * (crossprod(Xp) + ridge * diag(J))
  d <- 2 * crossprod(Xp, Yp)
  D <- (D + t(D)) / 2 + diag(1e-8, J) # 数値安定化

  Amat <- cbind(rep(1, J), diag(J)) # 列=制約, 先頭列が等式
  bvec <- c(1, rep(0, J))
  meq <- 1L

  sol <- tryCatch(
    quadprog::solve.QP(Dmat = D, dvec = d, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(sol)) {
    w <- rep(1 / J, J)
    conv <- FALSE
  } else {
    w <- pmax(sol$solution, 0)
    w <- w / sum(w)
    conv <- TRUE
  }

  gap <- .aggregate_tau(Y0, Yd, w)
  tau_hat <- mean(gap[ms$post])
  list(
    method = "SCM (classic)",
    w = as.numeric(w),
    tau_hat = tau_hat,
    gap = gap,
    converged = conv
  )
}

#---------------------------------
# 2) Ridge付き SCM（λ>0で安定化）
#---------------------------------
fit_scm_ridge <- function(data, lambda = 1e-2) {
  fit_scm_classic(data, ridge = lambda) # 上の関数を流用（目的関数にridgeを含める）
}

#---------------------------------
# 3) 単純 Difference-in-Differences（ドナーは等重み）
#---------------------------------
fit_did_simple <- function(data) {
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  T0 <- data$T0
  ms <- .get_masks(T, T0)

  donor_avg_pre <- colMeans(t(Yd[ms$pre, , drop = FALSE])) |> mean()
  donor_avg_post <- colMeans(t(Yd[ms$post, , drop = FALSE])) |> mean()

  tr_pre <- mean(Y0[ms$pre])
  tr_post <- mean(Y0[ms$post])

  tau_hat <- (tr_post - tr_pre) - (donor_avg_post - donor_avg_pre)

  # DiDは重みベクトルが不要なので、等重みを参考までに
  J <- ncol(Yd)
  w_eq <- rep(1 / J, J)
  gap <- .aggregate_tau(Y0, Yd, w_eq)

  list(
    method = "DiD (simple, donor unweighted mean)",
    w = w_eq,
    tau_hat = as.numeric(tau_hat),
    gap = gap,
    converged = TRUE
  )
}

#---------------------------------
# 4) 比較ラッパー（3手法を一括実行）
#---------------------------------
compare_methods <- function(data, lambda_ridge = 1e-2) {
  res1 <- fit_scm_classic(data)
  res2 <- fit_scm_ridge(data, lambda = lambda_ridge)
  res3 <- fit_did_simple(data)

  out <- list(
    classic_scm = res1,
    ridge_scm = res2,
    did_simple = res3
  )

  # 簡易サマリを表示
  cat("=== Comparison Summary ===\n")
  cat(sprintf(
    "[SCM   ] tau_hat = %.4f, converged = %s\n",
    res1$tau_hat,
    res1$converged
  ))
  cat(sprintf(
    "[Ridge ] tau_hat = %.4f, lambda = %.1e, converged = %s\n",
    res2$tau_hat,
    lambda_ridge,
    res2$converged
  ))
  cat(sprintf("[DiD   ] tau_hat = %.4f\n", res3$tau_hat))

  invisible(out)
}
