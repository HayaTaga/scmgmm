#========================
# 推論ユーティリティ
#========================

# 既存の calculate_hac_matrix() を単変量にラップ
.hac_var_post_mean <- function(
  gap,
  T0,
  kernel = "Bartlett",
  L = NULL,
  bw = NULL,
  prewhite = FALSE
) {
  gap <- as.numeric(gap)
  post <- seq_along(gap) > T0
  z <- matrix(gap[post], ncol = 1) # T_post × 1
  Omega <- calculate_hac_matrix(
    z,
    kernel = kernel,
    L = L,
    bw = bw,
    prewhite = prewhite
  ) # 1×1
  Tpost <- sum(post)
  as.numeric(Omega) / (Tpost^2) # Var( \bar g_post )
}

# ギャップと tau_hat を共通形式で
.scm_gap_and_tau <- function(data, w) {
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  gap <- as.numeric(Y0 - as.numeric(Yd %*% w))
  T0 <- data$T0
  tau_hat <- mean(gap[(T0 + 1):length(gap)])
  list(gap = gap, tau_hat = tau_hat)
}

# 事前・事後のMSPE（均二乗誤差）
.mspe_pre_post <- function(gap, T0) {
  pre <- seq_along(gap) <= T0
  post <- !pre
  c(pre = mean(gap[pre]^2), post = mean(gap[post]^2))
}

#========================
# 1) SCM: HAC CI（classic / ridge 共通）
#========================
infer_scm_hac_ci <- function(
  fit,
  data,
  alpha = 0.05,
  kernel = "Bartlett",
  L = NULL,
  bw = NULL,
  prewhite = FALSE
) {
  stopifnot(!is.null(fit$w))
  gt <- .scm_gap_and_tau(data, fit$w)
  var_tau <- .hac_var_post_mean(gt$gap, data$T0, kernel, L, bw, prewhite)
  se_tau <- sqrt(var_tau)
  z <- qnorm(1 - alpha / 2)
  ci <- c(gt$tau_hat - z * se_tau, gt$tau_hat + z * se_tau)
  list(
    method = "SCM-HAC (normal approx)",
    tau_hat = gt$tau_hat,
    se = se_tau,
    ci = ci,
    details = list(kernel = kernel, L = L, bw = bw, prewhite = prewhite)
  )
}

#========================
# 2) SCM: placebo-in-space（Abadie型, MSPE比, フィルタあり）
#========================
infer_scm_placebo_space <- function(
  data,
  estimator = c("classic", "ridge"),
  lambda_ridge = 1e-2,
  pre_mspe_ratio_max = 10, # フィルタ閾値（treatedに対する倍率）
  return_table = FALSE
) {
  estimator <- match.arg(estimator)
  # 処置ユニットの推定
  base_fit <- if (estimator == "classic") {
    fit_scm_classic(data)
  } else {
    fit_scm_ridge(data, lambda = lambda_ridge)
  }
  base_gap <- .scm_gap_and_tau(data, base_fit$w)$gap
  base_mspe <- .mspe_pre_post(base_gap, data$T0)
  base_ratio <- base_mspe["post"] / base_mspe["pre"]

  # 各ドナー i を擬似処置に（i を除いた残りをドナーに）
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  J <- ncol(Yd)
  T0 <- data$T0

  ratios <- rep(NA_real_, J)
  keep <- rep(FALSE, J)

  for (i in seq_len(J)) {
    # 擬似処置シリーズ：元の i 列を treated に、残り列を donor に
    Y0_i <- Yd[, i]
    Yd_i <- Yd[, -i, drop = FALSE]
    data_i <- make_scm_data(Y0_i, Yd_i, T0)

    fit_i <- if (estimator == "classic") {
      fit_scm_classic(data_i)
    } else {
      fit_scm_ridge(data_i, lambda = lambda_ridge)
    }
    gap_i <- .scm_gap_and_tau(data_i, fit_i$w)$gap
    mspe_i <- .mspe_pre_post(gap_i, T0)
    ratio_i <- mspe_i["post"] / mspe_i["pre"]

    ratios[i] <- ratio_i
    # 事前適合が悪すぎる placebo をフィルタ
    keep[i] <- (mspe_i["pre"] <= pre_mspe_ratio_max * base_mspe["pre"])
  }

  # フィルタ後の分布で one-sided p値（treated より大きい比がどれだけあるか）
  pool <- ratios[keep & is.finite(ratios)]
  pval <- if (length(pool)) {
    (sum(pool >= base_ratio) + 1) / (length(pool) + 1)
  } else {
    NA_real_
  }

  out <- list(
    method = sprintf("SCM placebo-in-space (%s)", estimator),
    treated_ratio = as.numeric(base_ratio),
    placebo_ratios = pool,
    p_value = pval,
    pre_filter_kept = sum(keep),
    pre_filter_total = J
  )
  if (return_table) {
    out$raw <- data.frame(idx = which(keep), ratio = pool)
  }
  out
}

#========================
# 3) SCM: placebo-in-time（簡易版, MSPE比）
#========================
infer_scm_placebo_time <- function(
  data,
  estimator = c("classic", "ridge"),
  lambda_ridge = 1e-2,
  placebo_T0 = NULL # 例:  floor(T0*0.6):(T0-5)
) {
  estimator <- match.arg(estimator)
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  T0 <- data$T0

  # デフォルトの疑似時点：事前の後半（十分な事後長を残す）
  if (is.null(placebo_T0)) {
    placebo_T0 <- seq(max(5, floor(T0 * 0.6)), T0 - 5)
  }

  # 実データの基準（本物の T0）
  base_fit <- if (estimator == "classic") {
    fit_scm_classic(data)
  } else {
    fit_scm_ridge(data, lambda = lambda_ridge)
  }
  base_gap <- .scm_gap_and_tau(data, base_fit$w)$gap
  base_ratio <- .mspe_pre_post(base_gap, T0)["post"] /
    .mspe_pre_post(base_gap, T0)["pre"]

  ratios <- numeric(length(placebo_T0))
  for (k in seq_along(placebo_T0)) {
    t0_k <- placebo_T0[k]
    data_k <- make_scm_data(Y0, Yd, t0_k) # 処置時点だけ入れ替え
    fit_k <- if (estimator == "classic") {
      fit_scm_classic(data_k)
    } else {
      fit_scm_ridge(data_k, lambda = lambda_ridge)
    }
    gap_k <- .scm_gap_and_tau(data_k, fit_k$w)$gap
    mspe_k <- .mspe_pre_post(gap_k, t0_k)
    ratios[k] <- mspe_k["post"] / mspe_k["pre"]
  }

  pval <- (sum(ratios >= base_ratio) + 1) / (length(ratios) + 1)

  list(
    method = sprintf("SCM placebo-in-time (%s)", estimator),
    treated_ratio = as.numeric(base_ratio),
    placebo_ratios = ratios,
    placebo_T0 = placebo_T0,
    p_value = pval
  )
}

#========================
# 4) DiD: Newey–West(HAC) CI
#========================
infer_did_hac_ci <- function(data, alpha = 0.05, L = NULL) {
  # fit_did_simple() を想定
  fit <- fit_did_simple(data)
  # DiD 推定量を「時点ごとの差の差」系列の平均として捉えて HAC
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T0 <- data$T0
  pre <- seq_len(nrow(Yd)) <= T0
  post <- !pre
  donor_mean <- rowMeans(Yd)
  dd_series <- (Y0 - donor_mean) # 各時点の treated - control 平均
  tau_hat <- mean(dd_series[post])

  # 事後平均の HAC 分散
  var_tau <- .hac_var_post_mean(
    dd_series,
    T0,
    kernel = "Bartlett",
    L = L,
    prewhite = FALSE
  )
  se_tau <- sqrt(var_tau)
  z <- qnorm(1 - alpha / 2)
  ci <- c(tau_hat - z * se_tau, tau_hat + z * se_tau)

  list(
    method = "DiD (HAC)",
    tau_hat = tau_hat,
    se = se_tau,
    ci = ci,
    details = list(kernel = "Bartlett", L = L)
  )
}
