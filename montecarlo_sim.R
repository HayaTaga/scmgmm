# =========================================================
# 1) DGP: 因子構造 + トレンド + AR(1) 誤差 + 介入
# =========================================================
dgp_panel <- function(
  T,
  T0,
  J,
  r = 1,
  rho = 0,
  sigma = 0.2,
  tau = 1,
  trend_sd = 0.02,
  load_sd = 1,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # 因子とロード
  Ft <- matrix(rnorm(T * r), T, r)
  Gamma <- matrix(rnorm((J + 1) * r, sd = load_sd), J + 1, r) # treated + J donors
  beta <- rnorm(J + 1, sd = trend_sd) # 個別トレンド
  mu <- rnorm(J + 1, sd = 0.5)

  # AR(1) 誤差
  U <- matrix(0, T, J + 1)
  eps <- matrix(rnorm(T * (J + 1), sd = sigma), T, J + 1)
  for (t in 2:T) {
    U[t, ] <- rho * U[t - 1, ] + eps[t, ]
  }
  U[1, ] <- eps[1, ]

  # 系列生成（treated=1列目）
  tvec <- 1:T
  Y <- outer(tvec, beta) +
    matrix(mu, T, J + 1, byrow = TRUE) +
    Ft %*% t(Gamma) +
    U
  # 介入効果
  Y[(T0 + 1):T, 1] <- Y[(T0 + 1):T, 1] + tau

  # 出力（SCMのI/Oに合わせる）
  list(
    Y0 = as.numeric(Y[, 1]),
    Ydonor = as.matrix(Y[, -1, drop = FALSE]),
    T0 = as.integer(T0)
  )
}

# =========================================================
# 2) 1回の実験を走らせるラッパ（各手法のCI/推定を集計）
#    依拠: 既存の関数群
# =========================================================
run_once <- function(
  T,
  J,
  r,
  rho,
  sigma,
  tau,
  seed = NULL,
  grid_tau = seq(max(0, tau - 2), tau + 2, by = 0.1),
  B_boot = 999
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  T0 <- floor(0.7 * T)
  d <- dgp_panel(
    T,
    T0,
    J,
    r = r,
    rho = rho,
    sigma = sigma,
    tau = tau,
    seed = seed
  )
  data <- make_scm_data(d$Y0, d$Ydonor, T0)

  # ---- 提案法：CNS（乗数ブートCI）----
  t0 <- proc.time()[3]
  res_cns <- try(
    find_ci_scm(data, grid_tau, B = B_boot, seed = seed),
    silent = TRUE
  )
  t1 <- proc.time()[3]
  ci_cns <- if (inherits(res_cns, "try-error")) c(NA, NA) else res_cns$ci
  cover_cns <- !any(is.na(ci_cns)) && (ci_cns[1] <= tau && tau <= ci_cns[2])
  len_cns <- if (any(is.na(ci_cns))) NA_real_ else diff(ci_cns)

  # ---- SCM-HAC（Bartlett）----
  f_sc <- fit_scm_classic(data)
  inf_sc_b <- infer_scm_hac_ci(f_sc, data, kernel = "Bartlett")
  ci_sc_b <- inf_sc_b$ci
  cover_sc_b <- (ci_sc_b[1] <= tau && tau <= ci_sc_b[2])
  len_sc_b <- diff(ci_sc_b)

  # ---- SCM-HAC（QS+prewhite）----
  inf_sc_qs <- infer_scm_hac_ci(f_sc, data, kernel = "QS", prewhite = TRUE)
  ci_sc_qs <- inf_sc_qs$ci
  cover_sc_qs <- (ci_sc_qs[1] <= tau && tau <= ci_sc_qs[2])
  len_sc_qs <- diff(ci_sc_qs)

  # ---- Ridge–SCM–HAC（QS+prewhite）----
  f_rg <- fit_scm_ridge(data, lambda = 1e-2)
  inf_rg <- infer_scm_hac_ci(f_rg, data, kernel = "QS", prewhite = TRUE)
  ci_rg <- inf_rg$ci
  cover_rg <- (ci_rg[1] <= tau && tau <= ci_rg[2])
  len_rg <- diff(ci_rg)

  # ---- DiD–HAC（Bartlett）----
  inf_did <- infer_did_hac_ci(data, alpha = 0.05, L = NULL)
  ci_did <- inf_did$ci
  cover_did <- (ci_did[1] <= tau && tau <= ci_did[2])
  len_did <- diff(ci_did)

  data.frame(
    T = T,
    J = J,
    r = r,
    rho = rho,
    sigma = sigma,
    tau = tau,
    cover_CNS = as.numeric(cover_cns),
    len_CNS = len_cns,
    time_CNS = t1 - t0,
    cover_SCM_B = as.numeric(cover_sc_b),
    len_SCM_B = len_sc_b,
    cover_SCM_QS = as.numeric(cover_sc_qs),
    len_SCM_QS = len_sc_qs,
    cover_Ridge_QS = as.numeric(cover_rg),
    len_Ridge_QS = len_rg,
    cover_DiD = as.numeric(cover_did),
    len_DiD = len_did
  )
}

# =========================================================
# 3) モンテカルロ（パラメータ格子×反復）
# =========================================================
run_mc <- function(
  nrep = 200,
  grid = expand.grid(
    T = c(60, 100),
    J = c(10, 25),
    r = c(1, 2),
    rho = c(0, 0.5),
    sigma = c(0.1, 0.2),
    tau = c(0, 1)
  ),
  B_boot = 999,
  seed0 = 123
) {
  out <- vector("list", nrow(grid) * nrep)
  k <- 1L
  for (i in seq_len(nrow(grid))) {
    gi <- grid[i, ]
    for (rep in seq_len(nrep)) {
      s <- seed0 + i * 10000 + rep
      out[[k]] <- run_once(
        gi$T,
        gi$J,
        gi$r,
        gi$rho,
        gi$sigma,
        gi$tau,
        seed = s,
        B_boot = B_boot
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

# =========================================================
# 4) 要約（被覆率・平均長・計算時間）
# =========================================================
summarise_mc <- function(df) {
  agg <- within(
    aggregate(
      . ~ T + J + r + rho + sigma + tau,
      data = df,
      FUN = function(x) {
        if (is.numeric(x)) c(mean = mean(x, na.rm = TRUE)) else NA
      }
    ),
    NULL
  )

  # ネストを平坦化
  unpack <- function(col) {
    sapply(col, function(x) if (is.null(x)) NA_real_ else x[1])
  }
  agg$cover_CNS <- unpack(agg$cover_CNS)
  agg$len_CNS <- unpack(agg$len_CNS)
  agg$time_CNS <- unpack(agg$time_CNS)
  agg$cover_SCM_B <- unpack(agg$cover_SCM_B)
  agg$len_SCM_B <- unpack(agg$len_SCM_B)
  agg$cover_SCM_QS <- unpack(agg$cover_SCM_QS)
  agg$len_SCM_QS <- unpack(agg$len_SCM_QS)
  agg$cover_Ridge_QS <- unpack(agg$cover_Ridge_QS)
  agg$len_Ridge_QS <- unpack(agg$len_Ridge_QS)
  agg$cover_DiD <- unpack(agg$cover_DiD)
  agg$len_DiD <- unpack(agg$len_DiD)
  agg
}
