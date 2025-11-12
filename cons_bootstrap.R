# SCM-GMM with Constrained Multiplier Bootstrap (CNS 2023 full version only)
# Author: ChatGPT
# Depends: stats, sandwich, quadprog, purrr (optional), parallel (optional)

# ---- Libraries ----
# Explicit imports (safe if script-sourced)
if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("Package 'sandwich' required.")
}
if (!requireNamespace("quadprog", quietly = TRUE)) {
  stop("Package 'quadprog' required.")
}

# -----------------------------
# Data structure expected
# -----------------------------
# data <- list(
#   Y0 = numeric T vector for treated unit,
#   Ydonor = T x J matrix for donor units (columns are donors),
#   X = T x K matrix of instruments (may include constant),
#   T0 = integer, last pre-treatment index (1..T-1)
# )

# -----------------------------
# 1) Moment functions
# -----------------------------
calculate_moments <- function(theta, data) {
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  J <- ncol(Yd)
  X <- if (!is.null(data$X)) {
    as.matrix(data$X)
  } else {
    matrix(1, nrow = T, ncol = 1)
  }
  T0 <- data$T0
  w <- theta[seq_len(J)]
  tau <- theta[J + 1]
  synth <- as.numeric(Yd %*% w)
  pre_mask <- seq_len(T) <= T0
  post_mask <- seq_len(T) > T0
  g_pre <- X[pre_mask, , drop = FALSE] *
    as.numeric(Y0[pre_mask] - synth[pre_mask])
  g_post <- X[post_mask, , drop = FALSE] *
    as.numeric((Y0[post_mask] - synth[post_mask]) - tau)
  g <- rbind(g_pre, g_post)
  colMeans(g)
}

calculate_moment_contributions <- function(theta, data) {
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  J <- ncol(Yd)
  X <- if (!is.null(data$X)) {
    as.matrix(data$X)
  } else {
    matrix(1, nrow = T, ncol = 1)
  }
  T0 <- data$T0
  w <- theta[seq_len(J)]
  tau <- theta[J + 1]
  synth <- as.numeric(Yd %*% w)
  pre_mask <- seq_len(T) <= T0
  post_mask <- seq_len(T) > T0
  g_pre <- X[pre_mask, , drop = FALSE] *
    as.numeric(Y0[pre_mask] - synth[pre_mask])
  g_post <- X[post_mask, , drop = FALSE] *
    as.numeric((Y0[post_mask] - synth[post_mask]) - tau)
  g <- rbind(g_pre, g_post)
  g
}

calculate_hac_matrix <- function(
  g_t_matrix,
  kernel = c("Bartlett", "QS"),
  L = NULL,
  bw = NULL,
  prewhite = FALSE
) {
  # HAC estimator for moment contributions (matrix input)
  # - kernel: "Bartlett" (Newey-West) or "QS" (Quadratic Spectral)
  # - L: lag truncation for Bartlett (if NULL, uses Andrews-style plug-in)
  # - bw: bandwidth for QS (if NULL, uses b = 1.3221 * T^(1/5))
  # - prewhite: VAR(1) prewhitening with recoloring (Andrews–Monahan)
  Z <- as.matrix(g_t_matrix)
  Z <- scale(Z, center = TRUE, scale = FALSE)
  T <- nrow(Z)
  m <- ncol(Z)
  kern <- match.arg(kernel)

  # Optional VAR(1) prewhitening
  if (prewhite && T >= 3) {
    Zlag <- Z[1:(T - 1), , drop = FALSE]
    Zlead <- Z[2:T, , drop = FALSE]
    XtX <- crossprod(Zlag)
    XtY <- crossprod(Zlag, Zlead)
    A <- tryCatch(solve(XtX, XtY), error = function(e) {
      matrix(0, nrow = m, ncol = m)
    })
    E <- Zlead - Zlag %*% A
    Omega_e <- calculate_hac_matrix(
      E,
      kernel = kern,
      L = L,
      bw = bw,
      prewhite = FALSE
    )
    IA <- diag(m) - A
    IA_t <- diag(m) - t(A)
    Omega <- solve(IA_t, Omega_e) %*% solve(IA)
    return(Omega)
  }

  if (kern == "Bartlett") {
    if (is.null(L)) {
      L <- max(0L, floor(4 * (T / 100)^(2 / 9)))
    }
    Gamma0 <- crossprod(Z) / T
    Omega <- Gamma0
    if (L >= 1) {
      for (k in 1:L) {
        w <- 1 - k / (L + 1)
        Zlead <- Z[(k + 1):T, , drop = FALSE]
        Zlag <- Z[1:(T - k), , drop = FALSE]
        Gk <- crossprod(Zlead, Zlag) / T
        Omega <- Omega + w * (Gk + t(Gk))
      }
    }
    return(Omega)
  } else {
    if (is.null(bw)) {
      bw <- 1.3221 * T^(1 / 5)
    }
    Gamma0 <- crossprod(Z) / T
    Omega <- Gamma0
    for (k in 1:(T - 1)) {
      x <- k / bw
      if (x == 0) {
        w <- 1
      } else {
        w <- 25 /
          (12 * pi^2 * x^2) *
          (sin(6 * pi * x) / (6 * pi * x) - cos(6 * pi * x))
      }
      if (abs(w) < .Machine$double.eps) {
        break
      }
      Zlead <- Z[(k + 1):T, , drop = FALSE]
      Zlag <- Z[1:(T - k), , drop = FALSE]
      Gk <- crossprod(Zlead, Zlag) / T
      Omega <- Omega + w * (Gk + t(Gk))
    }
    return(Omega)
  }
}

.simplex_from_free <- function(u) {
  u_pos <- pmax(u, 0)
  last <- pmax(1 - sum(u_pos), 0)
  w <- c(u_pos, last)
  if (sum(w) == 0) {
    w[length(w)] <- 1
  } else {
    w <- w / sum(w)
  }
  w
}

compute_test_stat <- function(lambda, data, W_matrix) {
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  J <- ncol(Yd)
  X <- if (!is.null(data$X)) {
    as.matrix(data$X)
  } else {
    matrix(1, nrow = T, ncol = 1)
  }
  K <- ncol(X)
  T0 <- data$T0

  # m 次元（= K）に合わせて W を検査・調整
  tmp_theta <- c(rep(1 / J, J), lambda)
  m_dim <- ncol(calculate_moment_contributions(tmp_theta, data))
  if (!is.matrix(W_matrix) || any(dim(W_matrix) != c(m_dim, m_dim))) {
    g_t0 <- calculate_moment_contributions(tmp_theta, data)
    Omega <- calculate_hac_matrix(g_t0, kernel = "QS", prewhite = TRUE)
    eig <- eigen(Omega, symmetric = TRUE)
    tol <- 1e-8
    eig$values[eig$values < tol] <- tol
    W_matrix <- eig$vectors %*%
      diag(1 / eig$values, nrow = length(eig$values)) %*%
      t(eig$vectors)
  }

  # ḡ(w) = c - B w を構成
  pre <- seq_len(T) <= T0
  post <- !pre
  # c = 平均_t [ X_t * (Y0_t - 1{post}*lambda) ]  （K×1）
  c_vec <- colMeans(rbind(
    X[pre, , drop = FALSE] * Y0[pre],
    X[post, , drop = FALSE] * (Y0[post] - lambda)
  ))

  # B の各列 j は 平均_t [ X_t * Yd_{t,j} ]  （K×J）
  Bmat <- sapply(seq_len(J), function(j) colMeans(X * Yd[, j]))
  if (is.null(dim(Bmat))) {
    Bmat <- matrix(Bmat, nrow = K, ncol = J)
  }

  # 目的関数：min_w (c - B w)' W (c - B w)
  # solve.QP 形式：min_w 1/2 w' D w - d' w
  D <- 2 * t(Bmat) %*% W_matrix %*% Bmat
  d <- 2 * t(Bmat) %*% W_matrix %*% c_vec
  D <- (D + t(D)) / 2 + diag(1e-8, J) # 数値安定化の微小リッジ

  # 制約：sum(w)=1（等式）, w>=0（不等式）
  # quadprog は列が各制約、先頭 meq 列が等式
  Amat <- cbind(rep(1, J), diag(J))
  bvec <- c(1, rep(0, J))
  meq <- 1L

  sol <- tryCatch(
    quadprog::solve.QP(Dmat = D, dvec = d, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(sol)) {
    # フォールバック：等分重み
    w_hat <- rep(1 / J, J)
  } else {
    w_hat <- pmax(sol$solution, 0)
    sw <- sum(w_hat)
    if (sw <= 0) w_hat[] <- 1 / J else w_hat <- w_hat / sw
  }

  # 検定統計量
  gbar <- c_vec - Bmat %*% w_hat
  I_n <- as.numeric(t(gbar) %*% W_matrix %*% gbar)

  list(
    I_n = I_n,
    w_hat = as.numeric(w_hat),
    conv = !is.null(sol),
    lambda = lambda,
    W_matrix = W_matrix
  )
}

# -----------------------------
# 2) Full constrained bootstrap per CNS (2023)
# -----------------------------

.jacobian_moments <- function(theta, data) {
  Y0 <- as.numeric(data$Y0)
  Yd <- as.matrix(data$Ydonor)
  T <- NROW(Yd)
  J <- NCOL(Yd)
  X <- if (!is.null(data$X)) {
    as.matrix(data$X)
  } else {
    matrix(1, nrow = T, ncol = 1)
  }
  K <- NCOL(X)
  T0 <- data$T0

  Dw <- sapply(seq_len(J), function(j) colMeans(X * as.numeric(-Yd[, j])))
  if (is.null(dim(Dw))) {
    Dw <- matrix(Dw, nrow = K, ncol = J, byrow = FALSE)
  } else {
    Dw <- matrix(Dw, nrow = K, ncol = J)
  }

  Dtau <- colMeans(rbind(
    matrix(0, nrow = T0, ncol = K),
    -X[(T0 + 1):T, , drop = FALSE]
  ))
  Dtau <- matrix(as.numeric(Dtau), nrow = K, ncol = 1)

  # Jmat: K × (J+1)
  Jmat <- cbind(Dw, Dtau)
  dim(Jmat) <- c(K, J + 1)
  return(Jmat)
}

.build_local_constraints <- function(w_hat, T, r_n) {
  J <- length(w_hat)
  thr <- r_n / sqrt(T)

  Aeq <- matrix(0, nrow = 1, ncol = J + 1)
  Aeq[1, 1:J] <- 1
  beq <- 0

  active_idx <- which(!is.na(w_hat) & w_hat <= thr)
  if (length(active_idx) == 0L) {
    return(list(Aeq = Aeq, beq = beq, Aineq = NULL, bineq = NULL))
  }

  Aineq <- matrix(0, nrow = length(active_idx), ncol = J + 1)
  for (k in seq_along(active_idx)) {
    Aineq[k, active_idx[k]] <- 1
  }
  bineq <- -sqrt(T) * w_hat[active_idx]

  list(Aeq = Aeq, beq = beq, Aineq = Aineq, bineq = bineq)
}

.solve_qp_boot <- function(
  Jmat,
  W,
  Gstar,
  Aeq,
  beq,
  Aineq = NULL,
  bineq = NULL,
  ridge = 1e-8
) {
  if (!is.matrix(Jmat)) {
    Jmat <- as.matrix(Jmat)
  }
  if (!is.matrix(W)) {
    W <- as.matrix(W)
  }
  m <- nrow(Jmat)
  p <- ncol(Jmat)

  # 1) W のブロードキャスト & 形の統一
  if (nrow(W) == 1L && ncol(W) == 1L) {
    W <- as.numeric(W[1, 1]) * diag(m)
  }
  if (any(dim(W) != c(m, m))) {
    stop(sprintf(
      "Dimension mismatch: W is %s x %s, expected %s x %s",
      dim(W)[1],
      dim(W)[2],
      m,
      m
    ))
  }

  # 2) Gstar を m×1 の列ベクトルへ強制
  Gstar <- matrix(as.numeric(Gstar), nrow = m, ncol = 1)

  # 3) QP の各成分（crossprodは行列次元が明確になるよう t( ) を使用）
  H <- 2 * t(Jmat) %*% W %*% Jmat + diag(ridge, p)
  f <- as.numeric(2 * t(Jmat) %*% W %*% Gstar)

  # 4) 制約
  Amat <- t(Aeq)
  bvec <- beq
  meq <- nrow(Aeq)
  if (!is.null(Aineq)) {
    Amat <- cbind(Amat, t(Aineq))
    bvec <- c(bvec, bineq)
  }

  sol <- tryCatch(
    quadprog::solve.QP(Dmat = H, dvec = f, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(sol)) {
    return(list(value = NA_real_, par = rep(NA_real_, p)))
  }

  h <- sol$solution
  resid <- Gstar - Jmat %*% h # m×1
  val <- as.numeric(t(resid) %*% W %*% resid)
  list(value = val, par = h)
}

.compute_boot_stats_full <- function(
  lambda,
  data,
  W_matrix,
  w_hat,
  B = 999,
  rn = NULL,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  theta_hat <- c(w_hat, lambda)
  G <- calculate_moment_contributions(theta_hat, data) # T × m
  Gc <- scale(G, center = TRUE, scale = FALSE)
  T <- nrow(Gc)

  Jmat <- .jacobian_moments(theta_hat, data) # m × p
  m <- nrow(Jmat)

  # --- 重要ポイント：W を m×m に合わせて再構築 ---
  if (!is.matrix(W_matrix) || any(dim(W_matrix) != c(m, m))) {
    Omega <- calculate_hac_matrix(G, kernel = "QS", prewhite = TRUE)
    eig <- eigen(Omega, symmetric = TRUE)
    tol <- 1e-8
    eig$values[eig$values < tol] <- tol
    W_matrix <- eig$vectors %*%
      diag(1 / eig$values, nrow = length(eig$values)) %*%
      t(eig$vectors)
  }
  # ---------------------------------------------------

  if (is.null(rn)) {
    rn <- max(1, sqrt(log(max(T, 10))))
  }
  cons <- .build_local_constraints(w_hat, T, rn)

  boot_stats <- numeric(B)
  for (b in seq_len(B)) {
    v <- rnorm(T)
    G_star <- as.numeric(colSums(Gc * v) / sqrt(T)) # m 次元
    qp <- .solve_qp_boot(
      Jmat,
      W_matrix,
      G_star,
      Aeq = cons$Aeq,
      beq = cons$beq,
      Aineq = cons$Aineq,
      bineq = cons$bineq
    )
    boot_stats[b] <- T * qp$value
  }
  boot_stats
}

compute_critical_value <- function(
  lambda_info,
  data,
  W_matrix,
  B = 999,
  rn = NULL,
  seed = NULL,
  alpha = 0.05
) {
  boot_stats <- .compute_boot_stats_full(
    lambda = lambda_info$lambda,
    data = data,
    W_matrix = W_matrix,
    w_hat = lambda_info$w_hat,
    B = B,
    rn = rn,
    seed = seed
  )
  list(
    c_n = stats::quantile(boot_stats, probs = 1 - alpha, na.rm = TRUE),
    boot_stats = boot_stats
  )
}

test_point_scm <- function(
  lambda,
  data,
  W_matrix,
  B = 999,
  rn = NULL,
  seed = NULL,
  alpha = 0.05
) {
  info <- compute_test_stat(lambda, data, W_matrix)
  # pass through possibly recomputed W to bootstrap
  crit <- compute_critical_value(
    info,
    data,
    info$W_matrix,
    B = B,
    rn = rn,
    seed = seed,
    alpha = alpha
  )
  I_n_scaled <- info$I_n * nrow(data$Ydonor)
  pval <- mean(crit$boot_stats >= I_n_scaled, na.rm = TRUE)
  list(
    p_value = pval,
    I_n = info$I_n,
    I_n_scaled = I_n_scaled,
    w_hat = info$w_hat,
    boot = crit
  )
}

.find_W_matrix <- function(data, w_init = NULL, tau_init = 0) {
  Yd <- as.matrix(data$Ydonor)
  J <- ncol(Yd)
  if (is.null(w_init)) {
    w_init <- rep(1 / J, J)
  }
  theta1 <- c(w_init, tau_init)
  g_t <- calculate_moment_contributions(theta1, data)
  Omega <- calculate_hac_matrix(g_t, kernel = "QS", prewhite = TRUE)
  eig <- eigen(Omega, symmetric = TRUE)
  tol <- 1e-8
  eig$values[eig$values < tol] <- tol
  Omega_inv <- eig$vectors %*%
    diag(1 / eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)
  Omega_inv
}

find_ci_scm <- function(
  data,
  tau_grid,
  B = 999,
  W_matrix = NULL,
  rn = NULL,
  seed = NULL,
  alpha = 0.05,
  parallel = FALSE
) {
  if (is.null(W_matrix)) {
    W_matrix <- .find_W_matrix(data)
  }
  runner <- function(tau) {
    test_point_scm(
      tau,
      data,
      W_matrix,
      B = B,
      rn = rn,
      seed = seed,
      alpha = alpha
    )$p_value
  }
  if (parallel) {
    p_values <- parallel::mclapply(
      tau_grid,
      runner,
      mc.cores = max(1L, parallel::detectCores() - 1L)
    )
    p_values <- unlist(p_values)
  } else {
    p_values <- vapply(tau_grid, runner, numeric(1))
  }
  accepted <- tau_grid[p_values > alpha]
  ci <- if (length(accepted)) range(accepted) else c(NA_real_, NA_real_)
  list(ci = ci, grid = tau_grid, p_values = p_values, alpha = alpha)
}

make_scm_data <- function(Y0, Ydonor, T0, X = NULL) {
  stopifnot(length(Y0) == nrow(Ydonor))
  if (is.null(X)) {
    X <- matrix(1, nrow = length(Y0), ncol = 1)
  }
  list(
    Y0 = as.numeric(Y0),
    Ydonor = as.matrix(Ydonor),
    X = as.matrix(X),
    T0 = as.integer(T0)
  )
}

# Example usage (pseudo)
set.seed(1)
T <- 100
T0 <- 70
J <- 10
Ydonor <- matrix(rnorm(T * J), T, J)
w_star <- runif(J)
w_star <- w_star / sum(w_star)
Y0 <- as.numeric(Ydonor %*% w_star + 0.2 * rnorm(T))
tau_true <- 1
Y0[(T0 + 1):T] <- Y0[(T0 + 1):T] + tau_true
data <- make_scm_data(Y0, Ydonor, T0)
tau_grid <- seq(0, 2, by = 0.1)
res <- find_ci_scm(data, tau_grid, B = 999)
print(res$ci)

source("scm_method.R")
set.seed(1)
T <- 100
T0 <- 70
J <- 10
Ydonor <- matrix(rnorm(T * J), T, J)
w_star <- runif(J)
w_star <- w_star / sum(w_star)
Y0 <- as.numeric(Ydonor %*% w_star + 0.2 * rnorm(T))
tau_true <- 1
Y0[(T0 + 1):T] <- Y0[(T0 + 1):T] + tau_true

data <- make_scm_data(Y0, Ydonor, T0)

cmp <- compare_methods(data, lambda_ridge = 1e-2)
# それぞれの推定値
cmp$classic_scm$tau_hat
cmp$ridge_scm$tau_hat
cmp$did_simple$tau_hat

source("comp_infer.R")
# 既に data がある前提
scm1 <- fit_scm_classic(data)
scm2 <- fit_scm_ridge(data, lambda = 1e-2)

# 1) HAC CI
infer_scm_hac_ci(scm1, data)
infer_scm_hac_ci(scm2, data, kernel = "QS", prewhite = TRUE)

# 2) placebo-in-space（preMSPEフィルタ10倍）
infer_scm_placebo_space(data, estimator = "classic", pre_mspe_ratio_max = 10)

# 3) placebo-in-time
infer_scm_placebo_time(data, estimator = "ridge")

# 4) DiD の HAC CI
infer_did_hac_ci(data)


source("montecarlo_sim.R")
set.seed(42)
grid_small <- expand.grid(
  T = c(60, 100),
  J = c(10),
  r = 1,
  rho = c(0, 0.5),
  sigma = 0.2,
  tau = c(0, 1)
)
res_small <- run_mc(nrep = 20, grid = grid_small, B_boot = 499, seed0 = 42)
tab_small <- summarise_mc(res_small)
print(
  tab_small[order(tab_small$tau, tab_small$rho, tab_small$T), ],
  row.names = FALSE
)

# 本番（例：各シナリオ 200反復、B=999）
# res <- run_mc(nrep = 200, grid = your_grid, B_boot = 999, seed0 = 2025)
# tab <- summarise_mc(res)
# write.csv(tab, "mc_summary.csv", row.names = FALSE)
