# =========================================================
# SCM-GMM (CNS) + AE(2023)-style Rn calibration & test inversion
# =========================================================
# Depends: stats, sandwich, quadprog, (optional) parallel, (optional) lpSolve

# ---- Libraries ----
if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("Package 'sandwich' required.")
}
if (!requireNamespace("quadprog", quietly = TRUE)) {
  stop("Package 'quadprog' required.")
}

# -----------------------------
# Utilities (AE-style)
# -----------------------------
sqrtm_sym <- function(M) {
  M <- as.matrix(M)
  E <- eigen((M + t(M)) / 2, symmetric = TRUE)
  vals <- pmax(E$values, 0)
  E$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(E$vectors)
}

# Rn calibration for moments (identity constraints on moments)
# We approximate AE's sup-stat calibration on studentized inequalities
# by using the moments' HAC covariance Omega and simulating max |U_i|
# where U ~ N(0, Omega). Return quantile pRn of sup-norm.
select_Rn_moments <- function(Omega, B = 500, pRn = 0.95, seed = 210) {
  set.seed(seed)
  m <- ncol(Omega)
  Sroot <- sqrtm_sym(Omega)
  draws <- replicate(B, {
    z <- rnorm(m)
    u <- Sroot %*% z
    max(abs(u)) # sup-norm (two-sided). If one-sided is desired, use max(u)
  })
  max(stats::quantile(draws, pRn, na.rm = TRUE), 0)
}

# -----------------------------
# Data structure helper
# -----------------------------
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
  pre <- seq_len(T) <= T0
  post <- !pre
  g_pre <- X[pre, , drop = FALSE] * as.numeric(Y0[pre] - synth[pre])
  g_post <- X[post, , drop = FALSE] * as.numeric((Y0[post] - synth[post]) - tau)
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
  pre <- seq_len(T) <= T0
  post <- !pre
  g_pre <- X[pre, , drop = FALSE] * as.numeric(Y0[pre] - synth[pre])
  g_post <- X[post, , drop = FALSE] * as.numeric((Y0[post] - synth[post]) - tau)
  rbind(g_pre, g_post) # T × K
}

calculate_hac_matrix <- function(
  g_t_matrix,
  kernel = c("Bartlett", "QS"),
  L = NULL,
  bw = NULL,
  prewhite = FALSE
) {
  Z <- as.matrix(g_t_matrix)
  Z <- scale(Z, center = TRUE, scale = FALSE)
  T <- nrow(Z)
  m <- ncol(Z)
  kern <- match.arg(kernel)

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
    return(solve(IA_t, Omega_e) %*% solve(IA))
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
      w <- if (x == 0) {
        1
      } else {
        25 /
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

# -----------------------------
# 2) Test statistic (outer problem) with QP for w (sum=1, w>=0)
# -----------------------------
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

  # Moment covariance -> W check
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

  # Build c and B: ḡ(w) = c - B w  (K×1, K×J)
  pre <- seq_len(T) <= T0
  post <- !pre
  c_vec <- colMeans(rbind(
    X[pre, , drop = FALSE] * Y0[pre],
    X[post, , drop = FALSE] * (Y0[post] - lambda)
  ))
  Bmat <- sapply(seq_len(J), function(j) colMeans(X * Yd[, j]))
  if (is.null(dim(Bmat))) {
    Bmat <- matrix(Bmat, nrow = K, ncol = J)
  }

  D <- 2 * t(Bmat) %*% W_matrix %*% Bmat
  d <- 2 * t(Bmat) %*% W_matrix %*% c_vec
  D <- (D + t(D)) / 2 + diag(1e-8, J)

  Amat <- cbind(rep(1, J), diag(J)) # eq first col, then inequalities
  bvec <- c(1, rep(0, J))
  meq <- 1L

  sol <- tryCatch(
    quadprog::solve.QP(Dmat = D, dvec = d, Amat = Amat, bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  w_hat <- if (is.null(sol)) {
    rep(1 / J, J)
  } else {
    ww <- pmax(sol$solution, 0)
    sw <- sum(ww)
    if (sw <= 0) rep(1 / J, J) else ww / sw
  }
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
# 3) Inner bootstrap QP (CNS)
# -----------------------------
.jacobian_moments <- function(theta, data) {
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

  Dw <- sapply(seq_len(J), function(j) colMeans(X * as.numeric(-Yd[, j])))
  if (is.null(dim(Dw))) {
    Dw <- matrix(Dw, nrow = K, ncol = J)
  } else {
    Dw <- matrix(Dw, nrow = K, ncol = J)
  }
  Dtau <- colMeans(rbind(
    matrix(0, nrow = T0, ncol = K),
    -X[(T0 + 1):T, , drop = FALSE]
  ))
  Dtau <- matrix(as.numeric(Dtau), nrow = K, ncol = 1)
  Jmat <- cbind(Dw, Dtau)
  dim(Jmat) <- c(K, J + 1)
  Jmat
}

.build_local_constraints <- function(w_hat, T, r_n) {
  # AE風: r_n を校正。閾値は r_n / sqrt(T)
  J <- length(w_hat)
  thr <- r_n / sqrt(T)
  Aeq <- matrix(0, nrow = 1, ncol = J + 1)
  Aeq[1, 1:J] <- 1
  beq <- 0
  active <- which(!is.na(w_hat) & w_hat <= thr)
  if (length(active) == 0L) {
    return(list(Aeq = Aeq, beq = beq, Aineq = NULL, bineq = NULL))
  }
  Aineq <- matrix(0, nrow = length(active), ncol = J + 1)
  for (k in seq_along(active)) {
    Aineq[k, active[k]] <- 1
  }
  bineq <- -sqrt(T) * w_hat[active] # local h_j >= -sqrt(T) * w_hat_j  (CNSのローカル近似)
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
  Jmat <- as.matrix(Jmat)
  W <- as.matrix(W)
  m <- nrow(Jmat)
  p <- ncol(Jmat)
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
  Gstar <- matrix(as.numeric(Gstar), nrow = m, ncol = 1)
  H <- 2 * t(Jmat) %*% W %*% Jmat + diag(ridge, p)
  f <- as.numeric(2 * t(Jmat) %*% W %*% Gstar)
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
  resid <- Gstar - Jmat %*% h
  val <- as.numeric(t(resid) %*% W %*% resid)
  list(value = val, par = h)
}

.compute_boot_stats_full <- function(
  lambda,
  data,
  W_matrix,
  w_hat,
  B = 999,
  r_n = NULL,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  theta_hat <- c(w_hat, lambda)
  G <- calculate_moment_contributions(theta_hat, data)
  Gc <- scale(G, center = TRUE, scale = FALSE)
  T <- nrow(Gc)
  Jmat <- .jacobian_moments(theta_hat, data) # m × p
  m <- nrow(Jmat)

  # W refresh to m×m if needed
  if (!is.matrix(W_matrix) || any(dim(W_matrix) != c(m, m))) {
    Omega <- calculate_hac_matrix(G, kernel = "QS", prewhite = TRUE)
    eig <- eigen(Omega, symmetric = TRUE)
    tol <- 1e-8
    eig$values[eig$values < tol] <- tol
    W_matrix <- eig$vectors %*%
      diag(1 / eig$values, nrow = length(eig$values)) %*%
      t(eig$vectors)
  }

  if (is.null(r_n)) {
    r_n <- max(1, sqrt(log(max(T, 10))))
  } # fallback if Rn not supplied
  cons <- .build_local_constraints(w_hat, T, r_n)

  boot_stats <- numeric(B)
  for (b in seq_len(B)) {
    v <- rnorm(T)
    G_star <- as.numeric(colSums(Gc * v) / sqrt(T))
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

# -----------------------------
# 4) Critical value via bootstrap & TestPoint (p-value)
# -----------------------------
.compute_Rn_from_data <- function(
  data,
  w_for_Rn = NULL,
  B = 500,
  pRn = 0.95,
  seed = 210
) {
  # Rn calibrates on moments' covariance Omega
  Yd <- as.matrix(data$Ydonor)
  T <- nrow(Yd)
  J <- ncol(Yd)
  if (is.null(w_for_Rn)) {
    w_for_Rn <- rep(1 / J, J)
  }
  theta <- c(w_for_Rn, 0)
  G <- calculate_moment_contributions(theta, data)
  Omega <- calculate_hac_matrix(G, kernel = "QS", prewhite = TRUE)
  select_Rn_moments(Omega, B = B, pRn = pRn, seed = seed)
}

compute_critical_value <- function(
  lambda_info,
  data,
  W_matrix,
  B = 999,
  r_n = NULL,
  seed = NULL,
  alpha = 0.05
) {
  boot_stats <- .compute_boot_stats_full(
    lambda = lambda_info$lambda,
    data = data,
    W_matrix = W_matrix,
    w_hat = lambda_info$w_hat,
    B = B,
    r_n = r_n,
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
  r_n = NULL,
  seed = NULL,
  alpha = 0.05
) {
  info <- compute_test_stat(lambda, data, W_matrix)
  crit <- compute_critical_value(
    info,
    data,
    info$W_matrix,
    B = B,
    r_n = r_n,
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

# -----------------------------
# 5) AE-style: grid + inversion with uniroot (find CI)
#     - pRn から Rn を校正して r_n に投入
# -----------------------------
find_ci_scm <- function(
  data,
  tau_grid,
  B = 999,
  W_matrix = NULL,
  r_n = NULL,
  pRn = 0.95,
  B_Rn = 500,
  seed = NULL,
  alpha = 0.05,
  parallel = FALSE,
  refine = TRUE,
  tol = 1e-3
) {
  if (is.null(W_matrix)) {
    W_matrix <- .find_W_matrix(data)
  }
  # AE-style Rn calibration (if r_n not supplied)
  if (is.null(r_n)) {
    Rn <- .compute_Rn_from_data(
      data,
      w_for_Rn = NULL,
      B = B_Rn,
      pRn = pRn,
      seed = if (is.null(seed)) 210 else seed
    )
    r_n <- max(1e-8, Rn) # feed into local constraints threshold
  }

  runner <- function(tau) {
    test_point_scm(
      tau,
      data,
      W_matrix,
      B = B,
      r_n = r_n,
      seed = seed,
      alpha = alpha
    )$p_value
  }
  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
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
  ci_grid <- if (length(accepted)) range(accepted) else c(NA_real_, NA_real_)
  if (!refine || any(is.na(ci_grid))) {
    return(list(
      ci = ci_grid,
      grid = tau_grid,
      p_values = p_values,
      alpha = alpha,
      r_n = r_n
    ))
  }

  # --- AEの fzero に相当：端点付近を uniroot で精密化 ---
  f_p <- function(x) runner(x) - (1 - (1 - alpha)) # = p(x) - alpha
  # 下端
  lower <- ci_grid[1]
  upper <- ci_grid[2]
  # bracket for lower endpoint
  iL <- max(which(tau_grid <= lower))
  iL <- max(1, iL - 1)
  aL <- tau_grid[iL]
  bL <- tau_grid[min(length(tau_grid), iL + 1)]
  if (is.na(f_p(aL)) || is.na(f_p(bL))) {
    ci_low <- lower
  } else {
    if (f_p(aL) * f_p(bL) > 0) {
      # 少し外側に広げる（AEの while に相当）
      step <- max(0.2 * (bL - aL), 1e-2)
      repeat {
        aL <- aL - step
        if (aL < min(tau_grid)) {
          break
        }
        if (!is.na(f_p(aL)) && f_p(aL) * f_p(bL) <= 0) break
      }
    }
    rootL <- try(
      uniroot(f_p, lower = min(aL, bL), upper = max(aL, bL), tol = tol),
      silent = TRUE
    )
    ci_low <- if (inherits(rootL, "try-error")) lower else rootL$root
  }
  # bracket for upper endpoint
  iU <- min(which(tau_grid >= upper))
  iU <- min(length(tau_grid) - 1, max(1, iU))
  aU <- tau_grid[iU]
  bU <- tau_grid[iU + 1]
  if (is.na(f_p(aU)) || is.na(f_p(bU))) {
    ci_up <- upper
  } else {
    if (f_p(aU) * f_p(bU) > 0) {
      step <- max(0.2 * (bU - aU), 1e-2)
      repeat {
        bU <- bU + step
        if (bU > max(tau_grid)) {
          break
        }
        if (!is.na(f_p(aU)) && f_p(aU) * f_p(bU) <= 0) break
      }
    }
    rootU <- try(
      uniroot(f_p, lower = min(aU, bU), upper = max(aU, bU), tol = tol),
      silent = TRUE
    )
    ci_up <- if (inherits(rootU, "try-error")) upper else rootU$root
  }

  list(
    ci = c(ci_low, ci_up),
    grid = tau_grid,
    p_values = p_values,
    alpha = alpha,
    r_n = r_n
  )
}

# -----------------------------
# 6) Weighting matrix from HAC (same as before)
# -----------------------------
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
  eig$vectors %*%
    diag(1 / eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)
}

# -----------------------------
# 7) Example (smoke test)
# -----------------------------
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
res <- find_ci_scm(data, tau_grid, B = 999, pRn = 0.95, B_Rn = 500, seed = 1)
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
