# ==============================================================================
# DynaPrism_instaprism.R
#
# InstaPrism-style (derandomized, fixed-point) engine for DynaPrism.
#
# WHAT THIS IS
# ------------
# BayesPrism/DynaPrism infer theta by Gibbs sampling: each iteration DRAWS
#   Z_g   ~ Multinomial(X_g, prob = phi_g o theta o exp(beta * niche(theta)))
#   theta ~ Dirichlet(colSums(Z) + alpha)
# and averages the post-burn-in draws. InstaPrism (Hu & Chikina 2024) replaces
# both draws with their conditional means, turning the sampler into a
# deterministic fixed-point iteration:
#   Z_g   <- X_g * B_g                                   (B_g = column-normalized prob)
#   theta <- (colSums(Z) + alpha) / sum(colSums(Z) + alpha)
# iterated to convergence. This removes the chain-length x per-gene sampling loop
# and vectorizes the read assignment: typically 1-2 orders of magnitude faster at
# effectively identical results FOR THE VANILLA MODEL. The niche factor
# exp(beta * niche(theta)) is carried through unchanged, using the current theta
# at each step exactly as sample.Z.theta_n.Omega() / sample.theta_n.Omega() do.
#
# ------------------------------------------------------------------------------
# IMPORTANT CAVEAT (read before trusting these numbers in the thesis)
# ------------------------------------------------------------------------------
# For vanilla BayesPrism the derandomization is exact IN EXPECTATION. DynaPrism
# adds a NONLINEAR, theta-dependent term exp(beta * (theta . Omega)). The Gibbs
# sampler averages exp(beta * niche(theta^(i))) over sampled theta^(i); this
# engine evaluates exp(beta * niche(theta*)) at the fixed-point mean. Because exp
# is convex, by Jensen's inequality E[exp(beta*niche)] >= exp(beta*niche(E)),
# so this is a MEAN-FIELD APPROXIMATION of DynaPrism, not an exact
# derandomization. The gap grows with beta and with the posterior spread of theta
# (larger for low-count / collinear states) -- exactly the beta ~ 5 regime you
# operate in. Validate ip-vs-Gibbs theta / Z / ExpSpe across a beta sweep before
# using ip numbers in a claim; use relax < 1 if the fixed point oscillates.
#
# SAMPLING CONTROL
# ----------------
# This engine uses the ip form of sampling.control (validated by
# valid.ip.control() in DynaPrism.R): {max.iter, n.cores, alpha, beta, tol,
# relax}. alpha, max.iter, tol, relax enter the sampling hash; beta rides the
# run-id suffix; n.cores is runtime. There is no seed (deterministic), no
# burn.in / thinning (no chain).
#
# SOURCE ORDER
#   source("DynaPrism.R"); source("Omega_tensor.R")   # sampling.control lives here
#   source("DynaPrism_instaprism.R")
# ==============================================================================


# ------------------------------------------------------------------------------
# Fixed-point kernel: joint Z + theta (initial pass, cell-state level)
# Drop-in analog of sample.Z.theta_n.Omega(). Deterministic; no RNG.
# ------------------------------------------------------------------------------
ip.Z.theta_n.Omega <- function(X_n, phi, alpha, Omega, beta, max.iter,
                               compute.elbo = FALSE, tol = 1e-5, relax = 1) {
  
  K_names <- rownames(phi); G_names <- colnames(phi)
  K_dim   <- length(K_names); G_dim <- length(G_names)
  
  phi_mat <- as.matrix(phi)                           # K x G
  X_vec   <- as.numeric(X_n[G_names])
  
  use_omega <- !is.null(Omega) && length(Omega) > 1L
  if (use_omega) {
    G_shared_names <- intersect(G_names, dimnames(Omega)$Gene)
    use_omega      <- length(G_shared_names) > 0L
  }
  if (use_omega) {
    G_shared_idx <- which(G_names %in% G_shared_names)
    Omega_flat   <- matrix(Omega[K_names, K_names, G_shared_names],
                           nrow = K_dim, ncol = K_dim * length(G_shared_names))
  }
  
  theta     <- rep(1 / K_dim, K_dim)
  converged <- FALSE
  delta     <- NA_real_
  
  for (i in seq_len(max.iter)) {
    prob.mat <- phi_mat * theta                       # [k,g] = phi[k,g] * theta[k]
    
    if (use_omega) {
      np_vec <- as.numeric(theta %*% Omega_flat)       # sum_s theta_s Omega[s,r,g]
      np_mat <- matrix(np_vec, nrow = K_dim, ncol = length(G_shared_names))  # [receiver, g]
      e      <- beta * np_mat
      cmax   <- apply(e, 2, max)                        # per-gene max: cancels under
      prob.mat[, G_shared_idx] <-                       # column norm, prevents overflow
        prob.mat[, G_shared_idx] * exp(sweep(e, 2, cmax, "-"))
    }
    
    csum <- colSums(prob.mat); csum[csum == 0] <- 1
    B     <- prob.mat / rep(csum, each = K_dim)         # K x G, each gene-column sums to 1
    Z_KxG <- B * rep(X_vec, each = K_dim)               # [k,g] = X_g * B[k,g]
    
    theta_raw <- rowSums(Z_KxG) + alpha                 # colSums(Z) over genes + Dirichlet prior
    theta_new <- theta_raw / sum(theta_raw)             # Dirichlet posterior MEAN
    if (relax < 1) theta_new <- (1 - relax) * theta + relax * theta_new
    
    delta <- max(abs(theta_new - theta))
    theta <- theta_new
    if (delta < tol) { converged <- TRUE; break }
  }
  
  if (!converged)
    warning(sprintf(paste0("ip.Z.theta_n.Omega: no convergence after %d iters ",
                           "(delta = %.3g). Raise max.iter or lower relax."),
                    max.iter, delta))
  
  gibbs.constant <- 0
  if (compute.elbo) {
    # Log joint at the fixed point (a point estimate, not a chain mean). NB: both
    # Z_KxG and prob_norm are K x G here; the Gibbs-kernel ELBO multiplies a G x K
    # Z by a K x G prob_norm, non-conformable unless G == K -- worth a check in
    # DynaPrism.R if you ever enable compute.elbo on the Gibbs path.
    prob_norm <- B; prob_norm[prob_norm == 0] <- .Machine$double.eps
    log_lik   <- sum(Z_KxG * log(prob_norm))
    log_prior <- sum((alpha - 1) * log(pmax(theta, .Machine$double.eps)))
    log_coef  <- -sum(lfactorial(Z_KxG))
    gibbs.constant <- log_lik + log_prior + log_coef
  }
  
  list(Z_n            = t(Z_KxG),                        # G x K, matches Gibbs shape
       theta_n        = theta,                           # length K, rownames(phi) order
       theta.cv_n     = rep(NA_real_, K_dim),            # no posterior spread for a point estimate
       gibbs.constant = gibbs.constant)
}


# ------------------------------------------------------------------------------
# Fixed-point kernel: theta only (final pass, cell-type level)
# Drop-in analog of sample.theta_n.Omega(). Omega_type NULL/empty -> pure
# InstaPrism update (non-Omega-modulated final pass).
# ------------------------------------------------------------------------------
ip.theta_n.Omega <- function(X_n, phi, alpha, Omega_type, beta, max.iter,
                             tol = 1e-5, relax = 1) {
  
  K_names <- rownames(phi); G_names <- colnames(phi)
  K_dim   <- length(K_names); G_dim <- length(G_names)
  
  phi_mat <- as.matrix(phi)
  X_vec   <- as.numeric(X_n[G_names])
  
  use_omega <- !is.null(Omega_type) && length(Omega_type) > 1L
  if (use_omega) {
    G_shared_names <- intersect(G_names, dimnames(Omega_type)$Gene)
    use_omega      <- length(G_shared_names) > 0L
  }
  if (use_omega) {
    G_shared_idx <- which(G_names %in% G_shared_names)
    Omega_flat   <- matrix(Omega_type[K_names, K_names, G_shared_names],
                           nrow = K_dim, ncol = K_dim * length(G_shared_names))
  }
  
  theta     <- rep(1 / K_dim, K_dim)
  converged <- FALSE
  delta     <- NA_real_
  
  for (i in seq_len(max.iter)) {
    prob.mat <- phi_mat * theta
    if (use_omega) {
      np_vec <- as.numeric(theta %*% Omega_flat)
      np_mat <- matrix(np_vec, nrow = K_dim, ncol = length(G_shared_names))
      e      <- beta * np_mat
      cmax   <- apply(e, 2, max)
      prob.mat[, G_shared_idx] <- prob.mat[, G_shared_idx] * exp(sweep(e, 2, cmax, "-"))
    }
    csum <- colSums(prob.mat); csum[csum == 0] <- 1
    B     <- prob.mat / rep(csum, each = K_dim)
    Z_KxG <- B * rep(X_vec, each = K_dim)
    
    theta_raw <- rowSums(Z_KxG) + alpha
    theta_new <- theta_raw / sum(theta_raw)
    if (relax < 1) theta_new <- (1 - relax) * theta + relax * theta_new
    
    delta <- max(abs(theta_new - theta))
    theta <- theta_new
    if (delta < tol) { converged <- TRUE; break }
  }
  
  if (!converged)
    warning(sprintf(paste0("ip.theta_n.Omega: no convergence after %d iters ",
                           "(delta = %.3g). Raise max.iter or lower relax."),
                    max.iter, delta))
  
  list(theta_n    = theta,
       theta.cv_n = rep(NA_real_, K_dim))
}


# ------------------------------------------------------------------------------
# Dispatchers  (mirror run.gibbs.refPhi.*.Omega; knobs read from sampling.control,
# no MCMC index machinery, no time estimate)
# ------------------------------------------------------------------------------

run.ip.refPhi.ini.Omega <- function(gibbsSampler.obj, compute.elbo = FALSE) {
  
  phi   <- gibbsSampler.obj@reference@phi
  X     <- gibbsSampler.obj@X
  Omega <- gibbsSampler.obj@Omega
  sc    <- gibbsSampler.obj@gibbs.control            # ip-form sampling.control
  alpha    <- sc$alpha; beta <- sc$beta; n.cores <- sc$n.cores
  max.iter <- sc$max.iter; tol <- sc$tol; relax <- sc$relax
  tmp.dir  <- tempdir(check = TRUE)
  
  cat("Start InstaPrism-style fixed point (initial, cell-state)...\n")
  
  cpu.fun <- function(n) {
    file.name.X_n <- file.path(tmp.dir, paste0("mixture_", n, ".rdata"))
    load(file.name.X_n)                              # -> X_n
    ip.Z.theta_n.Omega(X_n = X_n, phi = phi, alpha = alpha, Omega = Omega,
                       beta = beta, max.iter = max.iter,
                       compute.elbo = compute.elbo, tol = tol, relax = relax)
  }
  
  gibbs.list <- .omega_cluster_lapply(
    n_items    = nrow(X),
    cpu.fun    = cpu.fun,
    n.cores    = n.cores,
    data.vars  = c("phi", "alpha", "Omega", "beta", "max.iter",
                   "compute.elbo", "tol", "relax", "tmp.dir"),
    data.envir = environment(),
    fun.vars   = c("ip.Z.theta_n.Omega"))
  if (n.cores <= 1) cat("\n")
  
  elbo_per_sample <- if (compute.elbo) {
    setNames(vapply(gibbs.list, `[[`, numeric(1L), "gibbs.constant"), rownames(X))
  } else numeric(0)
  
  jointPost <- newJointPost(bulkID     = rownames(X),
                            geneID     = colnames(X),
                            cellType   = rownames(phi),
                            gibbs.list = gibbs.list)
  
  list(jointPost = jointPost, elbo_per_sample = elbo_per_sample)
}


run.ip.refPhi.final.Omega <- function(gibbsSampler.obj) {
  
  phi        <- gibbsSampler.obj@reference@phi
  X          <- gibbsSampler.obj@X
  Omega_type <- gibbsSampler.obj@Omega
  sc         <- gibbsSampler.obj@gibbs.control
  alpha    <- sc$alpha; beta <- sc$beta; n.cores <- sc$n.cores
  max.iter <- sc$max.iter; tol <- sc$tol; relax <- sc$relax
  tmp.dir  <- tempdir(check = TRUE)
  
  cat("Start InstaPrism-style fixed point (final, cell-type)...\n")
  
  cpu.fun <- function(n) {
    file.name.X_n <- file.path(tmp.dir, paste0("mixture_", n, ".rdata"))
    load(file.name.X_n)                              # -> X_n
    ip.theta_n.Omega(X_n = X_n, phi = phi, alpha = alpha,
                     Omega_type = Omega_type, beta = beta, max.iter = max.iter,
                     tol = tol, relax = relax)
  }
  
  gibbs.list <- .omega_cluster_lapply(
    n_items    = nrow(X),
    cpu.fun    = cpu.fun,
    n.cores    = n.cores,
    data.vars  = c("phi", "alpha", "Omega_type", "beta", "max.iter",
                   "tol", "relax", "tmp.dir"),
    data.envir = environment(),
    fun.vars   = c("ip.theta_n.Omega"))
  if (n.cores <= 1) cat("\n")
  
  newThetaPost(bulkID     = rownames(X),
               cellType   = rownames(phi),
               gibbs.list = gibbs.list)
}


#' InstaPrism-engine analog of run.gibbs.Omega (dispatch on the `final` flag).
run.ip.Omega <- function(gibbsSampler.obj, final, compute.elbo = FALSE) {
  if (!is(gibbsSampler.obj@reference, "refPhi"))
    stop("run.ip.Omega supports refPhi references only (refTumor uses the Gibbs path).")
  if (!final)
    run.ip.refPhi.ini.Omega(gibbsSampler.obj, compute.elbo = compute.elbo)
  else
    run.ip.refPhi.final.Omega(gibbsSampler.obj)
}


# ------------------------------------------------------------------------------
# Top-level runner: mirror of run.prism.Omega() using the fixed-point engine.
# Reuses .run_meta_Omega() for identity; `sampling.control` is the ip form. The
# "ip" method tag is the only identity marker distinguishing the engines.
# ------------------------------------------------------------------------------
run.prism.Omega.ip <- function(prism,
                               omega.control    = NULL,
                               beta,
                               method           = "ip",
                               mixture.name     = "mix",
                               reference.name   = "ref",
                               Omega            = NULL,
                               Omega_type       = NULL,
                               n.cores          = 1,
                               update.gibbs     = TRUE,
                               compute.elbo     = FALSE,
                               sampling.control = list(),
                               opt.control      = list(),
                               store.C.input    = c("value", "hash")) {
  
  store.C.input <- match.arg(store.C.input)
  
  # sampling.control travels in the 11th positional slot (received as
  # gibbs.control) and is validated as ip form because method = "ip".
  meta          <- .run_meta_Omega(prism, beta, method, mixture.name, reference.name,
                                   omega.control, Omega, Omega_type, update.gibbs,
                                   opt.control, sampling.control, n.cores, store.C.input)
  method           <- meta$method
  opt.control      <- meta$opt.control
  sampling.control <- meta$sampling.control
  omega.control    <- meta$omega.control
  store.C.input    <- meta$store.C.input
  omega.log        <- meta$omega.log
  
  # --- build tensors (unless pre-built ones were supplied) --------------------
  if (is.null(Omega) && !is.null(omega.control)) {
    message("Building state-level Omega from omega.control$state ...")
    Omega <- build_Omega_from_control(prism@phi_cellState@phi, omega.control$state)
  }
  if (update.gibbs && isTRUE(opt.control$update.omega) &&
      is.null(Omega_type) && !is.null(omega.control)) {
    message("Building type-level Omega from omega.control$type ...")
    Omega_type <- build_Omega_from_control(prism@phi_cellType@phi, omega.control$type)
  }
  if (is.null(Omega))
    stop("No state-level Omega available: supply `omega.control` or a pre-built `Omega`.")
  
  # --- write mixture slices to disk for memory-frugal workers -----------------
  tmp.dir <- tempdir(check = TRUE)
  for (n in seq_len(nrow(prism@mixture))) {
    X_n <- prism@mixture[n, ]
    save(X_n, file = file.path(tmp.dir, paste0("mixture_", n, ".rdata")))
  }
  on.exit(unlink(tmp.dir, recursive = TRUE), add = TRUE)
  
  make_control_param <- function() list(
    sampling.control = sampling.control,      # canonical key (ip form)
    opt.control      = opt.control,
    update.gibbs     = update.gibbs,
    method           = method,
    mixture.name     = mixture.name,
    reference.name   = reference.name,
    store.C.input    = store.C.input,
    omega.log        = omega.log
  )
  
  # --- initial fixed point: cell-state level, Omega-modulated -----------------
  gibbsSampler.ini.cs <- new("gibbsSamplerOmega",
                             reference     = prism@phi_cellState,
                             X             = prism@mixture,
                             Omega         = Omega,
                             gibbs.control = sampling.control)   # slot carries ip form
  
  ini_result       <- run.ip.Omega(gibbsSampler.ini.cs, final = FALSE,
                                   compute.elbo = compute.elbo)
  jointPost.ini.cs <- ini_result$jointPost
  elbo_per_sample  <- ini_result$elbo_per_sample
  
  jointPost.ini.ct <- mergeK(jointPost.obj = jointPost.ini.cs, map = prism@map)
  
  if (!update.gibbs) {
    return(new("BayesPrismOmega",
               prism                       = prism,
               elbo_per_sample             = elbo_per_sample,
               posterior.initial.cellState = jointPost.ini.cs,
               posterior.initial.cellType  = jointPost.ini.ct,
               control_param               = make_control_param()))
  }
  
  # --- MAP/MLE reference update (never Omega-modulated) -----------------------
  opt.control.base <- opt.control
  opt.control.base$update.omega <- NULL
  
  psi <- updateReference(Z           = jointPost.ini.ct@Z,
                         phi_prime   = prism@phi_cellType,
                         map         = prism@map,
                         key         = prism@key,
                         opt.control = opt.control.base)
  
  # --- final fixed point: cell-type level; Omega-modulated iff update.omega ---
  if (isTRUE(opt.control$update.omega) && !is.null(Omega_type)) {
    gibbsSampler.update <- new("gibbsSamplerOmega",
                               reference     = psi,
                               X             = prism@mixture,
                               Omega         = Omega_type,
                               gibbs.control = sampling.control)
  } else {
    # Pure InstaPrism final pass: sentinel Omega (length 1) => niche term skipped.
    gibbsSampler.update <- new("gibbsSamplerOmega",
                               reference     = psi,
                               X             = prism@mixture,
                               Omega         = array(0, dim = c(1L, 1L, 1L)),
                               gibbs.control = sampling.control)
  }
  theta_f <- run.ip.Omega(gibbsSampler.update, final = TRUE)
  
  new("BayesPrismOmega",
      prism                       = prism,
      elbo_per_sample             = elbo_per_sample,
      posterior.initial.cellState = jointPost.ini.cs,
      posterior.initial.cellType  = jointPost.ini.ct,
      reference.update            = psi,
      posterior.theta_f           = theta_f,
      control_param               = make_control_param())
}


# ------------------------------------------------------------------------------
# Cache-aware wrapper (mirror of save.run.Omega; method defaults to "ip").
# ------------------------------------------------------------------------------
save.run.ip <- function(prism,
                        omega.control    = NULL,
                        beta,
                        method           = "ip",
                        mixture.name     = "mix",
                        reference.name   = "ref",
                        Omega            = NULL,
                        Omega_type       = NULL,
                        n.cores          = 1,
                        update.gibbs     = TRUE,
                        compute.elbo     = FALSE,
                        sampling.control = list(),
                        opt.control      = list(),
                        store.C.input    = c("value", "hash"),
                        cache.dir,
                        overwrite        = FALSE) {
  
  stopifnot(!missing(cache.dir), is.character(cache.dir), length(cache.dir) == 1L)
  store.C.input <- match.arg(store.C.input)
  
  run.id <- .run_meta_Omega(prism, beta, method, mixture.name, reference.name,
                            omega.control, Omega, Omega_type, update.gibbs,
                            opt.control, sampling.control, n.cores, store.C.input)$run.id
  
  if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
  cache.file <- file.path(cache.dir, paste0(run.id, ".rds"))
  
  if (file.exists(cache.file) && !overwrite) {
    message("Loading cached result: ", cache.file)
    return(readRDS(cache.file))
  }
  
  res <- run.prism.Omega.ip(prism = prism, omega.control = omega.control, beta = beta,
                            method = method, mixture.name = mixture.name,
                            reference.name = reference.name, Omega = Omega,
                            Omega_type = Omega_type, n.cores = n.cores,
                            update.gibbs = update.gibbs, compute.elbo = compute.elbo,
                            sampling.control = sampling.control, opt.control = opt.control,
                            store.C.input = store.C.input)
  message("Saving result: ", cache.file)
  saveRDS(res, cache.file)
  res
}