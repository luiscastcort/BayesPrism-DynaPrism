# ==============================================================================
# DynaPrism – InstaPrism-style fixed-point acceleration
#
# Key addition vs. previous version:
#   run.instaprism.ini.Omega now extracts elbo_per_sample from gibbs.list
#   BEFORE newJointPost aggregates it away, and returns
#   list(jointPost, elbo_per_sample).  run.instaprism.Omega unpacks this and
#   stores elbo_per_sample in result@elbo_per_sample.
# ==============================================================================


# ------------------------------------------------------------------------------
# Control validation
# ------------------------------------------------------------------------------

valid.instaprism.control.Omega <- function(control = list()) {
  ctrl <- list(n.iter = 100, tol = 1e-6, n.cores = 1, seed = 123,
               alpha = 1, beta = 0.4, compute.elbo = FALSE)
  namc    <- names(control)
  unknown <- namc[!(namc %in% names(ctrl))]
  if (length(unknown) > 0)
    stop("unknown names in instaprism.control: ", paste(unknown, collapse = ", "))
  ctrl[namc] <- control
  if (ctrl$alpha  < 0) stop("alpha must be non-negative")
  if (ctrl$beta   < 0) stop("beta must be non-negative")
  if (ctrl$n.iter < 1) stop("n.iter must be >= 1")
  if (ctrl$tol   <= 0) stop("tol must be positive")
  ctrl
}


# ------------------------------------------------------------------------------
# Core fixed-point update (per sample)
# ------------------------------------------------------------------------------

#' ELBO at convergence (compute.elbo = TRUE):
#'
#'   ELBO_n = Σ_g Σ_k E[Z_nkg] · log p̃_nkg        <- log likelihood
#'           + Σ_k (alpha-1) · log theta_nk          <- log Dirichlet prior
#'           - Σ_g Σ_k lgamma(E[Z_nkg] + 1)          <- combinatorial term
#'
#' Because Z_mat is continuous, lfactorial is replaced by lgamma(Z+1).
#' This gives a smooth, noise-free ELBO as a function of beta.
instaprism_fp_n <- function(X_n,
                            phi,
                            alpha,
                            Omega,
                            beta,
                            n.iter       = 100,
                            tol          = 1e-6,
                            return_Z     = TRUE,
                            compute.elbo = FALSE) {
  
  K_names        <- rownames(phi)
  G_names        <- colnames(phi)
  G_shared_names <- intersect(G_names, dimnames(Omega)[[3]])
  G_shared_idx   <- which(G_names %in% G_shared_names)
  
  K_dim <- length(K_names)
  
  phi_mat <- as.matrix(phi)
  X_vec   <- as.numeric(X_n[G_names])
  
  Omega_flat <- matrix(Omega[K_names, K_names, G_shared_names],
                       nrow = K_dim,
                       ncol = K_dim * length(G_shared_names))
  
  theta_n        <- rep(1/K_dim, K_dim)
  names(theta_n) <- K_names
  
  for (i in seq_len(n.iter)) {
    theta_old <- theta_n
    
    np_vec <- as.numeric(theta_n %*% Omega_flat)
    np_mat <- matrix(np_vec, nrow = K_dim, ncol = length(G_shared_names))
    
    prob_mat <- phi_mat * theta_n
    prob_mat[, G_shared_idx] <- prob_mat[, G_shared_idx] * exp(beta * np_mat)
    
    col_sums <- colSums(prob_mat)
    col_sums[col_sums == 0] <- 1
    prob_mat_norm <- sweep(prob_mat, 2, col_sums, "/")
    
    Z_mat   <- sweep(prob_mat_norm, 2, X_vec, "*")   # [K x G]
    Z_nk    <- rowSums(Z_mat)
    theta_n <- (Z_nk + alpha) / sum(Z_nk + alpha)
    
    if (sum(abs(theta_n - theta_old)) < tol) break
  }
  
  # ELBO at convergence (prob_mat_norm and Z_mat are from the final iteration)
  elbo <- 0
  if (compute.elbo) {
    prob_norm_safe <- prob_mat_norm
    prob_norm_safe[prob_norm_safe == 0] <- .Machine$double.eps
    
    log_lik   <- sum(Z_mat * log(prob_norm_safe))
    log_prior <- sum((alpha - 1) * log(pmax(theta_n, .Machine$double.eps)))
    log_coef  <- -sum(lgamma(Z_mat + 1))
    
    elbo <- log_lik + log_prior + log_coef
  }
  
  out <- list(theta_n        = theta_n,
              theta.cv_n     = rep(NA_real_, K_dim),
              gibbs.constant = elbo)
  
  if (return_Z) out$Z_n <- t(Z_mat)   # [G x K]
  
  out
}


# ------------------------------------------------------------------------------
# Internal: cross-platform parallel lapply
# ------------------------------------------------------------------------------

.parallel_lapply <- function(X_list, FUN, n.cores, .envir = parent.frame()) {
  if (n.cores == 1)
    return(lapply(X_list, FUN))
  if (.Platform$OS.type == "unix")
    return(parallel::mclapply(X_list, FUN, mc.cores = n.cores))
  # Windows / socket cluster fallback
  cl <- parallel::makeCluster(n.cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  vars <- ls(envir = .envir)
  parallel::clusterExport(cl, varlist = vars, envir = .envir)
  parallel::parLapply(cl, X_list, FUN)
}


# ------------------------------------------------------------------------------
# Initial pass dispatcher — returns list(jointPost, elbo_per_sample)
# ------------------------------------------------------------------------------

run.instaprism.ini.Omega <- function(gibbsSampler.obj,
                                     Omega,
                                     instaprism.control) {
  
  phi          <- gibbsSampler.obj@reference@phi
  X            <- gibbsSampler.obj@X
  alpha        <- instaprism.control$alpha
  beta         <- instaprism.control$beta
  n.iter       <- instaprism.control$n.iter
  tol          <- instaprism.control$tol
  n.cores      <- instaprism.control$n.cores
  compute.elbo <- instaprism.control$compute.elbo
  
  cat("Running InstaPrism fixed-point (initial pass, cell-state level)...\n")
  
  fp_fn <- instaprism_fp_n   # capture locally for Windows clusterExport
  
  worker <- function(n) {
    fp_fn(X_n          = X[n, ],
          phi          = phi,
          alpha        = alpha,
          Omega        = Omega,
          beta         = beta,
          n.iter       = n.iter,
          tol          = tol,
          return_Z     = TRUE,
          compute.elbo = compute.elbo)
  }
  
  gibbs.list <- .parallel_lapply(seq_len(nrow(X)), worker,
                                 n.cores = n.cores,
                                 .envir  = environment())
  
  # Extract per-sample ELBO BEFORE newJointPost aggregates it away
  elbo_per_sample <- if (compute.elbo) {
    setNames(vapply(gibbs.list, `[[`, numeric(1L), "gibbs.constant"),
             rownames(X))
  } else {
    numeric(0)
  }
  
  jointPost <- newJointPost(bulkID     = rownames(X),
                            geneID    = colnames(X),
                            cellType  = rownames(phi),
                            gibbs.list = gibbs.list)
  
  list(jointPost = jointPost, elbo_per_sample = elbo_per_sample)
}


# ------------------------------------------------------------------------------
# Final pass dispatcher
# ------------------------------------------------------------------------------

run.instaprism.final.Omega <- function(gibbsSampler.obj,
                                       Omega_type,
                                       instaprism.control) {
  
  phi          <- gibbsSampler.obj@reference@phi
  X            <- gibbsSampler.obj@X
  alpha        <- instaprism.control$alpha
  beta         <- instaprism.control$beta
  n.iter       <- instaprism.control$n.iter
  tol          <- instaprism.control$tol
  n.cores      <- instaprism.control$n.cores
  compute.elbo <- instaprism.control$compute.elbo
  
  if (!is.null(Omega_type)) {
    cat("Running InstaPrism fixed-point (final pass, cell-type level, with Omega_type)...\n")
    Omega_used <- Omega_type
  } else {
    cat("Running InstaPrism fixed-point (final pass, cell-type level, no niche pressure)...\n")
    K_names    <- rownames(phi)
    G_names    <- colnames(phi)
    Omega_used <- array(0, dim      = c(length(K_names), length(K_names), length(G_names)),
                        dimnames = list(K_names, K_names, G_names))
  }
  
  fp_fn <- instaprism_fp_n
  
  worker <- function(n) {
    fp_fn(X_n          = X[n, ],
          phi          = phi,
          alpha        = alpha,
          Omega        = Omega_used,
          beta         = beta,
          n.iter       = n.iter,
          tol          = tol,
          return_Z     = FALSE,
          compute.elbo = compute.elbo)
  }
  
  gibbs.list <- .parallel_lapply(seq_len(nrow(X)), worker,
                                 n.cores = n.cores,
                                 .envir  = environment())
  
  newThetaPost(bulkID     = rownames(X),
               cellType   = rownames(phi),
               gibbs.list = gibbs.list)
}


# ------------------------------------------------------------------------------
# Top-level function
# ------------------------------------------------------------------------------

run.instaprism.Omega <- function(prism,
                                 Omega,
                                 Omega_type          = NULL,
                                 n.cores             = 1,
                                 update.gibbs        = TRUE,
                                 instaprism.control  = list(),
                                 opt.control         = list()) {
  
  if (!"n.cores" %in% names(instaprism.control)) instaprism.control$n.cores <- n.cores
  if (!"n.cores" %in% names(opt.control))        opt.control$n.cores        <- n.cores
  stopifnot(is.logical(update.gibbs) & length(update.gibbs) == 1)
  
  opt.control        <- valid.opt.control(opt.control)
  instaprism.control <- valid.instaprism.control.Omega(instaprism.control)
  
  if (prism@phi_cellState@pseudo.min == 0)
    instaprism.control$alpha <- max(1, instaprism.control$alpha)
  
  gibbsSampler.ini.cs <- new("gibbsSamplerOmega",
                             reference     = prism@phi_cellState,
                             X             = prism@mixture,
                             Omega         = Omega,
                             gibbs.control = list(alpha = instaprism.control$alpha,
                                                  beta  = instaprism.control$beta))
  
  # Unpack list(jointPost, elbo_per_sample)
  ini_result       <- run.instaprism.ini.Omega(
    gibbsSampler.obj   = gibbsSampler.ini.cs,
    Omega              = Omega,
    instaprism.control = instaprism.control)
  jointPost.ini.cs  <- ini_result$jointPost
  elbo_per_sample   <- ini_result$elbo_per_sample
  
  jointPost.ini.ct <- mergeK(jointPost.obj = jointPost.ini.cs,
                             map           = prism@map)
  
  if (!update.gibbs) {
    return(new("BayesPrismOmega",
               prism                       = prism,
               Omega                       = Omega,
               Omega_type                  = Omega_type,
               elbo_per_sample             = elbo_per_sample,
               posterior.initial.cellState = jointPost.ini.cs,
               posterior.initial.cellType  = jointPost.ini.ct,
               control_param = list(instaprism.control = instaprism.control,
                                    opt.control        = opt.control,
                                    update.gibbs       = update.gibbs,
                                    method             = "instaprism")))
  }
  
  psi <- updateReference(Z           = jointPost.ini.ct@Z,
                         phi_prime   = prism@phi_cellType,
                         map         = prism@map,
                         key         = prism@key,
                         opt.control = opt.control)
  
  gibbsSampler.update <- new("gibbsSampler",
                             reference     = psi,
                             X             = prism@mixture,
                             gibbs.control = list(alpha = instaprism.control$alpha,
                                                  beta  = instaprism.control$beta))
  
  theta_f <- run.instaprism.final.Omega(
    gibbsSampler.obj   = gibbsSampler.update,
    Omega_type         = Omega_type,
    instaprism.control = instaprism.control)
  
  new("BayesPrismOmega",
      prism                       = prism,
      Omega                       = Omega,
      Omega_type                  = Omega_type,
      elbo_per_sample             = elbo_per_sample,
      posterior.initial.cellState = jointPost.ini.cs,
      posterior.initial.cellType  = jointPost.ini.ct,
      reference.update            = psi,
      posterior.theta_f           = theta_f,
      control_param = list(instaprism.control = instaprism.control,
                           opt.control        = opt.control,
                           update.gibbs       = update.gibbs,
                           method             = "instaprism"))
}