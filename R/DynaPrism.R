# ==============================================================================
# DynaPrism – core functions
#
# Key additions vs. previous version:
#   1. BayesPrismOmega gains elbo_per_sample slot ("numeric", named by bulk
#      sample ID).  Populated when compute.elbo = TRUE; numeric(0) otherwise.
#   2. sample.Z.theta_n.Omega: ELBO formula corrected to the full log joint
#      (log-likelihood + log-prior + combinatorial term), replacing the old
#      incomplete multinomial-coefficient-only formula.
#   3. run.prism.Omega: new compute.elbo = FALSE parameter; plumbed through
#      to the initial-pass sampler; per-sample ELBO extracted before
#      newJointPost discards it.
#   4. run.gibbs.Omega / run.gibbs.refPhi.ini.Omega: compute.elbo forwarded;
#      ini dispatcher now returns list(jointPost, elbo_per_sample).
# ==============================================================================


estimate.gibbs.time.Omega <- function(gibbsSampler.obj,
                                      final,
                                      chain.length = 50) {
  
  ref           <- gibbsSampler.obj@reference
  X             <- gibbsSampler.obj@X
  gibbs.control <- gibbsSampler.obj@gibbs.control
  
  ptm <- proc.time()
  
  if (!final) {
    sample.Z.theta_n.Omega(X_n       = X[1, ],
                           phi       = ref@phi,
                           alpha     = gibbs.control$alpha,
                           Omega     = gibbsSampler.obj@Omega,
                           beta      = gibbs.control$beta,
                           gibbs.idx = get.gibbs.idx(
                             list(chain.length = chain.length,
                                  burn.in      = chain.length * gibbs.control$burn.in /
                                    gibbs.control$chain.length,
                                  thinning     = gibbs.control$thinning)))
  } else {
    if (is(ref, "refPhi")) {
      sample.theta_n(X_n       = X[1, ],
                     phi       = ref@phi,
                     alpha     = gibbs.control$alpha,
                     gibbs.idx = get.gibbs.idx(
                       list(chain.length = chain.length,
                            burn.in      = chain.length * gibbs.control$burn.in /
                              gibbs.control$chain.length,
                            thinning     = gibbs.control$thinning)))
    }
    if (is(ref, "refTumor")) {
      phi_1      <- rbind(ref@psi_mal[1, ], ref@psi_env)
      nonzero.idx <- apply(phi_1, 2, max) > 0
      sample.theta_n(X_n       = X[1, nonzero.idx, drop = FALSE],
                     phi       = phi_1[, nonzero.idx, drop = FALSE],
                     alpha     = gibbs.control$alpha,
                     gibbs.idx = get.gibbs.idx(
                       list(chain.length = chain.length,
                            burn.in      = chain.length * gibbs.control$burn.in /
                              gibbs.control$chain.length,
                            thinning     = gibbs.control$thinning)))
    }
  }
  
  total.time     <- as.numeric((proc.time() - ptm)["elapsed"])
  estimated.time <- gibbs.control$chain.length / chain.length * total.time *
    ceiling(nrow(X) / gibbs.control$n.cores) * 2
  current.time   <- Sys.time()
  
  cat("Current time: ",              as.character(current.time),                  "\n")
  cat("Estimated time to complete: ", my_seconds_to_period(estimated.time),       "\n")
  cat("Estimated finishing time: ",   as.character(current.time + estimated.time), "\n")
}


# ------------------------------------------------------------------------------
# Classes
# ------------------------------------------------------------------------------

setClass("gibbsSamplerOmega",
         slots = c(
           reference     = "reference",
           X             = "matrix",
           Omega         = "array",
           gibbs.control = "list"
         ),
         prototype = list(
           reference     = NULL,
           X             = matrix(),
           gibbs.control = list(chain.length = NULL, burn.in = NULL, thinning = NULL,
                                n.cores = NULL, seed = NULL, alpha = NULL, beta = NULL)
         ),
         validity = function(object) {
           errors   <- character()
           ref.gene <- if (is(object@reference, "refPhi"))   colnames(object@reference@phi)
           else                                   colnames(object@reference@psi_mal)
           if (!identical(ref.gene, colnames(object@X)))
             errors <- c(errors, "Gene names do not match between reference and X")
           if (length(errors) == 0) TRUE else errors
         }
)


#' An S4 class to represent the output of run.prism.Omega
#'
#' @slot elbo_per_sample Named numeric vector (bulk sample → ELBO).
#'   Populated when compute.elbo = TRUE; numeric(0) otherwise.
setClass("BayesPrismOmega",
         slots = c(
           prism                       = "prism",
           Omega                       = "array",
           Omega_type                  = "ANY",
           elbo_per_sample             = "numeric",   # NEW: named, length N or 0
           posterior.initial.cellState = "jointPost",
           posterior.initial.cellType  = "jointPost",
           reference.update            = "reference",
           posterior.theta_f           = "theta_f",
           control_param               = "list"
         ),
         prototype = list(
           prism                       = new("prism"),
           Omega_type                  = NULL,
           elbo_per_sample             = numeric(0),  # numeric(0) = not computed
           posterior.initial.cellState = new("jointPost"),
           posterior.initial.cellType  = new("jointPost"),
           reference.update            = NULL,
           posterior.theta_f           = NULL,
           control_param = list(
             gibbs.control = list(chain.length = NULL, burn.in = NULL, thinning = NULL,
                                  n.cores = NULL, seed = NULL, alpha = NULL, beta = NULL),
             opt.control   = list(trace = NULL, maxit = NULL,
                                  optimizer = NA_character_, n.cores = NULL),
             update.gibbs  = logical()
           )
         )
)


# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

valid.gibbs.control.Omega <- function(control) {
  ctrl <- list(chain.length = 1000, burn.in = 500, thinning = 2,
               n.cores = 1, seed = 123, alpha = 1, beta = 0.4)
  namc <- names(control)
  if (!all(namc %in% names(ctrl)))
    stop("unknown names in gibbs.control: ",
         paste(namc[!(namc %in% names(ctrl))], collapse = ", "))
  ctrl[namc] <- control
  if (ctrl$alpha < 0) stop("alpha must be non-negative")
  if (ctrl$beta  < 0) stop("beta must be non-negative")
  ctrl
}


# ------------------------------------------------------------------------------
# Top-level run function
# ------------------------------------------------------------------------------

#' Run DynaPrism deconvolution (Gibbs sampler)
#'
#' @param compute.elbo logical; compute per-sample log joint (ELBO) during the
#'   initial Gibbs pass and store in result@elbo_per_sample.  Default FALSE.
#'   Adds ~10-20% to run time.
run.prism.Omega <- function(prism,
                            Omega,
                            Omega_type    = NULL,
                            n.cores       = 1,
                            update.gibbs  = TRUE,
                            compute.elbo  = FALSE,      # NEW
                            gibbs.control = list(),
                            opt.control   = list()) {
  
  if (!"n.cores" %in% names(gibbs.control)) gibbs.control$n.cores <- n.cores
  if (!"n.cores" %in% names(opt.control))   opt.control$n.cores   <- n.cores
  stopifnot(is.logical(update.gibbs) & length(update.gibbs) == 1)
  stopifnot(is.numeric(n.cores)      & length(n.cores)      == 1)
  
  opt.control   <- valid.opt.control(opt.control)
  gibbs.control <- valid.gibbs.control.Omega(gibbs.control)
  
  if (prism@phi_cellState@pseudo.min == 0)
    gibbs.control$alpha <- max(1, gibbs.control$alpha)
  
  # Write mixture slices to disk for parallel workers
  tmp.dir <- tempdir(check = TRUE)
  for (n in seq_len(nrow(prism@mixture))) {
    X_n <- prism@mixture[n, ]
    save(X_n, file = file.path(tmp.dir, paste0("mixture_", n, ".rdata")))
  }
  
  # Initial Gibbs: cell-state level, Omega-modulated
  gibbsSampler.ini.cs <- new("gibbsSamplerOmega",
                             reference     = prism@phi_cellState,
                             X             = prism@mixture,
                             Omega         = Omega,
                             gibbs.control = gibbs.control)
  
  # run.gibbs.Omega now returns list(jointPost, elbo_per_sample)
  ini_result       <- run.gibbs.Omega(gibbsSampler.ini.cs,
                                      Omega        = Omega,
                                      final        = FALSE,
                                      compute.elbo = compute.elbo)   # NEW
  jointPost.ini.cs  <- ini_result$jointPost
  elbo_per_sample   <- ini_result$elbo_per_sample
  
  jointPost.ini.ct <- mergeK(jointPost.obj = jointPost.ini.cs,
                             map           = prism@map)
  
  if (!update.gibbs) {
    unlink(tmp.dir, recursive = TRUE)
    return(new("BayesPrismOmega",
               prism                       = prism,
               Omega                       = Omega,
               Omega_type                  = Omega_type,
               elbo_per_sample             = elbo_per_sample,
               posterior.initial.cellState = jointPost.ini.cs,
               posterior.initial.cellType  = jointPost.ini.ct,
               control_param = list(gibbs.control = gibbs.control,
                                    opt.control   = opt.control,
                                    update.gibbs  = update.gibbs)))
  }
  
  # MAP/MLE reference update
  psi <- updateReference(Z           = jointPost.ini.ct@Z,
                         phi_prime   = prism@phi_cellType,
                         map         = prism@map,
                         key         = prism@key,
                         opt.control = opt.control)
  
  # Final Gibbs: cell-type level
  gibbsSampler.update <- new("gibbsSampler",
                             reference     = psi,
                             X             = prism@mixture,
                             gibbs.control = gibbs.control)
  
  if (!is.null(Omega_type)) {
    theta_f <- run.gibbs.refPhi.final.Omega(gibbsSampler.obj = gibbsSampler.update,
                                            Omega_type       = Omega_type,
                                            compute.elbo     = FALSE)
  } else {
    theta_f <- run.gibbs(gibbsSampler.update, final = TRUE)
  }
  
  unlink(tmp.dir, recursive = TRUE)
  new("BayesPrismOmega",
      prism                       = prism,
      Omega                       = Omega,
      Omega_type                  = Omega_type,
      elbo_per_sample             = elbo_per_sample,
      posterior.initial.cellState = jointPost.ini.cs,
      posterior.initial.cellType  = jointPost.ini.ct,
      reference.update            = psi,
      posterior.theta_f           = theta_f,
      control_param = list(gibbs.control = gibbs.control,
                           opt.control   = opt.control,
                           update.gibbs  = update.gibbs))
}


# ------------------------------------------------------------------------------
# Gibbs dispatch
# ------------------------------------------------------------------------------

#' Returns list(jointPost, elbo_per_sample) for the initial pass.
#' final=TRUE is dead code — kept for safety, throws a clear error.
run.gibbs.Omega <- function(gibbsSampler.obj,
                            Omega,
                            final,
                            if.estimate  = TRUE,
                            compute.elbo = FALSE) {   # NEW
  
  if (final) cat("Run Gibbs sampling using updated reference ...\n")
  else       cat("Run Gibbs sampling...\n")
  
  if (if.estimate)
    estimate.gibbs.time.Omega(gibbsSampler.obj = gibbsSampler.obj, final = final)
  
  if (is(gibbsSampler.obj@reference, "refPhi")) {
    if (!final)
      return(run.gibbs.refPhi.ini.Omega(gibbsSampler.obj = gibbsSampler.obj,
                                        Omega            = Omega,
                                        compute.elbo     = compute.elbo))
    else
      stop("run.gibbs.Omega with final=TRUE is not supported. ",
           "Call run.gibbs.refPhi.final.Omega() directly.")
  }
  
  if (is(gibbsSampler.obj@reference, "refTumor"))
    return(run.gibbs.refTumor(gibbsSampler.obj = gibbsSampler.obj))
}


# ------------------------------------------------------------------------------
# Initial Gibbs dispatcher — returns list(jointPost, elbo_per_sample)
# ------------------------------------------------------------------------------

run.gibbs.refPhi.ini.Omega <- function(gibbsSampler.obj,
                                       Omega,
                                       compute.elbo) {
  
  phi           <- gibbsSampler.obj@reference@phi
  X             <- gibbsSampler.obj@X
  gibbs.control <- gibbsSampler.obj@gibbs.control
  alpha         <- gibbs.control$alpha
  beta          <- gibbs.control$beta
  gibbs.idx     <- get.gibbs.idx(gibbs.control)
  seed          <- gibbs.control$seed
  
  cat("Start run...\n")
  
  if (gibbs.control$n.cores > 1) {
    sfInit(parallel = TRUE, cpus = gibbs.control$n.cores, type = "SOCK")
    tmp.dir <- tempdir(check = TRUE)
    sfExport("phi", "alpha", "Omega", "beta", "gibbs.idx", "seed",
             "compute.elbo", "sample.Z.theta_n.Omega", "rdirichlet", "tmp.dir")
    
    cpu.fun <- function(n) {
      if (!is.null(seed)) set.seed(seed)
      file.name.X_n <- file.path(tmp.dir, paste0("mixture_", n, ".rdata"))
      load(file.name.X_n)
      sample.Z.theta_n.Omega(X_n          = X_n,
                             phi          = phi,
                             alpha        = alpha,
                             Omega        = Omega,
                             beta         = beta,
                             gibbs.idx    = gibbs.idx,
                             compute.elbo = compute.elbo)
    }
    environment(cpu.fun) <- globalenv()
    gibbs.list <- sfLapply(seq_len(nrow(X)), cpu.fun)
    sfStop()
    
  } else {
    cpu.fun <- function(n) {
      if (!is.null(seed)) set.seed(seed)
      cat(n, " ")
      sample.Z.theta_n.Omega(X_n          = X[n, ],
                             phi          = phi,
                             alpha        = alpha,
                             Omega        = Omega,
                             beta         = beta,
                             gibbs.idx    = gibbs.idx,
                             compute.elbo = compute.elbo)
    }
    gibbs.list <- lapply(seq_len(nrow(X)), cpu.fun)
    cat("\n")
  }
  
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
# Final Gibbs dispatcher (cell-type level, Omega_type-modulated)
# ------------------------------------------------------------------------------

run.gibbs.refPhi.final.Omega <- function(gibbsSampler.obj,
                                         Omega_type,
                                         compute.elbo = FALSE) {
  
  phi           <- gibbsSampler.obj@reference@phi
  X             <- gibbsSampler.obj@X
  gibbs.control <- gibbsSampler.obj@gibbs.control
  alpha         <- gibbs.control$alpha
  beta          <- gibbs.control$beta
  gibbs.idx     <- get.gibbs.idx(gibbs.control)
  seed          <- gibbs.control$seed
  
  cat("Start run...\n")
  
  if (gibbs.control$n.cores > 1) {
    sfInit(parallel = TRUE, cpus = gibbs.control$n.cores, type = "SOCK")
    tmp.dir <- tempdir(check = TRUE)
    sfExport("phi", "alpha", "Omega_type", "beta", "gibbs.idx", "seed",
             "sample.theta_n.Omega", "rdirichlet", "tmp.dir")
    
    cpu.fun <- function(n) {
      if (!is.null(seed)) set.seed(seed)
      file.name.X_n <- file.path(tmp.dir, paste0("mixture_", n, ".rdata"))
      load(file.name.X_n)
      sample.theta_n.Omega(X_n        = X_n,
                           phi        = phi,
                           alpha      = alpha,
                           Omega_type = Omega_type,
                           beta       = beta,
                           gibbs.idx  = gibbs.idx)
    }
    environment(cpu.fun) <- globalenv()
    gibbs.list <- sfLapply(seq_len(nrow(X)), cpu.fun)
    sfStop()
    
  } else {
    cpu.fun <- function(n) {
      if (!is.null(seed)) set.seed(seed)
      cat(n, " ")
      sample.theta_n.Omega(X_n        = X[n, ],
                           phi        = phi,
                           alpha      = alpha,
                           Omega_type = Omega_type,
                           beta       = beta,
                           gibbs.idx  = gibbs.idx)
    }
    gibbs.list <- lapply(seq_len(nrow(X)), cpu.fun)
    cat("\n")
  }
  
  newThetaPost(bulkID     = rownames(X),
               cellType   = rownames(phi),
               gibbs.list = gibbs.list)
}


# ------------------------------------------------------------------------------
# Sampler: joint Z + theta  (initial pass, cell-state level)
# ------------------------------------------------------------------------------

#' ELBO formula: full log joint (dropping constant Σ_g log(X_ng!))
#'
#'   log P(Z_n^i, θ_n^i | X_n, φ, Ω, β) ≈
#'       Σ_k (α-1) log θ_nk^i             ← log Dirichlet prior
#'     + Σ_g Σ_k z_nkg^i log p̃_nkg^i     ← log multinomial likelihood
#'     - Σ_g Σ_k log(z_nkg^i !)           ← combinatorial term
#'
#' p̃_nkg^i is the normalised niche-modulated probability computed at the
#' START of iteration i (from θ^{i-1}).  Z_n^i was sampled from this
#' distribution.  θ_n^i is the Dirichlet draw from Z_n^i.  The one-step
#' lag between prob (from θ^{i-1}) and prior (from θ^i) is standard in
#' Gibbs ELBO estimation and is asymptotically consistent.
sample.Z.theta_n.Omega <- function(X_n,
                                   phi,
                                   alpha,
                                   Omega,
                                   beta,
                                   gibbs.idx,
                                   compute.elbo = FALSE) {
  
  K_names        <- rownames(phi)
  G_names        <- colnames(phi)
  G_shared_names <- intersect(G_names, dimnames(Omega)$Gene)
  G_shared_idx   <- which(G_names %in% G_shared_names)
  
  K_dim <- length(K_names)
  G_dim <- length(G_names)
  
  theta_n.i    <- rep(1/K_dim, K_dim)
  Z_n.i        <- array(NA, c(G_dim, K_dim))
  Z_n.sum      <- array(0,  c(G_dim, K_dim))
  theta_n.sum  <- rep(0, K_dim)
  theta_n2.sum <- rep(0, K_dim)
  elbo_sum     <- 0              # replaces old multinom.coef
  
  phi_mat <- as.matrix(phi)
  X_vec   <- as.numeric(X_n[G_names])
  
  Omega_flat <- matrix(Omega[K_names, K_names, G_shared_names],
                       nrow = K_dim,
                       ncol = K_dim * length(G_shared_names))
  
  for (i in seq_len(max(gibbs.idx))) {
    
    # Niche pressure and probability matrix (computed from θ^{i-1})
    niche_pressure_vec <- as.numeric(theta_n.i %*% Omega_flat)
    niche_pressure_mat <- matrix(niche_pressure_vec, nrow = K_dim,
                                 ncol = length(G_shared_names))
    
    prob.mat <- phi_mat * theta_n.i
    prob.mat[, G_shared_idx] <- prob.mat[, G_shared_idx] * exp(beta * niche_pressure_mat)
    
    # Sample Z | θ^{i-1}
    for (g in seq_len(G_dim))
      Z_n.i[g, ] <- rmultinom(n = 1, size = X_vec[g], prob = prob.mat[, g])
    
    # Sample θ^i | Z
    Z_nk.i    <- colSums(Z_n.i)
    theta_n.i <- rdirichlet(alpha = Z_nk.i + alpha)
    
    if (i %in% gibbs.idx) {
      Z_n.sum      <- Z_n.sum      + Z_n.i
      theta_n.sum  <- theta_n.sum  + theta_n.i
      theta_n2.sum <- theta_n2.sum + theta_n.i^2
      
      if (compute.elbo) {
        # Normalise prob.mat (from θ^{i-1}, consistent with the Z that was sampled)
        col_sums  <- colSums(prob.mat)
        col_sums[col_sums == 0] <- 1
        prob_norm <- sweep(prob.mat, 2, col_sums, "/")
        prob_norm[prob_norm == 0] <- .Machine$double.eps
        
        log_lik_i   <- sum(Z_n.i * log(prob_norm))                          # Σ z log p̃
        log_prior_i <- sum((alpha - 1) * log(pmax(theta_n.i,                # log Dir prior
                                                  .Machine$double.eps)))
        log_coef_i  <- -sum(lfactorial(Z_n.i))                              # -Σ log(z!)
        
        elbo_sum <- elbo_sum + log_lik_i + log_prior_i + log_coef_i
      }
    }
    if ((i %% 50) == 0) gc()
  }
  
  samples.size   <- length(gibbs.idx)
  Z_n            <- Z_n.sum      / samples.size
  theta_n        <- theta_n.sum  / samples.size
  theta.cv_n     <- sqrt(theta_n2.sum / samples.size - theta_n^2) / theta_n
  gibbs.constant <- elbo_sum / samples.size   # 0 when compute.elbo = FALSE
  
  list(Z_n            = Z_n,
       theta_n        = theta_n,
       theta.cv_n     = theta.cv_n,
       gibbs.constant = gibbs.constant)
}


# ------------------------------------------------------------------------------
# Sampler: theta only  (final pass, cell-type level)
# ------------------------------------------------------------------------------

sample.theta_n.Omega <- function(X_n,
                                 phi,
                                 alpha,
                                 Omega_type,
                                 beta,
                                 gibbs.idx) {
  
  K_names        <- rownames(phi)
  G_names        <- colnames(phi)
  G_shared_names <- intersect(G_names, dimnames(Omega_type)$Gene)
  G_shared_idx   <- which(G_names %in% G_shared_names)
  
  K_dim <- length(K_names)
  G_dim <- length(G_names)
  
  theta_n.i    <- rep(1/K_dim, K_dim)
  theta_n.sum  <- rep(0, K_dim)
  theta_n2.sum <- rep(0, K_dim)
  
  phi_mat <- as.matrix(phi)
  X_vec   <- as.numeric(X_n[G_names])
  
  Omega_flat <- matrix(Omega_type[K_names, K_names, G_shared_names],
                       nrow = K_dim,
                       ncol = K_dim * length(G_shared_names))
  
  for (i in seq_len(max(gibbs.idx))) {
    
    niche_pressure_vec <- as.numeric(theta_n.i %*% Omega_flat)
    niche_pressure_mat <- matrix(niche_pressure_vec, nrow = K_dim,
                                 ncol = length(G_shared_names))
    
    prob.mat <- phi_mat * theta_n.i
    prob.mat[, G_shared_idx] <- prob.mat[, G_shared_idx] * exp(beta * niche_pressure_mat)
    
    Z_n.i <- matrix(NA, G_dim, K_dim)
    for (g in seq_len(G_dim))
      Z_n.i[g, ] <- rmultinom(n = 1, size = X_vec[g], prob = prob.mat[, g])
    
    Z_nk.i    <- colSums(Z_n.i)
    theta_n.i <- rdirichlet(alpha = Z_nk.i + alpha)
    
    if (i %in% gibbs.idx) {
      theta_n.sum  <- theta_n.sum  + theta_n.i
      theta_n2.sum <- theta_n2.sum + theta_n.i^2
    }
    if ((i %% 50) == 0) gc()
  }
  
  samples.size <- length(gibbs.idx)
  theta_n      <- theta_n.sum  / samples.size
  theta.cv_n   <- sqrt(theta_n2.sum / samples.size - theta_n^2) / theta_n
  
  list(theta_n = theta_n, theta.cv_n = theta.cv_n)
}