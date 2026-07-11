# ==============================================================================
# DynaPrism - core functions
# ==============================================================================


# ------------------------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

.digest4 <- function(x) substr(digest::digest(x), 1L, 4L)


#' Construct a deterministic run id / filename stem.
#'
#' @param method one of "bp", "dp", "ip".
#' @param mixture.name,reference.name user-defined dataset tags.
#' @param h_opt,h_gibbs,h_omega 4-char hashes (opt+key, gibbs, omega inputs).
#' @param beta numeric scalar; rendered as e.g. 1.4 -> "b1p4", 0.7 -> "b0p7",
#'   5 -> "b5".
#' @return character scalar, e.g. "dp_blkHM_HMHM_eb68630dab8b_b1p4".
make_run_id <- function(method, mixture.name, reference.name,
                        h_opt, h_gibbs, h_omega, beta) {
  beta_tag <- paste0("b", gsub("\\.", "p", format(beta, trim = TRUE, scientific = FALSE)))
  paste0(method, "_", mixture.name, "_", reference.name, "_",
         h_opt, h_gibbs, h_omega, "_", beta_tag)
}


#' Normalise omega.control into nested state/type form.
#' Flat list -> treated as $state; $type defaults to $state when absent.
.normalize_omega_control <- function(omega.control) {
  if (is.null(omega.control)) return(NULL)
  nested <- length(names(omega.control)) > 0 &&
    all(names(omega.control) %in% c("state", "type"))
  if (!nested) omega.control <- list(state = omega.control)
  if (is.null(omega.control$type)) omega.control$type <- omega.control$state
  omega.control
}


#' Build an Omega tensor from a phi matrix and one build_Omega arg set.
#' `control` is a named list of build_Omega() arguments EXCLUDING phi
#' (ligand_target_matrix, C_input, lr_network, mask_threshold, masking_type).
build_Omega_from_control <- function(phi, control) {
  if (is.null(control))
    stop("build_Omega_from_control: `control` is NULL (no build_Omega arguments).")
  do.call(build_Omega, c(list(phi = phi), control))
}


# ------------------------------------------------------------------------------
# Omega fingerprint log  (hash preimage + experiment record + rebuild key)
# ------------------------------------------------------------------------------
# One per-resolution structure serves three jobs.  The C_input *hash* is always
# recorded; the C_input *value* is stored only when store.C.input == "value".
# The invariant priors (ligand_target_matrix, lr_network) are always hash-only.
# Hashes never depend on the stored C_input value, so storage cannot perturb
# identity: a "value" object and a "hash" object built from the same context
# yield the same run id.

.digest_full <- function(x) digest::digest(x)

# Stored log entry for one resolution's build_Omega arg set.
.omega_log_entry <- function(control, store.C.input = "value") {
  entry <- list(
    mask_threshold = control$mask_threshold,
    masking_type   = control$masking_type,
    ltm_hash       = .digest_full(control$ligand_target_matrix),
    lr_hash        = .digest_full(control$lr_network),
    C_input_hash   = .digest_full(control$C_input)
  )
  if (identical(store.C.input, "value")) entry$C_input <- control$C_input
  entry
}

# Hash preimage (spec) from a stored log entry, in canonical order, excluding
# the optional stored C_input value.  Run-time and recompute-time both route
# through this, so their hashes are byte-identical by construction.
.omega_spec_from_entry <- function(entry) {
  list(mask_threshold = entry$mask_threshold,
       masking_type   = entry$masking_type,
       ltm_hash       = entry$ltm_hash,
       lr_hash        = entry$lr_hash,
       C_input_hash   = entry$C_input_hash)
}

# Nested log (state/type) from a normalised omega.control, or tensor
# fingerprints when no control is available (pre-built Omega/Omega_type).
.omega_log_from_control <- function(omega.control, Omega, Omega_type,
                                    store.C.input = "value") {
  if (!is.null(omega.control)) {
    list(state = .omega_log_entry(omega.control$state, store.C.input),
         type  = .omega_log_entry(omega.control$type,  store.C.input))
  } else {
    list(Omega = .digest_full(Omega), Omega_type = .digest_full(Omega_type))
  }
}

# Hash preimage from a stored log: nested spec for control-built runs, or the
# tensor-fingerprint list as-is for pre-built runs.
.omega_spec_from_log <- function(omega.log) {
  if (length(names(omega.log)) && all(names(omega.log) %in% c("state", "type"))) {
    list(state = .omega_spec_from_entry(omega.log$state),
         type  = .omega_spec_from_entry(omega.log$type))
  } else {
    omega.log
  }
}


# ------------------------------------------------------------------------------
# Hashing core
# ------------------------------------------------------------------------------
# Single source of truth for the opt/gibbs/omega triplet.  Assumes opt.control
# and gibbs.control are already validated.  n.cores is excluded from opt and
# gibbs; beta is excluded from gibbs (it rides the readable run-id suffix);
# key and update.gibbs enter the opt hash.  omega.spec is the nested state/type
# spec (or the tensor-fingerprint list for pre-built runs).
.compute_hashes <- function(key, opt.control, gibbs.control, update.gibbs, omega.spec) {
  gibbs_for_hash <- gibbs.control; gibbs_for_hash$n.cores <- NULL; gibbs_for_hash$beta <- NULL
  opt_for_hash   <- opt.control;   opt_for_hash$n.cores   <- NULL
  c(opt   = .digest4(list(key = key, opt.control = opt_for_hash, update.gibbs = update.gibbs)),
    gibbs = .digest4(gibbs_for_hash),
    omega = .digest4(omega.spec))
}


# ------------------------------------------------------------------------------
# Control validators
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


#' Extend base valid.opt.control() with the update.omega toggle.
#' The base validator rejects unknown names, so update.omega is peeled off
#' before delegation and re-attached afterwards.
valid.opt.control.Omega <- function(control) {
  update.omega <- control$update.omega %||% FALSE
  control$update.omega <- NULL
  ctrl <- valid.opt.control(control)          # base BayesPrism validator
  stopifnot(is.logical(update.omega), length(update.omega) == 1L)
  ctrl$update.omega <- update.omega
  ctrl
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
           ref.gene <- if (is(object@reference, "refPhi")) colnames(object@reference@phi)
           else                                 colnames(object@reference@psi_mal)
           if (!identical(ref.gene, colnames(object@X)))
             errors <- c(errors, "Gene names do not match between reference and X")
           if (length(errors) == 0) TRUE else errors
         }
)


#' An S4 class to represent the output of run.prism.Omega
#'
#' @slot elbo_per_sample Named numeric vector (bulk sample -> ELBO).
#'   Populated when compute.elbo = TRUE; numeric(0) otherwise.
#' @slot control_param carries gibbs.control, opt.control, update.gibbs, method,
#'   dataset names, store.C.input, and the omega.log fingerprint record.  The
#'   hashes and run id are NOT stored; recompute them with run.hashes(obj) /
#'   run.id(obj).  The Omega tensors are not stored either; rebuild them with
#'   rebuild_Omega(obj, ...).
setClass("BayesPrismOmega",
         slots = c(
           prism                       = "prism",
           elbo_per_sample             = "numeric",
           posterior.initial.cellState = "jointPost",
           posterior.initial.cellType  = "jointPost",
           reference.update            = "reference",
           posterior.theta_f           = "theta_f",
           control_param               = "list"
         ),
         prototype = list(
           prism                       = new("prism"),
           elbo_per_sample             = numeric(0),
           posterior.initial.cellState = new("jointPost"),
           posterior.initial.cellType  = new("jointPost"),
           reference.update            = NULL,
           posterior.theta_f           = NULL,
           control_param = list(
             gibbs.control = list(chain.length = NULL, burn.in = NULL, thinning = NULL,
                                  n.cores = NULL, seed = NULL, alpha = NULL, beta = NULL),
             opt.control   = list(trace = NULL, maxit = NULL,
                                  optimizer = NA_character_, n.cores = NULL,
                                  update.omega = FALSE),
             update.gibbs  = logical(),
             method        = NA_character_,
             mixture.name  = NA_character_,
             reference.name= NA_character_,
             store.C.input = "value",
             omega.log     = list()
           )
         )
)


# ------------------------------------------------------------------------------
# Parallel helper (base `parallel`, PSOCK)
# ------------------------------------------------------------------------------

#' Run cpu.fun over seq_len(n_items), serially or on a PSOCK cluster.
#'
#' Data objects are exported from `data.envir`; loose helper functions
#' (e.g. the samplers, rdirichlet) are exported from the global environment.
#' clusterExport assigns into each worker's global env, so cpu.fun's environment
#' is set to globalenv() to resolve free variables against the exported copies.
.omega_cluster_lapply <- function(n_items, cpu.fun, n.cores,
                                  data.vars, data.envir, fun.vars) {
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = data.vars, envir = data.envir)
    if (length(fun.vars))
      parallel::clusterExport(cl, varlist = fun.vars, envir = globalenv())
    environment(cpu.fun) <- globalenv()
    parallel::parLapply(cl, seq_len(n_items), cpu.fun)
  } else {
    lapply(seq_len(n_items), cpu.fun)
  }
}


# ------------------------------------------------------------------------------
# Gibbs run-time estimator
# ------------------------------------------------------------------------------

estimate.gibbs.time.Omega <- function(gibbsSampler.obj,
                                      final,
                                      chain.length = 50) {
  
  ref           <- gibbsSampler.obj@reference
  X             <- gibbsSampler.obj@X
  gibbs.control <- gibbsSampler.obj@gibbs.control
  
  gi <- get.gibbs.idx(list(chain.length = chain.length,
                           burn.in      = chain.length * gibbs.control$burn.in /
                             gibbs.control$chain.length,
                           thinning     = gibbs.control$thinning))
  ptm <- proc.time()
  
  if (!final) {
    sample.Z.theta_n.Omega(X_n       = X[1, ],
                           phi       = ref@phi,
                           alpha     = gibbs.control$alpha,
                           Omega     = gibbsSampler.obj@Omega,
                           beta      = gibbs.control$beta,
                           gibbs.idx = gi)
  } else {
    # Final pass: use the Omega sampler when this object carries a tensor,
    # otherwise the base theta-only sampler (timing estimate only).
    if (is(ref, "refPhi")) {
      if (length(gibbsSampler.obj@Omega) > 1L) {
        sample.theta_n.Omega(X_n        = X[1, ],
                             phi        = ref@phi,
                             alpha      = gibbs.control$alpha,
                             Omega_type = gibbsSampler.obj@Omega,
                             beta       = gibbs.control$beta,
                             gibbs.idx  = gi)
      } else {
        sample.theta_n(X_n       = X[1, ],
                       phi       = ref@phi,
                       alpha     = gibbs.control$alpha,
                       gibbs.idx = gi)
      }
    }
  }
  
  total.time     <- as.numeric((proc.time() - ptm)["elapsed"])
  estimated.time <- gibbs.control$chain.length / chain.length * total.time *
    ceiling(nrow(X) / gibbs.control$n.cores) * 2
  current.time   <- Sys.time()
  
  cat("Current time: ",               as.character(current.time),                  "\n")
  cat("Estimated time to complete: ", my_seconds_to_period(estimated.time),        "\n")
  cat("Estimated finishing time: ",   as.character(current.time + estimated.time), "\n")
}


# ------------------------------------------------------------------------------
# Top-level run functions
# ------------------------------------------------------------------------------

# Internal: validate/normalise controls, build the omega fingerprint log, and
# compute the run meta (hashes + run id) through the shared hashing core.  Used
# by run.prism.Omega() and save.run.Omega(); the object accessors run.hashes()/
# run.id() rebuild the same identity from a stored object, so a standalone id
# can never drift from a stamped one.
.run_meta_Omega <- function(prism, beta, method,
                            mixture.name, reference.name,
                            omega.control, Omega, Omega_type,
                            update.gibbs, opt.control, gibbs.control, n.cores,
                            store.C.input = "value") {
  
  stopifnot(is.numeric(beta), length(beta) == 1L, beta >= 0)
  stopifnot(is.logical(update.gibbs), length(update.gibbs) == 1L)
  stopifnot(is.numeric(n.cores),      length(n.cores)      == 1L)
  method        <- match.arg(method, c("dp", "bp", "ip"))
  store.C.input <- match.arg(store.C.input, c("value", "hash"))
  
  if (!"n.cores" %in% names(gibbs.control)) gibbs.control$n.cores <- n.cores
  if (!"n.cores" %in% names(opt.control))   opt.control$n.cores   <- n.cores
  
  opt.control   <- valid.opt.control.Omega(opt.control)
  gibbs.control <- valid.gibbs.control.Omega(gibbs.control)
  gibbs.control$beta <- beta                       # top-level beta is authoritative
  
  omega.control <- .normalize_omega_control(omega.control)
  
  # Auto-demote label when no niche pressure is active anywhere.
  niche.active <- !is.null(omega.control) || !is.null(Omega) || !is.null(Omega_type)
  if (!niche.active) method <- "bp"
  
  # Fingerprint log (also the hash preimage and rebuild key).  Deriving the
  # spec from the log means the run-time hash and any later run.hashes(obj)
  # recompute traverse identical structure.
  omega.log  <- .omega_log_from_control(omega.control, Omega, Omega_type, store.C.input)
  omega.spec <- .omega_spec_from_log(omega.log)
  
  hashes <- .compute_hashes(prism@key, opt.control, gibbs.control,
                            update.gibbs, omega.spec)
  run.id <- make_run_id(method, mixture.name, reference.name,
                        hashes[["opt"]], hashes[["gibbs"]], hashes[["omega"]], beta)
  
  list(method        = method,
       opt.control   = opt.control,
       gibbs.control = gibbs.control,
       omega.control = omega.control,
       update.gibbs  = update.gibbs,
       store.C.input = store.C.input,
       omega.log     = omega.log,
       hashes        = hashes,
       run.id        = run.id)
}


#' Run DynaPrism deconvolution (niche-aware Gibbs sampler)
#'
#' Pure compute: runs the sampler and returns the BayesPrismOmega object,
#' stamped with its run id for provenance.  Does NOT read or write disk --
#' use save.run.Omega() for cache-aware persistence.
#'
#' @param prism a `prism` object (from BayesPrism::new.prism()).
#' @param omega.control named list of build_Omega() arguments (excluding phi).
#'   Nested form list(state = <args>, type = <args>) is supported; a flat list
#'   is taken as the state set; $type defaults to $state.  May be NULL if a
#'   pre-built `Omega` (and, when update.omega = TRUE, `Omega_type`) is supplied.
#' @param beta niche-pressure strength; authoritative over gibbs.control$beta and
#'   used for the run-id suffix.  Applied to both passes.
#' @param method run label for the filename: "dp" (niche-aware, default),
#'   "bp" (auto when no Omega is active) or "ip".
#' @param mixture.name,reference.name user-defined dataset tags for the filename.
#' @param Omega,Omega_type optional pre-built tensors; when supplied the
#'   corresponding build step is skipped.
#' @param update.gibbs run the final (reference-update) Gibbs pass. Default TRUE.
#' @param compute.elbo compute per-sample log joint during the initial pass.
#' @param opt.control list; gains `update.omega` (default FALSE) toggling whether
#'   the final pass is Omega-modulated.
#' @param store.C.input "value" (default) stores C_input in the object so it is
#'   self-describing and rebuildable without re-supplying C; "hash" stores only
#'   C_input's fingerprint (leaner for beta sweeps or heavy C_slg / C_tlg), in
#'   which case rebuild_Omega() requires C_input to be passed back.  Does not
#'   affect identity: the run id is the same either way.
run.prism.Omega <- function(prism,
                            omega.control  = NULL,
                            beta,
                            method         = "dp",
                            mixture.name   = "mix",
                            reference.name = "ref",
                            Omega          = NULL,
                            Omega_type     = NULL,
                            n.cores        = 1,
                            update.gibbs   = TRUE,
                            compute.elbo   = FALSE,
                            gibbs.control  = list(),
                            opt.control    = list(),
                            store.C.input  = c("value", "hash")) {
  
  store.C.input <- match.arg(store.C.input)
  
  meta          <- .run_meta_Omega(prism, beta, method, mixture.name, reference.name,
                                   omega.control, Omega, Omega_type, update.gibbs,
                                   opt.control, gibbs.control, n.cores, store.C.input)
  method        <- meta$method
  opt.control   <- meta$opt.control
  gibbs.control <- meta$gibbs.control
  omega.control <- meta$omega.control
  store.C.input <- meta$store.C.input
  omega.log     <- meta$omega.log
  
  if (prism@phi_cellState@pseudo.min == 0)
    gibbs.control$alpha <- max(1, gibbs.control$alpha)
  
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
    gibbs.control = gibbs.control,
    opt.control   = opt.control,
    update.gibbs  = update.gibbs,
    method        = method,
    mixture.name  = mixture.name,
    reference.name= reference.name,
    store.C.input = store.C.input,
    omega.log     = omega.log
  )
  
  # --- initial Gibbs: cell-state level, Omega-modulated -----------------------
  gibbsSampler.ini.cs <- new("gibbsSamplerOmega",
                             reference     = prism@phi_cellState,
                             X             = prism@mixture,
                             Omega         = Omega,
                             gibbs.control = gibbs.control)
  
  ini_result       <- run.gibbs.Omega(gibbsSampler.ini.cs,
                                      final        = FALSE,
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
  # updateReference() re-validates opt.control with the stock base validator,
  # which rejects unknown names, so strip our custom flag(s) before the call.
  opt.control.base <- opt.control
  opt.control.base$update.omega <- NULL
  
  psi <- updateReference(Z           = jointPost.ini.ct@Z,
                         phi_prime   = prism@phi_cellType,
                         map         = prism@map,
                         key         = prism@key,
                         opt.control = opt.control.base)
  
  # --- final Gibbs: cell-type level; Omega-modulated iff update.omega ---------
  if (isTRUE(opt.control$update.omega) && !is.null(Omega_type)) {
    gibbsSampler.update <- new("gibbsSamplerOmega",
                               reference     = psi,
                               X             = prism@mixture,
                               Omega         = Omega_type,
                               gibbs.control = gibbs.control)
    theta_f <- run.gibbs.Omega(gibbsSampler.update, final = TRUE, compute.elbo = FALSE)
  } else {
    gibbsSampler.update <- new("gibbsSampler",
                               reference     = psi,
                               X             = prism@mixture,
                               gibbs.control = gibbs.control)
    theta_f <- run.gibbs(gibbsSampler.update, final = TRUE)
  }
  
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
# Cache-aware wrapper: run and persist (or load if already present)
# ------------------------------------------------------------------------------

#' Run DynaPrism and save the result under its run id (or load the cache)
#'
#' Computes the run id; if "<cache.dir>/<run.id>.rds" already exists (and
#' overwrite = FALSE) loads and returns it, otherwise runs run.prism.Omega()
#' and saveRDS() the result there.  Accepts every run.prism.Omega() argument
#' plus cache.dir and overwrite.
#'
#' Note on identity: the run id (and hence the cache key) is independent of
#' compute.elbo and store.C.input -- neither changes theta/Z.  A cached result
#' therefore satisfies a later request regardless of those settings; use
#' overwrite = TRUE (or run.prism.Omega() directly) if you specifically need
#' the ELBO recomputed or C_input re-stored differently.
save.run.Omega <- function(prism,
                           omega.control  = NULL,
                           beta,
                           method         = "dp",
                           mixture.name   = "mix",
                           reference.name = "ref",
                           Omega          = NULL,
                           Omega_type     = NULL,
                           n.cores        = 1,
                           update.gibbs   = TRUE,
                           compute.elbo   = FALSE,
                           gibbs.control  = list(),
                           opt.control    = list(),
                           store.C.input  = c("value", "hash"),
                           cache.dir,
                           overwrite      = FALSE) {
  
  stopifnot(!missing(cache.dir), is.character(cache.dir), length(cache.dir) == 1L)
  store.C.input <- match.arg(store.C.input)
  
  run_args <- list(prism = prism, omega.control = omega.control, beta = beta,
                   method = method, mixture.name = mixture.name,
                   reference.name = reference.name, Omega = Omega,
                   Omega_type = Omega_type, n.cores = n.cores,
                   update.gibbs = update.gibbs, compute.elbo = compute.elbo,
                   gibbs.control = gibbs.control, opt.control = opt.control,
                   store.C.input = store.C.input)
  
  # Filename comes from the shared meta helper (same identity run.prism.Omega
  # computes internally), so the cache key matches what the run would stamp.
  run.id <- .run_meta_Omega(prism, beta, method, mixture.name, reference.name,
                            omega.control, Omega, Omega_type, update.gibbs,
                            opt.control, gibbs.control, n.cores, store.C.input)$run.id
  
  if (!dir.exists(cache.dir)) dir.create(cache.dir, recursive = TRUE)
  cache.file <- file.path(cache.dir, paste0(run.id, ".rds"))
  
  if (file.exists(cache.file) && !overwrite) {
    message("Loading cached result: ", cache.file)
    return(readRDS(cache.file))
  }
  
  res <- do.call(run.prism.Omega, run_args)
  message("Saving result: ", cache.file)
  saveRDS(res, cache.file)
  res
}


# ------------------------------------------------------------------------------
# Identity accessors and Omega rebuild
# ------------------------------------------------------------------------------

#' Compute the opt/gibbs/omega hash triplet from loose control lists
#'
#' Mirror of the identity used inside a run, without sampling.  Provide `key`
#' (the malignant cell-type key, prism@key) and `update.gibbs` to match a real
#' run; the defaults are for quick inspection.  Covers control-built runs; for a
#' pre-built-Omega run recompute from the object with run.hashes(obj) instead.
control.hashes <- function(gibbs.control = list(), opt.control = list(),
                           omega.control = NULL, key = NA, update.gibbs = TRUE) {
  if (!"n.cores" %in% names(opt.control))   opt.control$n.cores   <- 1
  if (!"n.cores" %in% names(gibbs.control)) gibbs.control$n.cores <- 1
  opt.control   <- valid.opt.control.Omega(opt.control)
  gibbs.control <- valid.gibbs.control.Omega(gibbs.control)
  
  omega.control <- .normalize_omega_control(omega.control)
  omega.log     <- .omega_log_from_control(omega.control, NULL, NULL, "hash")
  omega.spec    <- .omega_spec_from_log(omega.log)
  
  .compute_hashes(key, opt.control, gibbs.control, update.gibbs, omega.spec)
}


#' Recompute the opt/gibbs/omega hash triplet from a BayesPrismOmega object
#'
#' Uses only what the object stores (validated controls, prism@key, and the
#' omega.log fingerprints), so no hashes need to be persisted.
run.hashes <- function(obj) {
  stopifnot(is(obj, "BayesPrismOmega"))
  cp <- obj@control_param
  .compute_hashes(obj@prism@key, cp$opt.control, cp$gibbs.control,
                  cp$update.gibbs, .omega_spec_from_log(cp$omega.log))
}


#' Recompute the run id (filename stem) from a BayesPrismOmega object
#'
#' Equals the filename save.run.Omega() wrote, since both converge on the same
#' controls and omega spec.
run.id <- function(obj) {
  stopifnot(is(obj, "BayesPrismOmega"))
  cp <- obj@control_param
  h  <- run.hashes(obj)
  make_run_id(cp$method, cp$mixture.name, cp$reference.name,
              h[["opt"]], h[["gibbs"]], h[["omega"]], cp$gibbs.control$beta)
}


#' Rebuild a state- or type-level Omega tensor from a stored run
#'
#' Reconstructs Omega using the reference phi carried in obj@prism and the
#' stored omega.log.  The invariant priors are always re-supplied and verified
#' against the stored fingerprints; C_input is taken from the object when it was
#' stored (store.C.input == "value"), otherwise it must be re-supplied and is
#' verified against the stored fingerprint.
#'
#' @param which "state" (initial pass) or "type" (reference-update pass).
rebuild_Omega <- function(obj, ligand_target_matrix, lr_network,
                          C_input = NULL, which = c("state", "type")) {
  stopifnot(is(obj, "BayesPrismOmega"))
  which     <- match.arg(which)
  omega.log <- obj@control_param$omega.log
  
  if (!(length(names(omega.log)) && all(names(omega.log) %in% c("state", "type"))))
    stop("This run has no build recipe (Omega was supplied pre-built); nothing to rebuild.")
  entry <- omega.log[[which]]
  
  if (!identical(.digest_full(ligand_target_matrix), entry$ltm_hash))
    stop("ligand_target_matrix does not match the fingerprint stored for this run.")
  if (!identical(.digest_full(lr_network), entry$lr_hash))
    stop("lr_network does not match the fingerprint stored for this run.")
  
  C <- if (!is.null(entry$C_input)) {
    entry$C_input                                   # stored value
  } else {
    if (is.null(C_input))
      stop("C_input was not stored for this run (store.C.input = 'hash'); ",
           "re-supply it via the C_input argument.")
    if (!identical(.digest_full(C_input), entry$C_input_hash))
      stop("Supplied C_input does not match the fingerprint stored for this run.")
    C_input
  }
  
  control <- list(ligand_target_matrix = ligand_target_matrix,
                  C_input              = C,
                  lr_network           = lr_network,
                  mask_threshold       = entry$mask_threshold,
                  masking_type         = entry$masking_type)
  phi <- if (which == "state") obj@prism@phi_cellState@phi else obj@prism@phi_cellType@phi
  build_Omega_from_control(phi, control)
}


# ------------------------------------------------------------------------------
# Gibbs dispatch  (initial + final off the `final` flag, mirroring run.gibbs)
# ------------------------------------------------------------------------------

#' Initial pass returns list(jointPost, elbo_per_sample).
#' Final pass returns a thetaPost.  The modulating tensor is read from
#' gibbsSampler.obj@Omega (state-level for ini, type-level for final).
run.gibbs.Omega <- function(gibbsSampler.obj,
                            final,
                            if.estimate  = TRUE,
                            compute.elbo = FALSE) {
  
  if (final) cat("Run Gibbs sampling using updated reference ...\n")
  else       cat("Run Gibbs sampling...\n")
  
  if (if.estimate)
    estimate.gibbs.time.Omega(gibbsSampler.obj = gibbsSampler.obj, final = final)
  
  if (is(gibbsSampler.obj@reference, "refPhi")) {
    if (!final)
      return(run.gibbs.refPhi.ini.Omega(gibbsSampler.obj = gibbsSampler.obj,
                                        compute.elbo     = compute.elbo))
    else
      return(run.gibbs.refPhi.final.Omega(gibbsSampler.obj = gibbsSampler.obj,
                                          compute.elbo     = compute.elbo))
  }
  
  if (is(gibbsSampler.obj@reference, "refTumor"))
    return(run.gibbs.refTumor(gibbsSampler.obj = gibbsSampler.obj))
}


# ------------------------------------------------------------------------------
# Initial Gibbs dispatcher -> list(jointPost, elbo_per_sample)
# ------------------------------------------------------------------------------

run.gibbs.refPhi.ini.Omega <- function(gibbsSampler.obj,
                                       compute.elbo = FALSE) {
  
  phi           <- gibbsSampler.obj@reference@phi
  X             <- gibbsSampler.obj@X
  Omega         <- gibbsSampler.obj@Omega
  gibbs.control <- gibbsSampler.obj@gibbs.control
  alpha         <- gibbs.control$alpha
  beta          <- gibbs.control$beta
  gibbs.idx     <- get.gibbs.idx(gibbs.control)
  seed          <- gibbs.control$seed
  n.cores       <- gibbs.control$n.cores
  tmp.dir       <- tempdir(check = TRUE)
  
  cat("Start run...\n")
  
  cpu.fun <- function(n) {
    if (!is.null(seed)) set.seed(seed)
    file.name.X_n <- file.path(tmp.dir, paste0("mixture_", n, ".rdata"))
    load(file.name.X_n)                       # -> X_n
    sample.Z.theta_n.Omega(X_n          = X_n,
                           phi          = phi,
                           alpha        = alpha,
                           Omega        = Omega,
                           beta         = beta,
                           gibbs.idx    = gibbs.idx,
                           compute.elbo = compute.elbo)
  }
  
  gibbs.list <- .omega_cluster_lapply(
    n_items   = nrow(X),
    cpu.fun   = cpu.fun,
    n.cores   = n.cores,
    data.vars = c("phi", "alpha", "Omega", "beta", "gibbs.idx", "seed",
                  "compute.elbo", "tmp.dir"),
    data.envir = environment(),
    fun.vars   = c("sample.Z.theta_n.Omega", "rdirichlet"))
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


# ------------------------------------------------------------------------------
# Final Gibbs dispatcher (cell-type level, Omega-modulated) -> thetaPost
# ------------------------------------------------------------------------------

run.gibbs.refPhi.final.Omega <- function(gibbsSampler.obj,
                                         compute.elbo = FALSE) {
  
  phi           <- gibbsSampler.obj@reference@phi
  X             <- gibbsSampler.obj@X
  Omega_type    <- gibbsSampler.obj@Omega
  gibbs.control <- gibbsSampler.obj@gibbs.control
  alpha         <- gibbs.control$alpha
  beta          <- gibbs.control$beta
  gibbs.idx     <- get.gibbs.idx(gibbs.control)
  seed          <- gibbs.control$seed
  n.cores       <- gibbs.control$n.cores
  tmp.dir       <- tempdir(check = TRUE)
  
  cat("Start run...\n")
  
  cpu.fun <- function(n) {
    if (!is.null(seed)) set.seed(seed)
    file.name.X_n <- file.path(tmp.dir, paste0("mixture_", n, ".rdata"))
    load(file.name.X_n)                       # -> X_n
    sample.theta_n.Omega(X_n        = X_n,
                         phi        = phi,
                         alpha      = alpha,
                         Omega_type = Omega_type,
                         beta       = beta,
                         gibbs.idx  = gibbs.idx)
  }
  
  gibbs.list <- .omega_cluster_lapply(
    n_items   = nrow(X),
    cpu.fun   = cpu.fun,
    n.cores   = n.cores,
    data.vars = c("phi", "alpha", "Omega_type", "beta", "gibbs.idx", "seed", "tmp.dir"),
    data.envir = environment(),
    fun.vars   = c("sample.theta_n.Omega", "rdirichlet"))
  if (n.cores <= 1) cat("\n")
  
  newThetaPost(bulkID     = rownames(X),
               cellType   = rownames(phi),
               gibbs.list = gibbs.list)
}


# ------------------------------------------------------------------------------
# Sampler: joint Z + theta  (initial pass, cell-state level)  [unchanged]
# ------------------------------------------------------------------------------

#' ELBO formula: full log joint (dropping constant Sum_g log(X_ng!))
#'
#'   log P(Z_n^i, theta_n^i | X_n, phi, Omega, beta) ~=
#'       Sum_k (alpha-1) log theta_nk^i          <- log Dirichlet prior
#'     + Sum_g Sum_k z_nkg^i log p~_nkg^i         <- log multinomial likelihood
#'     - Sum_g Sum_k log(z_nkg^i !)               <- combinatorial term
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
  elbo_sum     <- 0
  
  phi_mat <- as.matrix(phi)
  X_vec   <- as.numeric(X_n[G_names])
  
  Omega_flat <- matrix(Omega[K_names, K_names, G_shared_names],
                       nrow = K_dim,
                       ncol = K_dim * length(G_shared_names))
  
  for (i in seq_len(max(gibbs.idx))) {
    
    niche_pressure_vec <- as.numeric(theta_n.i %*% Omega_flat)
    niche_pressure_mat <- matrix(niche_pressure_vec, nrow = K_dim,
                                 ncol = length(G_shared_names))
    
    prob.mat <- phi_mat * theta_n.i
    prob.mat[, G_shared_idx] <- prob.mat[, G_shared_idx] * exp(beta * niche_pressure_mat)
    
    for (g in seq_len(G_dim))
      Z_n.i[g, ] <- rmultinom(n = 1, size = X_vec[g], prob = prob.mat[, g])
    
    Z_nk.i    <- colSums(Z_n.i)
    theta_n.i <- rdirichlet(alpha = Z_nk.i + alpha)
    
    if (i %in% gibbs.idx) {
      Z_n.sum      <- Z_n.sum      + Z_n.i
      theta_n.sum  <- theta_n.sum  + theta_n.i
      theta_n2.sum <- theta_n2.sum + theta_n.i^2
      
      if (compute.elbo) {
        col_sums  <- colSums(prob.mat)
        col_sums[col_sums == 0] <- 1
        prob_norm <- sweep(prob.mat, 2, col_sums, "/")
        prob_norm[prob_norm == 0] <- .Machine$double.eps
        
        log_lik_i   <- sum(Z_n.i * log(prob_norm))
        log_prior_i <- sum((alpha - 1) * log(pmax(theta_n.i, .Machine$double.eps)))
        log_coef_i  <- -sum(lfactorial(Z_n.i))
        
        elbo_sum <- elbo_sum + log_lik_i + log_prior_i + log_coef_i
      }
    }
    if ((i %% 50) == 0) gc()
  }
  
  samples.size   <- length(gibbs.idx)
  Z_n            <- Z_n.sum      / samples.size
  theta_n        <- theta_n.sum  / samples.size
  theta.cv_n     <- sqrt(theta_n2.sum / samples.size - theta_n^2) / theta_n
  gibbs.constant <- elbo_sum / samples.size
  
  list(Z_n            = Z_n,
       theta_n        = theta_n,
       theta.cv_n     = theta.cv_n,
       gibbs.constant = gibbs.constant)
}


# ------------------------------------------------------------------------------
# Sampler: theta only  (final pass, cell-type level)  [unchanged]
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