avg_per_label <- function(sce, labels, assay = "logcounts"){
  stopifnot(assay %in% names(assays(sce)))
  label_names <- unique(as.vector(labels))
  avg_expr <- sapply(label_names, function(cl){
    # Subset to the cell type
    sub_sce <- sce[, labels == cl]
    # Calculate the mean of selected assay across cells for each gene
    rowMeans(as.matrix(assay(sub_sce, assay)))
  })
  return(t(avg_expr))
}



aggregate_stratified_counts <- function(counts_mat, sample_names, label_names,
                                        include_pseudobulk = FALSE, 
                                        target_sum = 1e6, log_base = exp(1)) {
  # Aggregate counts by Sample_State pairs
  group_id <- paste0(sample_names, ":::", label_names)
  pseudobulk <- rowsum(as(counts_mat, "matrix"), group_id) # [Sample_State x Genes]
  
  # Get abundance weights (number of cells per Sample_State)
  ncells <- table(group_id)
  
  # Get shifted logarithm 
  lib_sizes <- rowSums(pseudobulk)
  lib_sizes[lib_sizes == 0] <- 1
  
  expr_shifted_log <- log(sweep(pseudobulk, 1, lib_sizes, "/") * target_sum + 1, base = log_base)
  expr_shifted_log[is.na(expr_shifted_log)] <- 0
  
  meta <- data.frame(
    id     = rownames(expr_shifted_log),
    Sample = sub(":::.*", "", rownames(expr_shifted_log)),
    Label  = sub(".*:::", "", rownames(expr_shifted_log)),
    ncells = as.numeric(ncells[rownames(expr_shifted_log)]),
    stringsAsFactors = FALSE
  )
  if (include_pseudobulk){
    return(list(logexpr = t(expr_shifted_log), pseudoexp = t(pseudobulk), meta = meta))
  } else {
    return(list(logexpr = t(expr_shifted_log), meta = meta))
  }
}



compute_C_tensor <- function(expr_mat, meta, lig_genes, tar_genes, min_samples = 3) {
  samples <- unique(meta$Sample)
  all_labels <- unique(meta$Label) 
  
  lig_genes <- intersect(lig_genes, rownames(expr_mat))
  tar_genes <- intersect(tar_genes, rownames(expr_mat))
  
  C_tsr <- array(0, dim = c(length(all_labels), length(lig_genes), length(tar_genes)),
                 dimnames = list(Receiver = all_labels, Ligand = lig_genes, Gene = tar_genes))
  
  for (receiver in all_labels) {
    # Check if receiver exists in enough samples
    present_samples <- unique(meta$Sample[meta$Label == receiver])
    n_present <- length(present_samples)
    
    if (n_present < min_samples) {
      message(sprintf("Skipping %s: Only found in %d samples (minimum required: %d)", 
                      receiver, n_present, min_samples))
      next
    }
    
    message(sprintf("Computing Context for %s (%d samples)...", receiver, n_present))
    
    eff_niche_lig_mat <- matrix(NA, nrow=length(samples), ncol=length(lig_genes), dimnames=list(samples, lig_genes))
    rec_tar_mat <- matrix(NA, nrow=length(samples), ncol=length(tar_genes), dimnames=list(samples, tar_genes))
    
    for(smp in samples) {
      # Identify the Niche for this sample
      smp_niche_meta <- meta[meta$Sample == smp & meta$Label != receiver, ]
      
      if(nrow(smp_niche_meta) > 0) {
        w <- smp_niche_meta$ncells
        smp_lig_expr <- expr_mat[lig_genes, smp_niche_meta$id, drop=FALSE]
        eff_niche_lig_mat[smp, ] <- as.numeric(smp_lig_expr %*% w) / sum(w)
      }
      
      # Identify the Receiver for this sample
      rec_id <- meta$id[meta$Sample == smp & meta$Label == receiver]
      if(length(rec_id) > 0) {
        rec_tar_mat[smp, ] <- expr_mat[tar_genes, rec_id]
      }
    }
    
    # Only keep genes/ligands that vary across the samples where the receiver exists
    v_lig <- apply(eff_niche_lig_mat, 2, function(x) sd(x, na.rm=TRUE))
    v_tar <- apply(rec_tar_mat, 2, function(x) sd(x, na.rm=TRUE))
    
    valid_lig_idx <- which(!is.na(v_lig) & v_lig > 0)
    valid_tar_idx <- which(!is.na(v_tar) & v_tar > 0)
    
    if(length(valid_lig_idx) > 0 && length(valid_tar_idx) > 0) {
      # Spearman Correlation
      cor_block <- suppressWarnings(cor(eff_niche_lig_mat[, valid_lig_idx, drop=FALSE], 
                                        rec_tar_mat[, valid_tar_idx, drop=FALSE], 
                                        method = "spearman", 
                                        use = "pairwise.complete.obs"))
      
      C_tsr[receiver, names(valid_lig_idx), names(valid_tar_idx)] <- cor_block
    }
  }
  return(C_tsr)
}



get_trans_correlations_parallel <- function(expr_mat, meta, lig_genes, tar_genes, n_cores = 6, rho_threshold = 0.1) {
  
  lig_genes <- intersect(lig_genes, rownames(expr_mat))
  tar_genes <- intersect(tar_genes, rownames(expr_mat))
  
  abundance_matrix <- as.matrix(unclass(table(meta$Sample, meta$Label)))
  samples <- rownames(abundance_matrix)
  all_labels <- colnames(abundance_matrix)
  
  message(sprintf("Initializing snowfall with %d cores...", n_cores))
  sfInit(parallel = TRUE, cpus = n_cores)
  
  sfLibrary(stats)
  sfLibrary(data.table)
  sfExport("expr_mat", "meta", "lig_genes", "tar_genes", 
           "abundance_matrix", "samples", "all_labels", "rho_threshold")
  
  cor_worker <- function(rec_s) {
    niche_states <- all_labels[all_labels != rec_s]
    num_samples <- length(samples)
    
    eff_niche_lig_mat <- matrix(NA, nrow=num_samples, ncol=length(lig_genes), dimnames=list(samples, lig_genes))
    rec_tar_mat <- matrix(NA, nrow=num_samples, ncol=length(tar_genes), dimnames=list(samples, tar_genes))
    
    for(smp in samples){
      w_n_all <- abundance_matrix[smp, niche_states]
      active_senders <- names(w_n_all[w_n_all > 0])
      if(length(active_senders) > 0) {
        smp_meta <- meta[meta$Sample == smp & meta$Label %in% active_senders, ]
        if(nrow(smp_meta) > 0) {
          w_vals <- w_n_all[smp_meta$Label]
          eff_niche_lig_mat[smp, ] <- as.numeric(expr_mat[lig_genes, smp_meta$id, drop=FALSE] %*% w_vals) / sum(w_n_all)
        }
      }
      rec_id <- meta$id[meta$Sample == smp & meta$Label == rec_s]
      if(length(rec_id) > 0) rec_tar_mat[smp, ] <- expr_mat[tar_genes, rec_id]
    }
    
    v_lig <- apply(eff_niche_lig_mat, 2, function(x) sd(x, na.rm=TRUE))
    v_tar <- apply(rec_tar_mat, 2, function(x) sd(x, na.rm=TRUE))
    
    valid_ligs <- names(v_lig[!is.na(v_lig) & v_lig > 0])
    valid_tars <- names(v_tar[!is.na(v_tar) & v_tar > 0])
    
    if(length(valid_ligs) < 2 || length(valid_tars) < 2) return(NULL)
    
    # Compute Correlation Matrix
    rho_mtx <- cor(eff_niche_lig_mat[, valid_ligs, drop=FALSE], 
                   rec_tar_mat[, valid_tars, drop=FALSE], 
                   method = "spearman", 
                   use = "pairwise.complete.obs")
    
    # Identify significant indices only
    hit_idx <- which(abs(rho_mtx) > rho_threshold, arr.ind = TRUE)
    if(nrow(hit_idx) == 0) return(NULL)
    
    # Extract only the values we need
    sig_rhos <- rho_mtx[hit_idx]
    
    # Calculate P-values only for the hits
    n_eff <- nrow(eff_niche_lig_mat)
    t_vals <- sig_rhos * sqrt((n_eff - 2) / (1 - sig_rhos^2))
    sig_pvals <- 2 * pt(-abs(t_vals), df = n_eff - 2)
    
    # Build the data.table
    res <- data.table(
      Receiver = rec_s,
      Ligand   = rownames(rho_mtx)[hit_idx[, 1]],
      Gene     = colnames(rho_mtx)[hit_idx[, 2]],
      Rho      = sig_rhos,
      Pval     = sig_pvals
    )
    
    # Garbage collect within worker
    rm(rho_mtx, hit_idx, t_vals)
    return(res)
  }
  
  message("Running trans-correlations...")
  results_list <-  
    sfStop()
  
  all_corrs_dt <- rbindlist(results_list)
  
  return(all_corrs_dt)
}



compute_stouffer_consensus <- function(cor_dt) {
  message("Calculating Stouffer consensus and consistency statistics...")
  
  # Ensure the input is a data.table
  if (!is.data.table(cor_dt)) setDT(cor_dt)
  
  # Group by Ligand and Target (Gene) to calculate statistics across all receivers
  consensus_dt <- cor_dt[, {
    # Calculate individual Z-scores (Stabilized)
    # Using lower.tail = FALSE and 1e-15 to prevent Infinity
    zs <- qnorm(pmax(Pval, 1e-15) / 2, lower.tail = FALSE) * sign(Rho)
    zs[!is.finite(zs)] <- 0
    
    # Count observations
    n_val <- .N
    
    # Consistency calculation
    pos_count <- sum(Rho > 0, na.rm = TRUE)
    neg_count <- sum(Rho < 0, na.rm = TRUE)
    consist_val <- max(pos_count, neg_count) / n_val
    
    # Return the aggregated row
    .(
      n_labels    = n_val,
      StoufferZ   = sum(zs, na.rm = TRUE) / sqrt(n_val),
      Consistency = consist_val,
      AvgRho      = mean(Rho, na.rm = TRUE)
    )
  }, by = .(Ligand, Target = Gene)] # Renaming Gene to Target per your request
  
  # Sort by absolute Z-score for easy inspection
  setorder(consensus_dt, -abs(StoufferZ))
  
  return(consensus_dt)
}



get_standardized_intensity <- function(phi_mtx, target_sum = 1e4, log_base = exp(1)) {
  # Scaling to 10k library size approximates actual molecule counts per cell
  # Using log_base = exp(1) for consistency with the Gibbs Sampler exponential
  return(log(phi_mtx * target_sum + 1, base = log_base))
}



build_Omega <- function(phi, ligand_target_matrix, C_tsr, lr_network, 
                        mask_threshold = 0.1, masking_type = "soft"){
  
  C_tsr[is.na(C_tsr)] <- 0
  
  # Standardize the Reference Intensity by converting the linear probabilities (phi) to log-intensities
  avg_expr <- get_standardized_intensity(phi)
  
  label_names <- rownames(avg_expr)
  match.arg(masking_type, c("hard", "soft"))
  
  # Identify common ligands and target genes
  common_ligands <- intersect(colnames(ligand_target_matrix), colnames(avg_expr))
  common_genes   <- intersect(rownames(ligand_target_matrix), colnames(avg_expr))
  
  # L = [Sender x Ligand] (Log-intensities)
  L_matrix <- as.matrix(avg_expr[, common_ligands])
  
  # W = [Ligand x TargetGene] (NicheNet Priors)
  W_matrix <- t(as.matrix(ligand_target_matrix[common_genes, common_ligands]))
  
  # Build Receptor Responsiveness Mask
  relevant_lr = lr_network %>% 
    filter(from %in% common_ligands) %>%
    filter(to %in% colnames(avg_expr))
  
  R_matrix <- matrix(0, nrow = length(label_names), ncol = length(common_ligands),
                     dimnames = list(label_names, common_ligands))
  
  for(s in label_names) {
    for(lig in common_ligands) {
      receptors = relevant_lr %>% filter(from == lig) %>% pull(to)
      
      if(length(receptors) > 0) {
        total_phi_receptor <- get_standardized_intensity(sum(phi[s, receptors]))
        
        if(total_phi_receptor > mask_threshold) {
          if(masking_type == "hard"){
            R_matrix[s, lig] <- 1
            } else {
              R_matrix[s, lig] <- total_phi_receptor
            }
        } 
      }
    }
  }
  
  # Assemble the Omega Tensor [Sender x Receiver x Gene]
  Omega <- array(0, dim = c(length(label_names), length(label_names), length(common_genes)),
                 dimnames = list(Sender = label_names, Receiver = label_names, Gene = common_genes))
  
  for(s in label_names) {
    # Get the Context slice for this specific receiver [Ligands x Genes]
    C_slice <- C_tsr[s, common_ligands, common_genes]
    
    # Weight the NicheNet prior by the observed context
    signed_W_slice <- W_matrix * C_slice 
    
    for(k in label_names) {
      if(k == s) next 
      
      # Masked Ligand Input: (Sender Expression * Receiver Sensitivity)
      masked_ligand_input <- L_matrix[k, ] * R_matrix[s, ]
      
      # Final Potential calculation via dot product
      # [1 x Ligands] %*% [Ligands x Genes]
      Omega[k, s, ] <- as.numeric(masked_ligand_input %*% signed_W_slice)
    }
  }
  
  # Final Normalization (Zero-Preserving)
  valid_vals <- Omega[Omega != 0]
  if(length(valid_vals) > 0) {Omega <- Omega / sd(valid_vals)}
  
  return(Omega)
}