sc_exp_tensor <- function(assay_mtx, sample_vec, cell_vec, fun = "sum"){
  # 1. Validation Logic
  if (ncol(assay_mtx) != length(sample_vec) | ncol(assay_mtx) != length(cell_vec)) {
    stop("Error: Metadata vector length doesn't match number of assay columns.")
  }
  
  # 2. Map Cells to Groups (Sample_State)
  group_id <- paste(sample_vec, cell_vec, sep = "_")
  group_map <- model.matrix(~ 0 + factor(group_id))
  colnames(group_map) <- levels(factor(group_id))
  
  # 3. Aggregate Expression (Sum)
  # (Genes x Cells) %*% (Cells x Groups) = (Genes x Groups)
  exp_mtx <- assay_mtx %*% group_map
  
  # 4. Handle "Mean" Logic
  if (fun == "mean") {
    # colSums of the binary map gives the number of cells per group
    group_cell_counts <- colSums(group_map)
    # Divide each gene's sum by the number of cells in that group
    # sweep() is efficient for matrix-vector division
    exp_mtx <- sweep(exp_mtx, 2, group_cell_counts, "/")
  }
  
  # 5. Tensor Initialization
  samples <- unique(sample_vec)
  genes   <- rownames(assay_mtx)
  cells   <- unique(cell_vec)
  
  exp_tensor <- array(0, 
                      dim = c(length(samples), length(genes), length(cells)),
                      dimnames = list(samples, genes, cells))
  
  # 6. Populate the Tensor
  groups <- colnames(exp_mtx)
  for (c in cells) {
    for (n in samples) {
      gr <- paste(n, c, sep = "_")
      if (gr %in% groups) {
        # Note: Fixed the variable names to match exp_tensor and exp_mtx
        exp_tensor[n, , c] <- exp_mtx[, gr]
      }
    }
  }
  
  return(exp_tensor)
}



# Takes a fraction estimation matrix of #samples x #states or #sample x #type and a ground truth matrix
eval_frac <- function(inf_frac_mtx, true_frac_mtx, metrics = c("MAE", "SCorr", "CCorr", "MAECorr")){
  # Align samples and states
  shared_samples <- intersect(rownames(true_frac_mtx), rownames(inf_frac_mtx))
  shared_states  <- intersect(colnames(true_frac_mtx), colnames(inf_frac_mtx))
  
  # Subset to aligned data
  E <- inf_frac_mtx[shared_samples, shared_states]
  G <- true_frac_mtx[shared_samples, shared_states]
  
  # Initialize results vector
  results <- rep(NA, length(metrics))
  names(results) <- metrics
  
  # --- Mean Absolute Error (Standard) ---
  if ("MAE" %in% metrics){
    results["MAE"] <- mean(abs(E - G))
  }
  
  # --- Sample-wise Spearman Correlation (SCorr) ---
  if ("SCorr" %in% metrics){
    sample_corrs <- sapply(1:nrow(E), function(i) {
      cor(E[i,], G[i,], method = "spearman", use = "complete.obs")
    })
    results["SCorr"] <- mean(sample_corrs, na.rm = TRUE)
  }
  
  # --- Cell-type-wise Correlation (CCorr) ---
  if ("CCorr" %in% metrics){
    type_corrs <- sapply(1:ncol(E), function(j) {
      cor(E[,j], G[,j], method = "pearson", use = "complete.obs")
    })
    results["CCorr"] <- mean(type_corrs, na.rm = TRUE)
  }
  
  # --- MAE of Correlation Matrices (MAECorr) ---
  if ("MAECorr" %in% metrics) {
    # Calculate sample-pairwise Pearson correlation matrices
    # We transpose (t) because cor() correlates columns (samples)
    P_E <- cor(t(E), method = "pearson", use = "pairwise.complete.obs")
    P_G <- cor(t(G), method = "pearson", use = "pairwise.complete.obs")
    
    # Calculate the mean of absolute differences between the two n x n matrices
    results["MAECorr"] <- mean(abs(P_E - P_G), na.rm = TRUE)
  }
  
  return(results)
}

eval_CTSE <- function(inf_exp_tsr, true_exp_tsr, metrics = c("ExpSCorr", "CvSpe", "ExpSpe", "ExpMAE", "ExpRMSE"), markers_list = NULL){
  
  # 1. FORCED ALIGNMENT - Ensures genes and states are identical across tensors
  shared_samples <- intersect(dimnames(true_exp_tsr)[[1]], dimnames(inf_exp_tsr)[[1]])
  shared_genes   <- intersect(dimnames(true_exp_tsr)[[2]], dimnames(inf_exp_tsr)[[2]])
  shared_cells   <- intersect(dimnames(true_exp_tsr)[[3]], dimnames(inf_exp_tsr)[[3]])
  
  message(paste("Aligned Data: Samples:", length(shared_samples), 
                "| Genes:", length(shared_genes), 
                "| States:", length(shared_cells)))
  
  # Subset and CAST to ensure we are working with standard arrays
  E <- inf_exp_tsr[shared_samples, shared_genes, shared_cells]
  G <- true_exp_tsr[shared_samples, shared_genes, shared_cells]
  
  G_dim <- length(shared_genes)
  K_dim <- length(shared_cells)
  results_list <- list()
  
  # --- 1. ExpSCorr ---
  if ("ExpSCorr" %in% metrics){
    sexp_matrix <- matrix(NA, nrow = G_dim, ncol = K_dim, dimnames = list(shared_genes, shared_cells))
    for(k in shared_cells){
      E_k <- E[,,k]; G_k <- G[,,k]
      valid_samples <- which(rowSums(G_k) > 0)
      if(length(valid_samples) < 2) next
      
      # Use a flat sapply to avoid deep nesting
      sexp_matrix[, k] <- suppressWarnings(sapply(1:G_dim, function(g){
        ve <- E_k[valid_samples, g]; vg <- G_k[valid_samples, g]
        if(sd(ve) == 0 || sd(vg) == 0) return(NA)
        return(cor(ve, vg, method = "spearman"))
      }))
    }
    results_list[["ExpSCorr_Matrix"]]  <- sexp_matrix
    results_list[["ExpSCorr_Summary"]] <- colMeans(sexp_matrix, na.rm = TRUE)
  }
  
  # --- 2. CvSpe (RE-OPTIMIZED FOR STACK SAFETY) ---
  if ("CvSpe" %in% metrics){
    cvspe_matrix <- matrix(NA, nrow = G_dim, ncol = K_dim, dimnames = list(shared_genes, shared_cells))
    
    for(g in 1:G_dim){
      # Simplify: Get matrix for gene g
      E_g <- E[, g, ] 
      
      # Step-by-step to avoid stack buildup
      cor_mtx <- suppressWarnings(cor(E_g, method = "pearson", use = "pairwise.complete.obs"))
      
      # Check if correlation was possible
      if(!all(is.na(cor_mtx))){
        # Count valid non-diagonal entries
        n_comps <- rowSums(!is.na(cor_mtx)) - 1
        sum_c   <- rowSums(cor_mtx, na.rm = TRUE) - 1
        cvspe_matrix[g, ] <- ifelse(n_comps > 0, sum_c / n_comps, NA)
      }
      
      # EVERY 1000 GENES: Clear the stack/memory
      if(g %% 1000 == 0) gc(verbose = FALSE)
    }
    
    results_list[["CvSpe_Matrix"]] <- cvspe_matrix
    
    if(!is.null(markers_list)){
      results_list[["CvSpe_Summary"]] <- sapply(shared_cells, function(k) {
        m <- intersect(markers_list[[k]], shared_genes)
        mean(cvspe_matrix[m, k], na.rm = TRUE)
      })
    }
  }
  
  # --- 3. ExpSpe (CCC) ---
  if ("ExpSpe" %in% metrics){
    # Use colMeans to prevent large apply overhead
    Z_bar_inf  <- apply(E, c(2, 3), mean, na.rm = TRUE)
    Z_bar_true <- apply(G, c(2, 3), mean, na.rm = TRUE)
    
    inf_log <- log1p(Z_bar_inf); true_log <- log1p(Z_bar_true)
    
    results_list[["ExpSpe_Vector"]] <- sapply(1:G_dim, function(g){
      x <- true_log[g, ]; y <- inf_log[g, ]
      if(sd(x, na.rm=T)==0 || sd(y, na.rm=T)==0) return(NA)
      rho <- cor(x, y, use="complete.obs")
      s2x <- var(x, na.rm=T); s2y <- var(y, na.rm=T)
      mx <- mean(x, na.rm=T); my <- mean(y, na.rm=T)
      return((2 * rho * sqrt(s2x) * sqrt(s2y)) / (s2x + s2y + (mx - my)^2))
    })
    results_list[["ExpSpe_Summary"]] <- mean(results_list[["ExpSpe_Vector"]], na.rm = TRUE)
  }
  
  # --- 4. ExpMAE & ExpRMSE ---
  if (any(c("ExpMAE", "ExpRMSE") %in% metrics)){
    sum_E <- apply(E, c(1,3), sum) + 1e-9
    sum_G <- apply(G, c(1,3), sum) + 1e-9
    norm_E <- sweep(E, c(1,3), sum_E, "/") * 1e6
    norm_G <- sweep(G, c(1,3), sum_G, "/") * 1e6
    diff <- norm_E - norm_G
    
    if("ExpMAE" %in% metrics) results_list[["ExpMAE_Summary"]]  <- colMeans(apply(abs(diff), c(2,3), mean, na.rm=T))
    if("ExpRMSE" %in% metrics) results_list[["ExpRMSE_Summary"]] <- colMeans(sqrt(apply(diff^2, c(2,3), mean, na.rm=T)))
  }
  
  return(results_list)
}