avg_per_label <- function(sce, labels, assay = "logcounts"){
  stopifnot(assay %in% names(assays(sce)))
  cell_labels <- unique(as.vector(labels))
  avg_expr <- sapply(cell_labels, function(cl){
    # Subset to the cell type
    sub_sce <- sce[, labels == cl]
    # Calculate the mean of selected assay across cells for each gene
    rowMeans(as.matrix(assay(sub_sce, assay)))
  })
  return(t(avg_expr))
}



build_Omega <- function(avg_expr, ligand_target_matrix, lr_network, mask_threshold = 0.1, masking_type = "hard"){
  cell_labels <- rownames(avg_expr)
  match.arg(masking_type, c("hard", "soft"))
  
  # Identify common ligands and target genes
  common_ligands <- intersect(colnames(ligand_target_matrix), colnames(avg_expr))
  common_genes <- intersect(rownames(ligand_target_matrix), colnames(avg_expr))
  
  # L = [Sender x Ligand]
  L_matrix <- as.matrix(avg_expr[, common_ligands])
  # W = [Ligand x TargetGene]
  W_matrix <- t(as.matrix(ligand_target_matrix[common_genes, common_ligands]))
  
  
  # Get receptors for the ligands we are using
  relevant_lr = lr_network %>% 
    filter(from %in% common_ligands) %>%
    filter(to %in% colnames(avg_expr)) # Ensure receptor exists in your data
  
  
  # Define a function to check if a cell type is 'responsive' to a ligand
  # We consider a cell type responsive if it expresses ANY receptor for that ligand
  # above a certain threshold (e.g., mean logcounts > 0.1)
  R_matrix <- matrix(0, 
                     nrow = length(cell_labels), 
                     ncol = length(common_ligands),
                     dimnames = list(cell_labels, common_ligands))
  
  for(s in cell_labels) {
    for(lig in common_ligands) {
      # Get all receptors for this ligand
      receptors = relevant_lr %>% filter(from == lig) %>% pull(to)
      
      # Check expression of these receptors in this cell type
      if(length(receptors) > 0) {
        rec_expr = avg_expr[s, receptors]
        if(any(rec_expr > mask_threshold)) {
          if(masking_type == "hard"){
            R_matrix[s, lig] <- 1 # Hard masking: 1 if responsive, 0 if not
          } else if(masking_type == "soft"){
            R_matrix[s, lig] <- sum(rec_expr) # Soft masking: use receptor level
          }
        }
      }
    }
  }
  
  # Initialize the 3D Tensor
  num_senders <- length(cell_labels)
  num_receivers <- length(cell_labels)
  num_genes <- length(common_genes)
  
  Omega <- array(0, 
                 dim = c(num_senders, num_receivers, num_genes),
                 dimnames = list(Sender = cell_labels, 
                                 Receiver = cell_labels, 
                                 Gene = common_genes))
  
  # Fill the tensor
  for(s in cell_labels) { # For each Receiver
    for(k in cell_labels) { # For each Sender
      if(k == s) next 
      
      # Logic: Ligands from k * Mask for s * NicheNet Weights
      # we use element-wise multiplication for the Ligand-Mask part
      masked_ligand_potential = L_matrix[k, ] * R_matrix[s, ]
      
      # Now multiply by the NicheNet Target weights
      # [1 x Ligands] %*% [Ligands x Genes]
      Omega[k, s, ] <- masked_ligand_potential %*% W_matrix
    }
  }
  
  # Normalize 
  valid_indices <- which(Omega != 0)
  top_val <- quantile(abs(Omega), 0.95)
  
  Omega[valid_indices] <- (Omega[valid_indices]) / top_val
  
  return(Omega)
}