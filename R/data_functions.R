# Get a SingleCellExperiment object from a set of assay.mtx, genes.tsv and cells.tsv file
mtx_to_sce <- function(assay_path, genes_path, cells_path, assay_type = "counts"){
  metric <- readMM(assay_path)
  genes <-  fread(genes_path, header = FALSE)$V1
  cells <- fread(cells_path, header = TRUE)
  barcodes <- cells[[1]]
  
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = as.data.frame(cells)
  )
  
  return(sce)
}