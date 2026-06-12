library(dplyr)
library(Seurat)
library(SeuratObject)

#' get_expressed_genes_mod: Get expressed genes mod
#'
#' @return Gets expressed genes mod
#' @export
get_expressed_genes_mod <- function (ident, seurat_obj, pct = 0.1, assay_oi = NULL)
{
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  if (!"RNA" %in% names(seurat_obj@assays)) {
    if ("Spatial" %in% names(seurat_obj@assays)) {
      if (class(seurat_obj@assays$Spatial@data) != "matrix" &
          class(seurat_obj@assays$Spatial@data) != "dgCMatrix") {
        warning("Spatial Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$Spatial@data' for default or 'seurat_obj@assays$SCT@data' for when the single-cell transform pipeline was applied")
      }
      if (sum(dim(seurat_obj@assays$Spatial@data)) == 0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$Spatial@data'")
      }
    }
  }
  else {
    # if (class(seurat_obj@assays$RNA@data) != "matrix" & class(seurat_obj@assays$RNA@data) !=
    #     "dgCMatrix") {
    #   warning("Seurat object should contain a matrix of normalized expression data. Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data or seurat_obj@assays$SCT@data for when the single-cell transform pipeline was applied")
    # }
    if ("integrated" %in% names(seurat_obj@assays)) {
      if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$integrated@data)) ==
          0)
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$integrated@data' for integrated data")
    }
    else if ("SCT" %in% names(seurat_obj@assays)) {
      if (sum(dim(seurat_obj@assays$RNA@data)) == 0 & sum(dim(seurat_obj@assays$SCT@data)) ==
          0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data' for default or 'seurat_obj@assays$SCT@data' for data corrected via SCT")
      }
    }
    else {
      if (sum(dim(seurat_obj@assays$RNA@data)) == 0) {
        stop("Seurat object should contain normalized expression data (numeric matrix). Check 'seurat_obj@assays$RNA@data'")
      }
    }
  }
  if (sum(ident %in% unique(Idents(seurat_obj))) != length(ident)) {
    stop("One or more provided cell clusters is not part of the 'Idents' of your Seurat object")
  }
  if (!is.null(assay_oi)) {
    if (!assay_oi %in% Seurat::Assays(seurat_obj)) {
      stop("assay_oi should be an assay of your Seurat object")
    }
  }
  cells_oi = Idents(seurat_obj) %>% .[Idents(seurat_obj) %in%
                                        ident] %>% names()
  if (!is.null(assay_oi)) {
    cells_oi_in_matrix = intersect(colnames(seurat_obj[[assay_oi]]@data),
                                   cells_oi)
    exprs_mat = seurat_obj[[assay_oi]]@data %>% .[, cells_oi_in_matrix]
  }
  else {
    if ("integrated" %in% names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat integration workflow. The expressed genes are now defined based on the integrated slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$integrated@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$integrated@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$integrated@data %>%
        .[, cells_oi_in_matrix]
    }
    else if ("SCT" %in% names(seurat_obj@assays) & !"Spatial" %in%
             names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat single-cell transform workflow. The expressed genes are defined based on the SCT slot. You can change this via the assay_oi parameter of the get_expressed_genes() functions. Recommended assays: RNA or SCT")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$SCT@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$SCT@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
    }
    else if ("Spatial" %in% names(seurat_obj@assays) & !"SCT" %in%
             names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat spatial object. The expressed genes are defined based on the Spatial slot. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! ;-) )")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$Spatial@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$Spatial@data %>% .[,
                                                       cells_oi_in_matrix]
    }
    else if ("Spatial" %in% names(seurat_obj@assays) & "SCT" %in%
             names(seurat_obj@assays)) {
      warning("Seurat object is result from the Seurat spatial object, followed by the SCT workflow. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on spot-based data (because you want to look at cell-cell interactions, and not at spot-spot interactions! The expressed genes are defined based on the SCT slot, but this can be changed via the assay_oi parameter.")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$SCT@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$Spatial@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$SCT@data %>% .[, cells_oi_in_matrix]
    }
    else {
      if (sum(cells_oi %in% colnames(seurat_obj@assays$RNA@data)) ==
          0)
        stop("None of the cells are in colnames of 'seurat_obj@assays$RNA@data'. The expression matrix should contain cells in columns and genes in rows.")
      cells_oi_in_matrix = intersect(colnames(seurat_obj@assays$RNA@data),
                                     cells_oi)
      if (length(cells_oi_in_matrix) != length(cells_oi))
        stop("Not all cells of interest are in your expression matrix (seurat_obj@assays$RNA@data). Please check that the expression matrix contains cells in columns and genes in rows.")
      exprs_mat = seurat_obj@assays$RNA@data %>% .[, cells_oi_in_matrix]
    }
  }
  n_cells_oi_in_matrix = length(cells_oi_in_matrix)
  if (n_cells_oi_in_matrix < 5000) {
    genes = exprs_mat %>% apply(1, function(x) {
      sum(x > 0)/n_cells_oi_in_matrix
    }) %>% .[. >= pct] %>% names()
  }
  else {
    splits = split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
    genes = splits %>% lapply(function(genes_indices, exprs,
                                       pct, n_cells_oi_in_matrix) {
      begin_i = genes_indices[1]
      end_i = genes_indices[length(genes_indices)]
      exprs = exprs[begin_i:end_i, , drop = FALSE]
      genes = exprs %>% apply(1, function(x) {
        sum(x > 0)/n_cells_oi_in_matrix
      }) %>% .[. >= pct] %>% names()
    }, exprs_mat, pct, n_cells_oi_in_matrix) %>% unlist() %>%
      unname()
  }
  return(genes)
}

#' remove_reciprocalRows: Remove reciprocal rows
#'
#' @return Removes reciprocal rows
#' @export
remove_reciprocalRows <- function(data){
  data$V1V2 <- paste0(data$Gene1,"-",data$Gene2)
  data$V2V1 <- paste0(data$Gene2,"-",data$Gene1)
  for(i in 1:nrow(data)){
    x <- data[i,]
    if(x$Gene1 %in% x$Gene2){
      next
    }
    if(x$V1V2 %in% data$V2V1){
      data <- data[-i,]
    }
    else{
      next
    }
  }

  return(data)
}
