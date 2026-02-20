
#' Pre-process gene expression data
#'
#' A function to pre-process the data to get it ready for the BayesClint algorithm.
#' @param cnts list of gene count matrices across tissue sections, of dimensions number of genes by number of cells. Make sure both the rows and columns are labelled with gene names and cell names.
#' @param info list of spatial coordinates across tissue sections, of dimensions number of cells by number of coordinates.
#' @param K number of spatial domains (for Potts model calculations)
#' @param cutoff_sample only retain quality cells with greater than cutoff_sample counts across genes
#' @param cutoff_feature only retain quality genes with at least cutoff_feature percent of cells with nonzero counts
#' @param doBatchCorrect Logical. Whether to perform batch effect adjustment with Seurat v3 (Stuart et al., 2019) to align expression data from different tissue sections.
#' @export
#' @examples
#' BayesClint_preprocess()
#' TODO: complete the example above, once you make and document the starmap mpfc dataset and/or the basic simulation scenario dataset
BayesClint_preprocess <- function(cnts, info, K, cutoff_sample = 100, cutoff_feature = 0.1,
                                  doBatchCorrect = F, nc_quant_list = NULL) {

  Np <- length(cnts)

  section_idx <- factor(rep(1:Np, times = sapply(cnts, ncol)))

  # Quality Control ----
  pre_tdataList <- cnts # recall that tdataList generally refers to genes by cells, whereas dataList refers to cells by genes

  # subsetting out the low-quality cells.
  cells2remove <- c() # append elements into this vector, if not in it already
  for (m in 1:Np) {
    section <- pre_tdataList[[m]]
    c2remove <- which(apply(section, 2, sum) < cutoff_sample)
    cells2remove <- base::append(cells2remove, c2remove)
  }
  # eliminate duplicates
  cells2remove <- unique(cells2remove)

  if (length(cells2remove) > 0) {
    pre_tdataList <- lapply(pre_tdataList, function(section) {
      section <- section[, -cells2remove]
    })
    info <- lapply(info, function(section) {
      section <- section[-cells2remove, ]
    })
  }

  # subsetting out the low-quality genes
  genes2remove <- c() # append elements into this vector, if not in it already
  for (m in 1:Np) {
    section <- pre_tdataList[[m]]
    g2remove <- which(apply(section, 1, function(x) mean(x == 0)) > (1 - cutoff_feature))
    genes2remove <- base::append(genes2remove, g2remove)
  }
  # eliminate duplicates
  genes2remove <- unique(genes2remove)

  if (length(genes2remove) > 0) {
    pre_tdataList <- lapply(pre_tdataList, function(section) {
      section <- section[-genes2remove, ]
    })
  }

  cnts <- pre_tdataList
  # ----



  # Batch effect removal, normalization, and scaling ----
  seuList <- list()
  for (i in 1:Np) {
    suppressWarnings(seuList[[i]] <- SeuratObject::CreateSeuratObject(counts = cnts[[i]], meta.data = info[[i]], project = paste0("section", i)))
  }

  for (i in 1:length(seuList)) {
    seuList[[i]] <- Seurat::NormalizeData(seuList[[i]], verbose = FALSE)
    seuList[[i]] <- Seurat::FindVariableFeatures(seuList[[i]], selection.method = "vst",
                                                 nfeatures = 2000, verbose = FALSE)
  }

  if (doBatchCorrect) {
    # find anchors
    anchors <- Seurat::FindIntegrationAnchors(object.list = seuList)
    # integrate data
    suppressWarnings(merged_obj <- Seurat::IntegrateData(anchorset = anchors))
    # switch to integrated assay. The variable features of this assay are automatically
    # set during IntegrateData
    Seurat::DefaultAssay(merged_obj) <- "integrated"
  } else {
    merged_obj <- merge(x = seuList[[1]], y = list(seuList[[2]], seuList[[3]]))
  }

  # Run the standard workflow for visualization and clustering
  merged_obj <- Seurat::ScaleData(merged_obj, split.by = section_idx, verbose = FALSE)
  merged_obj <- Seurat::RunPCA(merged_obj, npcs = 30, verbose = FALSE)
  # ----



  # Spatial pre-processing ----
  xy <- lapply(info, function(info.i){
    as.matrix(info.i[, 1:2])
  }) # a list of spatial coordinates matrices

  W_list <- list()
  for (m in 1:Np) {

    knn_obj <- spdep::knearneigh(xy[[m]], k = 4) # at least 4 neighbors
    nb_obj <- spdep::knn2nb(knn_obj, sym = T)
    W_list[[m]] <- (spdep::nb2mat(nb_obj) != 0)*1

  }

  Wtriplet_list <- vector(mode = "list", length = Np)
  Wbegfin_list <- vector(mode = "list", length = Np)
  edges_list <- vector(mode = "list", length = Np)

  for (m in 1:Np) {
    W <- W_list[[m]]
    W.quants <- common.Wcheckformat(W)
    Wtriplet_list[[m]] <- W.quants$W.triplet
    Wbegfin_list[[m]] <- W.quants$W.begfin

    edge_mat <- W.quants$W.triplet[, 1:2]
    edge_mat_sorted <- t(apply(edge_mat, 1, sort))
    edges_list[[m]] <- edge_mat_sorted[!duplicated(edge_mat_sorted), ]

  }



  # skip below chunk if you've loaded in nc_quant_list already
  if (is.null(nc_quant_list)) {

    nc_quant_list <- vector(mode = "list", length = Np)

    for (i in 1:Np) {
      nc_quant_list[[i]] <- getNC_quantities(subbetas = seq(0, 4, by = 0.01), nvertex = nrow(xy[[i]]), ncolor = K, edges = edges_list[[i]], n = 1200, burn = 200)
    }

  }

  subbetas = seq(0, 4, by = 0.01)

  # Calculate pth-order polynomial interpolations
  p_order <- 10
  M_vec <- 1:Np

  interp_coef_mat <- matrix(NA, nrow = p_order + 1, ncol = Np)
  colnames(interp_coef_mat) <- paste0("M", M_vec)
  row.names(interp_coef_mat) <- paste0("coeff", 0:p_order)

  BETA <- matrix(subbetas, length(subbetas), p_order, byrow = F)
  pow <- matrix(1:p_order, length(subbetas), p_order, byrow = T)
  X <- BETA^pow

  for (m in 1:Np) {

    coeff <- lm(nc_quant_list[[m]][, 1] ~ X)$coef
    interp_coef_mat[, m] <- coeff

  }



  obj_list <- Seurat::SplitObject(merged_obj, split.by = "orig.ident")

  dataList <- lapply(obj_list, function(x) t(as.matrix(x[["integrated"]]$data)))



  return(list(Wtriplet_list = Wtriplet_list, Wbegfin_list = Wbegfin_list, interp_coef_mat = interp_coef_mat, dataList = dataList, nc_quant_list = nc_quant_list))

}



#' Pre-process gene expression data
#'
#' A function to pre-process the data to get it ready for the BayesClint algorithm.
#' @param preprocess_results output from BASS_preprocess, containing
#' @export
#' @examples
#' BayesClint_run()
BayesClint_run <- function(preprocess_results) {



}



# BayesClint_postprocess



#' Pre-process gene expression data
#'
#' A function to pre-process the data to get it ready for the BayesClint algorithm.
#' @param fdr false discovery rate
#' @@export
#' @examples
#' hello()
# BayesClint_performance
