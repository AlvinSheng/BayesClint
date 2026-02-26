
#' Pre-process gene expression data by doing quality control, batch effect removal (if desired), log-normalization (using the counts-per-10K approach as done in Seurat by default), and scaling
#'
#' A function to pre-process the data to get it ready for the BayesClint algorithm.
#' @param cnts list of gene count matrices across tissue sections, of dimensions number of genes by number of cells. Make sure both the rows and columns are labelled with gene names and cell names.
#' @param info list of spatial coordinates across tissue sections, of dimensions number of cells by number of coordinates.
#' @param cutoff_sample only retain quality cells with greater than cutoff_sample counts across genes
#' @param cutoff_feature only retain quality genes with at least cutoff_feature percent of cells with nonzero counts
#' @param doBatchCorrect Logical. Whether to perform batch effect adjustment with Seurat v3 (Stuart et al., 2019) to align expression data from different tissue sections.
#' @returns a list of pre-processed gene expression matrices, one for each tissue section.
#' @export
#' @examples
#' BayesClint_preprocess()
#' TODO: complete the example above, once you make and document the starmap mpfc dataset and/or the basic simulation scenario dataset
BayesClint_preprocess <- function(cnts, info, cutoff_sample = 100, cutoff_feature = 0.1,
                                  doBatchCorrect = F, doScaling = T) {

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



  obj_list <- Seurat::SplitObject(merged_obj, split.by = "orig.ident")

  if (doScaling) {
    # I use the "scale.data" slot, corresponding to the z-scored/variance-stabilized data (done separately by tissue section), rather than "counts" (un-normalized raw counts) or "data" (normalized counts)
    # it is equivalent to using the scale(., center = T, scale = T) on "data"
    dataList <- lapply(obj_list, function(x) t(as.matrix(x[["integrated"]]$scale.data)))
  } else {
    # I use the "data" slot, corresponding to the normalized counts, rather than "counts" (un-normalized raw counts) or "scale.data" (z-scored/variance-stabilized data)
    dataList <- lapply(obj_list, function(x) t(as.matrix(x[["integrated"]]$scale.data)))
  }



  return(list(info = info, dataList = dataList))

}



#' Running the BayesClint MCMC algorithm
#'
#' A function to pre-process the data to get it ready for the BayesClint algorithm.
#' @param dataList pre-processed dataset, like the one outputted by BayesClint_preprocess
#' @param info list of spatial coordinates across tissue sections, of dimensions number of cells by number of coordinates. This should be the info returned by the BayesClint_preprocess's quality control.
#' @param C number of cell types
#' @param K number of spatial domains
#' @param r number of components in the factor model
#' @param nbrsample number of MCMC iterations retained
#' @param burnin number of initial MCMC iterations discarded as burn-in
#' @param probvarsel prior probability for an entry in the factor loadings matrix to be nonzero
#' @param nc_quant_list the quantities used to compute the normalizing constant of the Potts model within the MCMC.
#' The nc_quant_list object will already be calculated within the BayesClint_run function;
#' however, since this step takes a long time, one can import nc_quant_list (precalculated by an earlier run of BayesClint_run)
#' through this argument to skip this computationally intensive step. However, if the nc_quant_list argument is set to NULL,
#' BayesClint_run will compute it and return it as one of the outputs.
#' @param chainNbr MCMC chain number. The default is 1 for one MCMC number 1. If you want to run N multiple MCMC chains, make a loop in R as for(i in 1:N) BayesClint_run(..., chainNbr=i)
#' @returns the MCMC results from the BayesClint algorithm, and optionally the intermediary data file nc_quant_list.
#' @export
#' @examples
#' BayesClint_run()
BayesClint_run <- function(dataList, info, C, K, r, nbrsample, burnin, probvarsel, nc_quant_list = NULL, chainNbr = 1) {

  Np <- length(dataList)

  P <- unique(sapply(dataList, ncol))
  if (length(P) != 1) {
    stop("Each element of dataList must have a number of columns equal to the number of features.")
  }

  return_nc_quant_list <- FALSE
  if (is.null(nc_quant_list)) {
    return_nc_quant_list <- TRUE
  }



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
  # ----



  # Running the MCMC ----
  IndicVar <- 0;

  n <- NULL;
  for (i in 1:Np){
    n[i]=nrow(dataList[[i]])
  }
  datasets=do.call("rbind", dataList)

  if (nbrsample<=20){
    stop("Please specify a larger number of MCMC iterations")
  }

  if (is.null(r)){ # could rbind the datasets and do one PCA truncation instead
    mysvd=lapply(1:Np, function(i)  svd(dataList[[i]]))
    mysumsvd=lapply(1:Np, function(i) cumsum(mysvd[[i]]$d)/max(cumsum(mysvd[[i]]$d))*100)
    KMax=max(unlist(lapply(1:Np, function(i) min(which(mysumsvd[[i]]>=80, arr.ind = TRUE))))) #chooses maximum from Np cumulative proportions
    r=min(KMax+1,10)
  }



  ptm <- proc.time()[3]

  EstIMat1 <- vector(mode = "list", length = Np)
  EstIMat2 <- vector(mode = "list", length = Np)
  for (i in 1:Np) {
    EstIMat1[[i]] <- sample(1:C, size = n[i], replace = T)
    EstIMat2[[i]] <- sample(1:K, size = n[i], replace = T)
  }

  results <- mcmcfn(n = as.integer(n), P = as.integer(P), r = as.integer(r), Np = as.integer(Np), datasets = datasets, IndVar = as.integer(IndicVar), nbrsample = as.integer(nbrsample),
                    burninsample = as.integer(burnin), CompoSelMean = as.double(rep(0,r)), VarSelMean = as.double(rep(0,r*P)), VarSelMeanGlobal = as.double(rep(0,P)),
                    priorcompsel = c(1, 1), probvarsel = as.double(probvarsel),
                    C = C, EstMat1 = matrix(rnorm(r*C), nrow = r, ncol = C),
                    EstMat2 = diag(r),
                    EstIMat1 = EstIMat1,
                    K = K, EstIMat2 = EstIMat2,
                    EstMat3 = matrix(runif(C*K), nrow = C, ncol = K),
                    Wtriplet_list = Wtriplet_list, Wbegfin_list = Wbegfin_list,
                    beta = rep(1, Np), interp_coef_mat = interp_coef_mat, chainNbr = chainNbr, streamline = T);

  timing <- proc.time()[3] - ptm

  results$timing <- timing
  # ----



  if (return_nc_quant_list) {
    return(list(results = results, nc_quant_list = nc_quant_list))
  } else {
    return(list(results = results))
  }

}



#' Post-process the posterior sampling results
#'
#' Post-process the posterior sampling results to address the label switching issue associated with
#' the sampling of cell type labels and spatial domain labels based on the ECR-1 algorithm, estimate
#' the cell type labels and spatial domain labels as the mode of all their posterior samples, and
#' estimate the cell type composition in each spatial domain based on the final estimates of cell type and
#' spatial domain labels.
#' @param zeta2mcmc MCMC samples of the cell type labels
#' @param kappa2mcmc MCMC samples of the spatial domain labels
#' @param mu2mcmc MCMC samples of the mu parameter
#' @param C number of cell types
#' @param K number of spatial domains
#' @returns the cell type and spatial domain cluster labels, that have been corrected for label switching.
#' Also returns the permutations returned by the label.switching::label.switching() function, which can be used to permute other parameters
#' depending on the cluster labels, like mu, which is returned in this function as well (if mu2mcmc is provided).
#' Also returns the cell type compositions calculated on the basis of the corrected cluster labels.
#' @export
#' @examples
#' BayesClint_postprocess()
BayesClint_postprocess <- function(zeta2mcmc, kappa2mcmc,
                                   mu2mcmc = NULL,
                                   C, K) {

  n <- sapply(zeta2mcmc, nrow)
  Np <- length(zeta2mcmc)

  zeta_z <- t(do.call("rbind", zeta2mcmc))
  zeta_K <- C
  zeta_ls <- label.switching::label.switching(method = "ECR-ITERATIVE-1", z = zeta_z, K = zeta_K)

  kappa_z <- t(do.call("rbind", kappa2mcmc))
  kappa_K <- K
  kappa_ls <- label.switching::label.switching(method = "ECR-ITERATIVE-1", z = kappa_z, K = kappa_K)

  # estimate the cell type compositions theta on the basis of the permuted zeta and kappa
  theta_est <- prop.table(table(c(zeta_ls$clusters),
                                c(kappa_ls$clusters)),
                          margin = 2)
  theta_est[is.nan(theta_est)] <- 0 # corresponding to spatial domains with no elements

  # permute mu, if mu2mcmc object provided
  if (!is.null(mu2mcmc)) {
    permute_res <- label.switching::permute.mcmc(
      aperm(mu2mcmc, c(3, 2, 1)),
      zeta_ls$permutations$`ECR-ITERATIVE-1`)

    mu2mcmc <- aperm(permute_res$output, c(3, 2, 1)) # reversing the dimension reordering
  }

  # split the corrected labels back into their tissues
  zeta_est <- split(zeta_ls$clusters, rep(1:Np, times = n))
  kappa_est <- split(kappa_ls$clusters, rep(1:Np, times = n))

  return(list(zeta_est = zeta_est, kappa_est = kappa_est,
              mu2mcmc = mu2mcmc, theta_est = theta_est))

}


