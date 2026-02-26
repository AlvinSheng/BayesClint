
# function to get integration quantities needed for beta

getNC_quantities <- function(subbetas = seq(0, 4, by = 0.01), nvertex, ncolor, edges,
                             n = 1200, burn = 200) {

  EUs <- rep(0, length(subbetas))
  SUs <- rep(0, length(subbetas))
  for (i in 1:length(subbetas)) {
    colors <- PottsUtils::SW(n, nvertex, ncolor, edges, beta = subbetas[i])
    EUs[i] <- mean(apply(colors[, -(1:burn)], 2, function(x) sum(x[edges[, 1]] == x[edges[, 2]])))
    SUs[i] <- sd(apply(colors[, -(1:burn)], 2, function(x) sum(x[edges[, 1]] == x[edges[, 2]]))) / sqrt(n - burn)
  }

  return(cbind(EUs, SUs))

}



# adapted from duncanplee CARBayes GitHub
common.Wcheckformat <- function(W)
{
  #### Check W is a matrix of the correct dimension
  if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
  n <- nrow(W)
  if(ncol(W)!= n) stop("W is not a square matrix.", call.=FALSE)


  #### Check validity of inputed W matrix
  if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
  if(min(W)<0) stop("W has negative elements.", call.=FALSE)
  if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
  if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)


  #### Create the triplet form
  ids <- which(W > 0, arr.ind = T)
  W.triplet <- cbind(ids, W[ids])
  W.triplet <- W.triplet[ ,c(2,1,3)]

  #W.triplet <- c(NA, NA, NA)
  #for(i in 1:n)
  #{
  #    for(j in 1:n)
  #    {
  #        if(W[i,j]>0)
  #        {
  #            W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))
  #        }else{}
  #    }
  #}
  #W.triplet <- W.triplet[-1, ]
  n.triplet <- nrow(W.triplet)
  W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
  n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)


  #### Create the start and finish points for W updating
  W.begfin <- cbind(c(1, cumsum(n.neighbours[-n])+1), cumsum(n.neighbours))
  #W.begfin <- array(NA, c(n, 2))
  #temp <- 1
  #for(i in 1:n)
  #{
  #    W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
  #    temp <- temp + n.neighbours[i]
  #}


  #### Return the critical quantities
  results <- list(W=W, W.triplet=W.triplet, n.triplet=n.triplet, W.triplet.sum=W.triplet.sum, n.neighbours=n.neighbours, W.begfin=W.begfin, n=n)
  return(results)
}



# For Model Diagnostics

# Diagnosing estimation performance by examining the MCMC draws
# ggplot boxplots and/or error bars
#
# param mcmc_df matrix of dimension number of parameters x number of MCMC iterations
# To get informative x-labels, do something like the following: mcmc_df <- get the first slice of ups2mcmc, then: row.names(mcmc_df) <- paste0("ups2_", 1:7)
# param truth the true values of the estimated parameters, if available.
mcmc_boxplot_cred_ints <- function(mcmc_df, truth = NULL, boxplot = T, errorbar = T, main = "MCMC boxplots and credible intervals") {

  # Transform the matrix data.
  # t(mcmc_df) produces a matrix where each column is a Parameter.
  df <- melt(t(mcmc_df))
  names(df) <- c("Observation", "Parameter", "Value")
  df$Parameter <- as.factor(df$Parameter)  # Ensure Parameter is a factor

  if (!is.null(truth)) {
    # Create a data frame for the red points (one per Parameter)
    df_points <- data.frame(Parameter = levels(df$Parameter),
                            Truth = truth)
  } else {
    df_points <- data.frame(Parameter = levels(df$Parameter))
  }

  # Compute the lower and upper quantiles for each Parameter
  df_error <- df %>%
    group_by(Parameter) %>%
    summarise(qlower = quantile(Value, probs = 0.005),
              qupper = quantile(Value, probs = 0.995))

  # Plot using ggplot2
  p <- ggplot() +
    labs(title = main,
         x = "Parameter", y = "Value") + coord_cartesian(ylim = c(min(df_error$qlower), max(df_error$qupper)))

  if (boxplot) {
    p <- p + geom_boxplot(data = df, aes(x = Parameter, y = Value), outlier.shape = NA, size = 1)
  }

  if (errorbar) {
    p <- p + geom_errorbar(data = df_error,
                           aes(x = Parameter, ymin = qlower, ymax = qupper),
                           color = "magenta", width = 0.2)
  }

  if (!is.null(truth)) {
    p <- p + geom_point(data = df_points, aes(x = 1:nrow(mcmc_df), y = Truth),
                        color = "red", size = 3) # adjust size if needed
  }

  return(p)

}



# Plotting a ggplot that displays the bayesian p-value calculated from the posterior predictive distribution, for all data points.
#
# Arguments:
# data_rbind: a matrix containing all of the observed values, where the data are rbinded across the tissue sections
# Xppd_array: an array containing the PPD of each datapoint, in which there are 1000s of MCMC draws for each datapoint
data_points_bayesian_pval <- function(data_rbind, Xppd_array) {

  # flatten both data_rbind and Xppd_array along each genetic feature
  data_rbind_flatten <- c(data_rbind)
  Xppd_array_flatten <- do.call("rbind", apply(Xppd_array, 2, function(mat) mat, simplify = F))

  # calculate the Bayesian p-values
  bayesian_pvals <- rep(NA, length(data_rbind_flatten))
  for (i in 1:length(bayesian_pvals)) {
    bayesian_pvals[i] <- mean(Xppd_array_flatten[i, ] > data_rbind_flatten[i])
  }

  bayesian_pvals_df <- data.frame(var_idx = 1:length(bayesian_pvals), bayesian_pvals = bayesian_pvals)

  p <- ggplot(bayesian_pvals_df, aes(x = var_idx, y = bayesian_pvals)) +
    geom_point()

  return(list(bayesian_pvals = bayesian_pvals, pval_plot = p))

}



# For the calculation of the threshold corresponding to a given BFDR, I step the threshold down from
# 1 to 0 where the step size is the smallest difference among the ordered PPIs, to make the
# feature selections as powerful as possible while still being below the FDR (which will likely be 0.05)

# Calculate Bayesian false discovery rate (BFDR)
bfdr_precise <- function(PPI, alpha){

  # find the smallest step size to use
  PPI_ordered <- sort(PPI)
  PPI_diffs <- PPI_ordered - c(NA, PPI_ordered[-length(PPI_ordered)])
  seqby <- min(PPI_diffs[PPI_diffs > 0], na.rm = T)

  for (c in seq(1,0,by=-seqby)) {

    BFDR <- sum((1 - PPI)*((1 - PPI) < c))/sum((1 - PPI) < c)

    if (BFDR < alpha){
      return(c)
      stop;
    }
  }

}



# bfdr <- function(PPI, alpha){
#   for (c in seq(1,0.1,by=-0.1)) {
#
#     BFDR <- sum((1 - PPI)*((1 - PPI) < c))/sum((1 - PPI) < c)
#
#     if (BFDR < alpha){
#       return(c)
#       stop;
#     }
#   }
# }



cal_bfdr <- Vectorize(function(PPI, threshold){
  return(sum((1 - PPI)*((1 - PPI) < threshold))/sum((1 - PPI) < threshold))
}, vectorize.args = "threshold")



# function that globally permutes the labels to best match a ground truth
label_truth_permute <- function(est_labels, true_labels) {

  # table() handles numeric vs character automatically
  tab <- table(est_labels, true_labels)

  # Find the optimal column index for each row index
  # solve_map[1] = index of the column that best matches row 1
  solve_map <- solve_LSAP(tab, maximum = TRUE)

  # Get the actual character names from the column headers
  # colnames(tab) contains ("A", "B", "C"...)
  best_char_labels <- colnames(tab)[solve_map]

  # Map the original estimated labels to these character labels
  # We use match() to ensure we handle the mapping correctly
  final_labels <- best_char_labels[match(est_labels, rownames(tab))]

  return(final_labels)

}



# function that globally permutes the cell-type and spatial-domain labels to best match the truth, so the pi estimates would be meaningful
# even though only BayesClint and BASS does joint multi-scale clustering among all of the comparison methods, you technically can still use this function for all of the clustering methods. You use these methods for both scales of clustering anyways, although not jointly.
zeta_kappa_truth_permute <- function(concat_zeta, true_zeta, concat_kappa, true_kappa) {

  # for zeta
  zeta_perm <- label_truth_permute(concat_zeta, true_zeta)
  zeta_perm <- factor(zeta_perm, levels = levels(true_zeta))

  # for kappa
  kappa_perm <- label_truth_permute(concat_kappa, true_kappa)
  kappa_perm <- factor(kappa_perm, levels = levels(true_kappa))

  # estimate theta on the basis of the permuted zeta and kappa
  est_theta <- prop.table(table(zeta_perm, kappa_perm), margin = 2)
  est_theta[is.nan(est_theta)] <- 0 # corresponding to spatial domains with no elements

  return(list(zeta_perm = zeta_perm, kappa_perm = kappa_perm, est_theta = est_theta))

}


