
# Function to plot cluster labels with real data spatial point pattern superimposed
# Takes in the S_dim parameter instead of assuming that it's something like x_dim = 0:13, y_dim = 0:10.
# Would automatically adjust for non-unit square tiles
plot_clust_spp = function(spp, cluster_lab, S_dim, main = "", tile_x_incr = 0, tile_y_incr = 0, truncate = F, levels = NULL, types, type_labels = NULL){

  in_mat <- cluster_lab

  # getting the width and height of the tiles
  tile_width <- S_dim$x_dim[2] - S_dim$x_dim[1]
  tile_height <- S_dim$y_dim[2] - S_dim$y_dim[1]

  # getting the centers
  X <- S_dim$x_dim[-1] - tile_width/2 + tile_x_incr
  Y <- S_dim$y_dim[-1] - tile_height/2 + tile_y_incr

  X_range <- c(S_dim$x_dim[1], S_dim$x_dim[length(S_dim$x_dim)])
  Y_range <- c(S_dim$y_dim[1], S_dim$y_dim[length(S_dim$y_dim)])

  df <- expand.grid(X=X, Y=Y)

  if (is.null(levels)) {
    df$Z <- as.character(t(in_mat))
  } else {
    df$Z <- factor(t(in_mat), levels = levels)
  }

  df2 <- data.frame(x = spp$x, y = spp$y, Type = spp$marks$zeta)
  if (truncate) {
    df2 <- df2[df2$x > X_range[1] & df2$x < tail(X_range, 1) &
                 df2$y > Y_range[1] & df2$y < tail(Y_range, 1), ]
  }

  if (!is.null(type_labels)) {
    df2$Type <- factor(df2$Type, levels = types)

    num_types <- length(type_labels)
    col_vec <- gg_color_hue(num_types - 1)

    col_vec <- c(col_vec, "gray")
  }

  ggplot(df) +
    geom_tile(aes(X, Y, fill= Z), show.legend = T) +
    guides(fill=guide_legend(override.aes=list(shape=NA))) +
    scale_fill_grey(name = "Spatial Domain", drop = FALSE) +
    ggtitle(main) +
    geom_point(data = df2, mapping = aes(x, y, col = Type), size = 2, show.legend = T) +
    labs(color = "Cell Type") +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    # scale_color_manual(labels = type_labels, values = col_vec, drop = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme_classic()

}



# Map spatial domain label to cell types in that domain
# 1 -> (1, 2, 3), 2 -> (2, 3, 4)
# 3 -> (3, 4, 1), 4 -> (4, 1, 2)
map_k2c <- function(k)
{
  case_when(
    k == 1 ~ c(1, 2, 3),
    k == 2 ~ c(2, 3, 4),
    k == 3 ~ c(3, 4, 1),
    k == 4 ~ c(4, 1, 2)
  )
}



# Simulation Version 3: Still working with rectangular tissue sections, but generating spatial point patterns for the main simulation study now
#
# Arguments (Inputs)
# kappa_grid_list: labelled grids for the tissues sections, indicating the boundaries of the spatial domains
# cell_intens: number of cells per grid region in kappa_grid (element of kappa_grid_list)
# C: number of cell types
# K: number of spatial domains
# theta: cell-type proportions. Will be used to generate cell types, only for the factor model simulation. Otherwise, the splatter package would generate the cell types.
#
# Output:
# list of SPPs where each point is a cell with two marks, the cell type and spatial domain label
simulation_ver3_spps <- function(kappa_grid_list, cell_intens, C, K, theta = NULL) {

  Np <- length(kappa_grid_list)

  spp_list <- vector(mode = "list", length = Np)

  for (m in 1:Np) {

    kappa_grid <- kappa_grid_list[[m]]

    sub_spp_list <- vector(mode = "list", length = K)

    for (k in 1:K) { # generating sub-SPPs for each spatial domain label

      if (k %in% c(kappa_grid)) {
        w <- owin(mask = (kappa_grid == k))
        k_spp <- runifpoint(cell_intens * sum(kappa_grid == k),
                            win = w)

        if (is.null(theta)) {
          mark_df <- data.frame(zeta = rep(NA, length = npoints(k_spp)),
                                kappa = rep(k, length = npoints(k_spp)))
        } else {
          mark_df <- data.frame(zeta = sample(1:C, size = npoints(k_spp), replace = T, prob = theta[, k]),
                                kappa = rep(k, length = npoints(k_spp)))
        }
        mark_df$zeta <- factor(mark_df$zeta, levels = 1:C)
        mark_df$kappa <- factor(mark_df$kappa, levels = 1:K)
        marks(k_spp) <- mark_df
      } else { # for the case when the spatial domain label doesn't exist within tissue section
        next # leaving sub_spp_list[[k]] at NULL
      }

      sub_spp_list[[k]] <- k_spp

    }

    attr(sub_spp_list, "class") <- "ppplist"

    # superimposing the SPPs in sub_spp_list
    spp_list[[m]] <- superimpose(... = sub_spp_list, W = owin(xrange = c(0, ncol(kappa_grid)), yrange = c(0, nrow(kappa_grid))))

  }

  return(spp_list)

}



# Simulation Version 3: Generating the genetic expressions for each cell in the SPPs
#
# Arguments (Inputs)
# spp_list: list of SPPs where each point is a cell with two marks, the cell type and spatial domain label
# A: (sparse) loadings matrix
# P: number of genetic features
# mu: r x C matrix of cell means
# Sigma: r x r variance-covariance matrix for the factors
# EstUps2: the variances of the independent errors, for each tissue section and each gene
# tgints: Np x P matrix of tissue section-gene-specific intercepts
#
# Output:
# Two lists, where the first list has the genetic expression counts and the second list has the meta-information for each cell (spatial coordinates and cell-type/spatial-domain labels)
simulation_ver3_exprs <- function(spp_list, A, P, mu, Sigma,
                                  EstUps2 = matrix(1, nrow = length(spp_list), ncol = P),
                                  tgints = matrix(0, nrow = length(spp_list), ncol = P)) {

  n <- sapply(spp_list, npoints)

  Np <- length(spp_list)

  r <- nrow(mu)

  zeta <- lapply(spp_list, function(x) x$marks$zeta)

  # Generating factors U
  U <- vector(mode = "list", length = Np)
  for (m in 1:Np) {
    U[[m]] <- matrix(NA, nrow = length(zeta[[m]]), ncol = r)
    for (i in 1:n[m]) {
      U[[m]][i, ] <- mvrnorm(n = 1, mu = mu[, zeta[[m]][i] ], Sigma = Sigma)
    }
  }

  # Error matrix E randomly generated from N(0, diag(\upsilon_{m1}, ..., \upsilon_{mP}))
  E <- vector(mode = "list", length = Np)
  X <- vector(mode = "list", length = Np)
  for (m in 1:Np) {
    E[[m]] <- mvrnorm(n[m], rep(0, P), Sigma = diag(EstUps2[m, ]))
    X[[m]] <- U[[m]] %*% A + E[[m]] + matrix(rep(tgints[m, ], n[m]), nrow = n[m], ncol = P, byrow = T)
  }



  Xcounts <- vector(mode = "list", length = Np)

  # Now, do Poisson counts based on the generated log gene expression
  for (m in 1:Np) {
    Xcounts[[m]] <- matrix(NA, nrow = n[m], ncol = P)
    for (i in 1:n[m]) {
      for (j in 1:P) {
        Xcounts[[m]][i, j] <- rpois(1, lambda = exp(X[[m]][i, j]))
      }
    }
  }


  # convert the results into a form similar to starmap_cnts and starmap_info
  simver_cnts <- vector(mode = "list", length = Np)
  simver_info <- vector(mode = "list", length = Np)
  for (m in 1:Np) {
    simver_cnts[[m]] <- t(Xcounts[[m]])
    rownames(simver_cnts[[m]]) <- paste0("gene_", seq_len(nrow(simver_cnts[[m]])))
    colnames(simver_cnts[[m]]) <- paste0("cell_", seq_len(ncol(simver_cnts[[m]])))

    simver_info[[m]] <- data.frame(x = spp_list[[m]]$x, y = spp_list[[m]]$y,
                                   zeta = spp_list[[m]]$marks$zeta, kappa = spp_list[[m]]$marks$kappa)
    rownames(simver_info[[m]]) <- paste0("cell_", seq_len(ncol(simver_cnts[[m]])))
  }

  return(list(simver_cnts = simver_cnts, simver_info = simver_info, U = U))

}



# Simulation Version 4: Generating the log-normalized genetic expressions for each cell in the SPPs
#
# Arguments (Inputs)
# spp_list: list of SPPs where each point is a cell with two marks, the cell type and spatial domain label
# A: (sparse) loadings matrix
# P: number of genetic features
# mu: r x C matrix of cell means
# Sigma: r x r variance-covariance matrix for the factors
# EstUps2: the variances of the independent errors, for each tissue section and each gene
#
# Output:
# Two lists, where the first list has the genetic expression counts and the second list has the meta-information for each cell (spatial coordinates and cell-type/spatial-domain labels)
simulation_ver4_exprs <- function(spp_list, A, P, mu, Sigma,
                                  EstUps2 = matrix(1, nrow = length(spp_list), ncol = P)) {

  n <- sapply(spp_list, npoints)

  Np <- length(spp_list)

  r <- nrow(mu)

  zeta <- lapply(spp_list, function(x) x$marks$zeta)

  # Generating factors U
  U <- vector(mode = "list", length = Np)
  for (m in 1:Np) {
    U[[m]] <- matrix(NA, nrow = length(zeta[[m]]), ncol = r)
    for (i in 1:n[m]) {
      U[[m]][i, ] <- mvrnorm(n = 1, mu = mu[, zeta[[m]][i] ], Sigma = Sigma)
    }
  }

  # Error matrix E randomly generated from N(0, diag(\upsilon_{m1}, ..., \upsilon_{mP}))
  E <- vector(mode = "list", length = Np)
  X <- vector(mode = "list", length = Np)
  for (m in 1:Np) {
    E[[m]] <- mvrnorm(n[m], rep(0, P), Sigma = diag(EstUps2[m, ]))
    X[[m]] <- U[[m]] %*% A + E[[m]]
  }



  # convert the results into a form similar to starmap_cnts and starmap_info
  simver_expr <- vector(mode = "list", length = Np)
  simver_info <- vector(mode = "list", length = Np)
  for (m in 1:Np) {
    simver_expr[[m]] <- t(X[[m]])
    rownames(simver_expr[[m]]) <- paste0("gene_", seq_len(nrow(simver_expr[[m]])))
    colnames(simver_expr[[m]]) <- paste0("cell_", seq_len(ncol(simver_expr[[m]])))

    simver_info[[m]] <- data.frame(x = spp_list[[m]]$x, y = spp_list[[m]]$y,
                                   zeta = spp_list[[m]]$marks$zeta, kappa = spp_list[[m]]$marks$kappa)
    rownames(simver_info[[m]]) <- paste0("cell_", seq_len(ncol(simver_expr[[m]])))
  }

  return(list(simver_expr = simver_expr, simver_info = simver_info, U = U))

}



# LaTeX formatting functions, initially from main_simulation_results_oct7.Rmd and updated for december results

# function to collect all the results from BayesClint and competing feature selection methods, given a metric
res_collect_metric1 <- function(metric_name, model_names, sel_type = "active") {

  sim_results_mean <- matrix(NA, nrow = 16, ncol = length(model_names))
  sim_results_sd <- matrix(NA, nrow = 16, ncol = length(model_names))

  ctr <- 1 # for indexing the factor combos
  for (Np in c(1, 3)) {
    for (theta_type in c("arb", "reg")) {
      for (P in c(200, 1000)) {
        for (num_deg in c(40, 80)) {
          factor_string <- paste0("Np", Np, "_theta", theta_type, "_P", P, "_numdeg", num_deg)

          model_results <- readRDS(here("modeling_files", paste0("simver", simver, "_", factor_string), "monte_carlo_results.rds"))

          load(here("modeling_files", paste0("simver", simver, "_", factor_string), "comparison_methods_results_batch1.RData"))
          load(here("modeling_files", paste0("simver", simver, "_", factor_string), "comparison_methods_results_batch2.RData")) # just for the SPC results

          if (sel_type == "active") {

            if (str_extract(metric_name, "(?<=_).*") == "acc") {
              model_results_list <- list(model_results, spc_results, hvgs_results,
                                         corpair_results)
            } else if (str_extract(metric_name, "(?<=_).*") == "auc") {
              model_results_list <- list(model_results, hvgs_results, corpair_results)
            }

          } else {

            if (str_extract(metric_name, "(?<=_).*") == "acc") {
              model_results_list <- list(model_results, spc_results, sparkx_results,
                                         deseq_results)
            } else if (str_extract(metric_name, "(?<=_).*") == "auc") {
              model_results_list <- list(model_results, sparkx_results,
                                         deseq_results)
            }

          }

          names(model_results_list) <- model_names

          metric_names <- c("zeta_ARI",
                            "kappa_ARI",
                            paste0("zeta_ARI_single", 1:Np),
                            paste0("kappa_ARI_single", 1:Np),
                            "active_auc",
                            "deg_auc",
                            "active_acc",
                            "deg_acc",
                            "model_timing")

          for (model_str in model_names) {
            # for (idx in 1:num_MCs) {
            # for (metric_name in main_metric_names) {

            sim_results_mean[ctr, model_names == model_str] <- mean(model_results_list[[model_str]][, metric_names == metric_name])
            sim_results_sd[ctr, model_names == model_str] <- sd(model_results_list[[model_str]][, metric_names == metric_name])

            # }
            # }
          }

          ctr <- ctr + 1

        }
      }
    }
  }

  return(list(sim_results_mean = sim_results_mean, sim_results_sd = sim_results_sd))

}



# function to collect all the results from BayesClint and competing clustering methods, given a metric
res_collect_metric2 <- function(metric_name, model_names) {

  sim_results_mean <- matrix(NA, nrow = 16, ncol = length(model_names))
  sim_results_sd <- matrix(NA, nrow = 16, ncol = length(model_names))

  ctr <- 1 # for indexing the factor combos
  for (Np in c(1, 3)) {
    for (theta_type in c("arb", "reg")) {
      for (P in c(200, 1000)) {
        for (num_deg in c(40, 80)) {
          factor_string <- paste0("Np", Np, "_theta", theta_type, "_P", P, "_numdeg", num_deg)

          model_results <- readRDS(here("modeling_files", paste0("simver", simver, "_", factor_string), "monte_carlo_results.rds"))

          load(here("modeling_files", paste0("simver", simver, "_", factor_string), "comparison_methods_results_batch2.RData"))
          load(here("modeling_files", paste0("simver", simver, "_", factor_string), "comparison_methods_results_batch1.RData")) # just for the BASS results

          model_results_list <- list(model_results, spc_kmeans_results,
                                     spc_mclust_results, precast_results,
                                     bass_results)
          names(model_results_list) <- model_names

          metric_names <- c("zeta_ARI",
                            "kappa_ARI",
                            paste0("zeta_ARI_single", 1:Np),
                            paste0("kappa_ARI_single", 1:Np),
                            "active_auc",
                            "deg_auc",
                            "active_acc",
                            "deg_acc",
                            "model_timing")

          for (model_str in model_names) {
            # for (idx in 1:num_MCs) {
            # for (metric_name in main_metric_names) {

            sim_results_mean[ctr, model_names == model_str] <- mean(model_results_list[[model_str]][, metric_names == metric_name])
            sim_results_sd[ctr, model_names == model_str] <- sd(model_results_list[[model_str]][, metric_names == metric_name])

            # }
            # }
          }

          ctr <- ctr + 1

        }
      }
    }
  }

  return(list(sim_results_mean = sim_results_mean, sim_results_sd = sim_results_sd))

}



# function to format the sim_results_mean and sim_results_sd in a nice way
latex_table_format <- function(sim_results_mean, sim_results_sd, digits, mult) {

  rounded_mean_mat <- round(sim_results_mean * mult, digits = digits)
  rounded_sd_mat <- round(sim_results_sd * mult, digits = digits)

  # 1. Format both matrices to 1 decimal place as strings
  fmt_mat1 <- apply(rounded_mean_mat, 2, function(x) sprintf(paste0("%.", digits ,"f"), x))
  fmt_mat2 <- apply(rounded_sd_mat, 2, function(x) sprintf(paste0("%.", digits ,"f"), x))

  # 2. Combine the formatted strings element-wise: "Value_Table1 (Value_Table2)"
  combined_mat <- paste0(fmt_mat1, " (", fmt_mat2, ")")
  dim(combined_mat) <- dim(rounded_mean_mat) # Restore matrix dimensions

  best_model_idxs <- apply(sim_results_mean, 1, which.max)
  for (i in 1:nrow(combined_mat)) {
    combined_mat[i, best_model_idxs[i]] <- paste0("\\textbf{", rounded_mean_mat[i, best_model_idxs[i]], "} (", rounded_sd_mat[i, best_model_idxs[i]], ")")
  }

  # 3. Collapse the columns for each row using " & "
  # The number 1 in apply means apply the function to rows
  output_strings <- apply(combined_mat, 1, paste, collapse = " & ")

  # 4. preface each string with the factor label, and the cursive numeral if appropriate
  ctr <- 1
  for (Np in c(1, 3)) {
    for (theta_type in c("irr", "reg")) {
      for (P in c(200, 1000)) {
        for (num_deg in c(40, 80)) {
          factor_string <- paste0(Np, " & ", theta_type, " & ", P, " & ", num_deg, " & ")
          output_strings[ctr] <- paste0(factor_string, output_strings[ctr])

          # if ((ctr %% 4) == 1) {
          #   output_strings[ctr] <- paste0("\\textit{", ctr %/% 4 + 1, "} & ", output_strings[ctr])
          # } else {
          #   output_strings[ctr] <- paste0(" & ", output_strings[ctr])
          # }

          ctr <- ctr + 1
        }
      }
    }
  }

  # # adding \hdashline
  # output_strings <- base::append(output_strings, "\\hdashline", after = 4)
  # output_strings <- base::append(output_strings, "\\hdashline", after = 9)
  # output_strings <- base::append(output_strings, "\\hdashline", after = 14)

  # print(output_strings)

  # return(output_strings)

  cat(paste0(paste(output_strings, collapse = " \\\\ \n"), " \\\\"))

}


