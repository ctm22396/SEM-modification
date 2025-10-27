library(tidyverse)
library(lavaan)
source("./R/extract/helpers.R")
source("./R/extract/extract_matrices.R")
source("./R/extract/extract_bases_correlations.R")

partial_mat <- function(mat, chosen, remain) {
  if (length(remain) == 0) {
    part_mat <- mat[chosen, chosen, drop = FALSE]
    
    return(part_mat)
  }
  
  if (length(remain) == 1) {
    part_mat <- mat[chosen, chosen, drop = FALSE] -
      tcrossprod(mat[chosen, remain, drop = FALSE]) / mat[remain, remain]
    
    return(part_mat)
  }
  
  part_mat <-  mat[chosen, chosen, drop = FALSE] -
    mat[chosen, remain, drop = FALSE] %*%
    pracma::pinv(mat[remain, remain, drop = FALSE]) %*%
    mat[remain, chosen, drop = FALSE]
  
  return(part_mat)
}

safe_varimax <- function (x, normalize = TRUE, eps = 1e-05) {
  if (ncol(x) <= 1) {
    z <- x
    TT <- diag(ncol(x))
    dimnames(z) <- dimnames(x)
    class(z) <- "loadings"
    return(list(loadings = z, rotmat = TT))
  } else {
    return(varimax(x, normalize, eps))
  }
}

hclust_reorder <- function(merge, height) {
  # Stack to hold the nodes to visit
  stack <- c(nrow(merge))
  
  # Vector to hold the final ordered leaves
  order_vector <- integer(0)
  
  # While the stack is not empty
  while (length(stack) > 0) {
    # Pop the last element from the stack
    current_node_idx <- stack[length(stack)]
    stack <- stack[-length(stack)]
    
    # If it's a leaf node (negative index)
    if (current_node_idx < 0) {
      order_vector <- c(order_vector, -current_node_idx)
      next
    }
    
    # It's an internal node, get its children
    sub_cluster1_idx <- merge[current_node_idx, 1]
    sub_cluster2_idx <- merge[current_node_idx, 2]
    
    # Get max height of children (for leaf, height is 0)
    height1 <- if (sub_cluster1_idx < 0) 0 else height[sub_cluster1_idx]
    height2 <- if (sub_cluster2_idx < 0) 0 else height[sub_cluster2_idx]
    
    # Push the tighter cluster first (LIFO stack means it gets popped first)
    if (height1 > height2) {
      stack <- c(stack, sub_cluster2_idx, sub_cluster1_idx)
    } else {
      stack <- c(stack, sub_cluster1_idx, sub_cluster2_idx)
    }
  }
  return(rev(order_vector))
}

cut_cluster <- function(out, n_clusts) {
  centering <- out$centering
  weighting <- out$weighting
  R_final <- out$R
  
  clusters <- split(seq_len(ncol(R_final)), cutree(out, k = n_clusts))
  clust_ind <- map(clusters, \(x) rowSums(diag(ncol(L))[, x, drop = FALSE])) %>%
    do.call(what = cbind)
  
  if (weighting) {
    # weighted_clust_ind <- clust_ind %*% sqrt(solve(crossprod(clust_ind)))
    Rsq_clust <- scale(R_final^2, center = centering, scale = FALSE) %*% clust_ind
    weight_matrix <- diag(1 / sqrt(colSums(Rsq_clust^2)))
    weighted_clust_ind <- clust_ind %*% weight_matrix
  } else {
    weighted_clust_ind <- clust_ind
  }
  
  R <- clustered_varimax(R_final, clust_ind, centering = centering)$loadings
  R_sum <- R^2 %*% clust_ind
  
  min_proj <- min(apply(R_sum, 1, max))
  clust_members <- apply(R_sum, 2, \(x) names(x[x >= min_proj]), simplify = FALSE)
  
  clust_order <- order(order(out$order) %*% clust_ind / colSums(clust_ind))
  
  cluster_info <- list(
    clust_members = clust_members,
    clusters = clusters,
    clust_ind = clust_ind,
    R_nclusts = R,
    R_final = R_final,
    R_sum = R_sum,
    min_proj_nclusts = min_proj,
    min_proj_final = rev(out$height)[n_clusts - 1],
    order = clust_order
  )
  
  cluster_info
}

hierarchical_MaxVar_clustering <- function(L) {
  p <- nrow(L)
  k <- ncol(L)
  bias <- (p - 1) / p
  cov_mat <- cov(L^2) * bias
  R <- L^2
  
  # mimic hclust behavior
  clust_dex <- -seq_len(k)
  clust_merge <- matrix(NA, nrow = k - 1, ncol = 2)
  
  diff_crit <- numeric(k - 1)
  varimax_crit <- numeric(k)
  min_diag_crit <- numeric(k - 1)
  for (i in -(clust_dex[-1] + 1)) {
    varimax_crit[i] <- sum(diag(cov_mat))
    min_diag_crit[i] <- min(apply(R, 1, max))
    
    uptri_dex <- which(upper.tri(cov_mat), arr.ind = TRUE)
    uptri_cov_mat <- cov_mat[uptri_dex]
    best_dex <- uptri_dex[which.max(uptri_cov_mat), ]
    
    clust_ind <- diag(k - i + 1)
    clust_ind[, best_dex[1]] <- rowSums(clust_ind[, best_dex])
    clust_ind <- clust_ind[, -best_dex[2], drop = FALSE]
    
    R <- R %*% clust_ind
    diff_crit[i] <- 2*cov_mat[best_dex, best_dex][1, 2]
    cov_mat <- crossprod(clust_ind, cov_mat * bias) %*% clust_ind
    
    clust_merge[i, ] <- clust_dex[best_dex]
    clust_dex[best_dex[1]] <- i
    clust_dex <- clust_dex[-best_dex[2]]
  }
  
  hclust_data <- list(
    merge = clust_merge,
    height = min_diag_crit,
    method = "single",
    dist.method = "crossprod_norm_sum",
    diff_crit = diff_crit,
    varimax_crit = varimax_crit,
    order = hclust_reorder(clust_merge, min_diag_crit),
    labels = colnames(L),
    call = match.call()
  )
  structure(hclust_data, class = "hclust")
}

clustered_varimax_rect <- function(x, I, Q, centering = FALSE, normalize = FALSE, eps = 1e-08) {
  nc <- ncol(as.matrix(x))
  if (nc < 2) {
    z <- x
    TT <- diag(nc)
    dimnames(z) <- dimnames(x)
    class(z) <- "loadings"
    return(list(loadings = z, rotmat = TT))
  }
  
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  
  p <- nrow(x)
  TT <- diag(nc)
  if (centering) {
    H <- diag(p) - matrix(1/p, p, p)
  } else {
    H <- diag(p)
  }
  II <- tcrossprod(I)
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT %*% Q
    
    B <- crossprod(x, (H %*% z^2 %*% II) * z) %*% t(Q)
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    
    if (d < dpast * (1 + eps)) 
      break
  }
  z <- x %*% TT
  
  if (normalize) {
    z <- z * sc
  }
  
  dimnames(z) <- dimnames(x)
  class(z) <- "loadings"
  list(loadings = z, rotmat = TT)
}

clustered_varimax <- function(x, I, centering = FALSE, normalize = FALSE, eps = 1e-08) {
  nc <- ncol(as.matrix(x))
  if (nc < 2) {
    z <- x
    TT <- diag(nc)
    dimnames(z) <- dimnames(x)
    class(z) <- "loadings"
    return(list(loadings = z, rotmat = TT))
  }
  
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  
  p <- nrow(x)
  TT <- diag(nc)
  if (centering) {
    H <- diag(p) - matrix(1/p, p, p)
  } else {
    H <- diag(p)
  }
  II <- tcrossprod(I)
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT
    
    B <- crossprod(x, (H %*% z^2 %*% II) * z)
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    
    if (d < dpast * (1 + eps)) 
      break
  }
  z <- x %*% TT
  
  if (normalize) {
    z <- z * sc
  }
  
  dimnames(z) <- dimnames(x)
  class(z) <- "loadings"
  list(loadings = z, rotmat = TT)
}

# TODO: Writeup
# - important to try to pick a cut point where min_diag (height) > 0.5 so
#   row assignment is unique
# yes centering, no weighting, but its so weird with samp_select
# no centering, yes weighting, is good for samp_select, though... very weird for full matrix though
# - also looks really good for clusters_full with 8 or 9 n_clusts
# - not sure about sqrt or not to sqrt... pretty sure sqrt is right as it ends
#   up weighting it so variance is penalized by number of components
# - however the non-sqrt nclusts = 7 for clusters_full looks nice too
# TODO
# TODO wait cov2cor uptri_cov_mat looks GREAT with centering = FALSE
# TODO
# - especially with the clusters_full mat. Its like... so good.
# - just look at how the within error-covs are at the front of the pack...
# - it's still fine with the full matrix
# - but the fact that I get 11 clusters with the same min_proj really
#   goes to show how good just a little filtering
# TODO STUDY: looks like centering produces the same cluster structure for 
#      clusters_full, n_clusts = 9
# - however, the rotation yields much lower min_diag
# TODO: Consider switching to the BD rot after min_diag passes 0.5
hierarchical_MaxVar_clustering_SimpRot <- function(L, centering = FALSE, weighting = TRUE) {
  p <- nrow(L)
  k <- ncol(L)
  bias <- (p - 1) / p
  R <- L
  diag_mat <- diag(ncol(L))
  
  # mimic hclust behavior
  clust_dex <- -seq_len(k)
  clust_merge <- matrix(NA, nrow = k - 1, ncol = 2)
  
  clusters <- as.list(-clust_dex)
  clust_ind <- map(clusters, \(x) rowSums(diag_mat[, x, drop = FALSE])) %>%
    do.call(what = cbind)
  
  diff_crit <- numeric(k - 1)
  varimax_crit <- numeric(k)
  min_diag_crit <- numeric(k - 1)
  for (i in -(clust_dex[-1] + 1)) {
    if (weighting) {
      # weighted_clust_ind <- clust_ind %*% sqrt(solve(crossprod(clust_ind)))
      Rsq_clust <- scale(R^2, center = centering, scale = FALSE) %*% clust_ind
      weight_matrix <- diag(1 / sqrt(colSums(Rsq_clust^2)))
      weighted_clust_ind <- clust_ind %*% weight_matrix
    } else {
      weighted_clust_ind <- clust_ind
    }
    
    R <- clustered_varimax(R, weighted_clust_ind, centering = centering)$loadings
    # no more weighting because weigthing depends on R before rotation
    R_sum <- R^2 %*% clust_ind
    
    if (centering) {
      cov_mat <- cov(R_sum) * bias
    } else {
      cov_mat <- crossprod(R_sum) / p
    }
    vars <- diag(cov_mat)
    
    varimax_crit[i] <- mean(vars)
    min_diag_crit[i] <- min(apply(R^2 %*% clust_ind, 1, max))
    
    uptri_dex <- which(upper.tri(cov_mat), arr.ind = TRUE)
    uptri_cov_mat <- cov_mat[uptri_dex]
    drop_vars <- apply(uptri_dex, 1, \(x) sum(vars[-x]))
    
    row_vars <- vars[uptri_dex[, 1]]
    col_vars <- vars[uptri_dex[, 2]]
    
    # should be some way to handle this differently with weighting
    if (weighting) {
      # TODO old_method 
      # clust_sizes <- colSums(clust_ind)
      # row_sizes <- clust_sizes[uptri_dex[, 1]]
      # col_sizes <- clust_sizes[uptri_dex[, 2]]
      # comb_vars <- row_vars*row_sizes^2 + col_vars*col_sizes^2 + 2*uptri_cov_mat*row_sizes*col_sizes
      # comb_vars <- comb_vars / (row_sizes + col_sizes)^2
      # new_totals <- 2*uptri_cov_mat # *row_sizes*col_sizes / (row_sizes + col_sizes)^2
      
      # cov2cor method
      new_totals <- uptri_cov_mat/sqrt(row_vars*col_vars)
    } else {
      comb_vars <- row_vars + col_vars + 2*uptri_cov_mat
      new_totals <- (drop_vars + comb_vars) / (length(clusters) - 1)
    }
    

    # new_totals <- (drop_vars + comb_vars) / (length(clusters) - 1)
    best_dex <- uptri_dex[which.max(new_totals), ]
    diff_crit[i] <- max(new_totals) - ifelse(weighting, 0, varimax_crit[i])
    
    clusters[[best_dex[1]]] <- unlist(clusters[best_dex])
    clusters <- clusters[-best_dex[2]]
    clust_ind <- map(clusters, \(x) rowSums(diag_mat[, x, drop = FALSE])) %>%
      do.call(what = cbind)
    
    clust_merge[i, ] <- clust_dex[best_dex]
    clust_dex[best_dex[1]] <- i
    clust_dex <- clust_dex[-best_dex[2]]
  }
  
  varimax_crit[k] <- 0
  
  hclust_data <- list(
    merge = clust_merge,
    height = min_diag_crit,
    method = "single",
    dist.method = "crossprod_norm_sum",
    diff_crit = diff_crit,
    varimax_crit = varimax_crit,
    order = hclust_reorder(clust_merge, min_diag_crit),
    labels = colnames(L),
    R = R,
    centering = centering,
    weighting = weighting,
    call = match.call()
  )
  structure(hclust_data, class = "hclust")
}

samp_mats <- read_rds("./data/lvShipley_writeup_samp_mats.rds")
samp_ortho_std <- samp_mats$samp_ortho_std
# a <- unlist(samp_mats$clusters_full)
a <- samp_mats$samp_select
# a <- -1
# a <- colnames(samp_ortho_std)
mat <- samp_ortho_std[a, a]
L <- compute_basis(mat)

out <- hierarchical_MaxVar_clustering_SimpRot(L, centering = FALSE, weighting = TRUE)
plot(out)
res <- cut_cluster(out, sum(out$height > 0.75) + 1)
