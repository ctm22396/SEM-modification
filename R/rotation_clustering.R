library(tidyverse)
library(lavaan)
source("./R/extract/helpers.R")
source("./R/extract/extract_matrices.R")
source("./R/extract/extract_bases_correlations.R")

transitive_cov2cor <- function(x) {
  A <- cov2cor(crossprod(x))
  for (i in seq_len(ncol(x) - 1)) {
    A <- cov2cor(crossprod(A))
  }
  return(cov2cor(A))
}

transitive_closure <- function(x) {
  A <- crossprod(x) > eps
  diag(A) <- TRUE
  for (i in seq_len(ncol(x) - 1)) {
    A <- crossprod(A) > eps
    diag(A) <- TRUE
  }
  return(A > eps)
}

bitransitive_closure <- function(x) {
  A <- transitive_closure(x)
  B <- transitive_closure(t(x))
  return((B %*% x %*% A) > eps)
}

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

clustered_varimax <- function(x, I = diag(ncol(x)), centering = FALSE, normalize = FALSE, eps = 1e-08) {
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
  II <- tcrossprod(tcrossprod(I) %*% MASS::ginv(tcrossprod(I)))
  if (centering) {
    H <- diag(p) - matrix(1/p, p, p)
  } else {
    H <- diag(p)
  }
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT
    
    B <- crossprod(x, (H %*% z^2 %*% II) * z)
    # B <- crossprod(x, z^3 %*% diag(rowSums(I)))
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

clustered_parsimax <- function(x, I = diag(ncol(x)), centering = FALSE, normalize = FALSE, eps = 0) {
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
  q <- ncol(x)
  kappa <- 0 # (q - 1)/(p + q - 2)
  M <- (1 - kappa)*diag(p) + kappa*matrix(1, p, p)
  # INI <- tcrossprod(I %*% N, I)
  for (i in 1L:1000L) {
    z <- x %*% TT
    B <- crossprod(x, (M %*% z^2 %*% II) * z)
    # B <- crossprod(x, (H %*% z^2 %*% II) * z)
    # Technically equivalent to centering = FALSE
    # B <- crossprod(x, z^3 %*% diag(rowSums(I)))
    sB <- La.svd(B)
    # sB <- La.svd(TT + 0.2*B)
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


hierarchical_MaxVar_clustering_ParsiRot <- function(L, rotation = TRUE,
                                                    centering = FALSE,
                                                    weighting = TRUE) {
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
  corFnorm_crit <- numeric(k)
  min_diag_crit <- numeric(k - 1)
  for (i in -(clust_dex[-1] + 1)) {
    q <- ncol(clust_ind)
    kappa <- (q - 1)/(p + q - 2)
    M <- (1 - kappa)*diag(p) + kappa*matrix(1, p, p)
    if (weighting) {
      # weighted_clust_ind <- clust_ind %*% sqrt(solve(crossprod(clust_ind)))
      R_sum <- R^2 %*% clust_ind
      if (centering) {
        cov_mat <- crossprod(R_sum) - kappa * tcrossprod(colMeans(R_sum)) * p
        # cov_mat <- cov(R_sum) * bias
      } else {
        cov_mat <- crossprod(R_sum)
      }
      weight_matrix <- diag(1 / sqrt(diag(cov_mat)))
      weighted_clust_ind <- clust_ind %*% weight_matrix
      # weighted_clust_ind <- clust_ind
    } else {
      weighted_clust_ind <- clust_ind
    }
    
    if (rotation) {
      R <- clustered_parsimax(R, weighted_clust_ind, centering = centering)$loadings
    }
    # no more weighting because weigthing depends on R before rotation
    R_sum <- R^2 %*% clust_ind
    
    if (centering) {
      cov_mat <- (crossprod(R_sum) - kappa * tcrossprod(colMeans(R_sum)) * p) / p
      # cov_mat <- cov(R_sum) * bias
    } else {
      cov_mat <- crossprod(R_sum) / p
    }
    vars <- diag(cov_mat)
    
    clust_sizes <- colSums(clust_ind)
    varimax_crit[i] <- sum(clust_sizes * vars) / sum(clust_sizes)
    corFnorm_crit[i] <- 2 - sum(cov2cor(cov_mat)^2) / ncol(cov_mat)
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
    corFnorm_crit = corFnorm_crit,
    order = hclust_reorder(clust_merge, min_diag_crit),
    labels = colnames(L),
    R = R,
    centering = centering,
    weighting = weighting,
    call = match.call()
  )
  structure(hclust_data, class = "hclust")
}

clustered_cubemax <- function(x, I = diag(ncol(x)), normalize = FALSE, eps = 1e-08, ...) {
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
  q <- ncol(I)
  TT <- diag(nc)
  
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT
    
    B <- crossprod(x, z * I)
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
hierarchical_MaxVar_clustering_SimpRot <- function(L, rotation = TRUE,
                                                   centering = FALSE,
                                                   weighting = TRUE) {
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
  corFnorm_crit <- numeric(k)
  min_diag_crit <- numeric(k - 1)
  for (i in -(clust_dex[-1] + 1)) {
    # if (weighting) {
    #   # weighted_clust_ind <- clust_ind %*% sqrt(solve(crossprod(clust_ind)))
    #   # Rsq_clust <- scale(R^2, center = TRUE, scale = FALSE) %*% clust_ind
    #   # weight_matrix <- diag(1 / sqrt(colSums(Rsq_clust^2)))
    #   # weighted_clust_ind <- clust_ind %*% weight_matrix
    #   weighted_clust_ind <- clust_ind
    # } else {
    #   weighted_clust_ind <- clust_ind
    # }
    
    if (rotation) {
      R <- clustered_varimax(R, clust_ind, centering = TRUE)$loadings
    }
    # no more weighting because weigthing depends on R before rotation
    R_sum <- R^2 %*% clust_ind
    
    if (centering) {
      # q <- ncol(clust_ind)
      # kappa <- (q - 1)/(p + q - 2)
      # cov_mat <- crossprod(R_sum) - kappa * tcrossprod(colMeans(R_sum)) * p
      cov_mat <- cov(R_sum) * bias
      # R_weight <- t(t(R_sum) - colMeans(R_sum) / colSums(clust_ind))
      # cov_mat <- crossprod(R_weight) / p
    } else {
      cov_mat <- crossprod(R_sum) / p
    }
    vars <- diag(cov_mat)
    
    clust_sizes <- colSums(clust_ind)
    varimax_crit[i] <- sum(clust_sizes * vars) / sum(clust_sizes)
    corFnorm_crit[i] <- 2 - sum(cov2cor(cov_mat)^2) / ncol(cov_mat)
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
    corFnorm_crit = corFnorm_crit,
    order = hclust_reorder(clust_merge, min_diag_crit),
    labels = colnames(L),
    R = R,
    centering = centering,
    weighting = weighting,
    call = match.call()
  )
  structure(hclust_data, class = "hclust")
}

# TODO: Best yet, this clustered avg_projmax + angle merging 
avg_projmax <- function(x, R, C, normalize = FALSE, centering = FALSE, eps = 1e-08) {
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
  q <- ncol(R)
  TT <- diag(nc)
  
  if (centering) {
    H <- diag(p) - matrix(1/p, p, p)
    # J <- diag(q) - matrix(1/q, q, q)
  } else {
    H <- diag(p)
    # J <- diag(q)
  }
  
  M <- H %*% tcrossprod(R, C)
  # U <- R %*% tcrossprod(J, R)
  # V <- tcrossprod(C)
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT
    
    B <- crossprod(x, (z * M))
    # B <- crossprod(x, (U %*% z^2 %*% V) * z)
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

project_onto_Stiefel <- function(x) {
  sX <- La.svd(x)
  sX$u %*% sX$vt
}

clustered_minmax <- function(x, I = x*0 + 1, normalize = FALSE, eps = 1e-08, ...) {
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
  q <- ncol(I)
  TT <- diag(nc)
  d <- 0
  Z <- matrix(0, p, p)
  step <- 1 / nc
  for (i in 1L:1000L) {
    z <- x %*% TT
    q <- diag(tcrossprod(z * I))
    k <- which.min(q)
    P_k <- tcrossprod(seq_len(p) == k)
    B <- crossprod(x, P_k %*% z * I)
    TT_past <- TT
    sB <- La.svd(TT + step*B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    print(min(q))
    B_past <- B
    B <- crossprod(x, P_k %*% (z %*% TT) * I)
    # step <- get_step(TT - TT_past, B - B_past)
    # print(step)
    if (abs(step) < eps)
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

hierarchical_MaxVar_clustering_avgMax <- function(L, centering = FALSE, rotation = TRUE) {
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
  corFnorm_crit <- numeric(k)
  min_diag_crit <- numeric(k - 1)
  for (i in -(clust_dex[-1] + 1)) {
    # Ensure each cluster has variable associated with it.
    # Eventually capture any that exceed 0.5 (unique association for rows)
    # C <- apply(R^2 %*% clust_ind, 2, \(x) x >= max(x)) | (R^2 %*% clust_ind > 0.5)
    # Alternatively, only rotate when unique association achieved for all rows
    C <- (R^2 %*% clust_ind) > 0.5
    if (rotation & all(rowSums(C) > 0)) {
      R <- avg_projmax(R, C, clust_ind, centering = centering)$loadings
      
      # Rotate within clusters based of C and clust_ind
      # mask <- tcrossprod(C, clust_ind)
      # R <- R %*% clustered_varimax(R * mask, centering = FALSE, normalize = TRUE, eps = eps)$rotmat
    }
    R_sum <- R^2 %*% clust_ind
    
    if (centering) {
      # q <- ncol(clust_ind)
      # # kappa <- 1 / p # varimax
      # kappa <- q / (q + p - 2) # parsimax
      # # kappa <- q / p / 2 # equamax
      # R_kappa <- t(t(R_sum) - kappa * colMeans(R_sum)) -
      #   (1 - kappa) * rowSums(R_sum) / k +
      #   kappa * (1 - kappa) * mean(R_sum)
      # cov_mat <- crossprod(R_kappa) / p
      cov_mat <- cov(R_sum) * bias
      # R_weight <- t(t(R_sum) - colMeans(R_sum) / colSums(clust_ind))
      # cov_mat <- crossprod(R_weight) / p
    } else {
      cov_mat <- crossprod(R_sum) / p
    }
    vars <- diag(cov_mat)
    
    clust_sizes <- colSums(clust_ind)
    varimax_crit[i] <- sum(clust_sizes * vars) / sum(clust_sizes)
    corFnorm_crit[i] <- mean(diag(crossprod(C, R^2) %*% clust_ind / colSums(C)))
    # 2 - sum(cov2cor(cov_mat)^2) / ncol(cov_mat)
    min_diag_crit[i] <- min(apply(R^2 %*% clust_ind, 1, max))
    
    uptri_dex <- which(upper.tri(cov_mat), arr.ind = TRUE)
    uptri_cov_mat <- cov_mat[uptri_dex]
    drop_vars <- apply(uptri_dex, 1, \(x) sum(vars[-x]))
    
    row_vars <- vars[uptri_dex[, 1]]
    col_vars <- vars[uptri_dex[, 2]]
    
    # min_diag criteria
    # comb_dex <- combn(length(clusters), 2)
    # if (length(clusters) > 2) {
    #   drop_maxs <- apply(comb_dex, 2, \(x) apply(R^2 %*% clust_ind[, -x], 1, max))
    #   comb_maxs <- apply(comb_dex, 2, \(x) R^2 %*% rowSums(clust_ind[, x]))
    #   new_maxs <- pmax(drop_maxs, comb_maxs)
    #   new_min <- apply(new_maxs, 2, min)
    #   best_dex <- comb_dex[, which.max(new_min)]
    #   diff_crit[i] <- max(new_min)
    # } else {
    #   best_dex <- comb_dex[, 1]
    #   diff_crit[i] <- 1
    # }
    
    new_totals <- uptri_cov_mat/sqrt(row_vars*col_vars)
    best_dex <- uptri_dex[which.max(new_totals), ]
    diff_crit[i] <- max(new_totals)
    
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
    dist.method = "avg_projmax",
    diff_crit = diff_crit,
    varimax_crit = varimax_crit,
    corFnorm_crit = corFnorm_crit,
    order = hclust_reorder(clust_merge, min_diag_crit),
    labels = colnames(L),
    R = R,
    centering = FALSE,
    weighting = FALSE,
    call = match.call()
  )
  structure(hclust_data, class = "hclust")
}

# Non variance version
# hierarchical_MaxVar_clustering_avgMax <- function(L, centering = FALSE, rotation = TRUE) {
#   p <- nrow(L)
#   k <- ncol(L)
#   bias <- (p - 1) / p
#   R <- L
#   diag_mat <- diag(ncol(L))
#   
#   # mimic hclust behavior
#   clust_dex <- -seq_len(k)
#   clust_merge <- matrix(NA, nrow = k - 1, ncol = 2)
#   
#   clusters <- as.list(-clust_dex)
#   clust_ind <- map(clusters, \(x) rowSums(diag_mat[, x, drop = FALSE])) %>%
#     do.call(what = cbind)
#   
#   diff_crit <- numeric(k - 1)
#   varimax_crit <- numeric(k)
#   corFnorm_crit <- numeric(k)
#   min_diag_crit <- numeric(k - 1)
#   for (i in -(clust_dex[-1] + 1)) {
#     C <- apply(R^2 %*% clust_ind, 1, \(x) x >= max(x)) %>% t()
#     if (rotation) {
#       R <- avg_projmax(R, C, clust_ind, centering = centering)$loadings
#     }
#     cov_mat <- crossprod(C, scale(R^2, center = centering, scale = FALSE)) %*% clust_ind
#     
#     vars <- diag(cov_mat)
#     
#     clust_sizes <- colSums(clust_ind)
#     varimax_crit[i] <- sum(clust_sizes * vars) / sum(clust_sizes)
#     corFnorm_crit[i] <- mean(diag(crossprod(C, scale(R^2, center = centering, scale = FALSE)) %*% clust_ind / colSums(C)))
#     # 2 - sum(cov2cor(cov_mat)^2) / ncol(cov_mat)
#     min_diag_crit[i] <- min(apply(R^2 %*% clust_ind, 1, max))
#     
#     uptri_dex <- which(upper.tri(cov_mat), arr.ind = TRUE)
#     uptri_cov_mat <- cov_mat[uptri_dex] + t(cov_mat)[uptri_dex]
#     
#     row_vars <- vars[uptri_dex[, 1]]
#     col_vars <- vars[uptri_dex[, 2]]
#     
#     new_totals <- uptri_cov_mat#/sqrt(row_vars*col_vars)
#     best_dex <- uptri_dex[which.max(new_totals), ]
#     diff_crit[i] <- max(new_totals)
#     
#     clusters[[best_dex[1]]] <- unlist(clusters[best_dex])
#     clusters <- clusters[-best_dex[2]]
#     clust_ind <- map(clusters, \(x) rowSums(diag_mat[, x, drop = FALSE])) %>%
#       do.call(what = cbind)
#     
#     clust_merge[i, ] <- clust_dex[best_dex]
#     clust_dex[best_dex[1]] <- i
#     clust_dex <- clust_dex[-best_dex[2]]
#   }
#   
#   varimax_crit[k] <- 0
#   
#   hclust_data <- list(
#     merge = clust_merge,
#     height = min_diag_crit,
#     method = "single",
#     dist.method = "avg_projmax",
#     diff_crit = diff_crit,
#     varimax_crit = varimax_crit,
#     corFnorm_crit = corFnorm_crit,
#     order = hclust_reorder(clust_merge, min_diag_crit),
#     labels = colnames(L),
#     R = R,
#     centering = FALSE,
#     weighting = FALSE,
#     call = match.call()
#   )
#   structure(hclust_data, class = "hclust")
# }


bdd <- function(x, cluster_members) {
  L <- compute_basis(x)
  nc <- ncol(L)
  
  # Step 1: no col indicator
  J <- map(cluster_members, \(membs) rownames(L) %in% membs) %>%
    do.call(what = cbind)
  L_BD1 <- compute_basis(tcrossprod(L) * tcrossprod(J))[, seq_len(nc)]
  I1 <- crossprod(L_BD1^2, J) > eps
  
  Q1 <- reduce(La.svd(crossprod(L, L_BD1))[-1], `%*%`)
  sign_fix1 <- c(rep(1, nc - 1), det(Q1))
  Q1 <- sweep(Q1, 2, sign_fix1, "*")
  
  R_BD1 <- L %*% Q1
  
  # Step 2: with col indicator
  L_BD2 <- compute_basis(cov2cor(tcrossprod(R_BD1 * tcrossprod(J, I1))))
  I2 <- crossprod(L_BD2^2, J) > eps
  
  Q2 <- reduce(La.svd(crossprod(R_BD1, L_BD2))[-1], `%*%`)
  sign_fix2 <- c(rep(1, nc - 1), det(Q2))
  Q2 <- sweep(Q2, 2, sign_fix2, "*")
  R_BD2 <- R_BD1 %*% Q2
  
  list(R = R_BD2, rotmat = Q1 %*% Q2, clust_ind = I2)
}


# Sandbox -----------------------------------------------------------------

# mat <- cov2cor(read_rds("./data/JAMO_cov_mat.rds")) # WORKS WAY BETTER WITH CENTERING

# TODO: WHY CENTERING SO WEIRD FOR SAMP_SELECT


samp_mats <- read_rds("./data/lvShipley_writeup_samp_mats.rds")
samp_ortho_std <- samp_mats$samp_ortho_std
a <- unlist(samp_mats$clusters_full)
# a <- samp_mats$samp_select
# a <- samp_mats$samp_select[-c(15:17, 19, 20, 21, 23, 26, 27)]
# a <- -1
# a <- colnames(samp_ortho_std)
mat <- samp_ortho_std[a, a]
L <- compute_basis(mat)
R <- clustered_varimax(L, eps = 0)$loadings

L2 <- compute_basis(MASS::ginv(mat) + 0*mat)
R2 <- varimax(L2, normalize = FALSE, eps = eps)$loadings
A <- bitransitive_closure(R2^2 > 0.023) # works for clusters_full

# out <- hierarchical_MaxVar_clustering_avgMax(R, centering = FALSE, rotation = TRUE)
out <- hierarchical_MaxVar_clustering_SimpRot(R, rotation = TRUE, centering = FALSE, weighting = TRUE)
plot(out, hang = -1)
res <- cut_cluster(out, 10)
res$clust_members[res$order]
plot_correlations_clustered(abs(mat), res$clust_members[res$order])
# map2(res$clust_members, res$clusters, \(x, y) res$R_final[x, y, drop = FALSE])[orig_order] %>%
#   map(tcrossprod) %>% map(\(x) x^2) %>% map(plot_correlations)

# test that mutually disjoint
map(seq_along(res$clust_members), \(i) unlist(res$clust_members[-i])) %>%
  map2(res$clust_members, intersect) %>%
  unlist()

resid_proj <- samp_ortho_std[1, a] %*% MASS::ginv(mat) %*% res$R_final
R_full <- rbind(resids = as.vector(resid_proj), res$R_final)

1 - pbeta(resid_proj^2 %*% res$clust_ind, colSums(res$clust_ind)/2, (ncol(R) - colSums(res$clust_ind))/2)[res$order]

# sig clusters when a = -1; avgMax no centering; n_clusts = 10
# maps to alph/eps and gam/delt
# I guess kinda stinks bc missing beta/delt
res$clust_members[res$order][c(2, 10)] 
