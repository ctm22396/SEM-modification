clustered_varimax_pooling <- function(x, I, normalize = FALSE, eps = 1e-08) {
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
  G <- matrix(1/p, p, p)
  H <- diag(p) - G
  I_pinv <- tcrossprod(solve(crossprod(I)), I)
  P <- I %*% I_pinv
  Q <- diag(nc) - P
  D <- diag(colSums(I_pinv))
  d <- 0
  # I_pinv %*% diag(crossprod(R^2 %*% Q, diag(p) - H) %*% R^2 %*% Q + crossprod(R^2, H) %*% R^2) / p
  # OR: I_pinv %*% (colMeans((H %*% z^2)^2) + colMeans((G %*% z^2 %*% Q)^2))
  # OR: rowMeans(tcrossprod(I_pinv, (H %*% z^2)^2 + (G %*% z^2 %*% Q)^2))
  for (i in 1L:1000L) {
    z <- x %*% TT
    B <- crossprod(x, ((H %*% z^2) + G %*% z^2 %*% Q) %*% D * z)
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

compute_pooled_means <- function(x, I) {
  stopifnot(ncol(x) == nrow(I))
  
  nc <- ncol(x)
  p <- nrow(x)
  
  I_pinv <- tcrossprod(solve(crossprod(I)), I)
  P <- I %*% I_pinv
  
  rowMeans(tcrossprod(I_pinv, x^2 %*% P))
}

compute_pooled_variances <- function(x, I) {
  stopifnot(ncol(x) == nrow(I))
  
  nc <- ncol(x)
  p <- nrow(x)
  
  G <- matrix(1/p, p, p)
  H <- diag(p) - G
  I_pinv <- tcrossprod(solve(crossprod(I)), I)
  P <- I %*% I_pinv
  Q <- diag(nc) - P
  
  rowMeans(tcrossprod(I_pinv, (H %*% x^2)^2 + (G %*% x^2 %*% Q)^2))
}

pairwise_pooling <- function(mus, sigs) {
  stopifnot(length(mus) == 2 & length(sigs) == 2)
  
  # law of total variance: variance of means + mean of variances.
  # population variance of two values is difference squared divided by 4.
  diff(mus)^2/4 + mean(sigs)
}

hierarchical_MaxVar_clustering_PoolRot <- function(L) {
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
    R <- clustered_varimax_pooling(R, clust_ind)$loadings
    means <- compute_pooled_means(R, clust_ind)
    vars <- compute_pooled_variances(R, clust_ind)
    
    varimax_crit[i] <- mean(vars)
    min_diag_crit[i] <- min(apply(R^2 %*% clust_ind, 1, max))
    
    comb_dex <- combn(length(clusters), 2)
    drop_vars <- apply(comb_dex, 2, \(x) sum(vars[-x]))
    comb_vars <- apply(comb_dex, 2, \(x) pairwise_pooling(means[x], vars[x]))
    new_totals <- (comb_vars + drop_vars) / (length(clusters) - 1)
    best_dex <- comb_dex[, which.max(new_totals)]
    diff_crit[i] <- max(new_totals) - varimax_crit[i]
    
    # min_diag criteria
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

    
    clusters[[best_dex[1]]] <- unlist(clusters[best_dex])
    clusters <- clusters[-best_dex[2]]
    clust_ind <- map(clusters, \(x) rowSums(diag_mat[, x, drop = FALSE])) %>%
      do.call(what = cbind)
    
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
    R = R,
    call = match.call()
  )
  structure(hclust_data, class = "hclust")
}
