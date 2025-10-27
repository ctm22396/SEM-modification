library(lavaan)
library(Matrix)
# library(tidyverse)
library(purrr)
library(magrittr)

source("./extract/extract_matrices.R")
source("./extract/extract_bases_correlations.R")
source("./extract/helpers.R")

# start with full matrix M
# choose some column and place into a separate subset
# minimize some criterion related to the overlap of the subsets
# eventually, some condition triggers pruning of the columns to improve the paritioning


# start with full matrix M
# choose some subset of columns to focus on
# attempt to partition these columns
## criterion: pseudo frobenius inner product
## compute criterion between each subset
## normalize
## merge two with largest criterion
# at some point, prune membership based on (conditional) residual projection

pseudo_to_hclust <- function(fit) {
  clust_fit <- list()
  if (any(map_dbl(fit[[5]], length) > 1)) {
    clust_fit$labels <- seq_along(fit[[5]])
    clust_fit$order <- map(fit[[5]], match, table = unlist(fit[[1]])) %>%
      map_dbl(min) %>%
      order()
  } else {
    clust_fit$labels <- unlist(fit[[5]])
    clust_fit$order <- order(match(clust_fit$labels, unlist(fit[[1]])))
  }
  
  clust_fit$merge <- fit[[4]]
  
  clust_fit$height <- round(1 - fit[[3]], 3)
  
  class(clust_fit) <- "hclust"
  
  return(clust_fit)
}

BCvarimax <- function (x, U, C, normalize = TRUE, eps = 1e-05) {
  nc <- ncol(x)
  
  if (nc < 2) 
    return(x)
  
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  
  p <- nrow(x)
  TT <- diag(nc)
  
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT
    
    BV <- t(x) %*% (z^3 - z %*% diag(drop(rep(1, p) %*% z^2))/p)
    BC <- crossprod(x, tcrossprod(U, C) * z) / 2
    B <- BV + BC
    
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if (d < dpast * (1 + eps)) 
      break
  }
  
  z <- x %*% TT
  dimnames(z) <- dimnames(x)
  
  if (normalize) 
    z <- z * sc
  
  class(z) <- "loadings"
  list(loadings = z, rotmat = TT)
}


Cvarimax <- function (x, clusters, normalize = TRUE, eps = 1e-05) {
  params <- unlist(clusters)
  x <- x[params, ]
  
  U <- map(clusters, ~ which(params %in% .x)) %>%
    map(~ seq_along(params) %in% .x) %>%
    map(as.numeric) %>%
    do.call(what = cbind)
  rownames(U) <- params
  
  
  nc <- ncol(x)
  
  if (nc < 2) 
    return(x)
  
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  
  p <- nrow(x)
  P_U <- diag(p) - U %*% solve(crossprod(U)) %*% t(U)
  
  TT <- diag(nc)
  d <- 0
  for (i in 1L:1000L) {
    z <- x %*% TT
    
    BV <- t(x) %*% (z^3 - z %*% diag(drop(rep(1, p) %*% z^2))/p)
    BC <- t(x) %*% ((P_U %*% z^2) * z)
    B <- BV + BC
    
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if (d < dpast * (1 + eps)) 
      break
  }
  z <- x %*% TT
  dimnames(z) <- dimnames(x)
  
  ## FIXME the rowSums(A) > 1 is placing a lot of misfit outside of the clusters
  
  w <- z[, colSums(z^2) >= 0.5]
  dex <- which(colSums(z^2) >= 0.5)
  A <- solve(crossprod(w^2)) %*% t(w^2) %*% U^2
  A <- round(A)
  A[A != 1] <- 0
  A[rowSums(A) > 1, ] <- 0
  A_order <- map(seq_along(clusters), ~ dex[which(A[, .x] == 1)]) %>%
    map2(clusters, ~ .x[order(-diag(var(z[.y, .x]^2)))]) %>%
    unlist() %>%
    c(setdiff(seq_len(nc), .))
  
  A <- solve(crossprod(z^2)) %*% t(z^2) %*% U
  A[round(A) != 1] <- 0
  A <- round(A)
  A[rowSums(A) > 1, ] <- 0
  A_order <- map(seq_along(clusters), ~ which(A[, .x] == 1)) %>%
    map2(clusters, ~ .x[order(-diag(var(z[.y, .x]^2)))]) %>%
    unlist() %>%
    c(setdiff(seq_len(nc), .))
  
  if (normalize) 
    z <- z * sc
  
  class(z) <- "loadings"
  list(loadings = z, rotmat = TT, A = A, U = U)
}

function(R, U) {
  total_mat <- crossprod(U, R^2)
  perc_mat <- sweep(total_mat, 1, diag(crossprod(U)), "/")
  
  clust_dexes <- map(clusters, ~ numeric(0))
  for (i in seq_len(ncol(perc_mat))) {
    dex <- arrayInd(which.max(perc_mat), dim(perc_mat))
    
  }
  
}


Cvarimax2 <- function (x, clusters, normalize = TRUE, eps = 1e-05) {
  params <- unlist(clusters)
  x <- x[params, ]
  U <- map(clusters, ~ which(params %in% .x)) %>%
    map(~ seq_along(params) %in% .x) %>%
    map(as.numeric) %>%
    do.call(what = cbind)
  rownames(U) <- params
  
  
  nc <- ncol(x)
  
  if (nc < 2) 
    return(x)
  
  if (normalize) {
    sc <- sqrt(drop(apply(x, 1L, function(x) sum(x^2))))
    x <- x/sc
  }
  
  p <- nrow(x)
  P_U <- diag(p) - U %*% solve(crossprod(U)) %*% t(U)
  
  TT <- diag(nc)
  d <- 0
  for (i in 1L:1000L) {
    z <- P_U %*% x %*% TT
    
    BV <- t(x) %*% (z^3 - z %*% diag(drop(rep(1, p) %*% z^2))/p)
    B <- BV
    
    sB <- La.svd(B)
    TT <- sB$u %*% sB$vt
    dpast <- d
    d <- sum(sB$d)
    if (d < dpast * (1 + eps)) 
      break
  }
  z <- x %*% TT
  dimnames(z) <- dimnames(x)
  
  A <- solve(crossprod(z^2)) %*% t(z^2) %*% U
  A[round(A) != 1] <- 0
  A <- round(A)
  A[rowSums(A) > 1, ] <- 0
  A_order <- map(seq_along(clusters), ~ which(A[, .x] == 1)) %>%
    map2(clusters, ~ .x[order(-diag(var(z[.y, .x]^2)))]) %>%
    unlist() %>%
    c(setdiff(seq_len(nc), .))
  
  if (normalize) 
    z <- z * sc
  
  class(z) <- "loadings"
  list(loadings = z[, A_order], rotmat = TT[, A_order])
}

# TODO: method where we don't cov2cor and take the min of max-inner_prod

#TODO: ADDED ADDITIVE PENALTY FOR CLUSTER SIZE

# NOTE: This is equivalent to the correlation when x and y are length 1
# Even cooler, this is the correlation when x and y and both a set of vectors
# that span a single dimension each.
inner_prod <- function(x, y, mat) {
  lx <- length(x)
  ly <- length(y)

  sum(La.svd(mat[x, y], nu = 0, nv = 0)$d) / sqrt(lx * ly)
}

# Alternative, measures off-diagonal block norm
# inner_prod <- function(x, y, mat) {
#   # xy <- union(x, y)
#   # if (any(mat[xy, xy]^2 < sqrt(.Machine$double.eps))) {
#   #   return(0)
#   # }
# 
#   # x_norm <- norm(mat[x, x, drop = FALSE], type = "F")
#   # y_norm <- norm(mat[y, y, drop = FALSE], type = "F")
#   # half_dist_norm <- norm(mat[x, y, drop = FALSE], type = "F")
# 
#   # sqrt(max(half_dist_norm^2 / (length(x) * length(y)), 0))
#   median(mat[x, y]^2) # / sqrt(median(mat[x, x]^2) * median(mat[y, y]^2))
#   # median((MASS::ginv(mat[x, x]) %*% mat[x, y] %*% MASS::ginv(mat[y, y]))^2) /
#   #   sqrt(median(mat[x, x]^2) * median(mat[y, y]^2))
# }

# inner_prod <- function(x, y, mat) {
#   # xy <- union(x, y)
#   # if (any(mat[xy, xy]^2 < sqrt(.Machine$double.eps))) {
#   #   return(0)
#   # }
#   
#   x_norm <- norm(mat[x, x, drop = FALSE], type = "F")
#   y_norm <- norm(mat[y, y, drop = FALSE], type = "F")
#   half_dist_norm <- norm(mat[x, y, drop = FALSE], type = "F")
#   
#   sqrt(max(half_dist_norm^2 / (x_norm * y_norm), 0))
# }

# # Measures sum of squared off-diag correlations
# inner_prod <- function(x, y, mat) {
#   # if (any(mat[x, y]^2 < sqrt(.Machine$double.eps))) {
#   #   return(-Inf)
#   # }
#   
#   full_norm <- norm(mat, type = "F")
#   x_norm <- norm(mat[x, x, drop = FALSE], type = "F")
#   y_norm <- norm(mat[y, y, drop = FALSE], type = "F")
#   half_dist_norm <- norm(mat[x, y, drop = FALSE], type = "F")
#   
#   1 - (x_norm + y_norm - 2*half_dist_norm) / full_norm
#   
#   # sqrt(max(half_dist_norm^2 / (x_norm * y_norm), 0))
# }

# inner_prod_dist <- function(x, y, mat) {
#   lx <- length(x)
#   ly <- length(y)
#   
#   sqrt(max(lx*ly - sum(La.svd(mat[x, y], nu = 0, nv = 0)$d), 0))
# }

# Slightly different, but way faster, Gives the right answer as well
# inner_prod <- function(x, y, mat) {
#   d <- diag(mat)
#   lx <- length(x)
#   ly <- length(y)
#   # lx <- sum(d[x])
#   # ly <- sum(d[y])
# 
#   norm(mat[x, y, drop = FALSE], type = "F") / sqrt(lx * ly)
# }
# 
# Lt <- t(compute_basis(samp_ortho_std[-1, -1]))
# compute_cancor <- function(x, y) {
#   A <- Lt[, x, drop = FALSE]
#   B <- Lt[, y, drop = FALSE]
#   
#   # Orthonormalize the input matrices (if not already orthonormal)
#   Q_A <- qr.Q(qr(A)) # Orthonormal basis for subspace A
#   Q_B <- qr.Q(qr(B)) # Orthonormal basis for subspace B
#   
#   # Compute the singular values of Q_A^T Q_B
#   cors <- svd(crossprod(Q_A, Q_B), nu = 0, nv = 0)$d
#   
#   return(cors)
# }

# Grassman based measures (this one is Fubini-Study)
# inner_prod3 <- function(x, y, mat) {
#   cors <- compute_cancor(x, y)
#   1 - acos(min(prod(cors), 1)) / (pi / 2)
# }

# inner_prod <- function(x, y, mat) {
#   tol <- sqrt(.Machine$double.eps)
#   mat_xy <- partial_mat(mat, x, y)
#   mat_yx <- partial_mat(mat, y, x)
# 
#   if (min(c(diag(mat_xy), diag(mat_yx))) < tol) {
#     1
#   } else {
#     median(mat[x, y]^2)
#     # max(mean(mat_xy^2) / mean(mat[x, x]^2),
#     #     mean(mat_yx^2) / mean(mat[y, y]^2))
#   }
# }

# inner_prod <- function(x, y, mat) {
#   X <- L[x, , drop = FALSE]
#   Ux <- svd(X)$v
#   
#   Y <- L[y, , drop = FALSE]
#   Uy <- svd(Y)$v
#   r <- min(rankM(X), rankM(Y))
#   sqrt(mean(svd(crossprod(Ux, Uy))$d[seq_len(r)]))
# }

# inner_prod <- function(x, y, mat) {
#   if (length(setdiff(union(x, y), intersect(x, y))) == 0) {
#     return(1)
#   }
#   mat_xy <- mat[x, y] %*% MASS::ginv(mat[y, y]) %*% mat[y, x]
#   mat_yx <- mat[y, x] %*% MASS::ginv(mat[x, x]) %*% mat[x, y]
#   (sum(diag(mat_xy)) / length(x) + sum(diag(mat_yx)) / length(y)) / 2
#   
# }


# inner_prod <- function(x, y, mat) {
#   xy <- union(x, y)
#   mat_x <- mat[xy, x] %*% MASS::ginv(mat[x, x]) %*% mat[x, xy]
#   mat_y <- mat[xy, y] %*% MASS::ginv(mat[y, y]) %*% mat[y, xy]
#   mat_xy <- mat[xy, xy]
# 
#   sqrt(max(1 - mean(diag(2*mat_xy - mat_x - mat_y)), 0))
# }

# X <- t(L)
# inner_prod <- function(x, y, mat) {
#   X_x_y <- X[, y, drop = FALSE] %*% MASS::ginv(crossprod(X[, y, drop = FALSE])) %*% crossprod(X[, y, drop = FALSE], X[, x, drop = FALSE])
#   X_y_x <- X[, x, drop = FALSE] %*% MASS::ginv(crossprod(X[, x, drop = FALSE])) %*% crossprod(X[, x, drop = FALSE], X[, y, drop = FALSE])
#   
#   dimnames(X_x_y) <- list(rownames(X), x)
#   dimnames(X_y_x) <- list(rownames(X), y)
#   
#   A <- named_svd(X_x_y)$u
#   B <- named_svd(X_y_x)$u
#   
#   if (ncol(A) < ncol(B)) {
#     tmp <- A
#     A <- B
#     B <- tmp
#   }
#   B <- B - A %*% (t(A) %*% B)
#   theta <- asin(min(1, max(svd(B)$d)))
#   return(1 - theta / (pi/2))
# }

# inner_prod <- function(x, y, mat) {
#   X_x_y <- X[, y] %*% MASS::ginv(crossprod(X[, y])) %*% crossprod(X[, y], X[, x])
#   X_y_x <- X[, x] %*% MASS::ginv(crossprod(X[, x])) %*% crossprod(X[, x], X[, y])
#   max(diag(crossprod(cbind(X_x_y, X_y_x)))) / length(c(x, y))
# }

clust_cov <- function(clusts, mat) {
  cov_mat <- matrix(NA, nrow = length(clusts), ncol = length(clusts))
  for (i in seq_along(clusts)) {
    for (j in seq_len(i)) {
      cov_mat[i, j] <- cov_mat[j, i] <- inner_prod(clusts[[i]], clusts[[j]], mat = mat)
    }
  }
  return(cov_mat)
}

res_mat <- function(clusters, mat) {
  res_mat <- matrix(0, nrow = length(clusters), ncol = length(clusters))
  for (i in seq_along(clusters)) {
    chosen <- unlist(clusters[i])
    remaining <- unlist(clusters[-i])
    part_mat <- mat[c("resids", chosen), c("resids", chosen)] -
      mat[c("resids", chosen), remaining] %*%
      MASS::ginv(mat[remaining, remaining]) %*%
      mat[remaining, c("resids", chosen)]
    
    res_mat[i, i] <- part_mat[1, -1] %*% MASS::ginv(part_mat[-1, -1]) %*% part_mat[-1, 1]
    
  }
  for (i in seq_along(clusters)) {
    for (j in seq_len(i - 1)) {
      chosen <- unlist(clusters[c(i, j)])
      remaining <- unlist(clusters[-c(i, j)])
      part_mat <- mat[c("resids", chosen), c("resids", chosen)] -
        mat[c("resids", chosen), remaining] %*%
        MASS::ginv(mat[remaining, remaining]) %*%
        mat[remaining, c("resids", chosen)]
      total_ij <- part_mat[1, -1] %*% MASS::ginv(part_mat[-1, -1]) %*% part_mat[-1, 1]
      
      res_mat[i, j] <- res_mat[j, i] <- (total_ij - res_mat[i, i] - res_mat[j, j]) / 2
    }
  }
  return(res_mat)
}

safe_names <- function(x) {
  if (length(x) > 0) {
    names(x)
  } else {
    character(0)
  }
}

find_spanner_candidates <- function(clusters, mat, thresh = 0.05, drop = TRUE) {
  remainers <- map(seq_along(clusters), ~ unlist(clusters[-.x]))
  parted <- map2(clusters, remainers, partial_mat, mat = mat) %>%
    map(diag)
  
  candidates <- parted %>%
    map(~ .x < thresh) %>%
    map(which) %>%
    map(safe_names)
  
  if (drop) {
    if (all(map2_lgl(clusters, candidates, setequal))) {
      stop("both equal...")
    }
    # only remove the ones from the non-equal cluster
    candidates <- map2(clusters, candidates, setdiff) %>%
      map2_lgl(clusters, ~ length(.x) < 2 & length(.y) > 1) %>%
      modify_if(candidates, ., ~ character(0))
  }
  return(candidates)
}


refine <- function(clusters, mat, thresh = 0.05) {
  spans <- unlist(find_spanner_candidates(clusters, mat, thresh, drop = FALSE))
  split_clusters <- c(map(clusters, setdiff, spans), as.list(spans)) %>%
    compact()
  
  return(split_clusters)
}

get_clusters <- function(fit, mat, k) {
  o <- unlist(fit[[5]])
  state <- fit %>%
    pluck(2) %>%
    get_state(k, fit[[5]], mat)
  
  clusters <- state %>%
    pluck(1)
  
  clusters <- clusters[order(map_dbl(clusters, \(x) min(match(x, o))))] %>%
    map(~ .x[order(match(.x, o))])
  
  return(clusters)
}

refined_state <- function(mat, clusts, k, n_reps = 5) {
  for (i in seq_len(n_reps)) {
    fit <- pseudo_frobenius(mat, clusts)
    state <- get_state(fit[[2]], k, clusts, mat)
    clusts <- refine(state[[1]], mat)
  }
  fit <- pseudo_frobenius(mat, clusts)
  state <- get_state(fit[[2]], k, clusts, mat)
  return(state)
}

term_pseudo_frobenius <- function(mat, clusts, thresh = 0.05) {
  N <- length(clusts) - 2
  merge <- matrix(NA, nrow = N, ncol = 2)
  crit <- rep(NA, N)
  cov_mat <- clust_cov(clusts, mat = mat)
  
  M <- N + 2 - rankM(mat)
  for (i in seq_len(N)) {
    cor_mat <- cov2cor(cov_mat)
    crit[i] <- max(cor_mat[lower.tri(cor_mat)])
    
    merge[i, ] <- which(cor_mat - diag(diag(cor_mat)) == crit[i], arr.ind = TRUE)[1, ]
    dex <- merge[i, 2]
    clusts[[dex]] <- unlist(clusts[merge[i, ]])
    clusts <- clusts[-merge[i, 1]]
    
    if (i >= M) {
      candidates <- find_spanner_candidates(clusts, mat, thresh, drop = FALSE)
      if (all(map_dbl(candidates, length) == 0)) {
        break
      }
    }
    
    cov_mat <- cov_mat[-merge[i, 1], -merge[i, 1], drop = FALSE]
    cov_mat[dex, ] <- cov_mat[, dex] <- map_dbl(clusts, inner_prod, clusts[[dex]], mat = mat)
  }
  return(list(clusts, merge, crit))
}

# Ordering such that always remove the higher index cluster
# pseudo_frobenius <- function(mat, clusts) {
#   N <- length(clusts) - 2
#   merge <- matrix(NA, nrow = N, ncol = 2)
#   crit <- rep(NA, N)
#   cov_mat <- clust_cov(clusts, mat = mat)
#   for (i in seq_len(N)) {
#     # cor_mat <- cov_mat
#     
#     raw_crit <- max(cov_mat[lower.tri(cov_mat)])
#     
#     # add penalty for cluster size
#     # idk this does really weird things to 1, 1 cluster joins
#     lens <- map_dbl(clusts, length)
#     size_mat <- 0 # 1 / sqrt(tcrossprod(lens))
#     cor_mat <- cov_mat + size_mat * raw_crit
#     
#     crit[i] <- max(cor_mat[lower.tri(cor_mat)])
#     
#     merge[i, ] <- which(cor_mat - diag(diag(cor_mat)) == crit[i], arr.ind = TRUE)[1, ]
#     dex <- merge[i, 2]
#     clusts[[dex]] <- unlist(clusts[rev(merge[i, ])])
#     clusts <- clusts[-merge[i, 1]]
#     
#     cov_mat <- cov_mat[-merge[i, 1], -merge[i, 1], drop = FALSE]
#     cov_mat[dex, ] <- cov_mat[, dex] <- map_dbl(clusts, inner_prod, clusts[[dex]], mat = mat)
#   }
#   return(list(clusts, merge, crit))
# }

# Ordering such that always remove the cluster farther away from middle
pseudo_frobenius <- function(mat, clusts) {
  orig_clusts <- clusts

  N <- length(clusts) - 1
  merge <- matrix(NA, nrow = N, ncol = 2)

  # mimic hclust behavior
  clust_dex <- -seq_along(clusts)
  clust_merge <- matrix(NA, nrow = N, ncol = 2)

  crit <- rep(NA, N)
  cov_mat <- clust_cov(clusts, mat = mat)
  for (i in seq_len(N)) {
    # set diag and non rank-def clusters to NA
    diag(cov_mat) <- NA

    # Find highest correlation
    crit[i] <- max(cov_mat, na.rm = TRUE)
    best_dex <- arrayInd(which.max(cov_mat), dim(cov_mat))

    # Order merging such that we merge toward the "center" cluster
    n_clusts <- length(clusts)
    dex_dist <- best_dex - (n_clusts + 1) / 2
    dex_order <- order(dex_dist^2, decreasing = FALSE)
    ordered_best_dex <- best_dex[1, ]

    clust_merge[i, ] <- clust_dex[ordered_best_dex]
    merge[i, ] <- ordered_best_dex

    # Order merging such that we merge toward the "center" cluster
    best_dex <- best_dex[order(dex_dist, decreasing = FALSE)]
    clusts[[ordered_best_dex[2]]] <- unlist(clusts[best_dex])

    cov_mat[-ordered_best_dex[1], ordered_best_dex[2]] <-
      map_dbl(clusts[-ordered_best_dex[1]], inner_prod, clusts[[ordered_best_dex[2]]], mat = mat)
    cov_mat[ordered_best_dex[2], -ordered_best_dex[1]] <- cov_mat[-ordered_best_dex[1], ordered_best_dex[2]]

    clusts <- clusts[-ordered_best_dex[1]]
    cov_mat <- cov_mat[-ordered_best_dex[1], -ordered_best_dex[1], drop = FALSE]

    clust_dex[ordered_best_dex[2]] <- i
    clust_dex <- clust_dex[-ordered_best_dex[1]]
  }

  crit[is.infinite(crit)] <- 0

  return(list(clusts, merge, crit, clust_merge, orig_clusts))
}

pseudo_frobenius_freeze <- function(mat, clusts, skips = character(0)) {
  orig_clusts <- clusts
  
  N <- length(clusts) - 1
  merge <- matrix(NA, nrow = N, ncol = 2)
  
  # mimic hclust behavior
  clust_dex <- -seq_along(clusts)
  clust_merge <- matrix(NA, nrow = N, ncol = 2)
  
  crit <- rep(NA, N)
  cov_mat <- clust_cov(clusts, mat = mat)
  for (i in seq_len(N)) {
    skip_dex <- which(map_lgl(clusts, \(x) any(x %in% skips)))
    if (length(skip_dex) < length(clusts)) {
      cov_mat[skip_dex, skip_dex] <- NA
    }
    
    # set diag and non rank-def clusters to NA
    diag(cov_mat) <- NA

    # Find highest correlation
    crit[i] <- max(cov_mat, na.rm = TRUE)
    best_dex <- arrayInd(which.max(cov_mat), dim(cov_mat))
    
    # Order merging such that we merge toward the "center" cluster
    n_clusts <- length(clusts)
    dex_dist <- best_dex - (n_clusts + 1) / 2
    dex_order <- order(dex_dist^2, decreasing = FALSE)
    ordered_best_dex <- best_dex[dex_order]
    
    clust_merge[i, ] <- clust_dex[ordered_best_dex]
    merge[i, ] <- ordered_best_dex
    
    # Order merging such that we merge toward the "center" cluster
    best_dex <- best_dex[order(dex_dist, decreasing = FALSE)]
    clusts[[ordered_best_dex[2]]] <- unlist(clusts[best_dex])
    
    cov_mat[-ordered_best_dex[1], ordered_best_dex[2]] <-
      map_dbl(clusts[-ordered_best_dex[1]], inner_prod, clusts[[ordered_best_dex[2]]], mat = mat)
    cov_mat[ordered_best_dex[2], -ordered_best_dex[1]] <- cov_mat[-ordered_best_dex[1], ordered_best_dex[2]]
    
    clusts <- clusts[-ordered_best_dex[1]]
    cov_mat <- cov_mat[-ordered_best_dex[1], -ordered_best_dex[1], drop = FALSE]
    
    clust_dex[ordered_best_dex[2]] <- i
    clust_dex <- clust_dex[-ordered_best_dex[1]]
  }
  
  crit[is.infinite(crit)] <- 0
  
  return(list(clusts, merge, crit, clust_merge, orig_clusts))
}

build_queue <- function(cor_mat, merge_thresh = 0.5*max(cor_mat[lower.tri(cor_mat)]),
                        tol = sqrt(.Machine$double.eps)) {
  cor_dims <- dim(cor_mat)
  cor_dexes <- seq_len(prod(cor_dims))
  I <- arrayInd(cor_dexes, cor_dims, useNames = TRUE)
  I <- I[I[, 1] < I[, 2], , drop = FALSE]
  
  cor_table <- cbind(I, cor_mat[I])
  cor_table <- cor_table[cor_table[, 3] >= tol, , drop = FALSE]
  
  cor_order <- order(cor_table[, 3], decreasing = TRUE)
  cor_table <- cor_table[cor_order, , drop = FALSE]
  
  queue_table <- matrix(NA, nrow = 0, ncol = 2)
  while(length(cor_table) > 0) {
    focus <- cor_table[1, 1:2]
    queue_table <- rbind(queue_table, focus)

    focus_mask <- cor_table[, 1] %in% focus | cor_table[, 2] %in% focus
    cor_mask <- focus_mask & cor_table[, 3] > tol
    
    bad_indices <- as.vector(cor_table[cor_mask, 1:2])
    drop_mask <- apply(cor_table, 1, function(x) any(x %in% bad_indices))
    
    if (sum(!drop_mask) == 0) {
      break
    }
    
    cor_table <- cor_table[!drop_mask, , drop = FALSE]
  }
  
  return(queue_table[order(queue_table[, 2], decreasing = TRUE), , drop = FALSE])
}


pseudo_frobenius_multi <- function(mat, clusts) {
  N <- length(clusts) - 1
  merge <- matrix(NA, nrow = N, ncol = 3)
  crit <- rep(NA, N)
  cor_mat <- clust_cov(clusts, mat = mat)
  queue_table <- build_queue(cor_mat)
  
  k <- 1
  for (i in seq_len(N)) {
    if (nrow(queue_table) == 0) {
      queue_table <- build_queue(cor_mat)
      k <- k + 1
    }
    
    row <- queue_table[1, 1]
    col <- queue_table[1, 2]
    queue_table <- queue_table[-1, , drop = FALSE]
    merge[i, ] <- c(row, col, k)
    
    crit[i] <- cor_mat[row, col]
    
    clusts[[row]] <- unlist(clusts[c(row, col)])
    cor_mat[row, ] <- cor_mat[, row] <- map_dbl(clusts, inner_prod, clusts[[row]], mat = mat)
    
    clusts <- clusts[-col]
    cor_mat <- cor_mat[-col, -col, drop = FALSE]
  }
  return(list(clusts, merge, crit))
}

pseudo_frobenius_fixed <- function(mat, clusts, remainder) {
  N <- length(remainder)
  merge <- matrix(NA, nrow = N, ncol = 2)
  crit <- rep(NA, N)
  cov_mat <- clust_cov(c(clusts, as.list(remainder)), mat = mat)
  for (i in seq_len(N)) {
    cor_mat <- cov2cor(cov_mat)
    crit[i] <- max(cor_mat[seq_along(clusts), -seq_along(clusts)])
    
    merge[i, ] <- which(cor_mat[seq_along(clusts), -seq_along(clusts), drop = FALSE] == crit[i], arr.ind = TRUE)[1, ]
    clust_dex <- merge[i, 1]
    rem_dex <- merge[i, 2]
    clusts[[clust_dex]] <- c(clusts[[clust_dex]], remainder[rem_dex])
    remainder <- remainder[-rem_dex]
    
    cov_mat <- cov_mat[-(rem_dex + length(clusts)), -(rem_dex + length(clusts)), drop = FALSE]
    cov_mat[clust_dex, ] <- cov_mat[, clust_dex] <- map_dbl(c(clusts, as.list(remainder)), inner_prod, clusts[[clust_dex]], mat = mat)
  }
  return(list(clusts, merge, crit))
}

pseudo_frobenius_rep <- function(mat, clusts, k, n_reps = 5, factor = 2, cov_mat = NULL) {
  if (is.null(cov_mat)) {
    cov_mat <- clust_cov(clusts, mat = mat)
  }
  
  while(length(clusts) > k) {
    cor_mat <- cov2cor(cov_mat)
    crit <- max(cor_mat[lower.tri(cor_mat)])
    
    merge <- which(cor_mat - diag(diag(cor_mat)) == crit, arr.ind = TRUE)[1, ]
    dex <- merge[2]
    clusts[[dex]] <- unlist(clusts[rev(merge)])
    clusts <- clusts[-merge[1]]
    
    cov_mat <- cov_mat[-merge[1], -merge[1], drop = FALSE]
    cov_mat[dex, ] <- cov_mat[, dex] <- map_dbl(clusts, inner_prod, clusts[[dex]], mat = mat)
    
    if (length(clusts) < factor*k) {
      for (i in seq_len(n_reps)) {
        split_clusts <- refine(clusts, mat, thresh = 0.15)
        clusts <- pseudo_frobenius_rep(mat, split_clusts,
                                       k = length(clusts),
                                       n_reps = 0,
                                       factor = 0)
      }
    }
  }
  return(clusts)
}

find_max_excl <- function(cor_mat, thresh = 0.05) {
  off_diag_mat <- cor_mat - diag(diag(cor_mat))
  diffs <- apply(off_diag_mat[, -(1:2)], 2, function(x) -diff(sort(x^2, decreasing = TRUE))[1])
  best <- which.max(diffs) + 2
  match <- which.max(off_diag_mat[, best]^2)
  if (max(diffs) < thresh) {
    return(NULL)
  } else {
    return(sort(c(match, best), decreasing = TRUE))
  }
}

pseudo_frobenius_excl <- function(mat, clusts, error, thresh = 0.05) {
  N <- length(error)
  merge <- matrix(NA, nrow = N, ncol = 2)
  crit <- rep(NA, N)
  clusts <- c(clusts, as.list(error))
  cov_mat <- clust_cov(clusts, mat = mat)
  for (i in seq_len(N)) {
    cor_mat <- cov2cor(cov_mat)
    out <- find_max_excl(cor_mat, thresh)
    if (is.null(out)) break
    merge[i, ] <- out
    crit[i] <- cor_mat[out[1], out[2]]
    
    dex <- merge[i, 2]
    clusts[[dex]] <- unlist(clusts[rev(merge[i, ])])
    clusts <- clusts[-merge[i, 1]]
    
    cov_mat <- cov_mat[-merge[i, 1], -merge[i, 1], drop = FALSE]
    cov_mat[dex, ] <- cov_mat[, dex] <- map_dbl(clusts, inner_prod, clusts[[dex]], mat = mat)
  }
  return(list(clusts, merge, crit))
}

get_state <- function(merge, k, params, matrix = mat) {
  clusts <- as.list(params)
  N <- length(clusts) - k
  for (i in seq_len(N)) {
    clusts[[merge[i, 2]]] <- unlist(clusts[rev(merge[i, ])])
    clusts <- clusts[-merge[i, 1]]
  }
  cov_mat <- clust_cov(clusts, mat = matrix)
  cor_mat <- cov2cor(cov_mat)
  return(list(clusts, cov_mat, cor_mat))
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



# plot(crit)
# map(clusts, ~ c("resids", .x)) %>%
#   map(~ info_mat[.x, .x]) %>%
#   map_dbl(~ .x[1, -1] %*% MASS::ginv(.x[-1, -1]) %*% .x[-1, 1])
# 
# USV <- named_svd(info_mat)
# L <- (USV$u %*% sqrt(USV$d))[, diag(USV$d) > sqrt(.Machine$double.eps)]
# mod_clusts <- clusts
# for (i in seq_along(clusts)) {
#   rotmat <- varimax(L[unlist(clusts[-i]), ])$rotmat
#   R <- L %*% rotmat
#   crit <- cov2cor(tcrossprod(R[, -order(R["resids", ]^2, decreasing = TRUE)[seq_along(clusts[-1])][-1]]))[clusts[[i]], "resids"]^2
#   mod_clusts[[i]] <- names(which(crit > 0.05))
# }
# 
# u_clusts <- as.list(unlist(mod_clusts))
# N <- length(u_clusts) - 1
# merge <- matrix(NA, nrow = N, ncol = 2)
# crit <- rep(NA, N)
# for (i in seq_len(N)) {
#   cov_mat <- clust_cov(u_clusts, mat = info_mat)
#   cor_mat <- cov2cor(cov_mat)
#   crit[i] <- max(cor_mat[lower.tri(cor_mat)])
#   if (crit[i] < sqrt(0.05)) break 
#   merge[i, ] <- which(cor_mat - diag(diag(cor_mat)) == crit[i], arr.ind = TRUE)[1, ]
#   u_clusts[[merge[i, 2]]] <- unlist(u_clusts[merge[i, ]])
#   u_clusts <- u_clusts[-merge[i, 1]]
# }
# 
# plot_signed_communalities(info_mat[unlist(u_clusts), unlist(u_clusts)],
#                           info_mat[1, unlist(u_clusts)], 1, thresh = 0)
# 
# max_names <- map_chr(clusts, ~ .x[which.max(info_mat["resids", .x]^2)])
# cor_mat <- cov2cor(clust_cov(clusts, mat = mat))
# dimnames(cor_mat) <- rerun(2, max_names)
# 
# plot_signed_communalities(cor_mat, info_mat["resids", max_names], 1, sort_chisq = TRUE)
# ---------------------------------
# Associate (potentially overlapping) subsets of constraint vectors to each of
# the rotated singular vectors


# ---------------------------------
# Iteratively select subset of vectors associated with the residual
# Take one: just keep repeating SVD, rotation, truncation on residual projection
