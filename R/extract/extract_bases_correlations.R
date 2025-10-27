# Different bases ----------------------------------------------------------
compute_basis <- function(X) {
  USV <- named_svd(X)
  (USV$u %*% sqrt(USV$d))[, diag(USV$d) > sqrt(.Machine$double.eps), drop = FALSE]
}

compute_rotmat <- function(matrix) {
  USV <- named_svd(matrix)
  pos_mask <- diag(USV$d) > sqrt(.Machine$double.eps)
  
  rotmat <- tcrossprod(USV$u[, pos_mask], USV$v[, pos_mask])
  
  return(rotmat)
}

closest_orthogonal_basis <- function(matrix) {
  USV <- named_svd(matrix)
  pos_mask <- diag(USV$d) > sqrt(.Machine$double.eps)
  
  S <- sqrt(USV$d)[pos_mask, pos_mask, drop = FALSE]
  V <- USV$u[, pos_mask]
  V_std <- V %*% solve(S, t(V))
  
  return(V_std)
}

sparse_basis <- function(matrix) {
  USV <- named_svd(matrix)
  pos_mask <- diag(USV$d) > sqrt(.Machine$double.eps)
  
  S <- sqrt(USV$d)[pos_mask, pos_mask, drop = FALSE]
  V <- USV$u[, pos_mask, drop = FALSE]
  
  if (ncol(V) > 1) {
    rotmat <- varimax(V %*% S, eps = 1e-8)$rotmat
    mag <- colSums((V %*% S %*% rotmat)^2)
    
    R_std <- (V %*% solve(S) %*% rotmat)[, order(mag, decreasing = TRUE)]
  } else {
    R_std <- V %*% solve(S)
  }
  colnames(R_std) <- paste0("R", seq_len(ncol(R_std)))
  
  return(R_std)
}

# Compute dual vectors -----------------------------------------------------

# correlation with moment vectors
cor_params_moments <- function(mats) {
  jac <- mats$jac[[1]]
  W <- mats$W[[1]]
  
  vecs <- cbind(jac, diag(ncol(W)))
  colnames(vecs) <- c(colnames(jac), colnames(W))
  vecs_cov <- crossprod(vecs, W) %*% vecs
  
  mask <- diag(vecs_cov) > sqrt(.Machine$double.eps)
  vecs_cor <- cov2cor(vecs_cov)
  vecs_cor[!mask, ] <- vecs_cor[, !mask] <- 0
  
  vecs_cor <- vecs_cor[seq_len(ncol(jac)), -seq_len(ncol(jac))]
  
  return(vecs_cor)
}

# correlation with moment vectors
cor_params_moments_keepall <- function(mats) {
  jac <- mats$jac[[1]]
  W <- mats$W[[1]]
  
  vecs <- cbind(jac, diag(ncol(W)))
  colnames(vecs) <- c(colnames(jac), colnames(W))
  vecs_cov <- crossprod(vecs, W) %*% vecs
  
  mask <- diag(vecs_cov) > sqrt(.Machine$double.eps)
  vecs_cor <- cov2cor(vecs_cov)
  vecs_cor[!mask, ] <- vecs_cor[, !mask] <- 0
  
  return(vecs_cor)
}

# full subspace closest orthogonal basis of the moments
cor_params_moments_COB <- function(mats) {
  jac <- mats$jac[[1]]
  W <- mats$W[[1]]
  
  V_std <- closest_orthogonal_basis(cov2cor(W))
  colnames(V_std) <- paste0("V", seq_len(ncol(V_std)))
  
  vecs <- cbind(jac, diag(1 / sqrt(diag(W))) %*% V_std)
  
  vecs_cor <- cov2cor(crossprod(vecs, W) %*% vecs)[colnames(jac), colnames(V_std)]
  colnames(vecs_cor) <- paste0("V*", colnames(W))
  
  return(vecs_cor)
}

cor_params_moments_sparse <- function(mats) {
  jac <- mats$jac[[1]]
  W <- mats$W[[1]]
  
  R_std <- sparse_basis(cov2cor(W))
  vecs <- cbind(jac, diag(1 / sqrt(diag(W))) %*% R_std)
  vecs_cor <- cov2cor(crossprod(vecs, W) %*% vecs)[colnames(jac), colnames(R_std)]
  
  return(vecs_cor)
}

cor_params_all_sparse <- function(mats) {
  info <- mats$info
  nvcov <- mats$nvcov
  pars <- colnames(nvcov)
  
  info_std <- sweep(info, 1, sqrt(diag(info)), '/') %>%
    sweep(2, sqrt(diag(info)), '/')
  
  R_std <- sparse_basis(info_std)
  vecs_cor <- info_std %*% R_std
  
  return(vecs_cor)
}

cor_params_all_sparse_free <- function(mats) {
  info <- mats$info
  nvcov <- mats$nvcov
  pars <- colnames(nvcov)
  
  info_free <- info[, pars] %*% nvcov %*% info[pars, ]
  info_free_std <- sweep(info_free, 1, sqrt(diag(info)), '/') %>%
    sweep(2, sqrt(diag(info)), '/')
  
  R_std <- sparse_basis(info_free_std)
  vecs_cor <- info_free_std %*% R_std
  
  return(vecs_cor)
}

cor_params_free_sparse <- function(mats) {
  info <- mats$info
  free <- mats$names$free
  info_std <- cov2cor(info)
  info_std_free <- info_std[free, free]
  
  R_std <- sparse_basis(info_std_free)
  vecs_cor <- info_std[, free] %*% R_std
  
  return(vecs_cor)
}

cor_params_free_COB <- function(mats) {
  info <- mats$info
  free <- mats$names$free
  
  info_std <- cov2cor(info)
  free_info_std <- info_std[free, free]
  
  V_std <- closest_orthogonal_basis(free_info_std)
  colnames(V_std) <- paste0("V*", free)
  vecs_cor <- info_std[, free] %*% V_std
  
  return(vecs_cor)
}




# Only vectors in the free param subspace ----------------------------------

cor_params_moments_free <- function(mats) {
  jac <- mats$jac[[1]]
  W <- mats$W[[1]]
  nvcov <- mats$nvcov
  pars <- colnames(nvcov)
  
  W_free <- W %*% jac[, pars] %*% tcrossprod(nvcov, jac[, pars]) %*% W
  
  mats_new <- mats
  mats_new$W <- list(W_free)
  vecs_cor <- cor_params_moments(mats_new)
  
  return(vecs_cor)
}

cor_params_moments_COB_free <- function(mats) {
  jac <- mats$jac[[1]]
  W <- mats$W[[1]]
  nvcov <- mats$nvcov
  pars <- colnames(nvcov)
  
  W_free <- W %*% jac[, pars] %*% tcrossprod(nvcov, jac[, pars]) %*% W
  mask <- diag(W_free) > sqrt(.Machine$double.eps)
  
  mats_new <- mats
  mats_new$W <- list(W_free[mask, mask])
  mats_new$jac <- list(jac[mask, ])
  vecs_cor <- cor_params_moments_COB(mats_new)
  
  return(vecs_cor)
}



# Only vectors in the constrained subspace ---------------------------------

cor_params_moments_ortho <- function(mats) {
  jac <- mats$jac[[1]]
  W <- mats$W[[1]]
  nvcov <- mats$nvcov
  pars <- colnames(nvcov)
  
  W_free <- W %*% jac[, pars] %*% tcrossprod(nvcov, jac[, pars]) %*% W
  W_ortho <- W - W_free
  mask <- diag(W_ortho) > sqrt(.Machine$double.eps)
  
  mats_new <- mats
  mats_new$W <- list(W_ortho[mask, mask])
  mats_new$jac <- list(jac[mask, ])
  vecs_cor <- cor_params_moments(mats_new)
  
  return(vecs_cor)
}

compute_cancor <- function(matrix, xdex, ydex) {
  xydex <- c(xdex, ydex)
  USV <- named_svd(matrix[xydex, xydex])
  mask <- diag(USV$d) > sqrt(.Machine$double.eps)
  
  L <- tcrossprod(sqrt(USV$d), USV$u)[mask, ]
  X <- L[, xdex]
  Y <- L[, ydex]
  
  rho <- sqrt(sum(cancor(X, Y, FALSE, FALSE)$cor^2))
  
  return(rho)
}

find_max_cancor <- function(clusters, matrix) {
  idx <- seq_along(clusters)
  best_dex <- c(1, 2)
  best_cor <- 0
  
  for (i in idx) {
    xdex <- clusters[[i]]
    for (j in idx[seq_len(i - 1)]) {
      ydex <- clusters[[j]]
      cor <- compute_cancor(matrix, xdex, ydex)
      if (cor > best_cor) {
        best_dex <- c(j, i)
        best_cor <- cor
      }
    }
  }
  
  return(list(dex = best_dex, cor = best_cor))
}

glom_clusters <- function(clusters, matrix) {
  cancor_out <- find_max_cancor(clusters, matrix)
  
  glom_dex <- cancor_out$dex
  clusters[[glom_dex[1]]] <- unlist(clusters[glom_dex])
  clusters <- clusters[-glom_dex[2]]
  
  return(list(clusters = clusters, cor = cancor_out$cor))
}

cancor_clustering <- function(matrix) {
  clusters <- as.list(colnames(matrix))
  mem_matrix <- matrix(seq_len(ncol(matrix)), nrow = ncol(matrix), ncol = 1)
  mem_cors <- 1
  rownames(mem_matrix) <- colnames(matrix)
  while (length(clusters) > 1) {
    glom_out <- glom_clusters(clusters, matrix)
    clusters <- glom_out$clusters
    glom_cor <- glom_out$cor
    
    mem_cors <- c(mem_cors, glom_cor)
    
    clust_map <- imap(clusters, ~ set_names(rep_along(.x, .y), .x))
    mem_matrix <- cbind(mem_matrix, unlist(clust_map)[rownames(mem_matrix)])
  }
  return(list(membership = mem_matrix, cor = mem_cors))
}
