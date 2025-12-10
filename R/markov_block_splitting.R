# CITATION: Fritzsche, David & Mehrmann, Volker & Szyld, Daniel & Virnik, Elena.
#           (2008). An SVD approach to identifying metastable states of Markov
#           chains. Electronic Transactions on Numerical Analysis. 29. 46-69. 
# LINK TO PDF: https://opus4.kobv.de/opus4-matheon/frontdoor/deliver/index/docId/375/file/4034_meta_preprint.pdf
# IDK IT DOESN'T WORK WITH SAMP_SELECT


# helper functions
rotate_by_block <- function(R, blocks) {
  R_rot <- map(
    blocks, \(dex) 
    clustered_varimax(R[, dex, drop = FALSE],
                      centering = FALSE, normalize = FALSE,
                      eps = eps)
    ) %>%
    map("loadings") %>%
    do.call(what = cbind)
  
  R_rot[, order(unlist(blocks)), drop = FALSE]
}

# rotate_by_block <- function(R, blocks) {
#   R_rot <- map(
#     blocks, \(dex) {
#       Z <- as_unit_vec(R[, dex, drop = FALSE]^2) %*% compute_basis(cov2cor(crossprod(R[, dex, drop = FALSE]^2)))
#       Q <- Z^2 / rowSums(Z^2)
#       Qmax(R[, dex, drop = FALSE], Q)
#     }
#   ) %>%
#     map("loadings") %>%
#     do.call(what = cbind)
#   
#   R_rot[, order(unlist(blocks)), drop = FALSE]
# }


is_mutually_disjoint <- function(set_family) {
  intersections <- map(seq_along(set_family), \(i) unlist(set_family[-i])) %>%
    map2(set_family, intersect) %>%
    unlist()
  
  length(intersections) == 0
}

# norm functions
v_norm <- function(x, v) {
  sum(rowSums(x) * v) / sum(v)
}

# inner product functions
v_prod <- function(x, sk, sl, v = rep_along(c(sk, sl), 1)) {
  v_norm(x[sk, sl, drop = FALSE], v)
}


block_ident_step <- function(x, dex = seq_len(ncol(x)), v = rep(1, ncol(x)), thresh = 0.5) {
  # compute second left singular vector
  singular2 <- eigen(x, symmetric = FALSE)$vectors[, 2]
  
  # sort the vector and permute B index
  perm <- order(singular2)
  x_perm <- x[perm, perm, drop = FALSE]
  perm_dex <- dex[perm]
  
  
  # identify blocks based on positive and negative values in vector
  split_mask <- singular2[perm] >= 0
  
  # compute norm
  ind1 <- which(split_mask)
  ind2 <- which(!split_mask)
  
  diag_norms <- map_dbl(list(ind1, ind2), \(i) v_prod(x_perm, i, i, v[i]))
  
  # if both surpass the thresh then split, otherwise x cannot be reduced further
  if (all(diag_norms > thresh)) {
    split_dex <- map(list(ind1, ind2), \(i) perm_dex[i])
  } else {
    split_dex <- list(perm_dex)
  }
  
  list(dex = split_dex, norm = diag_norms)
}

R_init <- clustered_varimax(L, centering = FALSE, eps = eps)$loadings
# U <- L %*% project_onto_Stiefel(crossprod(L, MASS::ginv(mat))) # Not really sure how to think about this method... but it works
# R_init <- U
R <- R_init
# Note the B's here are stochastic (rowSums = 1)
B <- crossprod(R^2) / colSums(R^2) # Dimension-based clustering
# B <- matpow(B, 3)
# B <- tcrossprod(as_unit_vec(R)^2, R^2)
# B <- U^2 # If using project_onto_Stiefel method
p <- eigen(t(B), symmetric = FALSE)$vector[, 1]
# p <- rep_along(p, 1)
pi <- p / sum(p)
thresh <- 0.825
# target_nblocks <- 10

blocks <- list(seq_len(ncol(B)))
norms <- 1
skip_ind <- rep_along(blocks, FALSE)

# Stop if number of blocks did not change. Irreducible.
while (!all(skip_ind)) {
  temp_blocks <- list()
  temp_norms <- numeric(0)
  temp_skip_ind <- logical(0)
  for (i in seq_along(blocks)) {
    dex <- blocks[[i]]
    
    # Skip work if already determined irreducible
    if (skip_ind[i]) {
      temp_skip_ind <- c(temp_skip_ind, TRUE)
      temp_norms <- c(temp_norms, norms[i])
      temp_blocks <- c(temp_blocks, list(dex))
      next
    }
    
    split_res <- block_ident_step(B[dex, dex, drop = FALSE], dex, pi[dex], thresh)
    split_blocks <- split_res$dex
    
    # Set to skip if block was irreducible or contains a singleton block 
    if (length(split_blocks) == 1) {
      temp_skip_ind <- c(temp_skip_ind, TRUE)
      temp_norms <- c(temp_norms, norms[i])
    } else {
      length_skip_ind <- map_dbl(split_blocks, length) == 1
      temp_skip_ind <- c(temp_skip_ind, length_skip_ind)
      temp_norms <- c(temp_norms, split_res$norm)
    }
    
    temp_blocks <- c(temp_blocks, split_blocks)
  }
  
  # test that the blocks are still mutually disjoint
  assertthat::assert_that(is_mutually_disjoint(temp_blocks))
  blocks <- temp_blocks
  skip_ind <- temp_skip_ind
  
  # Rotate before splitting
  I <- map(blocks, \(dex) seq_len(ncol(R)) %in% dex) %>% do.call(what = cbind)
  R <- clustered_varimax(R, I, centering = FALSE, eps = eps)$loadings
  R <- rotate_by_block(R, blocks)
  B <- crossprod(R^2) / colSums(R^2) # for dimension based clustering
  
  norms <- temp_norms
}

R <- clustered_varimax(R, I, centering = FALSE, eps = eps)$loadings
map(blocks, \(dex) rowSums(R[, dex, drop = FALSE]^2)) %>% do.call(what = cbind) %>% plot_correlations()

R_sum <- map(blocks, \(x) rowSums(R[, x, drop = FALSE]^2)) %>% do.call(what = cbind)
sort(apply(R_sum, 1, max))

clust_members <- apply(R_sum > 0.5, 2, \(x) names(which(x)), simplify = FALSE)
clust_members

alls <- map(clust_members, \(x) map_lgl(all_groups, \(y) all(y %in% x)))
anys <- map(clust_members, \(x) map_lgl(all_groups, \(y) any(y %in% x)))
map2(alls, anys, xor) %>% map_lgl(any)

norms
map_dbl(blocks, \(dex) v_prod(B, dex, dex, pi[dex]))
