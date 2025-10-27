# TODO: alter pcalg and skeleton to skip conditioning sets that nulls one of the
# subsets

# - prevent it from exploring additional conditioning sets after one subset has
#   been found to null

get_incomp_pairs <- function(incomp_matrix) {
  incomp_matrix_lower <- incomp_matrix & lower.tri(incomp_matrix)
  incomp_pair_indices <- which(incomp_matrix_lower, arr.ind = TRUE)
  
  return(incomp_pair_indices)
}

compute_min_conditional_partial_diag <- function(pair, downsets, upsets, elements) {
  subsets <- elements[pair]
  upset <- reduce(upsets[pair], union) %>% unlist()
  downset <- reduce(downsets[pair], union) %>% unlist()
  condset <- setdiff(downset, upset)
  
  conditional_partial_mats <- map(rev(subsets), union, condset) %>%
    map2(subsets, ., partial_mat, mat = mat)
  
  min_conditional_diag <- map(conditional_partial_mats, diag) %>% map_dbl(min)
  
  return(min_conditional_diag)
}

conditional_partial_diag <- function(x, y, elements, mat) {
  chosen <- elements[[x]]
  remain <- unlist(elements[y])
  
  # Jitter the matrix very slightly to brute-force past dgesdd error
  done <- FALSE
  while (!done) {
    try({
      jittered_mat <- mat + runif(1, 0, 1e-12)*mat
      conditional_diag <- diag(partial_mat(jittered_mat, chosen, remain))
      done <- TRUE
    })
  }
  min_value <- min(conditional_diag)
  
  return(min_value)
}

transitive_closure <- function(covers_matrix) {
  p <- ncol(covers_matrix)
  closure_result <- matpow::matpow(covers_matrix, p, squaring = TRUE)
  
  closed_matrix <- closure_result$prod1 > 0
  
  return(closed_matrix)
}

resolve_eq_classes <- function(new_covers_matrix, covers_matrix, elements) {
  diff_covers_matrix <- new_covers_matrix & !(covers_matrix) | diag(diag(covers_matrix))
  
  # Let x and y be rows of the matrix. We can count the number of covers they
  # match on with the product x %*% t(y). Then to ensure equality, we compare
  # this match count to the total number of covers in x and y.
  match_count_matrix <- tcrossprod(diff_covers_matrix)
  total_covers <- diag(match_count_matrix)
  
  counts_match_x <- sweep(match_count_matrix, 1, total_covers, "==")
  counts_match_y <- sweep(match_count_matrix, 2, total_covers, "==")
  
  diff_eqclasses_matrix <- counts_match_x & counts_match_y
  
  
  # Now just for new covers matrix in its entirety
  match_count_matrix <- tcrossprod(new_covers_matrix)
  total_covers <- diag(match_count_matrix)
  
  counts_match_x <- sweep(match_count_matrix, 1, total_covers, "==")
  counts_match_y <- sweep(match_count_matrix, 2, total_covers, "==")
  
  eqclasses_matrix <- counts_match_x & counts_match_y | diff_eqclasses_matrix
  eqclasses_matrix <- eqclasses_matrix & !diag(ncol(covers_matrix))
  eqclasses <- which(eqclasses_matrix & upper.tri(eqclasses_matrix), arr.ind = TRUE)
  eqclasses <- eqclasses[order(eqclasses[, 2], decreasing = TRUE), , drop = FALSE]
  
  # merge equivalence classes
  for (i in seq_len(nrow(eqclasses))) {
    eq <- eqclasses[i, ]
    elements[eq[1]] <- list(union(elements[[eq[1]]], elements[[eq[2]]]))
    new_covers_matrix[eq[1], ] <- colSums(new_covers_matrix[eq, ]) > 0
  }
  merged <- unique(eqclasses[, 2])
  kept <- setdiff(seq_along(elements), merged)
  elements <- elements[kept]
  new_covers_matrix <- new_covers_matrix[kept, kept]
  
  return(list(covers_matrix = new_covers_matrix, elements = elements))
}

# Logic: 
# - A set is covered by another if given some conditioning set, it has a null on
#   the diagonal controlling for the other set
# - The conditioning set is derived from the downset and some kth order subset of
#   the incomparables to both sets being tested.
# - First, determine all first-order coverings.
# - Then determine all second-order coverings, choosing two elements from the
#   incomparables set. 
# - I guess both elements would be a "cover" if they zero it out, in that case
# - Check to make sure that if two sets are incomparable, the upset of one must 
#   be disjoint from the downset of the other.

# Let's pick a smallish subset
# elements <- c(AD_cluster, BD_cluster[-c(5, 6, 8)], CD_cluster)

# Now we'll get started implementing

# TODO:
# 1. Something wrong with equivalence classes (merging AE and CD even if 
#    covers_matrix is not TRUE??)
# 2. TO FIGURE OUT: I don't believe that it should be possible for AE to cover CD
#    or vice versa using only elements that are incomparable to them
# 3. Remove all of the testing AE, CD stuff.
# 4. Remove the downsets bit... I don't think that makes much sense tbh.
# 5. Consider adding upsets to the check? E.g., upsets shouldn't have zeros
#    either. Also none of the conditioning sets should have any elements of an 
#    upset in them. 
# 6. Remember the upsets aren't being updated rn.
# 7. What about min partial covariance of 1 as the indicator of covering.

library(pcalg)
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(pcalg))

mat <- samp_ortho_std
elements <- all_groups
tol <- sqrt(.Machine$double.eps)

# cover matrix indicating if row has zero partial controlling for col and the
# union of the downsets for both row and col
n_els <- length(elements)
partials_matrix <- diag(length(elements))
covers_matrix <- partials_matrix > tol

tol <- sqrt(.Machine$double.eps)
downsets_index <- map(elements, \(x) numeric(0))
downsets <- map(elements, \(x) NULL)
upsets <- imap(elements, \(x, y) y)
covers_matrix <- diag(length(elements)) > 0

for (ord in seq_along(elements)) {
  changed <- TRUE
  step <- 1
  while (changed) {
    changed <- FALSE
    cat("Order:", ord, "Step:", step, "\n")
    # TODO: Evaluate if this is better than the covers_matrix implementation
    # commented out below
    # le_matrix <- transitive_closure(covers_matrix)
    # incomp_matrix <- !(le_matrix | t(le_matrix))
    incomp_matrix <- !(covers_matrix | t(covers_matrix))
    uncovered <- which(colSums(covers_matrix) == 1)
    
    # Parallel iteration
    output <- foreach(i = seq_len(ncol(incomp_matrix))) %dopar% {
      # incomps <- uncovered[which(incomp_matrix[i, uncovered])]
      incomps <- which(incomp_matrix[i, ])
      n_incomps <- length(incomps)
    
      # If not enough incomparables to fill subset, go to next.
      if (n_incomps < ord) {
        # returning is equivalent to next in dopar
        # next()
        return(list(covers_matrix = covers_matrix, changed = FALSE))
      }
      
      # Start looping through until no more subsets
      subset_gen <- list(nextSet = seq_len(ord), wasLast = FALSE)
      while (!subset_gen$wasLast) {
        subset_dex <- subset_gen$nextSet
        incomp_subset <- incomps[subset_dex]
        subset <- unlist(elements[incomp_subset])
        rank <- rankM(mat[subset, subset, drop = FALSE])
        
        diag_value <- conditional_partial_diag(i, incomp_subset, elements, mat)
        if (diag_value < tol) {
          covers_matrix[i, incomp_subset] <- TRUE
          changed <- TRUE
        }
        
        # Choose next subset of incomps of size 'ord'
        subset_gen <- getNextSet(n_incomps, ord, subset_dex)
      }
      
      return(list(covers_matrix = covers_matrix, changed = changed))
    }
    
    output_combined <- list(covers_matrix = covers_matrix, changed = changed) %>%
      list() %>% c(output) %>% transpose() %>%
      map(reduce, `|`)
    new_covers_matrix <- output_combined$covers_matrix
    diff_matrix <- new_covers_matrix & !(covers_matrix)
    changed <- any(diff_matrix)
    
    # Combine elements if they are mutually nulling
    eq_classes_resolved <- resolve_eq_classes(new_covers_matrix, covers_matrix, elements)
    elements <- eq_classes_resolved$elements
    covers_matrix <- eq_classes_resolved$covers_matrix
    dims <- map_dbl(elements, \(x) rankM(mat[x, x, drop = FALSE]))
    
    # lt_matrix <- covers_matrix & t(!covers_matrix) # gets rid of mutual covers
    # downsets_index <- apply(lt_matrix, 2, which, simplify = FALSE)
    # downsets <- map(downsets_index, \(x) elements[x]) %>% map(unlist)
    # downsets <- map(elements, \(x) NULL)
    step <- step + 1
  }
}

stopCluster(cl)

# for (i in seq_along(upsets)) {
#   covers_matrix[i, upsets[[i]]] <- TRUE
# }
# covers_matrix <- covers_matrix & t(!covers_matrix) # gets rid of mutual covers



m <- lt_matrix; colnames(m) <- rownames(m) <- map_chr(elements, first)
g <- graph_from_adjacency_matrix(m, mode = "directed", diag = TRUE)
df <- get.data.frame(g)
r <- endorelation(
  domain = lapply(unique(unlist(df[c("from", "to")])), sets::as.set),
  graph = df[c("from", "to")]
)
rmat <- relation_incidence(transitive_reduction(r))
mattr <- m[row.names(rmat), colnames(rmat)] * rmat
gtr <- graph_from_adjacency_matrix(mattr, mode = "directed", diag = FALSE)

plot(gtr, layout = layout_with_sugiyama)


# Old Stuff ----------------------------------------------------------------



while (TRUE) {
  le_matrix <- matpow::matpow(covers_matrix, ncol(covers_matrix), squaring = TRUE)$prod1 > 0
  lt_matrix <- le_matrix & t(!le_matrix) # gets rid of mutual covers
  incomp_matrix <- !(le_matrix | t(le_matrix))
  
  # plot_correlations(lt_matrix, sort = FALSE)
  
  downsets <- apply(lt_matrix, 2, which, simplify = FALSE) %>%
    map(\(x) elements[x]) %>% map(unlist)
  
  upsets <- apply(le_matrix, 1, which, simplify = FALSE) %>%
    map(\(x) elements[x]) %>% map(unlist)
  
  incomp_pair_indices <- get_incomp_pairs(incomp_matrix)
  
  min_partial_diags <- apply(incomp_pair_indices, 1,
                             compute_min_conditional_partial_diag,
                             downsets, upsets, elements) %>% t()
  
  covers_matrix[incomp_pair_indices[, 1:2]] <- min_partial_diags[, 1] < tol
  covers_matrix[incomp_pair_indices[, 2:1]] <- min_partial_diags[, 2] < tol
  
  if (all(min_partial_diags > tol)) break
  print("hi")
  
  # find equivalence classes
  eqclasses_matrix <- covers_matrix & t(covers_matrix)
  eqclasses_matrix <- matpow::matpow(eqclasses_matrix, ncol(eqclasses_matrix), squaring = TRUE)$prod1 > 0
  eqclasses_matrix <- eqclasses_matrix & !diag(ncol(covers_matrix))
  eqclasses <- which(eqclasses_matrix & upper.tri(eqclasses_matrix), arr.ind = TRUE)
  eqclasses <- eqclasses[order(eqclasses[, 2], decreasing = TRUE), ]
  
  # merge equivalence classes
  for (i in seq_len(nrow(eqclasses))) {
    eq <- eqclasses[i, ]
    elements[[eq[1]]] <- union(elements[[eq[1]]], elements[[eq[2]]])
    covers_matrix[eq[1], ] <- colSums(covers_matrix[eq, ]) > 0
  }
  merged <- unique(eqclasses[, 2])
  kept <- setdiff(seq_along(elements), merged)
  elements <- elements[kept]
  covers_matrix <- covers_matrix[kept, kept]
}

# Problem, this thing terminates super early. There are 12 sets that have no
# lequal relations whatsoever. Perhaps I need to proceed in phases.

# - Iterate until no additional cover relations
# - Then do something with the sets that have cover relations and those that do not
# - Perhaps add all elements with cover relations into downsets, so long as they 
#   do not cover (or greater than) the generator of the downset.

# - phase 2: add structural relations
while (TRUE) {
  le_matrix <- matpow::matpow(covers_matrix, ncol(covers_matrix), squaring = TRUE)$prod1 > 0
  lt_matrix <- le_matrix & t(!le_matrix) # gets rid of mutual covers
  incomp_matrix <- !(le_matrix | t(le_matrix))
  
  # plot_correlations(lt_matrix, sort = FALSE)
  
  downsets <- apply(lt_matrix, 2, which, simplify = FALSE) %>%
    map(\(x) elements[x]) %>% map(unlist) %>%
    map(union, unlist(structural_relations)) %>%
    map(union, unlist(cross_loadings)) %>%
    map(union, unlist(error_covs))
  
  upsets <- apply(le_matrix, 1, which, simplify = FALSE) %>%
    map(\(x) elements[x]) %>% map(unlist)
  
  incomp_pair_indices <- get_incomp_pairs(incomp_matrix)
  
  min_partial_diags <- apply(incomp_pair_indices, 1,
                             compute_min_conditional_partial_diag,
                             downsets, upsets, elements) %>% t()
  
  covers_matrix[incomp_pair_indices[, 1:2]] <- min_partial_diags[, 1] < tol
  covers_matrix[incomp_pair_indices[, 2:1]] <- min_partial_diags[, 2] < tol
  
  if (all(min_partial_diags > tol)) break
  print("hi")
  
  # find equivalence classes
  eqclasses_matrix <- covers_matrix & t(covers_matrix)
  eqclasses_matrix <- matpow::matpow(eqclasses_matrix, ncol(eqclasses_matrix), squaring = TRUE)$prod1 > 0
  eqclasses_matrix <- eqclasses_matrix & !diag(ncol(covers_matrix))
  eqclasses <- which(eqclasses_matrix & upper.tri(eqclasses_matrix), arr.ind = TRUE)
  eqclasses <- eqclasses[order(eqclasses[, 2], decreasing = TRUE), ]
  
  # merge equivalence classes
  for (i in seq_len(nrow(eqclasses))) {
    eq <- eqclasses[i, ]
    elements[[eq[1]]] <- union(elements[[eq[1]]], elements[[eq[2]]])
    covers_matrix[eq[1], ] <- colSums(covers_matrix[eq, ]) > 0
  }
  merged <- unique(eqclasses[, 2])
  kept <- setdiff(seq_along(elements), merged)
  elements <- elements[kept]
  covers_matrix <- covers_matrix[kept, kept]
}
# - phase 3: loadings/within to downset, setdiff upset

# - phase 4: add between to downset, setdiff upset