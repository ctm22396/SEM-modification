my_is.corr.matrix <- function(R) {
  tol <- sqrt(.Machine$double.eps)
  if (!is.matrix(R)) {
    return(FALSE)
  }
  if (max(abs(R)) > 1 + tol) {
    return(FALSE)
  }
  if (nrow(R) != ncol(R)) {
    return(FALSE)
  }
  if (!isSymmetric(R, tol = tol)) {
    return(FALSE)
  }
  if (all(abs(diag(R) - 1) > tol)) {
    return(FALSE)
  }
  eigenvalues <- eigen(R, only.values = TRUE)$values
  if (any(eigenvalues < -tol)) {
    return(FALSE)
  }
  return(TRUE)
}

my_hcsvd <- function(R, q = "Kaiser", linkage = "average", is.corr = TRUE,
                     max.iter, trace = TRUE) {
  if (is.corr && !my_is.corr.matrix(R)) {
    stop("R must be a correlation matrix. Set 'is.corr = FALSE' if you want to supply a data matrix")
  }
  if (!is.corr) {
    X <- R
    if (anyNA(X)) {
      stop("X contains missing value indicator (NA)")
    }
    R <- stats::cor(X)
  }
  if (!((q == "Kaiser") | (is.numeric(q) & (q > 0) & (q <=
    1)))) {
    stop(paste(q), " is an invalid argument for q")
  }
  LINKAGE <- c("average", "single", "RV")
  if (!(linkage %in% LINKAGE)) {
    stop(paste(linkage), " is an invalid linkage function")
  }
  if (missing(max.iter)) {
    max.iter <- 500
  }
  p <- ncol(R)
  if (length(colnames(R)) == 0 | length(rownames(R)) == 0) {
    dimnames(R) <- list(as.character(seq_len(p)), as.character(seq_len(p)))
  }
  if (length(unique(colnames(R))) != p) {
    stop("Variable names are not unique")
  }
  q.p <- c()
  dist.matrix <- matrix(0, p, p)
  labels <- colnames(R)
  dimnames(dist.matrix) <- list(labels, labels)
  merge <- matrix(0, p - 1, 2)
  height <- vector(length = p - 1)
  order <- labels
  sub.matrices <- list(colnames(R))
  cluster.count <- p - 2
  for (iter in 1:(p - 1)) {
    current.labels <- labels %in% sub.matrices[[1]]
    hcsvd.ht <- bdsvd:::hcsvd.ht(
      R = R[current.labels, current.labels],
      q = q, linkage = linkage, labels = labels[current.labels],
      max.iter = max.iter, trace = trace
    )
    q.p <- c(q.p, hcsvd.ht$q.p)
    cluster.rows <- labels %in% hcsvd.ht$cluster1
    cluster.cols <- labels %in% hcsvd.ht$cluster2
    dist.matrix[cluster.rows, cluster.cols] <- hcsvd.ht$cluster.distance
    dist.matrix[cluster.cols, cluster.rows] <- hcsvd.ht$cluster.distance
    sub.matrices <- sub.matrices[-1]
    height[p - iter] <- hcsvd.ht$cluster.distance
    for (i in 1:2) {
      cluster <- hcsvd.ht[[paste0("cluster", i)]]
      if (length(cluster) != 1) {
        merge[p - iter, i] <- cluster.count
        cluster.count <- cluster.count - 1
        sub.matrices <- append(sub.matrices, list(cluster))
      } else {
        merge[p - iter, i] <- -which(labels == cluster)
      }
    }
    order.cluster1 <- order %in% hcsvd.ht$cluster1
    order.cluster2 <- order %in% hcsvd.ht$cluster2
    order[order.cluster1 | order.cluster2] <- c(
      order[order.cluster2],
      order[order.cluster1]
    )
    cat(sprintf(
      "\rSplit %d out of %d (%.2f%%)           ",
      iter, p - 1, iter / (p - 1) * 100
    ))
  }
  ordered.height <- order(height)
  merge <- merge[ordered.height, ]
  height <- height[ordered.height]
  # Removed the reordering
  # not.changed <- matrix(TRUE, p - 1, 2)
  # for (i in seq_len(p - 1)) {
  #   change.idx <- which(merge == i)
  #   merge[merge == i & not.changed] <- which(ordered.height ==
  #     i)
  #   not.changed[change.idx] <- FALSE
  # }
  hclust <- list(merge = merge, height = height, order = match(
    order,
    labels
  ), labels = labels, method = linkage)
  class(hclust) <- "hclust"
  u.cor <- 1 - dist.matrix
  dist.matrix <- stats::as.dist(dist.matrix)
  attr(dist.matrix, "Size") <- p
  cat("\r======== FINISHED ========                    ")
  cat("\n")
  return(list(
    hclust = hclust, dist.matrix = dist.matrix,
    u.cor = u.cor, q.p = q.p
  ))
}
