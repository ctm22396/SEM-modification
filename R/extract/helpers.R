`%T>%` <- magrittr::`%T>%`

eps <- sqrt(.Machine$double.eps)

rankM <- function(x, tol = sqrt(.Machine$double.eps)) {
  sum(La.svd(x, 0, 0)$d > tol)
}

normat <- function(V) {
  mask <- diag(V) > 1e-5
  V[mask, mask] <- cov2cor(V[mask, mask])
  V
}

read_model <- function(file) {
  raw_text <- readr::read_file(file)
  parsed_text <- stringr::str_replace_all(raw_text, "\r\n", "\n")
  return(parsed_text)
}

named_svd <- function(x) {
  USV <- svd(x)
  dnames <- paste0('D', seq_along(USV$d))
  
  names <- list(
    d = rerun(2, dnames),
    u = list(rownames(x), dnames),
    v = list(colnames(x), dnames)
  )
  
  if (length(USV$d) <= 1) {
    USV$d <- as.matrix(USV$d)
  } else {
    USV$d <- diag(USV$d)
  }
  
  out <- map2(USV, names, `dimnames<-`)
  
  return(out)
}

named_symmetric_sqrt <- function(matrix) {
  matrix_sqrt <- lav_matrix_symmetric_sqrt(matrix)
  dimnames(matrix_sqrt) <- dimnames(matrix)
  
  return(matrix_sqrt)
}

as_unit_vec <- function(x) {
  onedim <- !is.matrix(x)
  if (onedim) {
    x <- as.matrix(x)
  }
  m <- ncol(x)
  u <- sweep(x, 2, sqrt(colSums(x^2)), '/')
  
  u[, seq_len(m), drop=onedim]
}

# Plotting Functions -------------------------------------------------------

load_order <- function(x) {
  o <- apply(x^2, 1, function(v) {
    c(which.max(v), -max(v))
  })
  o <- order(o[1, ], o[2, ])
  x[o, , drop = FALSE]
}

plot_correlations <- function(matrix, normalize = FALSE, sort = TRUE, text = FALSE, thresh = 0) {
  matrix <- matrix[rowSums(matrix^2) >= thresh, ]
  
  if (normalize) {
    val <- 1
    col_limits <- c(-1, 1)
  } else {
    val <- Inf
    a <- max(abs(range(matrix)))
    col_limits <- c(-a, a)
  }
  
  if (sort) {
    ordered_matrix <- load_order(matrix)  
  } else {
    ordered_matrix <- matrix  
  }
  
  # colnames(ordered_matrix) <- paste0("Z", seq_len(ncol(ordered_matrix)))
  
  long_matrix <- ordered_matrix %>%
    as.data.frame() %>%
    rownames_to_column("Parameter") %>%
    mutate(Parameter = fct_rev(as_factor(Parameter))) %>%
    pivot_longer(-Parameter, names_to = "Basis") %>%
    mutate(Basis = as_factor(Basis),
           value = pmin(pmax(value, -val), val))
  
  g <- ggplot(long_matrix, mapping = aes(x = Basis, y = Parameter, fill = value)) +
    geom_tile(color = "white")
  
  if (text) {
    g <- g + geom_text(mapping = aes(label = formatC(value, 2, 2, "f")))
  }
  g <- g +  
    scale_fill_gradientn(colors = c("dodgerblue1", "skyblue1", "white", "coral1", "brown1"),
                         guide = "none", limits = col_limits) +
    theme_minimal() +
    theme(panel.grid.major = element_blank())
  
  return(g)
}

plot_correlations_clustered <- function(mat, clusters) {
  nclusts <- length(clusters)
  total <- length(unlist(clusters))
  clust_sizes <- map_dbl(clusters, length)
  divisions <- c(0, cumsum(clust_sizes)) + 0.5
  
  horiz_data <- tibble(
    xstart = rep(divisions[-(nclusts + 1)], each = 2),
    xend = rep(divisions[-1], each = 2),
    ystart = total - c(divisions[1], rep(divisions[-c(1, nclusts + 1)], each = 2), divisions[(nclusts + 1)]) + 1,
    yend = total - c(divisions[1], rep(divisions[-c(1, nclusts + 1)], each = 2), divisions[(nclusts + 1)]) + 1
  )
  
  vert_data <- tibble(
    xstart = c(divisions[1], rep(divisions[-c(1, nclusts + 1)], each = 2), divisions[(nclusts + 1)]),
    xend = c(divisions[1], rep(divisions[-c(1, nclusts + 1)], each = 2), divisions[(nclusts + 1)]),
    ystart = total - rep(divisions[-(nclusts + 1)], each = 2) + 1,
    yend = total - rep(divisions[-1], each = 2) + 1
  )
  seg_data <- rbind(horiz_data, vert_data)
  
  plot_correlations(mat[unlist(clusters), unlist(clusters)], sort = FALSE) +
    geom_segment(mapping = aes(x = xstart, xend = xend, y = ystart, yend = yend), size = 1, inherit.aes = FALSE, data = seg_data) +
    geom_hline(yintercept = total - divisions + 1, size = 1, linetype = "dashed") +
    geom_vline(xintercept = divisions, size = 1, linetype = "dashed") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 10))
}

plot_signed_communalities <- function(communalities, resid, chisq, thresh = 0, sort = TRUE, sort_chisq = FALSE) {
  if (sort_chisq) {
    o <- names(sort((resid^2 * chisq)[resid^2 * chisq > thresh], decreasing = TRUE))
  } else {
    o <- names(resid)[resid^2 * chisq > thresh]
  }
  
  # flip the sign so it matches the resid sign
  communalities <- sweep(communalities, 2, sign(resid), FUN = '*')
  communalities <- pmax(pmin(communalities[, o], 1), -1)
  if (sort) {
    communalities <- communalities[rownames(load_order(communalities^2)), ]
  }
  
  chisq <- resid[o]^2 * chisq
  communalities <- rbind(Residual = abs(resid[o]), communalities)
  comm_tbl <- communalities %>%
    as.data.frame() %>%
    rownames_to_column("Modification") %>%
    mutate(Modification = fct_rev(as_factor(Modification)))
  
  g_comm <- pivot_longer(comm_tbl, -Modification, names_to = "Dimension") %>%
    mutate(Dimension = as_factor(Dimension)) %>%
    ggplot(mapping = aes(y = Modification, x = Dimension)) +
    geom_tile(mapping = aes(fill = value), color = "white", size = 0) +
    # geom_text(mapping = aes(label = formatC(value, 2, 2, "f"))) +
    scale_fill_gradientn(colors = c("dodgerblue1", "skyblue1", "white", "coral1", "brown1"),
                         guide = "none", limits = c(-1, 1)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank())
  
  chisq_tbl <- tibble(Dimension = colnames(comm_tbl)[-1], Chisq = chisq[o]) %>%
    mutate(Dimension = as_factor(Dimension))
  
  g_scree <- ggplot(chisq_tbl, mapping = aes(x = Dimension, y = Chisq)) +
    geom_path(group = 1) + geom_point() +
    ylim(0, NA) +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  g <- cowplot::plot_grid(g_scree, g_comm, nrow = 2, rel_heights = c(0.2, 0.8),
                          align = "v", axis = "lr")
  return(list(g, as.character(comm_tbl$Modification)))
}

source("./R/extract/unified_helper_funcs.R")
source("./R/extract/manual_grhess.R")

stdize_info <- function(info) {
  # get rid of effective zeros (floating point error)
  zeros <- diag(info) < sqrt(.Machine$double.eps)
  diag(info)[zeros] <- 0.1 # avoids warning from cov2cor
  std_info <- cov2cor(info)
  std_info[zeros, ] <- std_info[, zeros] <- NA
  return(std_info)
}

extract_info <- function(fit) {
  FIT <- extend_fit_object_all_free(fit)
  
  free_names = names(coef(fit))
  full_names = names(coef(FIT))
  fixed_names = setdiff(full_names, free_names)
  
  W <- lavInspect(fit, 'WLS.V')
  dimnames(W) <- rerun(2, inspect(FIT, 'delta.rownames'))
  jac <- lavInspect(FIT, 'delta') # changes in implied moments given unit changes in parameters
  resids <- vectorize_resids(fit) # change in implied moments to match sample moments
  all_vectors <- cbind(resids = resids, jac)
  
  full_info <- crossprod(all_vectors, W) %*% all_vectors
  std_full_info <- stdize_info(full_info)
  
  copla_info <- full_info[, free_names] %*%
    MASS::ginv(full_info[free_names, free_names]) %*%
    full_info[free_names, ]
  dimnames(copla_info) <- dimnames(full_info)
  std_copla_info <- stdize_info(copla_info)
  
  ortho_info <- full_info - copla_info
  std_ortho_info <- stdize_info(ortho_info)
  
  modinds <- ortho_info["resids", ]^2 / diag(ortho_info) * nobs(fit)
  modinds[is.na(diag(std_ortho_info))] <- 0
  
  info_list <- list(
    names = list(full = full_names, free = free_names, fixed = fixed_names),
    info = list(full = full_info, copla = copla_info, ortho = ortho_info),
    std_info = list(full = std_full_info, copla = std_copla_info, ortho = std_ortho_info),
    modinds = modinds,
    FIT = FIT,
    resids = list(vec = resids, mat = lavResiduals(fit, 'raw')$cov),
    W = W,
    jac = jac
  )
  
  return(info_list)
}

# use simstandard::fixed2free to get estimable model from pop model
