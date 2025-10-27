as_unit_vec <- function(x) {
  onedim <- !is.matrix(x)
  if (onedim) {
    x <- as.matrix(x)
  }
  m <- ncol(x)
  u <- sweep(x, 2, sqrt(colSums(x^2)), '/')
  
  u[, seq_len(m), drop=onedim]
}

named_vech <- function(x) {
  p <- ncol(x)
  mat_names <- rownames(x)
  
  vec <- lav_matrix_vech(x)
  col_idx <- lav_matrix_vech_col_idx(p)
  row_idx <- lav_matrix_vech_row_idx(p)
  names(vec) <- paste0(mat_names[col_idx], '~~', mat_names[row_idx])
  
  return(vec)
}

# plot_BluRd <- function(x, breaks = 12, ...) {
#   max_abs <- max(abs(x), na.rm = TRUE)
#   sym_breaks <- seq(-max_abs, max_abs, length.out = floor(breaks / 2) * 2)
#   palette_func <- partial(hcl.colors, palette = 'Blue - Red 2')
#   plot.matrix:::plot.matrix(x, breaks = sym_breaks, col = palette_func, ...)
# }

get_FIT <- function(object) {
  strict.exo <- object@Model@conditional.x
  
  FULL <- lavaan:::lav_partable_full(partable = object@ParTable,
                                     lavpta = object@pta,
                                     free = TRUE, start = TRUE,
                                     strict.exo = strict.exo)
  
  FULL$free <- rep(1L, nrow(FULL))
  FULL$user <- rep(10L, nrow(FULL))
  
  FIT <- lavaan:::lav_object_extended(object, add = FULL, all.free = TRUE)
  
  return(FIT)
}

# TODO: add support for just identified models too
# TODO: perhaps work with bases (standardize the jacobian columns)
#   values are standard normals
# TODO: add support for constrained variables (transformed unconstrained param space)

whats <- c('ordered', 'delta.rownames', 'nobs', 'data', 'sampstat', 'implied',
           'resid', 'WLS.V', 'Gamma', 'fit', 'est', 'delta', 'list')
whats_full <- c('delta', 'list')
extract_info <- function(object, names = whats) {
  info <- map(names, lavInspect, object = object) %>%
    set_names(names) %>%
    update_list(coef.vec = coef(object)) %>%
    update_list(coef.names = ~ names(coef.vec))
  
  if ('list' %in% names) {
    info$partable <- info$list
    info$list = NULL
  }
  
  return(info)
}

extend_jacobian_info <- function(info) {
  df <- info$fit['df'] # fit is a double vector for some reason
  if (df == 0) {
    ortho.names <- character(0)
  } else {
    ortho.names <- paste0('OR', seq_len(df))
  }
  ortho.names <- ortho.names
  basis.names <- c(info$coef.names, ortho.names)
  
  delta.basis <- qr(info$delta, tol = sqrt(.Machine$double.eps)) %>%
    qr.Q(complete = TRUE)
  dimnames(delta.basis) <- list(info$delta.rownames, basis.names)

  delta.ortho <- delta.basis[, ortho.names]
  colnames(delta.basis) <- paste0('B', seq_len(ncol(delta.basis)))
  
  delta.coorth <- cbind(info$delta, delta.ortho)
  
  out <- list(df = df, delta.ortho = delta.ortho, ortho.names = ortho.names,
              delta.coortho = delta.coorth, delta.basis = delta.basis,
              basis.names = basis.names)
  
  return(out)
}

build_info_object <- function(object) {
  info <- extract_info(object, names = whats) %>%
    c(extend_jacobian_info(.)) %>%
    modify_at(c('sampstat', 'implied', 'resid'),
              ~ c(.x$mean, named_vech(.x$cov))) %>%
    modify_at(c('WLS.V', "Gamma"), `dimnames<-`, rerun(2, .$delta.rownames))
  
  object.full <- get_FIT(object)
  info.full <- extract_info(object.full, names = whats_full) %>%
    set_names(paste0(names(.), '.full')) %>%
    update_list(coef.fixed.vec = ~ coef.vec.full[coef.names.full],
                coef.fixed.names = ~ setdiff(coef.names.full, info$coef.names)) %>%
    update_list(delta.fixed = ~ delta.full[, coef.fixed.names])
  
  out <- c(info, info.full)
  
  return(out)
}

# delta_theta/nu
#   linear vs. quadratic approx
find_ord2_root <- function(x, hess, grad, resid) {
  sum((map_dbl(hess, ~ x %*% .x %*% x) / 2 + grad %*% x - resid)^2)
}

delta_theta <- function(info, first_order = TRUE) {
  resid <- info$resid
  delta.coortho <- -info$delta.coortho
  
  if (first_order) {
    delta.coortho_inv <- solve(delta.coortho)
    approx <- (delta.coortho_inv %*% resid)[, 1]
  } else {
    # approximate hessian by cross product of jacobian
    approx_hess <- delta.coortho %>%
      array_tree(margin = 1) %>%
      map(tcrossprod)
    
    res <- optim(rep_along(resid, 0), find_ord2_root, method = "BFGS",
                 hess = approx_hess, grad = delta.coortho, resid = resid,
                 control = list(maxit=5000, abstol = 1e-15))
    if(res$convergence > 0) {
      warning("Second-order root finder did not converge.")
    }
    if(res$value > sqrt(.Machine$double.eps)) {
      warning("Second-order root finder did not reach 0.")
    }
    approx <- set_names(res$par, colnames(delta.coortho))
  }
  
  return(approx)
}

delta_theta_full <- function(info, first_order = TRUE) {
  resid <- info$resid
  delta.full <- info$delta.full
  
  if (first_order) {
    delta.full_inv <- MASS::ginv(delta.full)
    dimnames(delta.full_inv) <- rev(dimnames(delta.full))
    approx <- (delta.full_inv %*% resid)[, 1]
  } else {
    # approximate hessian by cross product of jacobian
    approx_hess <- delta.full %>%
      array_tree(margin = 1) %>%
      map(tcrossprod)
    
    res <- optim(rep_along(delta.full[1,], 0), find_ord2_root, method = "BFGS",
                 hess = approx_hess, grad = delta.full, resid = resid,
                 control = list(maxit=5000, abstol = 1e-15))
    if(res$convergence > 0) {
      warning("Second-order root finder did not converge.")
    }
    if(res$value > sqrt(.Machine$double.eps)) {
      warning("Second-order root finder did not converge.")
    }
    approx <- set_names(res$par, colnames(delta.full))
  }
  
  return(approx)
}

# delta theta ideas:
#  - compute gradient of params as if you already optimized in OR dimensions,
#  - then compute newton step you to the same place
#  - from some prelim tests, it gets

# projections

# covars

# LR, LM, test stats

# HANDLE CONSTRAINED PARAMS: 
#   - project theta (may contain eq.constrained params) into the unconstrained param space (using CON.jac)
#   - this can be done with the jac_thetas

# reparameterization (nu) form combinations of params (free)

# 
