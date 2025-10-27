library(purrr)
library(tidyr)

# USE Q[, 1:28] %*% backsolve(R, as.matrix(info$resid), 28, transpose = TRUE)

rankM <- Matrix::rankMatrix

named_svd <- function(x) {
  USV <- svd(x)
  dnames <- paste0('D', seq_along(USV$d))
  
  names <- list(
    d = rerun(2, dnames),
    u = list(rownames(x), dnames),
    v <- list(colnames(x), dnames)
  )
  
  USV$d <- diag(USV$d)
  out <- map2(USV, names, `dimnames<-`)
  
  return(out)
}

gsolve <- function(a, b = NULL, tol = .Machine$double.eps) {
  robj <- rankM(a)
  tol <- attr(robj, 'tol')
  if (robj < max(dim(a))) {
    S <- named_svd(a)
    nonzero <- diag(S$d) > tol
    diag(S$d)[nonzero] <- 1 / diag(S$d)[nonzero]
    a_inv <- S$u[, nonzero] %*% tcrossprod(S$d[nonzero, nonzero], S$v[, nonzero])
    if (is.null(b)) {
      return(a_inv)
    } else {
      return(t(a_inv) %*% b)
    }
  } else {
    if (is.null(b)) {
      return(solve(a, tol = tol))
    } else {
      return(solve(a, b, tol = tol))
    }
  }
}

mat_vec_funcs <- list(lambda = lavaan::lav_matrix_vec,
                      beta = lavaan::lav_matrix_vec,
                      psi = lavaan::lav_matrix_vechru,
                      theta = lavaan::lav_matrix_vechru)

build_name_mat <- function(x, sep) {
  dnames <- set_names(dimnames(x), c('row', 'col'))
  if(sep == '=~') {
    # fixes how loadings have the factor name first
    names(dnames) <- rev(names(dnames))
  }
  
  name_mat <- expand_grid(!!!rev(dnames)) %>%
    dplyr::relocate(row, .before = col) %>%
    apply(1, paste, collapse = sep) %>%
    matrix(nrow = nrow(x))
  
  return(name_mat)
}

vectorize_mat_names <- function(mats) {
  mats <- mats[names(mat_vec_funcs)]
  
  seps <- list(lambda = '=~', beta = '~', psi = '~~', theta = '~~')
  name_mats <- map2(mats, seps, build_name_mat)
  
  mat_names <- map2(mat_vec_funcs, name_mats, exec)
  
  return(mat_names)
}

vectorize_param_values <- function(mats) {
  mats <- mats[names(mat_vec_funcs)]
  vec_mats <- unlist(map2(mat_vec_funcs, mats, exec))
  
  mat_names <- unlist(vectorize_mat_names(mats))
  vec_mats <- set_names(vec_mats, mat_names)
  
  return(vec_mats)
}

vectorize_resids <- function(fit) {
  resids <- lavaan::lavResiduals(fit, 'raw')
  sig_resids <- lavaan::lav_matrix_vechru(resids$cov)
  
  return(c(resids$mean, sig_resids))
}

compute_jac_param_sigma <- function(mats) {
  lambda <- mats$lambda
  beta0 <- mats$beta
  psi <- mats$psi
  
  p <- nrow(lambda)
  q <- nrow(psi)
  
  K <- lavaan::lav_matrix_commutation(p, p)
  L <- lavaan::lav_matrix_duplication_ginv(p)
  D <- lavaan::lav_matrix_duplication(q)
  
  beta <- solve(diag(q) - beta0)
  BtPB <- beta %*% tcrossprod(psi, beta)
  LB <- lambda %*% beta
  LIK <- L %*% (diag(p^2) + K)
  
  J_lambda <- LIK %*% kronecker(lambda %*% BtPB, diag(p))
  J_beta <- LIK %*% kronecker(lambda %*% BtPB, LB)
  J_psi <- L %*% kronecker(LB, LB) %*% D
  J_theta <- diag(p * (p + 1) / 2)
  
  out <- cbind(J_lambda, J_beta, J_psi, J_theta)
  mat_names <- vectorize_mat_names(mats)
  dimnames(out) <- list(mat_names$theta, unlist(mat_names))
  
  return(out)
}

compute_hess_param_sigma <- function(mats) {
  lambda <- mats$lambda
  beta0 <- mats$beta
  psi <- mats$psi
  
  p <- nrow(lambda)
  q <- nrow(psi)
  
  D <- lavaan::lav_matrix_duplication(q)
  Kl <- lavaan::lav_matrix_commutation(p, q)
  Kb <- lavaan::lav_matrix_commutation(q, q)
  
  beta <- solve(diag(q) - beta0)
  BtPB <- beta %*% tcrossprod(psi, beta)
  LB <- lambda %*% beta
  
  # computes hessian matrix excluding theta (thetas all 0s) (neudecker 1991)
  mat_names <- vectorize_mat_names(mats)
  hess_names <- unlist(mat_names[c('lambda', 'beta', 'psi')])
  E <- diag(p)
  hess_sigmaij_param <- function(row, col) {
    Eij <- tcrossprod(E[, row, drop=FALSE], E[, col, drop=FALSE])
    Tij <- Eij + t(Eij)
    
    H11 <- kronecker(BtPB, Tij)
    H12 <- kronecker(beta, Tij %*% LB) %*% D
    H13 <- -(kronecker(BtPB, Tij %*% LB) +
               Kl %*% kronecker(Tij %*% lambda %*% BtPB, beta))
    
    H22 <- matrix(0, ncol(H12), ncol(H12))
    H23 <- -crossprod(D, kronecker(beta, crossprod(LB, Tij) %*% LB))
    
    H33 <- kronecker(BtPB %*% crossprod(lambda, Tij) %*% LB, t(beta)) %*% Kb +
      Kb %*% kronecker(crossprod(LB, Tij) %*% lambda %*% BtPB, beta) +
      kronecker(BtPB, crossprod(LB, Tij) %*% LB)
    
    H_s <- rbind(
      cbind(H11,       H12, H13),
      cbind(t(H12),    H22, H23),
      cbind(t(H13), t(H23), H33)
    )
    dimnames(H_s) <- rerun(2, hess_names)
    
    return(H_s)
  }
  
  indices <- expand_grid(row = seq_len(p), col = seq_len(p)) %>%
    dplyr::filter(row <= col) # vechru
  out <- map2(indices$row, indices$col, hess_sigmaij_param) %>%
    set_names(mat_names$theta)
  out <- do.call(abind::abind, c(out, list(along = 3)))
  
  return(out)
}

second_order_resid_objfunc <- function(x, resid, jac, hess) {
  p <- ncol(hess)
  x_hess <- x[seq_len(p)]
  
  hess_part1 <- apply(hess, 3, '%*%', y = x_hess)
  resid_hat <- crossprod(hess_part1, x_hess) / 2 + jac %*% x
  
  return(sum((resid_hat[, 1] - resid)^2))
}

grad_second_order_resid_objfunc <- function(x, resid, jac, hess) {
  p <- ncol(hess)
  d <- length(x) - p
  x_hess <- x[seq_len(p)]
  
  hess_part1 <- apply(hess, 3, '%*%', y = x_hess)
  resid_hat <- crossprod(hess_part1, x_hess) / 2 + jac %*% x
  diff <- resid_hat[, 1] - resid
  
  gr_hess_part <- cbind(t(hess_part1), matrix(0, nrow(jac), d))
  gr_part <- gr_hess_part + jac
  
  return(2 * colSums(apply(gr_part, 2, '*', diff)))
}

detect_bad_params <- function(free_parnames, full_parnames, std.lv = FALSE) {
  fac_labels <- full_parnames[str_detect(full_parnames, '=~')] %>%
    str_split('=~') %>%
    map_chr(1) %>% unique()
  
  reg_rev <-  free_parnames[!str_detect(free_parnames, '=~|~~')] %>%
    str_split('~') %>%
    map(rev) %>% map_chr(paste0, collapse='~')
  
  bad_params <- reg_rev
  
  if (std.lv) {
    var_pars <- full_parnames[str_detect(full_parnames, '~~')]
    fac_var_mask <- var_pars %>%
      str_split('~~') %>%
      map_lgl(~ all(.x %in% fac_labels))
    
    bad_params <- c(bad_params, var_pars[fac_var_mask])
  }
  
  return(full_parnames %in% bad_params)
}

delta_theta_full <- function(info, first_order = TRUE, allow_bad_params = FALSE) {
  resid <- info$resid
  delta.full <- info$delta.full
  
  if (allow_bad_params) {
    good_params <- colnames(delta.full)
  } else {
    bad_param_mask <- detect_bad_params(info$coef.names, info$coef.names.full,
                                        std.lv = TRUE)
    good_params <- colnames(delta.full)[!bad_param_mask]
  }
  
  delta.full <- delta.full[, good_params]
  
  if (first_order) {
    delta.full_inv <- MASS::ginv(delta.full)
    dimnames(delta.full_inv) <- rev(dimnames(delta.full))
    approx <- (delta.full_inv %*% resid)[, 1]
  } else {
    # approximate hessian by cross product of jacobian
    SJ <- compute_jac_param_sigma(info$est)[, good_params]
    SH <- compute_hess_param_sigma(info$est)
    SH_good_names <- dimnames(SH)[[1]][dimnames(SH)[[1]] %in% good_params]
    SH <- SH[SH_good_names, SH_good_names, rownames(SJ)]
    
    res <- optim(rep(0, ncol(SJ)), fn = second_order_resid_objfunc,
                 gr = grad_second_order_resid_objfunc,
                 resid = resid, jac = SJ, hess = SH,
                 method = "BFGS", control = list(maxit=500, abstol = 1e-12))
    if(res$convergence > 0) {
      warning("Second-order root finder did not converge.")
    }
    if(res$value > sqrt(.Machine$double.eps)) {
      warning("Second-order root finder did not converge.")
    }
    approx <- set_names(res$par, colnames(SJ))
    approx <- approx[colnames(delta.full)]
  }
  
  return(approx)
}


# Test Suite --------------------------------------------------------------

extend_fit_object_all_free <- function(object) {
  strict.exo <- object@Model@conditional.x
  
  FULL <- lavaan:::lav_partable_full(partable = object@ParTable,
                                     free = TRUE, start = TRUE,
                                     strict.exo = strict.exo)
  
  FULL$free <- rep(1L, nrow(FULL))
  FULL$user <- rep(10L, nrow(FULL))
  
  FIT <- lavaan:::lav_object_extended(object, add = FULL, all.free = TRUE)
  
  return(FIT)
}

test_params <- function(fit) {
  FIT <- extend_fit_object_all_free(fit)
  coef.FIT <- lavaan:::coef(FIT)
  
  mats <- lavaan::lavInspect(fit, "est")
  params <- vectorize_param_values(mats)
  
  diff <- params[names(coef.FIT)] - coef.FIT
  if (max(abs(diff)) > sqrt(.Machine$double.eps)) {
    warning("Est-mat vectorization and lavaan coef difference is nonzero.")
    return(diff)
  }
  "Cleared"
}

test_sigma_jac <- function(fit) {
  FIT <- extend_fit_object_all_free(fit)
  delta.FIT <- lavaan::lavInspect(FIT, 'delta')
  
  mats <- lavaan::lavInspect(fit, "est")
  SJ <- compute_jac_param_sigma(mats)
  
  diff <- SJ[rownames(delta.FIT), colnames(delta.FIT)] - delta.FIT
  if (max(abs(diff)) > sqrt(.Machine$double.eps)) {
    warning("Analytical grad and lavaan delta difference is nonzero.")
    return(diff)
  }
  "Cleared"
}

test_gr_resid_hat <- function(fit) {
  mats <- lavaan::lavInspect(fit, "est")
  fit_resids <- vectorize_resids(fit)
  SJ <- compute_jac_param_sigma(mats)
  SH <- compute_hess_param_sigma(mats)
  
  x <- rep(0, ncol(SJ))
  objfunc <- partial(second_order_resid_objfunc, resid = fit_resids,
                     jac = SJ, hess = SH)
  obj_val <- objfunc(x)
  
  gr_val <- grad_second_order_resid_objfunc(x, resid = fit_resids, jac = SJ,
                                            hess = SH)
  numeric_gr <- numericDeriv(quote(objfunc(x)), 'x')
  
  obj_diff <- obj_val - numeric_gr
  if (max(abs(obj_diff)) > 1e-5) {
    warning("Objective function and numericDeriv eval difference is nonzero.")
    return(obj_diff)
  }
  
  gr_diff <- gr_val - attr(numeric_gr, 'gradient')
  if (max(abs(gr_diff)) > 1e-5) {
    warning("Analytical grad and numericDeriv difference is nonzero.")
    return(gr_diff)
  }
  "Cleared"
}

._.model <- ' 
  f1 =~ start(1)*x1 + start(0.8)*x2 + start(1)*x3
  f2 =~ start(1.2)*x4 + start(0)*x5 + start(0.7)*x6
  f2 ~ start(0.2)*f1
  x1 ~~ start(1)*x1
  x2 ~~ start(1)*x2
  x3 ~~ start(1)*x3
  x4 ~~ start(1)*x4
  x5 ~~ start(1)*x5
  x6 ~~ start(1)*x6
  x2 ~~ start(-0.3)*x4
'
._.test_sim_data <- lavaan::simulateData(._.model, std.lv = TRUE, seed = 124)
._.test_fit <- lavaan::lavaan(._.model, data = ._.test_sim_data, std.lv = TRUE)
stopifnot(test_params(._.test_fit) == "Cleared")
stopifnot(test_sigma_jac(._.test_fit) == "Cleared")
stopifnot(test_gr_resid_hat(._.test_fit) == "Cleared")
