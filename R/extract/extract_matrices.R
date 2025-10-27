library(lavaan)

build_pop_fit_object <- function(fit) {
  ptable <- parTable(fit) %>%
    mutate(free = ifelse(ustart != 0, seq_along(free), 0))
  
  pop_fit <- sem(ptable, do.fit = FALSE)
  
  return(pop_fit)
}

add_sampstats_fit_object <- function(fit) {
  implied <- lavInspect(fit, 'implied', drop.list.single.group = FALSE) %>%
    transpose()
  implied$cov <- map(implied$cov, ~ 2 * .x)
  
  fit <- update(fit, sample.cov = implied$cov, sample.mean = implied$mean,
                sample.nobs = 2, do.fit = FALSE)
  return(fit)
}

build_full_fit_object <- function(fit) {
  FULL <- lavaan:::lav_partable_full(partable = fit@ParTable,
                                     free = TRUE, start = TRUE,
                                     strict.exo = FALSE)
  FULL$user <- rep(10L, nrow(FULL))
  FULL$free <- rep(1L, nrow(FULL))
  FIT <- lavaan:::lav_object_extended(fit, add = FULL, all.free = TRUE)
  
  return(FIT)
}

add_constraints_fit_object <- function(fit, FIT) {
  PTABLE <- parTable(FIT) %>%
    mutate(plabel = ifelse(user == 0, paste0(".f", id, "."), plabel),
           free = ifelse(user == 0 & free == 0, max(free) + cumsum(user == 0 & free == 0), free))
  
  con_ptable <- PTABLE %>%
    filter(startsWith(plabel, ".f")) %>%
    mutate(id = id, lhs = plabel, op = "==", rhs = as.character(start), plabel = NA,
           free = 0) %>%
    bind_rows(PTABLE, .)
  
  CON_FIT <- lavaan:::lav_object_extended(fit, add = con_ptable, all.free = TRUE)
  
  return(CON_FIT)
}

extract_matrices <- function(fit) {
  FIT <- build_full_fit_object(fit)
  
  ptable <- parTable(fit)
  PTABLE <- parTable(FIT)
  free_params <- lav_partable_labels(ptable)
  all_params <- lav_partable_labels(PTABLE)
  moment_names <- lavaan:::lav_object_inspect_delta_rownames(fit)
  
  all_names <- list(all = all_params, free = free_params, moment = moment_names)
  
  fit <- build_pop_fit_object(fit) %>%
    add_sampstats_fit_object()
  FIT@SampleStats <- fit@SampleStats
  
  JAC <- lavaan:::lav_object_inspect_delta(FIT, add.labels = TRUE, add.class = TRUE)
  W <- lavaan:::lav_object_inspect_wls_v(fit, add.labels = TRUE, add.class = TRUE) %>%
    unname() %>%
    map2(moment_names, ~ `dimnames<-`(.x, rerun(2, .y)))

  info_all <- lavaan:::lav_object_inspect_information(
    FIT, information = "expected", add.labels = TRUE, add.class = TRUE
  )
  
  # nVCOV <- lavaan:::lav_object_inspect_information(
  #   fit, information = "expected", add.labels = TRUE, add.class = TRUE,
  #   inverted = TRUE
  # )
  return(list(jac = JAC, W = W, info = info_all, names = all_names))
}

extract_con_jac <- function(fit) {
  fit <- build_pop_fit_object(fit) 
  
  ptable <- parTable(fit)
  free_params <- lav_partable_labels(ptable)
  
  print(length(coef(fit)))
  print(length(free_params))
  
  con_JAC <- fit@Model@ceq.JAC
  colnames(con_JAC) <- free_params
  rownames(con_JAC) <- map_chr(seq_len(nrow(con_JAC)), ~ paste0(".C", .x, "."))
  
  if (nrow(con_JAC) == 0) {
    null_JAC <- diag(length(free_params))
  } else {
    null_JAC <- fit@Model@eq.constraints.K
  }
  colnames(null_JAC) <- free_params
  rownames(null_JAC) <- map_chr(seq_len(nrow(null_JAC)), ~ paste0(".N", .x, "."))
  
  return(list(con = con_JAC, null = null_JAC))
}

extract_subspaces <- function(mats) {
  full_info <- mats$info
  nVCOV <- mats$nvcov 
  
  free_info <- full_info[, colnames(nVCOV)] %*% nVCOV %*% full_info[colnames(nVCOV), ]
  ortho_info <- full_info - free_info
  
  class(free_info) <- c("lavaan.matrix.symmetric", "matrix")
  class(ortho_info) <- c("lavaan.matrix.symmetric", "matrix")
  
  return(list(full = full_info, free = free_info, ortho = ortho_info))
}
