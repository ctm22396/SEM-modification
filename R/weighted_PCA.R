source("G:/My Drive/RESEARCH/Dissertation/Simulation/power_analysis/doublyNonCentralBetaPDF.R")
# Equivalent to average of modification indices divided by total LM
# The standardized statistic is, roughly, a strange average of t-distributed
# variables divided by the "degrees of freedom" K + 2, which puts it on the scale
# of cohen's D

# Based on stat from Kariya 1977 and King 1981

# Interpretation: 
# Essentially, it is a measure of how much closer the modindices (L %*% X) are
# to having a spherical distribution than they are to one with scatter matrix
# given by Sigma (LM covariance matrix)

# I believe this can also be made into a comparative test of two scatter matrices
# By dividing by colSums(colSums(L_other) %*% Y^2)

runeven_approx_error <- function(n, k, K, delta = K, power = 0.8, alpha = 0.05) {
  if (k > 0) {
    # Compute asymptotic unevenenness (pi) that achieves desired power
    pi <- est_prop(a = k, K = K, delta = delta, power = power)
    
    mu_x <- matrix(replicate(n, rnorm(k)), nrow = k)
    mu_u <- sweep(mu_x, 2, sqrt(colSums(mu_x^2)), "/")
    mu_y <- matrix(replicate(n, rnorm(K - k)), nrow = K - k)
    mu_v <- sweep(mu_y, 2, sqrt(colSums(mu_y^2)), "/")
    
    mu <- rbind(sqrt(pi) * mu_u, sqrt(1 - pi) * mu_v) * sqrt(delta)
  } else {
    mu_x <- matrix(replicate(n, rnorm(K)), nrow = K)
    mu_u <- sweep(mu_x, 2, sqrt(colSums(mu_x^2)), "/")
    
    mu <- mu_u * sqrt(delta)
  }
  
  return(mu)
}

targeted_rotation <- function(A, targets) {
  A_proj <- A[, targets] %*% MASS::ginv(A[targets, targets]) %*% A[targets, ]
  
  L <- compute_basis(A)
  L_targ <- compute_basis(A[targets, targets])
  L_comb <- cbind(A[, targets] %*% MASS::ginv(A[targets, targets]) %*% L_targ,
                  compute_basis(A - A_proj))
  
  
  inv_rot_targ <- crossprod(L, MASS::ginv(A)) %*% L_comb
  
  return(inv_rot_targ)
}

targeted_varimax_rotation <- function(A, targets) {
  A_proj <- A[, targets] %*% MASS::ginv(A[targets, targets]) %*% A[targets, ]
  
  L <- compute_basis(A)
  L_targ <- varimax(compute_basis(A[targets, targets]))$loadings
  L_comb <- cbind(A[, targets] %*% MASS::ginv(A[targets, targets]) %*% L_targ,
                  compute_basis(A - A_proj))
  
  
  inv_rot_targ <- crossprod(L, MASS::ginv(A)) %*% L_comb
  
  return(inv_rot_targ)
}

rtargeted_approx_error <- function(n, targets, Sigma, delta = rankM(Sigma), power = 0.8, alpha = 0.05) {
  k <- rankM(Sigma[, targets] %*% MASS::ginv(Sigma[targets, targets]) %*% Sigma[targets, ])
  K <- rankM(Sigma)
  mu_targ <- runeven_approx_error(n, k, K, delta, power, alpha)
  
  inv_rot_targ <- targeted_rotation(Sigma, targets)
  mu <- inv_rot_targ %*% mu_targ
  
  return(mu)
}

n <- 5e4

unmasked_mat <- samp_mats$std_info$ortho
unmasked_unstd_mat <- samp_mats$info$ortho
samp_o_mask <- !is.na(diag(unmasked_mat))
Q <- unmasked_mat[samp_o_mask, samp_o_mask]
Q_inv <- cov2cor(MASS::ginv(unmasked_unstd_mat[samp_o_mask, samp_o_mask][-1, -1])) +
  0*unmasked_unstd_mat[samp_o_mask, samp_o_mask][-1, -1]
q_obs <- Q[1, -1]

L <- compute_basis(Q[-1, -1])
L_inv <- compute_basis(Q_inv)
K <- ncol(L)

X <- replicate(n, rnorm(K))
Y <- sweep(X, 2, sqrt(colSums(X^2)), "/")

mu_null <- runeven_approx_error(n, 0, K)
X_null <- replicate(n, rnorm(K)) + mu_null
Y_null <- sweep(X_null, 2, sqrt(colSums(X_null^2)), "/")

k <- 1
mu_comp <- runeven_approx_error(n, k, K)
X_comp <- replicate(n, rnorm(K)) + mu_comp
Y_comp <- sweep(X_comp, 2, sqrt(colSums(X_comp^2)), "/")

mu_lcomp <- runeven_approx_error(n, k, K)
X_lcomp <- replicate(n, rnorm(K)) + rbind(mu_lcomp[-seq_len(k), ], mu_lcomp[seq_len(k), ])
Y_lcomp <- sweep(X_lcomp, 2, sqrt(colSums(X_lcomp^2)), "/")

# head(sort(((L^2) %*% colMeans(L^2))[, 1])) -- x2~~x7 is the least aligned constraint
# tail(sort(((L^2) %*% colMeans(L^2))[, 1])) -- x4~~x6 is the most aligned constraint
# targets <- tail(names(sort(((L^2) %*% colMeans(L^2))[, 1])))[-1][-2]
# targets <- "Gamma~~Delta"
# targets <- "Gamma~Alpha"
# targets <- "x4~~x6"
# targets <- c("Gamma~~Delta", "Gamma~Alpha", "b1~~e2")
targets <- misspecs
mu_mod <- rtargeted_approx_error(n, targets, Q[-1, -1])
X_mod <- replicate(n, rnorm(K)) + mu_mod
Y_mod <- sweep(X_mod, 2, sqrt(colSums(X_mod^2)), "/")


lambdas <- colMeans(L^2)
# V <- (colMeans(L^2) %*% (2*(diag(K)/K - 1/K^2) / (K + 2)) %*% colMeans(L^2))[1, 1]
V <- (mean(lambdas^2) - mean(lambdas)^2) * 2 / (K + 2)
stat_perf <- (colMeans((L %*% c(1, rep(0, K - 1)))^2) - 1/K) / sqrt(V) / sqrt(K + 2)
stat <- (colMeans((L %*% Y)^2) - 1/K) / sqrt(V) / sqrt(K + 2)
stat_null <- (colMeans((L %*% Y_null)^2) - 1/K) / sqrt(V) / sqrt(K + 2)
stat_comp <- (colMeans((L %*% Y_comp)^2) - 1/K) / sqrt(V) / sqrt(K + 2)
stat_lcomp <- (colMeans((L %*% Y_lcomp)^2) - 1/K) / sqrt(V) / sqrt(K + 2)
stat_mod <- (colMeans((L %*% Y_mod)^2) - 1/K) / sqrt(V) / sqrt(K + 2)

selected <- intersect(colnames(Q), con_select)
lambdas_sel <- colMeans(L[selected, ]^2)
# V_sel <- (colMeans(L[selected, ]^2) %*% (2*(diag(K)/K - 1/K^2) / (K + 2)) %*% colMeans(L[selected, ]^2))[1, 1]
V_sel <- (mean(lambdas_sel^2) - mean(lambdas_sel)^2) * 2 / (K + 2)
stat_sel <- (colMeans((L[selected, ] %*% Y)^2) - 1/K) / sqrt(V_sel) / sqrt(K + 2)
stat_msel <- (colMeans((L[selected, ] %*% Y_mod)^2) - 1/K) / sqrt(V_sel) / sqrt(K + 2)

# improve low comp power by selecting subset of components?
stat_c <- colSums((L[, -seq_len(K - 3*k)] %*% Y[-seq_len(K - 3*k), ])^2) 
stat_csel <- colSums((L[, -seq_len(K - 3*k)] %*% Y_lcomp[-seq_len(K - 3*k), ])^2)


# mean(stat_perf >= quantile(stat, 0.95))
mean(stat_null >= quantile(stat, 0.95))
mean(stat_comp >= quantile(stat, 0.95))
mean(stat_lcomp >= quantile(stat, 0.95))
mean(stat_mod >= quantile(stat, 0.95))

mean(stat_msel >= quantile(stat_sel, 0.95))
mean(stat_msel >= quantile(stat, 0.95))

mean(stat_csel >= quantile(stat_c, 0.95))
