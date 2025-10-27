factors <- c("Alpha", "Beta", "Gamma", "Delta", "Epsilon")
manifest <- c("a", "b", "c", "d", "e")

structural_relations <- expand_grid(x = factors, y = factors) %>%
  group_by(x, y) %>%
  summarize(value = list(map(c("→", "↔"), \(op) paste0(x, op, y)))) %>%
  ungroup() %>%
  pull(value) %>%
  unlist(recursive = FALSE) %>%
  map(intersect, unlist(clusters_full)) %>%
  compact()

cross_loadings <- expand_grid(x = factors, y = manifest) %>%
  group_by(x, y) %>%
  summarize(value = list(map_chr(1:3, \(d) paste0(x, "_", y, d)))) %>%
  ungroup() %>%
  pull(value) %>%
  map(intersect, unlist(clusters_full)) %>%
  compact()

error_covs_df <- expand_grid(x = manifest, y = manifest) %>%
  group_by(x, y) %>%
  summarize(value = list(unlist(map(1:3, \(d) map(1:3, \(e) paste0(x, d, "↔", y, e)))))) %>%
  ungroup()

within_error_covs <- error_covs_df %>%
  filter(x == y) %>%
  pull(value) %>%
  map(intersect, unlist(clusters_full)) %>%
  compact()

between_error_covs <- error_covs_df %>%
  filter(x != y) %>%
  pull(value) %>%
  map(intersect, unlist(clusters_full)) %>%
  compact()

all_groups <- c(structural_relations, cross_loadings, within_error_covs, between_error_covs)
clusters_facrel <- structural_relations
clusters_loading <- cross_loadings
clusters_within <- within_error_covs
clusters_between <- between_error_covs

all_groups_subset <- all_groups[-c(20, 28, 29, 30, 31, 32, 33, 36, 37, 42, 48, 49)]
# all_groups_subset <- all_groups[-(c(28, 29, 30, 31) + 4)]
# all_groups_subset <- as.list(colnames(mat)[-1])


fit <- pseudo_frobenius(samp_ortho_std, all_groups_subset)
state <- fit %>%
  pluck(2) %>%
  get_state(2, all_groups_subset, mat)
clusters <- state %>%
  pluck(1)
clusters <- clusters[order(map_dbl(clusters, \(x) min(match(x, o))))] %>%
  map(~ .x[order(match(.x, o))]) %>% 
  rev()
clusters <- c(clusters, list(setdiff(colnames(mat)[-1], unlist(clusters))))
acc_remaining <- c(list(character(0)), accumulate(clusters, c)[-length(clusters)])
mat_partials <- map(clusters, \(x) c("resids", x)) %>%
  map2(acc_remaining, partial_mat, mat = mat)

map(mat_partials, diag) %>% map_dbl(\(x) min(x[-1]))
map(mat_partials, \(x) x[-1, -1]) %>% map(svd, 0, 0) %>% map_dbl(\(x) sum(x$d > tol))
map_dbl(mat_partials, \(x) x[1, -1] %*% MASS::ginv(x[-1, -1]) %*% x[-1, 1])
map2(clusters, mat_partials, \(x, y) y[x, x])[-length(clusters)] %>%
  map(plot_correlations, sort = FALSE)

betweens <- error_covs[c(2:5, 7:9, 11:12, 14)]
cluster_betweens <- map(clusters, \(x) x[x %in% unlist(betweens)])
cluster_struct <- map2(clusters, cluster_betweens, setdiff)
# split_clusters <- map2(cluster_struct, cluster_betweens, list) %>%
#   unlist(recursive = FALSE)
split_clusters <- c(cluster_struct, cluster_betweens) %>% compact()
acc_remaining_split <- c(list(character(0)), accumulate(split_clusters, c)[-length(split_clusters)])
mat_partials_split <- map(split_clusters, \(x) c("resids", x)) %>%
  map2(acc_remaining_split, partial_mat, mat = mat)

map(mat_partials_split, diag) %>% map_dbl(\(x) min(x[-1]))
map(mat_partials_split, \(x) x[-1, -1]) %>% map(svd, 0, 0) %>% map_dbl(\(x) sum(x$d > tol))
map_dbl(mat_partials_split, \(x) x[1, -1] %*% MASS::ginv(x[-1, -1]) %*% x[-1, 1])
map2(split_clusters, mat_partials_split, \(x, y) y[x, x]) %>%
  map(plot_correlations, sort = FALSE)
# TODO: 2 main structural clusters, one residual cluster



# TODO: Summarize the grouping above

# TO SUMMARIZE PROCESS FOR DEFINING THIS ORDERING
# - I noticed the structural relations had 6df
# - So i wanted to define a 6 partition ordering
# - To do this, I found the structural relation that has nonzero partial w.r.t
# all of the other structural relations. This was Gamma-Delta.
# - I built up the cluster around Gamma-Delta by exploring the loading subsets
# that resulted in a zero partial, after controlling for them and all of the
# structural relations outside of the gamma-delta cluster
# - Then, I looked for error-cov subsets that resulted in a zero partial, after
# controlling for them and all of the structural relations and loadings outside
# of the gamma-delta
# - Then I repeated this until there was no subset that yielded zero partials.
# - Then I moved to the next structural relation Beta-Delta

# TODO: Ultimately introduce and summarize this approach as an alternative that
# you were able to find by trying enough ways to partition the model without zero
# partials

all_groups <- c(structural_relations, cross_loadings, error_covs)
all_flat <- unlist(all_groups)

AE_cluster <- c(structural_relations[c(3, 4)], 
                cross_loadings[c(4, 5)],
                error_covs[c(1)])
AE_cluster_flat <- unlist(AE_cluster)

BE_cluster <- c(structural_relations[c(9)], 
                cross_loadings[c(6, 8)],
                error_covs[c(0)])
BE_cluster_flat <- unlist(BE_cluster)

AC_cluster <- c(structural_relations[c(5, 6)], 
                cross_loadings[c(1, 2, 17, 18)],
                error_covs[c(6, 2)])
AC_cluster_flat <- unlist(AC_cluster)

AD_cluster <- c(structural_relations[c(1, 2, 10)], 
                cross_loadings[c(3)],
                error_covs[c(0)])
AD_cluster_flat <- unlist(AD_cluster)

BD_cluster <- c(structural_relations[c(7, 8, 11)], 
                cross_loadings[c(7, 9, 13)],
                error_covs[c(3, 4, 5)])
BD_cluster_flat <- unlist(BD_cluster)

# CD_cluster <- c(structural_relations[c(8, 9)], 
#                 cross_loadings[c(10, 11, 12, 14, 15, 16, 19, 20)],
#                 error_covs[c(10, 13, 15, 7, 8, 9, 11, 12, 14)])

CD_cluster <- c(structural_relations[c(12, 13, 14)], 
                cross_loadings[c(11, 19, 15, 16)],
                error_covs[c(11, 10, 13)])
CD_cluster_flat <- unlist(CD_cluster)

lower_cluster <- c(structural_relations[c(0)], 
                   cross_loadings[c(10, 12, 14, 20)],
                   error_covs[c(15, 7, 8, 9, 12, 14)])
lower_cluster_flat <- unlist(lower_cluster)

# Make lower-order cluster than CD 

# - Ultimate goal: want to keep the Gamma-Delta and Delta_c and Gamma_d in the 
# same level but push the other loadings and error-covs to a lower level.
# - Want to do this for the rest of the clusters too.
# - The idea is this is the necessary ordering. And i can cluster things however
# I want as long as I respect this ordering.
# - I might be able to push a-a and Beta_a from AE to AC

all_clusters <- list(AE_cluster_flat, BE_cluster_flat, AC_cluster_flat,
                     AD_cluster_flat, BD_cluster_flat, CD_cluster_flat)

all_clusters_flat <- unlist(all_clusters)

# TEMP REMOVE BELOW

# Focus as chosen
diag(partial_mat(mat, unlist(lower_cluster), all_clusters_flat)) %>% min()

# TEMP REMOVE ABOVE

partial_mat(mat, setdiff(AE_cluster_flat, unlist(error_covs)),
            setdiff(all_flat, all_clusters_flat)) %>%
  diag()

# Focus as chosen
map(structural_relations, setdiff, all_clusters_flat) %>%
  map(partial_mat, mat = mat, chosen = AE_cluster_flat) %>%
  map(diag) %>% map_dbl(min)

map(cross_loadings, setdiff, all_clusters_flat) %>%
  map(c, setdiff(unlist(structural_relations), all_clusters_flat)) %>%
  map(partial_mat, mat = mat, chosen = AE_cluster_flat) %>%
  map(diag) %>% map_dbl(min)

map(error_covs, setdiff, all_clusters_flat) %>%
  map(c, setdiff(unlist(structural_relations), all_clusters_flat)) %>%
  map(c, setdiff(unlist(cross_loadings), all_clusters_flat)) %>%
  map(partial_mat, mat = mat, chosen = AE_cluster_flat) %>%
  map(diag) %>% map_dbl(min)

# pairwised
types <- list(clusters_facrel, clusters_loading, clusters_between)

type_clusters <- map(all_clusters, \(x) map(types, map, \(y) intersect(y, x))) %>%
  map_depth(2, compact)

type_clusters_grouped <- map(all_clusters, \(x) map(types, map, \(y) intersect(y, x))) %>%
  map_depth(2, compact) %>%
  map_depth(2, unlist) %>%
  unlist(recursive = FALSE)

type_clusters_flat <- type_clusters %>%
  unlist(recursive = FALSE) %>%
  unlist(recursive = FALSE)

acc_remaining <- c(list(character(0)),
                   accumulate(type_clusters_grouped, c)[-length(type_clusters_grouped)])
partial_full <- map2(type_clusters_grouped, acc_remaining, partial_mat, mat = mat)
clust_mat <- do.call(partial_full, what = lav_matrix_bdiag)
rownames(clust_mat) <- colnames(clust_mat) <- unlist(type_clusters_grouped)

plot_correlations_clustered(clust_mat, type_clusters_grouped)

map(type_clusters_flat, setdiff, x = unlist(type_clusters_flat)) %>%
  map2_dbl(type_clusters_flat, ., \(x, y) max(mat[x, y]^2))

map(type_clusters_flat, setdiff, x = unlist(type_clusters_flat)) %>%
  map2_dbl(type_clusters_flat, ., \(x, y) max(clust_mat[x, y]^2))

acc_remaining_flat <- c(list(character(0)),
                        accumulate(type_clusters_flat, c)[-length(type_clusters_flat)])
partial_type_flat <- map(type_clusters_flat, \(x) c("resids", x)) %>%
  map2(acc_remaining_flat, partial_mat, mat = mat)

mat_type_flat <- map(partial_type_flat, \(x) x[, -1] %*% MASS::ginv(x[-1, -1]) %*% x[-1, ])

LM_type_flat <- map_dbl(mat_type_flat, \(x) x[1, 1])
x2 <- round(LM_type_flat * samp_mats$modinds[1], 2)
df <- map_dbl(mat_type_flat, \(x) sum(svd(x[-1, -1], nu = 0, nv = 0)$d > tol))
p <- round(pchisq(x2, df, lower.tail = FALSE), 2)
MI_type_flat <- map_dbl(mat_type_flat, \(x) max(x[1, -1]^2 / diag(x)[-1]))
mi <- round(MI_type_flat * samp_mats$modinds[1], 2)
p_mi <- round(pchisq(mi, 1, lower.tail = FALSE), 2)

names_flat <- map(type_clusters_flat, str_replace_all, "\\d+", "") %>%
  map_chr(first)

block_flat <- type_clusters %>%
  map(unlist, recursive = FALSE) %>%
  map2(., ., \(x, y) rep_along(x, first(unlist(y)))) %>%
  unlist()

results_table_flat <- tibble(name = names_flat, block = block_flat,
                             X2 = x2, df = df, p = p, MI = mi, p.MI = p_mi)

view(results_table_flat)

results_table_flat %>%
  filter(p < 0.05) %>%
  view()

results_table_flat %>%
  filter(p.MI < 0.05) %>%
  view()

