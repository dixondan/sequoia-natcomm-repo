cat("Checking R packages\n")
req_packages <- c(
  "coda", "matrixStats", "dplyr", "tidyr", "ggplot2",
  "bayesplot", "tidybayes", "parallel", "doParallel",
  "foreach", "rjags", "purrr", "tibble"
)

missing <- req_packages[!sapply(req_packages, require, character.only = TRUE)]

if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "), "\n",
    "Install them with install.packages() before running."
  )
}

data_intermediate <- "data/intermediate"
results_tables <- file.path("results", "tables")

dir.create(results_tables, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_intermediate, "model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv")
model_name <- "ensemble-bern_90_c-0_90_r-0_90_pb_recent"
iter_dir <- file.path(data_intermediate, "model_outputs", model_name, "iter_samples")

if (!file.exists(input_file)) {
  stop("Missing input data in data/intermediate/")
}

if (!dir.exists(iter_dir)) {
  stop("Missing model outputs in data/intermediate/model_outputs/")
}

rds_files <- list.files(
  iter_dir,
  pattern = "^samples_iter_[0-9]+\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No model output files found in iter_samples/")
}

cat("Reading and preprocessing tree data\n")
data_full <- read.csv(input_file, stringsAsFactors = FALSE)

cols_to_scale <- c(
  "rumple_0_30","cdensity_0_30","rumple_0_60","cdensity_0_60",
  "rumple_0_90","cdensity_0_90","rumple_0_120","cdensity_0_120",
  "ladder1","zq95",
  "solar_30m","solar_60m","solar_90m","solar_120m",
  "tpi_30m","tpi_60m","tpi_90m","tpi_120m",
  "slope_30m","slope_60m","slope_90m","slope_120m"
)

data_full_2 <- data_full %>%
  mutate(across(all_of(cols_to_scale), ~ as.numeric(scale(.)))) %>%
  mutate(
    ladder    = ladder1,
    cdensity  = cdensity_0_90,
    rumple    = rumple_0_90,
    slope     = slope_120m,
    tpi       = tpi_30m,
    solar     = solar_120m,
    upslope   = upslope_12,
    pb_recent = burned_recent,
    us_solar  = solar * upslope,
    us_slope  = slope * upslope
  )

cat("Reading pooled model outputs\n")

drop_log_vars <- function(chain) {
  keep <- !grepl("log", varnames(chain), ignore.case = TRUE)
  chain[, keep, drop = FALSE]
}

clean_chains <- lapply(rds_files, function(rds) {
  cat("Reading:", rds, "\n")
  raw_chains <- readRDS(rds)
  lapply(raw_chains, drop_log_vars)
})

flat_list <- unlist(clean_chains, recursive = FALSE)
big_post <- do.call(coda::mcmc.list, flat_list)

post_mat <- as.matrix(big_post)
rm(big_post, clean_chains, flat_list)
gc()

set.seed(2025)
n_sub <- 6000L
n_full <- nrow(post_mat)
keep <- sample.int(n_full, min(n_sub, n_full))
post_sub <- post_mat[keep, , drop = FALSE]
rm(post_mat)
gc()

cat("Preparing grove effects\n")
grove_ids <- sort(unique(data_full_2$grove_id))
eps_cols <- paste0("eps.star[", grove_ids, "]")

missing_eps <- setdiff(eps_cols, colnames(post_sub))
if (length(missing_eps) > 0) {
  stop("Missing eps.star columns in model outputs: ", paste(missing_eps, collapse = ", "))
}

eps_sub_full <- post_sub[, eps_cols, drop = FALSE]
grove_index <- match(data_full_2$grove_id, grove_ids)

cat("Building design matrices\n")
data2 <- data_full_2

X_base <- model.matrix(
  ~ zq95 + rumple + cdensity + ladder + pb_recent +
    upslope + slope + tpi + solar +
    us_solar + us_slope,
  data = data2
)

X_scen <- list(
  base_1  = X_base,
  nopb_2  = { X <- X_base; X[, "pb_recent"] <- 0; X },
  allpb_3 = { X <- X_base; X[, "pb_recent"] <- 1; X }
)

grp_ids <- c(as.character(grove_ids), "All")
idx_list <- split(seq_len(nrow(data2)), data2$grove_id)
idx_list$All <- seq_len(nrow(data2))

n_draws <- nrow(post_sub)
n_grps <- length(grp_ids)
n_scen <- length(X_scen)

counts <- array(
  NA_integer_,
  dim = c(n_draws, n_grps, n_scen),
  dimnames = list(NULL, grp_ids, names(X_scen))
)

params <- c(
  "(Intercept)" = "b_0.star",
  "zq95"        = "b_zq95",
  "rumple"      = "b_rumple",
  "cdensity"    = "b_cdensity",
  "ladder"      = "b_ladder",
  "pb_recent"   = "b_pb_recent",
  "upslope"     = "b_upslope",
  "slope"       = "b_slope",
  "tpi"         = "b_tpi",
  "solar"       = "b_solar",
  "us_solar"    = "b_us_solar",
  "us_slope"    = "b_us_slope"
)

missing_params <- setdiff(unname(params), colnames(post_sub))
if (length(missing_params) > 0) {
  stop("Missing fixed-effect parameters in model outputs: ", paste(missing_params, collapse = ", "))
}

cat("Running scenario simulations\n")
for (i in seq_len(n_draws)) {
  B <- as.numeric(post_sub[i, params])
  eps <- eps_sub_full[i, ]
  
  for (s in seq_len(n_scen)) {
    lin_full <- X_scen[[s]] %*% B + eps[grove_index]
    p <- plogis(lin_full)
    
    for (k in seq_len(n_grps)) {
      grp <- grp_ids[k]
      idx <- if (grp == "All") idx_list$All else idx_list[[grp]]
      counts[i, k, s] <- sum(rbinom(length(idx), 1, p[idx]))
    }
  }
}

summarize_scen <- function(x) {
  c(
    med = median(x),
    lo = quantile(x, 0.025),
    hi = quantile(x, 0.975)
  )
}

results <- purrr::map_dfr(seq_along(grp_ids), function(k) {
  stats_vec <- purrr::map(seq_len(n_scen), function(s) {
    summarize_scen(counts[, k, s])
  }) %>% purrr::flatten_dbl()
  
  names(stats_vec) <- unlist(
    lapply(names(X_scen), function(nm) paste0(nm, c("_med", "_lo", "_hi")))
  )
  
  tibble::tibble(
    grove = grp_ids[k],
    !!!as.list(stats_vec)
  )
})

cat("Writing scenario summary table\n")
write.csv(
  results,
  file.path(results_tables, "sequoia_survival_scenarios.csv"),
  row.names = FALSE
)

cat("Computing tree-level mean probabilities\n")
tree_probs <- array(
  NA_real_,
  dim = c(n_draws, nrow(data_full_2), length(X_scen)),
  dimnames = list(NULL, NULL, names(X_scen))
)

for (i in seq_len(n_draws)) {
  cat("Tree-level draw", i, "of", n_draws, "\n")
  B <- as.numeric(post_sub[i, params])
  eps <- eps_sub_full[i, ]
  
  for (s in names(X_scen)) {
    eta <- X_scen[[s]] %*% B + eps[grove_index]
    tree_probs[i, , s] <- plogis(eta)
  }
}

mean_prob <- apply(tree_probs, c(2, 3), mean)

data_full_2 <- data_full_2 %>%
  mutate(
    mortprob_base  = mean_prob[, "base_1"],
    mortprob_nopb  = mean_prob[, "nopb_2"],
    mortprob_allpb = mean_prob[, "allpb_3"],
    base_dead      = ifelse(mortprob_base >= 0.5, 1, 0),
    nopb_dead      = ifelse(mortprob_nopb >= 0.5, 1, 0),
    allpb_dead     = ifelse(mortprob_allpb >= 0.5, 1, 0)
  )

cat("Writing tree-level scenario outputs\n")
saveRDS(
  data_full_2,
  file.path(data_intermediate, "tree_level_scenarios.rds")
)

write.csv(
  data_full_2,
  file.path(results_tables, "tree_level_scenarios.csv"),
  row.names = FALSE
)

cat("-----scenarios successful------\n")

cat("\n----- Scenario Results -----\n\n")

print(
  results %>%
    mutate(across(-grove, ~ round(.x, 0)))
)

cat("\nInterpretation (All groves):\n")

all_row <- results %>% filter(grove == "All")

cat(sprintf(
  "Observed (base): %0.0f deaths [%0.0f, %0.0f]\n",
  all_row$base_1_med,
  all_row$base_1_lo,
  all_row$base_1_hi
))

cat(sprintf(
  "No prescribed burns: %0.0f deaths [%0.0f, %0.0f]\n",
  all_row$nopb_2_med,
  all_row$nopb_2_lo,
  all_row$nopb_2_hi
))

cat(sprintf(
  "All treated: %0.0f deaths [%0.0f, %0.0f]\n",
  all_row$allpb_3_med,
  all_row$allpb_3_lo,
  all_row$allpb_3_hi
))

cat("\nDifferences:\n")

cat(sprintf(
  "Burns prevented ~%0.0f deaths\n",
  all_row$nopb_2_med - all_row$base_1_med
))

cat(sprintf(
  "Universal treatment could have prevented ~%0.0f deaths\n",
  all_row$base_1_med - all_row$allpb_3_med
))