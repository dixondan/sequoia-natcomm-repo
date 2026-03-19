cat("Checking R packages\n")
req_packages <- c("rjags", "coda", "matrixStats", "dplyr", "doParallel", "parallel", "foreach")

missing <- req_packages[!sapply(req_packages, require, character.only = TRUE)]

if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "), "\n",
    "Install them with install.packages() before running."
  )
}

args <- commandArgs(trailingOnly = TRUE)

n_iter <- if (length(args) >= 1) as.numeric(args[1]) else 500
ncores <- if (length(args) >= 2) as.numeric(args[2]) else 4
test_mode <- if (length(args) >= 3) as.logical(args[3]) else FALSE

if (test_mode) {
  cat("TEST MODE ON\n")
  n_iter <- 2
  ncores <- 1
}

cat("Running simulations\n")
cat("Ensemble runs:", n_iter, "\n")
cat("MCMC iterations per run: 200\n")
cat("Cores:", ncores, "\n")
cat("Test mode:", test_mode, "\n")

data_intermediate <- "data/intermediate"

config_name <- "ensemble-bern_90_c-0_90_r-0_90_no_structure"

out_dir <- file.path(data_intermediate, "model_outputs", config_name)
iter_dir <- file.path(out_dir, "iter_samples")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(iter_dir, showWarnings = FALSE)

input_file <- file.path(data_intermediate, "model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv")
bern_file <- file.path(data_intermediate, "stage2-bern_matrix_gamma_6000.rds")
grid_dir <- file.path(data_intermediate, "grid_samples")
jags_model_file <- file.path(data_intermediate, "logistic_mixed_no_structure.jags")

if (!file.exists(input_file)) stop("Missing input data in data/intermediate/")
if (!file.exists(bern_file)) stop("Run gamma calibration first (bern_matrix missing)")
if (!dir.exists(grid_dir)) stop("Missing grid samples in data/intermediate/grid_samples/")

cat("Reading data\n")
data_full <- read.csv(input_file, stringsAsFactors = FALSE)
bern_matrix <- readRDS(bern_file)

if (ncol(bern_matrix) != nrow(data_full)) stop("Mismatch between bern_matrix and data_full")

data_full <- data_full %>%
  mutate(
    pb_recent = ifelse(!is.na(year_burnt) & year_burnt >= 2011 & year_burnt <= 2021, 1, 0)
  )

model_string <- "
model {
  for (i in 1:N) {
    y_obs[i] ~ dbern(p[i])
    logit(p[i]) <- mu[i]

    mu[i] <- b_0
          + b_zq95      * zq95[i]
          + b_pb_recent * pb_recent[i]
          + b_upslope   * upslope[i]
          + b_solar     * solar[i]
          + b_slope     * slope[i]
          + b_tpi       * tpi[i]
          + b_us_solar  * (upslope[i] * solar[i])
          + b_us_slope  * (upslope[i] * slope[i])
          + eps[grove_id[i]]

    log_lik[i] <- logdensity.bern(y_obs[i], p[i])
  }

  b_0         ~ dnorm(0, 0.01)
  b_zq95      ~ dnorm(0, 0.01)
  b_pb_recent ~ dnorm(0, 0.01)
  b_upslope   ~ dnorm(0, 0.01)
  b_solar     ~ dnorm(0, 0.01)
  b_slope     ~ dnorm(0, 0.01)
  b_tpi       ~ dnorm(0, 0.01)
  b_us_solar  ~ dnorm(0, 0.01)
  b_us_slope  ~ dnorm(0, 0.01)

  mean.eps <- mean(eps[])
  mean.eps.star <- mean(eps.star[])
  b_0.star <- b_0 + mean.eps

  tau.eps ~ dgamma(0.1, 0.1)
  sig.eps <- 1 / sqrt(tau.eps)

  for (j in 1:Ngroves) {
    eps[j] ~ dnorm(0, tau.eps)
    eps.star[j] <- eps[j] - mean.eps
  }
}
"

writeLines(model_string, jags_model_file)

cols_to_scale <- c(
  "rumple_0_30", "cdensity_0_30", "rumple_0_60", "cdensity_0_60",
  "rumple_0_90", "cdensity_0_90", "rumple_0_120", "cdensity_0_120",
  "ladder1", "zq95", "solar_30m", "solar_60m", "solar_90m", "solar_120m",
  "tpi_30m", "tpi_60m", "tpi_90m", "tpi_120m",
  "slope_30m", "slope_60m", "slope_90m", "slope_120m"
)

cat("Scaling covariates\n")
data_scaled <- data_full %>%
  mutate(across(all_of(cols_to_scale), ~ as.numeric(scale(.)))) %>%
  mutate(row_id = row_number())

set.seed(2025)
draw_idx <- sample.int(nrow(bern_matrix), size = n_iter)

adapt_iter <- if (test_mode) 100 else 2000
burn_iter <- if (test_mode) 100 else 1000

cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

res <- foreach::foreach(
  k = seq_len(n_iter),
  .packages = c("rjags", "dplyr")
) %dopar% {
  out_rds <- file.path(iter_dir, sprintf("samples_iter_%03d.rds", k))
  
  if (!file.exists(out_rds)) {
    draw_i <- draw_idx[k]
    
    sample_file <- file.path(grid_dir, sprintf("samples90m-iter%d.csv", k))
    if (!file.exists(sample_file)) stop("Missing: ", sample_file)
    
    sampled_ids <- read.csv(sample_file, stringsAsFactors = FALSE)$tree_id
    df_sub <- data_scaled %>% filter(tree_id %in% sampled_ids)
    
    if (nrow(df_sub) == 0) stop("No sampled trees matched for k = ", k)
    
    col_idx <- df_sub$row_id
    
    data_list <- list(
      N = nrow(df_sub),
      zq95 = df_sub$zq95,
      pb_recent = df_sub$pb_recent,
      upslope = df_sub$upslope_12,
      solar = df_sub$solar_120m,
      slope = df_sub$slope_120m,
      tpi = df_sub$tpi_30m,
      grove_id = df_sub$grove_id,
      Ngroves = length(unique(df_sub$grove_id)),
      y_obs = bern_matrix[draw_i, col_idx]
    )
    
    jm <- rjags::jags.model(
      jags_model_file,
      data_list,
      n.chains = 3,
      n.adapt = adapt_iter,
      quiet = TRUE
    )
    
    update(jm, burn_iter, progress.bar = "none")
    
    samp_k <- rjags::coda.samples(
      jm,
      variable.names = c(
        "eps.star",
        "b_0.star",
        "b_zq95",
        "b_pb_recent",
        "b_upslope",
        "b_tpi",
        "b_solar",
        "b_slope",
        "b_us_solar",
        "b_us_slope",
        "log_lik"
      ),
      n.iter = 200,
      progress.bar = "none"
    )
    
    saveRDS(samp_k, out_rds)
  }
  
  out_rds
}

parallel::stopCluster(cl)
gc()

write.csv(
  data.frame(file = unlist(res)),
  file.path(out_dir, "run_manifest.csv"),
  row.names = FALSE
)

cat("-----parallel simulations successful------\n")