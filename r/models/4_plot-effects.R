library(data.table)
library(matrixStats)
library(ggplot2)
library(bayesplot)
library(rjags)
library(coda)
library(dplyr)
library(tidyr)
library(tidybayes)
library(stringr)
library(here)

data_intermediate <- here("data", "intermediate")
results_figures <- here("results", "figures")
results_tables <- here("results", "tables")

model_name <- "ensemble-bern_90_c-0_90_r-0_90_pb_recent"
iter_dir <- file.path(data_intermediate, "model_outputs", model_name, "iter_samples")

if (!dir.exists(iter_dir)) {
  stop("iter_samples directory not found: ", iter_dir)
}

rds_files <- list.files(
  iter_dir,
  pattern = "^samples_iter_[0-9]+\\.rds$",
  full.names = TRUE
)

if (length(rds_files) == 0) {
  stop("No model output files found in iter_samples/")
}

input_file <- file.path(data_intermediate, "model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv")
if (!file.exists(input_file)) {
  stop("Missing input data in data/intermediate/")
}

dir.create(results_figures, showWarnings = FALSE, recursive = TRUE)
dir.create(results_tables, showWarnings = FALSE, recursive = TRUE)

data_full <- read.csv(input_file, stringsAsFactors = FALSE)

grove_lookup <- data_full %>%
  select(grove_id, grove) %>%
  distinct() %>%
  arrange(grove_id)

rand_names <- paste0("eps.star[", grove_lookup$grove_id, "]")
grove_labels <- setNames(grove_lookup$grove, rand_names)

fixed_pars <- c(
  "b_0.star","b_zq95","b_cdensity","b_rumple","b_ladder",
  "b_pb_recent","b_upslope","b_solar","b_slope",
  "b_us_solar","b_us_slope","b_tpi"
)

n_tot_fixed <- 0L
sum1_fixed <- setNames(numeric(length(fixed_pars)), fixed_pars)
sum2_fixed <- setNames(numeric(length(fixed_pars)), fixed_pars)
plot_list_fixed <- vector("list", length(fixed_pars))
names(plot_list_fixed) <- fixed_pars

n_tot_rand <- 0L
sum1_rand <- setNames(numeric(length(rand_names)), rand_names)
sum2_rand <- setNames(numeric(length(rand_names)), rand_names)
plot_list_rand <- vector("list", length(rand_names))
names(plot_list_rand) <- rand_names

for (rds in rds_files) {
  cat("Reading:", rds, "\n")
  chains <- readRDS(rds)
  mat <- as.matrix(chains)
  
  missing_fixed <- setdiff(fixed_pars, colnames(mat))
  missing_rand <- setdiff(rand_names, colnames(mat))
  
  if (length(missing_fixed) > 0) {
    stop("Missing fixed-effect parameters in saved model output: ", paste(missing_fixed, collapse = ", "))
  }
  if (length(missing_rand) > 0) {
    stop("Missing random-effect parameters in saved model output: ", paste(missing_rand, collapse = ", "))
  }
  
  fixed_mat <- mat[, fixed_pars, drop = FALSE]
  rand_mat <- mat[, rand_names, drop = FALSE]
  
  n_fixed <- nrow(fixed_mat)
  sum1_fixed <- sum1_fixed + colSums(fixed_mat)
  sum2_fixed <- sum2_fixed + colSums(fixed_mat^2)
  n_tot_fixed <- n_tot_fixed + n_fixed
  
  set.seed(0)
  idx_fixed <- sample(n_fixed, min(500, n_fixed))
  samp_fixed <- fixed_mat[idx_fixed, , drop = FALSE]
  for (p in fixed_pars) plot_list_fixed[[p]] <- c(plot_list_fixed[[p]], samp_fixed[, p])
  
  n_rand <- nrow(rand_mat)
  sum1_rand <- sum1_rand + colSums(rand_mat)
  sum2_rand <- sum2_rand + colSums(rand_mat^2)
  n_tot_rand <- n_tot_rand + n_rand
  
  set.seed(42)
  idx_rand <- sample(n_rand, min(500, n_rand))
  samp_rand <- rand_mat[idx_rand, , drop = FALSE]
  for (p in rand_names) plot_list_rand[[p]] <- c(plot_list_rand[[p]], samp_rand[, p])
  
  rm(mat, fixed_mat, rand_mat, samp_fixed, samp_rand)
  gc()
}

means_fixed <- sum1_fixed / n_tot_fixed
sds_fixed <- sqrt((sum2_fixed - n_tot_fixed * means_fixed^2) / (n_tot_fixed - 1))

stats_df <- data.frame(
  parameter = fixed_pars,
  mean = means_fixed,
  sd = sds_fixed,
  row.names = NULL,
  stringsAsFactors = FALSE
)

ci_mat_fixed <- t(sapply(fixed_pars, function(p) {
  qs <- quantile(
    plot_list_fixed[[p]],
    probs = c(0.025, 0.975, 0.05, 0.95, 0.075, 0.925, 0.10, 0.90)
  )
  c(
    lower95 = qs[1], upper95 = qs[2],
    lower90 = qs[3], upper90 = qs[4],
    lower85 = qs[5], upper85 = qs[6],
    lower80 = qs[7], upper80 = qs[8]
  )
}))

ci_df_fixed <- data.frame(
  parameter = rownames(ci_mat_fixed),
  ci_mat_fixed,
  row.names = NULL,
  stringsAsFactors = FALSE
)

full_stats <- merge(stats_df, ci_df_fixed, by = "parameter", sort = FALSE)

write.csv(
  full_stats,
  file.path(results_tables, "betas_stats_with_multi_CI.csv"),
  row.names = FALSE
)

mcmc_df_fixed <- bind_rows(
  lapply(names(plot_list_fixed), function(p) {
    data.frame(variable = p, samps = plot_list_fixed[[p]])
  }),
  .id = "ignored"
) %>% select(-ignored)

desired_order <- rev(c(
  "b_0.star","b_zq95","b_cdensity","b_rumple",
  "b_ladder","b_pb_recent","b_upslope","b_solar",
  "b_slope","b_us_solar","b_us_slope","b_tpi"
))

mcmc_long_fixed <- mcmc_df_fixed %>%
  mutate(variable = factor(variable, levels = desired_order))

years_since <- 10

new_labels <- c(
  "b_0.star" = "Intercept",
  "b_zq95" = "Sequoia height",
  "b_ladder" = "Ladder fuels",
  "b_tpi" = "Topographic \nPosition Index (TPI)",
  "b_solar" = "Heat Load \nIndex (HLI)",
  "b_upslope" = "Upslope burn",
  "b_pb_recent" = sprintf("Prescribed burn\n(last %d years)", years_since),
  "b_rumple" = "Height \nheterogeneity",
  "b_cdensity" = "Canopy Density",
  "b_us_solar" = "HLI * Upslope",
  "b_us_slope" = "Slope * Upslope",
  "b_slope" = "Slope"
)

effect_sizes_halfeye <- ggplot(mcmc_long_fixed, aes(x = samps, y = variable)) +
  stat_halfeye() +
  theme_bw(base_size = 12) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(
    x = "Effect size\nLog odds change in Pr(Sequoia mortality)\nfor a 1 standard deviation increase in covariate",
    y = NULL
  ) +
  scale_y_discrete(labels = new_labels) +
  scale_x_continuous(limits = c(-2.2, 1.1)) +
  theme(axis.title.x = element_text(size = 10))

ggsave(
  filename = file.path(results_figures, paste0(model_name, "--effects.png")),
  plot = effect_sizes_halfeye,
  width = 4.5,
  height = 4.8,
  units = "in"
)

print(effect_sizes_halfeye)

means_rand <- sum1_rand / n_tot_rand
sds_rand <- sqrt((sum2_rand - n_tot_rand * means_rand^2) / (n_tot_rand - 1))

stats_rand <- data.frame(
  parameter = rand_names,
  mean = means_rand,
  sd = sds_rand,
  row.names = NULL
)

ci_mat_rand <- t(sapply(rand_names, function(p) {
  quantile(plot_list_rand[[p]], probs = c(0.025, 0.975))
}))
colnames(ci_mat_rand) <- c("ci_lower", "ci_upper")

ci_df_rand <- data.frame(parameter = rownames(ci_mat_rand), ci_mat_rand, row.names = NULL)
stats_rand <- merge(stats_rand, ci_df_rand, by = "parameter", sort = FALSE)

write.csv(
  stats_rand,
  file.path(results_tables, "random_effects_CI.csv"),
  row.names = FALSE
)

mcmc_df_rand <- bind_rows(
  lapply(names(plot_list_rand), function(p) {
    data.frame(variable = p, samps = plot_list_rand[[p]])
  }),
  .id = "ignored"
) %>% select(-ignored)

mcmc_long_rand <- mcmc_df_rand %>%
  mutate(variable = factor(variable, levels = rev(rand_names)))

random_effects_plot <- ggplot(mcmc_long_rand, aes(x = samps, y = variable)) +
  stat_halfeye() +
  theme_bw(base_size = 12) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(
    x = "Random effect (log odds shift by grove)",
    y = NULL
  ) +
  scale_y_discrete(labels = grove_labels[rand_names]) +
  scale_x_continuous(limits = c(-4, 4))

ggsave(
  filename = file.path(results_figures, paste0(model_name, "--random-effects.png")),
  plot = random_effects_plot,
  width = 4.5,
  height = 6,
  units = "in"
)

print(random_effects_plot)

burn_row <- full_stats %>% filter(parameter == "b_pb_recent")

beta_est <- burn_row$mean
ci_lower <- burn_row$lower95
ci_upper <- burn_row$upper95

or_est <- exp(beta_est)
or_ci <- exp(c(ci_lower, ci_upper))

cat(sprintf(
  "Prescribed burn odds ratio = %.2f (95%% CI: [%.2f, %.2f])\n",
  or_est, or_ci[1], or_ci[2]
))

reduction <- (1 - or_est) * 100
reduction_ci <- c((1 - or_ci[2]), (1 - or_ci[1])) * 100

cat(sprintf(
  "Prescribed burns reduce mortality odds by %.0f%% (95%% CI: [%.0f%%, %.0f%%])\n",
  reduction, reduction_ci[1], reduction_ci[2]
))