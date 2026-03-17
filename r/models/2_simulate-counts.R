cat("Checking R packages\n")
req_packages <- c("dplyr", "tibble")

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

bern_file <- file.path(data_intermediate, "stage2-bern_matrix_gamma_6000.rds")
input_file <- file.path(data_intermediate, "model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv")

if (!file.exists(bern_file)) {
  stop("Missing data. Run r/models/0_calibrate-gamma.R first.")
}

if (!file.exists(input_file)) {
  stop("Missing data. Put model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv in data/intermediate/")
}

cat("Reading Bernoulli draws and tree data\n")
bern_matrix <- readRDS(bern_file)
dat <- read.csv(input_file, stringsAsFactors = FALSE)

n_draws <- nrow(bern_matrix)
n_trees <- ncol(bern_matrix)

if (n_trees != nrow(dat)) {
  stop("Mismatch: bern_matrix columns do not match number of rows in input data.")
}

cat("Summarizing total mortality\n")
total_dead <- rowSums(bern_matrix)
overall_med <- median(total_dead)
overall_ci <- quantile(total_dead, c(0.025, 0.975))

cat(sprintf(
  "Total dead out of %d trees:\n  median = %.0f\n  95%% CI = [%.0f, %.0f]\n\n",
  n_trees, overall_med, overall_ci[1], overall_ci[2]
))

cat("Summarizing mortality by grove\n")
grove_lookup <- dat %>%
  dplyr::select(grove_id, grove) %>%
  dplyr::distinct() %>%
  dplyr::mutate(grove_id = as.character(grove_id))

grove_ids <- grove_lookup$grove_id
n_groves <- length(grove_ids)

dead_by_grove <- matrix(
  0L,
  nrow = n_draws,
  ncol = n_groves,
  dimnames = list(NULL, grove_ids)
)

for (j in seq_along(grove_ids)) {
  this_id <- grove_ids[j]
  tree_cols <- which(as.character(dat$grove_id) == this_id)
  dead_by_grove[, j] <- rowSums(bern_matrix[, tree_cols, drop = FALSE])
}

posterior_by_id <- tibble::tibble(
  grove_id = grove_ids,
  median_dead = apply(dead_by_grove, 2, median),
  lower_95 = apply(dead_by_grove, 2, quantile, probs = 0.025),
  upper_95 = apply(dead_by_grove, 2, quantile, probs = 0.975)
) %>%
  dplyr::mutate(grove_id = as.character(grove_id))

df_by_grove <- posterior_by_id %>%
  dplyr::left_join(grove_lookup, by = "grove_id") %>%
  dplyr::relocate(grove, .after = grove_id) %>%
  tibble::add_row(
    grove_id = "ALL",
    grove = "All Groves",
    median_dead = overall_med,
    lower_95 = overall_ci[1],
    upper_95 = overall_ci[2]
  )

tree_counts <- dat %>%
  dplyr::mutate(grove_id = as.character(grove_id)) %>%
  dplyr::count(grove_id, name = "n_trees") %>%
  dplyr::bind_rows(
    tibble::tibble(grove_id = "ALL", n_trees = nrow(dat))
  )

df_by_grove_rates <- df_by_grove %>%
  dplyr::mutate(grove_id = as.character(grove_id)) %>%
  dplyr::left_join(tree_counts, by = "grove_id") %>%
  dplyr::mutate(
    median_rate = median_dead / n_trees,
    lower_rate = lower_95 / n_trees,
    upper_rate = upper_95 / n_trees
  )

cat("Writing outputs\n")
write.csv(
  df_by_grove,
  file.path(results_tables, "counts-table.csv"),
  row.names = FALSE
)

write.csv(
  df_by_grove_rates,
  file.path(results_tables, "counts-and-rates.csv"),
  row.names = FALSE
)

cat("-----simulate counts successful------\n")