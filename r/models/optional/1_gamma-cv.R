cat("Checking R packages\n")
req_packages <- c("rjags", "coda", "caret", "ggplot2")

missing <- req_packages[!sapply(req_packages, require, character.only = TRUE)]

if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "), "\n",
    "Install them with install.packages() before running."
  )
}

set.seed(2025)

data_intermediate <- "data/intermediate"
results_tables <- file.path("results", "tables", "optional")
results_figures <- file.path("results", "figures", "optional")

dir.create(results_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(results_figures, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_intermediate, "model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv")

if (!file.exists(input_file)) {
  stop("Missing data. Put model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv in data/intermediate/")
}

cat("Reading data\n")
dat <- read.csv(input_file, stringsAsFactors = FALSE)

## 1) Subset to labelled trees
lab <- which(!is.na(dat$dead_remap))
y_full <- dat$dead_remap[lab]

if (length(y_full) == 0) {
  stop("No labelled trees found in dead_remap.")
}

## 2) Extract and clamp raw scores
eps <- 1e-6
p_m0 <- pmin(pmax(dat$mortality_prob[lab], eps), 1 - eps)
p_s0 <- pmin(pmax(dat$survival_prob[lab], eps), 1 - eps)

## 3) Compute standardized logits
raw_m <- log(p_m0 / (1 - p_m0))
raw_s <- log(p_s0 / (1 - p_s0))

z_m <- as.numeric(scale(raw_m))
z_s <- as.numeric(scale(raw_s))

N <- length(y_full)

## 4) Create 5 stratified folds
folds <- caret::createFolds(y_full, k = 5, list = TRUE, returnTrain = FALSE)

## 5) Define gamma model in JAGS
model_txt <- "
model {
  alpha   ~ dnorm(0, 0.01)
  gamma_m ~ dnorm(0, 0.01)
  gamma_s ~ dnorm(0, 0.01)

  for (i in 1:N) {
    logit(pi[i]) <- alpha +
                    gamma_m * z_m[i] +
                    gamma_s * z_s[i]
    y[i] ~ dbern(pi[i])
  }
}
"

## 6) Helper functions
brier <- function(p, y) mean((p - y)^2)

metrics <- function(pred, truth) {
  tp <- sum(pred & truth)
  fp <- sum(pred & !truth)
  fn <- sum(!pred & truth)
  
  prec <- if (tp + fp == 0) NA else tp / (tp + fp)
  rec  <- if (tp + fn == 0) NA else tp / (tp + fn)
  f1   <- if (is.na(prec) || is.na(rec) || (prec + rec) == 0) NA else 2 * prec * rec / (prec + rec)
  
  c(F1 = f1, Precision = prec, Recall = rec)
}

summary_df <- function(mat, label) {
  data.frame(
    model = label,
    metric = colnames(mat),
    mean = colMeans(mat, na.rm = TRUE),
    sd = apply(mat, 2, sd, na.rm = TRUE),
    row.names = NULL
  )
}

## 7) Prepare results containers
met_names <- c("Brier", "F1", "Precision", "Recall")

cal_metrics <- matrix(
  NA,
  nrow = 5,
  ncol = 4,
  dimnames = list(paste0("Fold", 1:5), met_names)
)

unc_metrics <- cal_metrics

## 8) Cross-validation loop
n_draws <- 1000

cat("Running 5-fold cross-validation\n")
for (k in seq_along(folds)) {
  cat("Fold", k, "of", length(folds), "\n")
  
  test_idx <- folds[[k]]
  train_idx <- setdiff(seq_len(N), test_idx)
  
  y_cv <- y_full
  y_cv[test_idx] <- NA
  
  jm <- jags.model(
    textConnection(model_txt),
    data = list(
      N = N,
      z_m = z_m,
      z_s = z_s,
      y = y_cv
    ),
    n.chains = 3,
    n.adapt = 2000,
    quiet = TRUE
  )
  
  update(jm, 2000)
  
  samp <- as.matrix(
    coda.samples(
      jm,
      variable.names = c("alpha", "gamma_m", "gamma_s"),
      n.iter = n_draws
    )
  )
  
  b_cal <- b_unc <- numeric(n_draws)
  acc_cal <- acc_unc <- matrix(
    NA,
    nrow = n_draws,
    ncol = 3,
    dimnames = list(NULL, c("F1", "Precision", "Recall"))
  )
  
  truth <- y_full[test_idx] == 1
  p_raw_test <- p_m0[test_idx]
  
  for (d in seq_len(n_draws)) {
    a  <- samp[d, "alpha"]
    gm <- samp[d, "gamma_m"]
    gs <- samp[d, "gamma_s"]
    
    lp <- a + gm * z_m[test_idx] + gs * z_s[test_idx]
    p_cal <- plogis(lp)
    
    y_cal <- rbinom(length(p_cal), 1, p_cal) == 1
    y_unc <- rbinom(length(p_raw_test), 1, p_raw_test) == 1
    
    b_cal[d] <- brier(p_cal, truth)
    b_unc[d] <- brier(p_raw_test, truth)
    
    acc_cal[d, ] <- metrics(y_cal, truth)
    acc_unc[d, ] <- metrics(y_unc, truth)
  }
  
  cal_metrics[k, ] <- c(mean(b_cal), colMeans(acc_cal, na.rm = TRUE))
  unc_metrics[k, ] <- c(mean(b_unc), colMeans(acc_unc, na.rm = TRUE))
}

## 9) Summarize and write results
cal_summary <- summary_df(cal_metrics, "calibrated")
unc_summary <- summary_df(unc_metrics, "uncalibrated")
cv_summary <- rbind(cal_summary, unc_summary)

cat("\n--- Calibrated (gamma model) ---\n")
print(cal_summary, digits = 4)

cat("\n--- Raw (uncalibrated) ---\n")
print(unc_summary, digits = 4)

write.csv(
  cv_summary,
  file.path(results_tables, "cv_gamma_summary.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(fold = rownames(cal_metrics), cal_metrics, row.names = NULL),
  file.path(results_tables, "cv_gamma_fold_metrics_calibrated.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(fold = rownames(unc_metrics), unc_metrics, row.names = NULL),
  file.path(results_tables, "cv_gamma_fold_metrics_uncalibrated.csv"),
  row.names = FALSE
)

## 10) Simple comparison plot
plot_df <- cv_summary

g <- ggplot(plot_df, aes(x = metric, y = mean, fill = model)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.9),
    width = 0.2
  ) +
  labs(
    x = NULL,
    y = "Mean across folds",
    title = "5-fold CV: calibrated vs uncalibrated"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(results_figures, "cv_gamma_summary.png"),
  plot = g,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("-----cv gamma successful------\n")