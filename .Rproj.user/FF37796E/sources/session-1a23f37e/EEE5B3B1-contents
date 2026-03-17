cat("Checking R packages\n")
req_packages <- c("rjags", "coda", "lattice", "ggplot2")

missing <- req_packages[!sapply(req_packages, require, character.only = TRUE)]

if (length(missing) > 0) {
  stop(
    "Missing packages: ", paste(missing, collapse = ", "), "\n",
    "Install them with install.packages() before running."
  )
}

set.seed(2025)

data_intermediate <- "data/intermediate"
results_figures <- file.path("results", "figures", "calibration")

dir.create(data_intermediate, showWarnings = FALSE, recursive = TRUE)
dir.create(results_figures, showWarnings = FALSE, recursive = TRUE)

input_file <- file.path(data_intermediate, "model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv")

if (!file.exists(input_file)) {
  stop("Missing data. Put model-trees-Q12-sampled_75m-90m-120m-wgt-median.csv in data/intermediate/")
}

cat("Reading data\n")
dat <- read.csv(input_file, stringsAsFactors = FALSE)

if (nrow(dat) != length(unique(dat$tree_id))) {
  stop("Duplicate tree_id values detected in input data.")
}

y <- dat$dead_remap
p_m <- pmin(pmax(dat$mortality_prob, 1e-6), 1 - 1e-6)
p_s <- pmin(pmax(dat$survival_prob, 1e-6), 1 - 1e-6)

raw_m <- log(p_m / (1 - p_m))
raw_s <- log(p_s / (1 - p_s))

z_m <- as.numeric(scale(raw_m))
z_s <- as.numeric(scale(raw_s))
N <- length(y)

mu_raw_m <- mean(raw_m[is.finite(raw_m)], na.rm = TRUE)
sd_raw_m <- sd(raw_m[is.finite(raw_m)], na.rm = TRUE)
mu_raw_s <- mean(raw_s[is.finite(raw_s)], na.rm = TRUE)
sd_raw_s <- sd(raw_s[is.finite(raw_s)], na.rm = TRUE)

if (sd_raw_m <= 0 || !is.finite(sd_raw_m)) stop("Invalid sd_raw_m")
if (sd_raw_s <= 0 || !is.finite(sd_raw_s)) stop("Invalid sd_raw_s")

scalers <- list(
  mu_raw_m = mu_raw_m,
  sd_raw_m = sd_raw_m,
  mu_raw_s = mu_raw_s,
  sd_raw_s = sd_raw_s
)

saveRDS(scalers, file.path(data_intermediate, "calibration_scalers.rds"))

model_txt <- "
model {
  alpha   ~ dnorm(0, 1e-4)
  gamma_m ~ dnorm(0, 1e-4)
  gamma_s ~ dnorm(0, 1e-4)

  for (i in 1:N) {
    logit(pi[i]) <- alpha + gamma_m * z_m[i] + gamma_s * z_s[i]
    y[i] ~ dbern(pi[i])
  }
}
"

cat("Initializing JAGS model\n")
jm <- jags.model(
  textConnection(model_txt),
  data = list(N = N, z_m = z_m, z_s = z_s, y = y),
  n.chains = 3,
  n.adapt = 5000
)

cat("Running burn-in\n")
update(jm, 5000)

cat("Sampling posterior draws\n")
draws <- coda.samples(
  jm,
  variable.names = c("alpha", "gamma_m", "gamma_s"),
  n.iter = 2000
)

png(file.path(results_figures, "posterior_densities.png"), width = 1000, height = 1000, res = 120)
densityplot(
  draws,
  parms = c("alpha", "gamma_m", "gamma_s"),
  main = "Posterior densities",
  xlab = "Parameter value (logit scale)",
  scales = list(x = list(relation = "free", tick.number = 5))
)
dev.off()

print(gelman.diag(draws, multivariate = FALSE))
print(effectiveSize(draws))
print(summary(draws))

post_mat <- as.matrix(draws)
ndraw_total <- nrow(post_mat)

alpha_d <- post_mat[, "alpha"]
gamma_m_d <- post_mat[, "gamma_m"]
gamma_s_d <- post_mat[, "gamma_s"]

cat("Building posterior probability draws\n")
pi_draws <- sapply(seq_len(N), function(i) {
  plogis(alpha_d + gamma_m_d * z_m[i] + gamma_s_d * z_s[i])
})

cat("Simulating Bernoulli draws\n")
bern_matrix <- apply(pi_draws, 2, function(p) rbinom(ndraw_total, 1, p))

saveRDS(
  as.data.frame(post_mat[, c("alpha", "gamma_m", "gamma_s"), drop = FALSE]),
  file.path(data_intermediate, sprintf("stage1-param_draws_gamma_%d.rds", ndraw_total))
)

saveRDS(
  pi_draws,
  file.path(data_intermediate, sprintf("stage1-pi_draws_gamma_%d.rds", ndraw_total))
)

saveRDS(
  bern_matrix,
  file.path(data_intermediate, sprintf("stage2-bern_matrix_gamma_%d.rds", ndraw_total))
)

set.seed(2025)
example_idx <- sample.int(N, 9)
p_m_ex <- p_m[example_idx]

png(file.path(results_figures, "tree-examples.png"), width = 800, height = 800, res = 120)
old_par <- par(mfrow = c(3, 3), mar = c(4.2, 4.2, 2, 1))

for (j in seq_along(example_idx)) {
  idx <- example_idx[j]
  dens <- density(pi_draws[, idx], from = 0, to = 1, n = 512)
  plot(
    dens,
    main = paste0("Tree ", idx),
    xlab = "Calibrated probability",
    ylab = "Density",
    xlim = c(0, 1)
  )
  abline(v = p_m_ex[j], col = "red", lty = 2, lwd = 1)
  box()
}

dev.off()
par(old_par)

p_cal_mean <- colMeans(pi_draws)

df <- data.frame(
  raw = dat$mortality_prob,
  cal = p_cal_mean
)

g <- ggplot(df, aes(x = raw, y = cal)) +
  geom_point(alpha = 0.1, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Raw mortality probability", y = "Calibrated posterior mean") +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(results_figures, "prob-scatter.png"),
  plot = g,
  width = 5.5,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("-----calibration successful------\n")