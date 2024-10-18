# This script performs a grid search for optimal gbex hyperparameters.

library(tidyverse)
library(here)
library(ggpubr)
library(grf)
library(evd)
library(ismev)
library(future)
library(doFuture)
source("R/gbex_wrappers.R")
source("R/CI_plot_helpers_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/gbex/grid_search/"
parallel_strat <- "multisession" # "sequential", "multisession", "multicore"
n_workers <- 12 # (availableCores() - 1)#PARAM
err_handling <- "remove" # "stop", "remove", "pass"

seedR <- 0
seedGRF <- 1
seedgbex <- 1

# PARAM: Data
train_path <- "data_wrangled/imputations/X_forest.csv"
test_path <- "data_wrangled/X_test_sd.csv"
featdrop <- c("Y_feat", "WindDirection")
prop_valid <- 1 / 4

# PARAM: General
intermediate_method <- "grf" # qrn grf
interm_lvl <- 0.8
prob_lvls_predict <- c(0.9999)

# Params: QRN
interm_path <- "data/Switzerland/qrn_intermediate_quantile_best/"
par_qrn <- list(
  nb_fits = 3,
  rnn_type = "lstm",
  num_layers = 2,
  hidden_size = 256,
  p_drop = 0,
  L2_pen = 1e-6,
  seq_len = 10,
  learning_rate = 1e-3,
  n_epochs = 1e3,
  batch_size = 256,
  lr_decay = 0.4,
  patience_decay = 5,
  min_lr = 1e-5,
  patience_stop = 20,
  scale_features = TRUE,
  tol = 1e-4
)

# Params: GRF
num.trees <- 5e3
quantiles_fit <- c(0.1, 0.5, 0.9)
sample.fraction <- 0.5
mtry <- min(ceiling(sqrt(9) + 20), 9)
min.node.size <- 5
honesty <- TRUE
honesty.fraction <- 0.5
honesty.prune.leaves <- TRUE
alpha <- 0.05
imbalance.penalty <- 0

# Params: gbex
num_folds <- 5
gbex_params <- list(
  intermediate_q_feature = TRUE,
  scale_features = FALSE,
  Bmax = 3000,
  grid_lambda_ratio = c(5, 6, 7, 8, 9, 10),
  grid_depth = list(c(1, 0), c(1, 1), c(2, 1), c(2, 2), c(3, 1), c(3, 2), c(3, 3)),
  stratified = TRUE,
  lambda_scale = 0.01,
  min_leaf_size = NULL,
  sf = 0.75
)

## =============================== END PARAMETERS ===============================


start_time <- Sys.time()

check_directory(save_path, recursive = TRUE)
set.seed(seedR)

# Data
# Training data
dat <- read_csv(train_path) %>% select(-any_of(featdrop))
X_test <- read_csv(test_path) %>% as.matrix()

n_train_all <- nrow(dat)
n_test <- nrow(X_test)
n_valid <- floor(prop_valid * n_train_all)
n_train <- n_train_all - n_valid

shuffle_ind <- sample.int(n_train_all)
# safe_save_rds(shuffle_ind, "data_wrangled/other/shuffle_ind.rds")
shuffle_ind <- readRDS("data_wrangled/other/shuffle_ind.rds")

dat_sh <- dat[shuffle_ind, ]

X <- dat_sh %>%
  select(-Y) %>%
  as.matrix()
Y <- dat_sh[["Y"]]

X_train <- X[1:n_train, , drop = F]
y_train <- Y[1:n_train]
X_valid <- X[n_train + (1:n_valid), , drop = F]
y_valid <- Y[n_train + (1:n_valid)]

X_train_all <- X[1:(n_train + n_valid), , drop = F]
y_train_all <- Y[1:(n_train + n_valid)]



## ======= INTERMEDIATE QUANTILES =======

if (intermediate_method == "qrn") {
  stop("Not implemented yet.")
} else if (intermediate_method == "grf") {
  # GRF fit
  fit_grf <- quantile_forest(X_train, y_train,
    num.trees = num.trees, quantiles = quantiles_fit, sample.fraction = sample.fraction, mtry = mtry,
    min.node.size = min.node.size, honesty = honesty, honesty.fraction = honesty.fraction, honesty.prune.leaves = honesty.prune.leaves,
    alpha = alpha, imbalance.penalty = imbalance.penalty, seed = seedGRF
  )

  cat("\nElapsed time (GRF fit):\n")
  print(Sys.time() - start_time)

  # Out of bag quantiles prediction on X_train (for other QR ML methods, would need foldwise construction)
  intermediate_quantiles <- predict(fit_grf, newdata = NULL, quantiles = c(interm_lvl))$predictions
  valid_quantiles <- predict(fit_grf, newdata = X_valid, quantiles = c(interm_lvl))$predictions
  interm_quantiles_all <- rbind(intermediate_quantiles, valid_quantiles)
  # Predict intermediate quantiles on X_test with GRF
  pred_interm <- predict(fit_grf, newdata = X_test, quantiles = c(interm_lvl))$predictions
  pred_interm_fromall <- pred_interm
} else if (intermediate_method == "oracle") {
  stop("No ground truth for real data.")
}




## ======== GBEX CV =========

CV_results <- gbex_CV(
  X = X_train_all, y = y_train_all, intermediate_quantiles = interm_quantiles_all,
  interm_lvl = interm_lvl, intermediate_q_feature = gbex_params$intermediate_q_feature, scale_features = gbex_params$scale_features,
  num_folds = num_folds, Bmax = gbex_params$Bmax, grid_lambda_ratio = gbex_params$grid_lambda_ratio, grid_depth = gbex_params$grid_depth,
  stratified = gbex_params$stratified, lambda_scale = gbex_params$lambda_scale, min_leaf_size = gbex_params$min_leaf_size, sf = gbex_params$sf,
  parallel_strat = parallel_strat, n_workers = n_workers, seed = seedgbex, err_handling = err_handling
)

filename <- paste0(save_path, "results_gbex_ts_CV_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")

write_csv(CV_results$results_tibble, file = filename)

cat(paste0("\n Best params:\n", CV_results$best_params, "\n\n"))

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)
