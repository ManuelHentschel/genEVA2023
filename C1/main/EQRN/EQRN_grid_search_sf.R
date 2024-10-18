# This script performs a grid search for optimal EQRN hyperparameters.


library(tidyverse)
library(grf)
library(here)
library(ggpubr)
library(evd)
library(ismev)
library(torch)
library(future)
library(doFuture)

# devtools::install_github("opasche/EQRN")
library(EQRN)
source("R/CI_plot_helpers_loader.R")

## =============================== PARAMETERS ===============================
save_path <- "Results/EQRN_grid_search_sf/"

seedR <- 0
seedGRF <- 1
seedT <- seedR

# PARAM: Data
train_path <- "data_wrangled/imputations/X_forest.csv"
test_path <- "data_wrangled/X_test_sd.csv"
featdrop <- c("Y_feat", "WindDirection")
prop_valid <- 1 / 4

# PARAM: General
intermediate_method <- "grf" # qrn grf
interm_lvl <- 0.8
prob_lvls_predict <- c(0.9999) # q,CId,CIu

# PARAM: GRF
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

# PARAM: EQRN
params_list <- list(
  shape_fixed = TRUE,
  net_structure = list(c(5, 3, 3), c(20, 10), c(10, 5, 5), c(10, 10, 10), c(20, 10, 10), c(64, 64, 64), c(128, 128, 128), c(32, 32, 32, 32), c(256, 256, 256, 256)),
  hidden_fct = list(torch_tanh),
  p_drop = list(0),
  intermediate_q_feature = list(TRUE),
  learning_rate = list(1e-3, 1e-4),
  L2_pen = list(0, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3),
  shape_penalty = list(0),
  scale_features = TRUE,
  orthogonal_gpd = TRUE,
  n_epochs = 1000,
  batch_size = 256,
  lr_decay = 0.4,
  patience_decay = 5,
  min_lr = 1e-5,
  patience_stop = 20,
  nb_restarts = 3
)
hid_f_str <- "tanh"

# param_file <- "main/EQRN/parameters/grid_search/cosnorm2d_sigmoid_oracle.R"


parallel_strat <- "multisession" # "sequential", "multisession", "multicore"
n_workers <- 12 # (availableCores() - 1)#PARAM

## =============================== END PARAMETERS ===============================

# source(param_file)
check_directory(paste0(save_path, "networks/"), recursive = TRUE)
set.seed(seedR)

start_time <- Sys.time()

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


# GRF fit
fit_grf <- quantile_forest(X_train, y_train,
  num.trees = num.trees, quantiles = quantiles_fit, sample.fraction = sample.fraction, mtry = mtry,
  min.node.size = min.node.size, honesty = honesty, honesty.fraction = honesty.fraction, honesty.prune.leaves = honesty.prune.leaves,
  alpha = alpha, imbalance.penalty = imbalance.penalty, seed = seedGRF
)

cat("\nElapsed time (GRF fit):\n")
print(Sys.time() - start_time)


## ======= INTERMEDIATE QUANTILES =======

if (intermediate_method == "grf") {
  # Out of bag quantiles prediction on X_train (for other QR ML methods, would need foldwise construction)
  intermediate_quantiles <- predict(fit_grf, newdata = NULL, quantiles = c(interm_lvl))$predictions
  valid_quantiles <- predict(fit_grf, newdata = X_valid, quantiles = c(interm_lvl))$predictions
  interm_quantiles_all <- rbind(intermediate_quantiles, valid_quantiles)
  # Predict intermediate quantiles on X_test with GRF
  pred_interm <- predict(fit_grf, newdata = X_test, quantiles = c(interm_lvl))$predictions
  pred_interm_fromall <- pred_interm
} else if (intermediate_method == "qrn") {
  stop("Not implemented yet.")
}

## ======== TESTING (GRF and Ground truth) ========

# High quantile prediction with GRF
pred_grf_test <- predict(fit_grf, newdata = X_test, quantiles = prob_lvls_predict)$predictions

# UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = y_train, ntest = n_test)

# SEMI-CONDITIONAL predicted quantiles
pred_semicond <- predict_GPD_semiconditional(
  Y = y_train, interm_lvl = interm_lvl, thresh_quantiles = intermediate_quantiles,
  interm_quantiles_test = pred_interm, prob_lvls_predict = prob_lvls_predict
)

# Unconditional parameters and losses
uncond_losses_fixed <- unconditional_train_valid_GPD_loss(Y_train = y_train, interm_lvl = interm_lvl, Y_valid = y_valid)
uncond_losses_interm <- semiconditional_train_valid_GPD_loss(
  Y_train = y_train, Y_valid = y_valid,
  interm_quant_train = intermediate_quantiles,
  interm_quant_valid = valid_quantiles
)

## ======== EQRN grid FIT ========

grid_l <- purrr::cross(params_list, .filter = NULL)

`%fun%` <- set_doFuture_strategy(parallel_strat, n_workers = n_workers)

results_grid_fit <- foreach(params = grid_l, .errorhandling = "remove", .combine = rbind) %fun% {
  source("R/CI_plot_helpers_loader.R")

  # Naming standard for saved files
  stuct_str <- paste0(params$net_structure, "h", collapse = "_")
  if (is.character(params$hidden_fct)) {
    if (params$hidden_fct == "SSNN") {
      stuct_str <- paste0("sc", paste0(params$net_structure$scale, collapse = "_"), "_sh", paste0(params$net_structure$shape, collapse = "_"))
    }
  }
  lr_str <- str_replace(toString(params$learning_rate), "([.])", "p")
  params_string <- paste0(
    "eqrn_", stuct_str, "_", hid_f_str, "_", "n"[!params$intermediate_q_feature], "u_do", params$p_drop * 100,
    "_L2", str_replace(toString(params$L2_pen), "([.])", "p"), "_lr", lr_str
  )

  cat("======== Start: ", params_string, " ========\n")

  # Fit EQRN with intermediate OOB quantiles
  torch_manual_seed(seedT)
  fit_eqrn <- EQRN_fit_restart(X_train, y_train, intermediate_quantiles,
    interm_lvl = interm_lvl, number_fits = params$nb_restarts, shape_fixed = params$shape_fixed,
    net_structure = params$net_structure, hidden_fct = params$hidden_fct, p_drop = params$p_drop, intermediate_q_feature = params$intermediate_q_feature,
    learning_rate = params$learning_rate, L2_pen = params$L2_pen, shape_penalty = params$shape_penalty, scale_features = params$scale_features,
    n_epochs = params$n_epochs, batch_size = params$batch_size, X_valid = X_valid, y_valid = y_valid, quant_valid = valid_quantiles,
    lr_decay = params$lr_decay, patience_decay = params$patience_decay, min_lr = params$min_lr, patience_stop = params$patience_stop,
    orthogonal_gpd = params$orthogonal_gpd
  )
  EQRN_save(fit_eqrn, paste0(save_path, "networks/"), params_string)

  train_plot <- training_plot_eqrn(fit_eqrn, uncond_losses_interm)
  ggsave(paste0(params_string, "_training.png"),
    plot = train_plot, device = "png",
    path = save_path, width = 200, height = 150, units = "mm", dpi = 300
  )

  ## ======== PREDICTIONS ========

  # Final EQRN predictions on X_test
  pred_eqrn <- EQRN_predict(fit_eqrn, X_test, prob_lvls_predict, pred_interm, interm_lvl)
  pred_eqrn_val <- EQRN_predict(fit_eqrn, X_valid, prob_lvls_predict, valid_quantiles, interm_lvl)
  pred_eqrn_train <- EQRN_predict(fit_eqrn, X_train, prob_lvls_predict, intermediate_quantiles, interm_lvl)

  nb_prob_lvls_predict <- length(prob_lvls_predict)

  output <- c(
    shape_fixed = params$shape_fixed,
    net_structure = stuct_str,
    hidden_fct = hid_f_str,
    p_drop = params$p_drop,
    intermediate_q_feature = params$intermediate_q_feature,
    L2_pen = params$L2_pen,
    shape_penalty = params$shape_penalty,
    Train_loss = fit_eqrn$train_loss[length(fit_eqrn$train_loss)], Valid_loss = fit_eqrn$valid_loss[length(fit_eqrn$valid_loss)],
    learning_rate = params$learning_rate,
    n_epochs = params$n_epochs,
    lr_decay = params$lr_decay,
    patience_decay = params$patience_decay,
    min_lr = params$min_lr,
    patience_stop = params$patience_stop,
    orthogonal_gpd = params$orthogonal_gpd
  )
  RESULTS <- list(
    pred_eqrn = pred_eqrn, pred_grf_test = pred_grf_test, prob_lvls_predict = prob_lvls_predict,
    train_loss = fit_eqrn$train_loss, valid_loss = fit_eqrn$valid_loss,
    uncond_losses_interm = uncond_losses_interm, params_string = params_string,
    pred_eqrn_val = pred_eqrn_val, pred_eqrn_train = pred_eqrn_train, output = output
  )
  safe_save_rds(RESULTS, paste0(save_path, "Results_", params_string, ".rds"))
  cat("==== End: ", params_string, " ====\n")
  output
}
results_tibble <- tibble::as_tibble(results_grid_fit)
results_tibble <- bind_cols(n_train = n_train, p = ncol(X_train), interm_lvl = interm_lvl, results_tibble)

# Save grid-search results
filename <- paste0(save_path, "results_EQRN_grid_fit_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")

write_csv(results_tibble, file = filename)

end_doFuture_strategy()

end_time <- Sys.time()
cat("\nRun time:\n")
print(end_time - start_time)


# Sys.sleep(300)
# system('shutdown -t 30 -s')
