# This script produces the final quantile and CI predictions for Task C1.

library(tidyverse)
library(grf)
library(here)
library(evd)
library(ismev)
library(ggpubr)

# devtools::install_github("opasche/EQRN")
library(EQRN)

source("R/gbex_wrappers.R")
source("R/EGAM_wrappers.R")
source("R/CI_plot_helpers_loader.R")

## =============================== PARAMETERS ===============================

save_path <- "Results/EQRN_predictions_boot_transf/"

seedR <- 0
seedGRF <- 1
seedT <- seedR
force_refit <- FALSE

# PARAM: Data
train_path <- "data_wrangled/imputations/X_forest.csv"
test_path <- "data_wrangled/X_test_sd.csv"
featdrop <- c("Y_feat", "WindDirection")
prop_valid <- 1 / 4

# PARAM: General
intermediate_method <- "grf" # qrn grf
interm_lvl <- 0.8
prob_lvls_predict <- c(0.9999) # q,CId,CIu
confidence <- 0.5

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

# Params: EQRN
# path_eqrn <- paste0("Results/EQRN_grid_search_sf/networks/",intermediate_method,"/")#save_path
path_eqrn <- paste0("Results/EQRN_grid_search_sf/") # save_path
eqrn_params <- list(
  shape_fixed = TRUE, # TRUE FALSE
  net_structure = c(20, 10), # c(20,10) c(64,64,64)
  hidden_fct = torch_tanh,
  p_drop = 0,
  intermediate_q_feature = TRUE,
  learning_rate = 1e-4, # 1e-4 1e-4
  L2_pen = 1e-4, # 1e-4 5e-6
  shape_penalty = 0,
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

# Params: gbex
gbex_params <- list(
  B = 642,
  lambda = NULL,
  lambda_ratio = 5,
  lambda_scale = 0.01,
  depth = c(1, 0),
  sf = 0.75,
  intermediate_q_feature = TRUE,
  scale_features = FALSE
)

# Params: egam
egam_params <- list(
  model_shape = FALSE,
  intermediate_q_feature = gbex_params$intermediate_q_feature,
  scale_features = gbex_params$scale_features
)


# param_file <- "main/EQRN_iid/parameters/analysis/cosnorm_sigmoid_bestval.R"


parallel_strat <- "sequential" # "sequential", "multisession", "multicore"
n_workers <- 6


## =============================== END PARAMETERS ===============================


# source(param_file)
check_directory(save_path, recursive = TRUE)
check_directory(path_eqrn, recursive = TRUE)
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


# Naming standard for saved files
if (hid_f_str == "SSNN") {
  stuct_str <- paste0(
    "sc", paste0(eqrn_params$net_structure$scale, collapse = "_"), "_sh",
    paste0(eqrn_params$net_structure$shape, collapse = "_")
  )
} else {
  stuct_str <- paste0(eqrn_params$net_structure, "h", collapse = "_")
}
lr_str <- str_replace(toString(eqrn_params$learning_rate), "([.])", "p")
params_string <- paste0(
  "eqrn_", stuct_str, "_", hid_f_str, "_", "n"[!eqrn_params$intermediate_q_feature], "u_do", eqrn_params$p_drop * 100,
  "_L2", str_replace(toString(eqrn_params$L2_pen), "([.])", "p"), "_lr", lr_str
)

start_time <- Sys.time()


# GRF fit
if (!file.exists(paste0(save_path, "GRF_fit.rds")) | force_refit) {
  fit_grf <- quantile_forest(X_train, y_train,
    num.trees = num.trees, quantiles = quantiles_fit, sample.fraction = sample.fraction, mtry = mtry,
    min.node.size = min.node.size, honesty = honesty, honesty.fraction = honesty.fraction, honesty.prune.leaves = honesty.prune.leaves,
    alpha = alpha, imbalance.penalty = imbalance.penalty, seed = seedGRF
  )
  warning("Save not found, new intermediate GRF quantiles are fitted.")
  # safe_save_rds(fit_grf, paste0(path_data,"GRF_fit.rds"))
} else {
  fit_grf <- readRDS(paste0(path_data, "GRF_fit.rds"))
}
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

## ======== EQRN FIT ========

# Fit EQRN with intermediate OOB (or oracle) quantiles
if ((!dir.exists(paste0(path_eqrn, "networks/", params_string, "/"))) | force_refit) {
  warning("Save not found, new EQRN model is fitted.")
  torch_manual_seed(seedT)
  fit_eqrn <- EQRN_fit_restart(X_train, y_train, intermediate_quantiles,
    interm_lvl = interm_lvl, number_fits = eqrn_params$nb_restarts, shape_fixed = eqrn_params$shape_fixed,
    net_structure = eqrn_params$net_structure, hidden_fct = eqrn_params$hidden_fct, p_drop = eqrn_params$p_drop, intermediate_q_feature = eqrn_params$intermediate_q_feature,
    learning_rate = eqrn_params$learning_rate, L2_pen = eqrn_params$L2_pen, shape_penalty = eqrn_params$shape_penalty, scale_features = eqrn_params$scale_features,
    n_epochs = eqrn_params$n_epochs, batch_size = eqrn_params$batch_size, X_valid = X_valid, y_valid = y_valid, quant_valid = valid_quantiles,
    lr_decay = eqrn_params$lr_decay, patience_decay = eqrn_params$patience_decay, min_lr = eqrn_params$min_lr, patience_stop = eqrn_params$patience_stop,
    orthogonal_gpd = eqrn_params$orthogonal_gpd
  )
  # EQRN_save(fit_eqrn, paste0(path_eqrn, "networks/"), params_string)
} else {
  fit_eqrn <- EQRN_load(paste0(path_eqrn, "networks/"), params_string)
}
cat("\nElapsed time (EQRN fit):\n")
print(Sys.time() - start_time)


## ======== TESTING ========

# High quantile prediction with GRF
pred_grf_test <- predict(fit_grf, newdata = X_test, quantiles = prob_lvls_predict)$predictions

# UNCONDITIONAL predicted quantile(s) (Y quantile on X_train)
pred_unc <- predict_unconditional_quantiles(interm_lvl = interm_lvl, quantiles = prob_lvls_predict, Y = y_train_all, ntest = n_test)

# SEMI-CONDITIONAL predicted quantiles
pred_semicond <- predict_GPD_semiconditional(
  Y = y_train_all, interm_lvl = interm_lvl, thresh_quantiles = interm_quantiles_all,
  interm_quantiles_test = pred_interm, prob_lvls_predict = prob_lvls_predict
)

# Unconditional parameters and losses
uncond_losses_fixed <- unconditional_train_valid_GPD_loss(Y_train = y_train, interm_lvl = interm_lvl, Y_valid = y_valid)
uncond_losses_interm <- semiconditional_train_valid_GPD_loss(
  Y_train = y_train, Y_valid = y_valid,
  interm_quant_train = intermediate_quantiles,
  interm_quant_valid = valid_quantiles
)

## ======== GBEX and GPD GAM FITS =========
fit_gbex <- gbex_fit(
  X = X_train_all, y = y_train_all, intermediate_quantiles = interm_quantiles_all,
  interm_lvl = interm_lvl, intermediate_q_feature = gbex_params$intermediate_q_feature, scale_features = gbex_params$scale_features,
  B = gbex_params$B, lambda = gbex_params$lambda, lambda_ratio = gbex_params$lambda_ratio,
  lambda_scale = gbex_params$lambda_scale, depth = gbex_params$depth, sf = gbex_params$sf
)
pred_gbex <- gbex_predict(fit_gbex, X_test, to_predict = prob_lvls_predict, intermediate_quantiles = pred_interm, interm_lvl = interm_lvl)
pred_gbex_val <- gbex_predict(fit_gbex, X_valid, to_predict = prob_lvls_predict, intermediate_quantiles = valid_quantiles, interm_lvl = interm_lvl)
pred_gbex_train <- gbex_predict(fit_gbex, X_train, to_predict = prob_lvls_predict, intermediate_quantiles = intermediate_quantiles, interm_lvl = interm_lvl)

try({
  fit_egam <- fit_gpd_gam(
    X = X_train_all, y = y_train_all, intermediate_quantiles = interm_quantiles_all, interm_lvl = interm_lvl,
    model_shape = egam_params$model_shape,
    intermediate_q_feature = egam_params$intermediate_q_feature, scale_features = egam_params$scale_features
  )
  pred_egam <- predict_gpd_gam(fit_egam, X_test,
    to_predict = prob_lvls_predict,
    intermediate_quantiles = pred_interm, interm_lvl = interm_lvl
  )
})

## ===========================

train_plot <- training_plot_eqrn(fit_eqrn, uncond_losses_interm, show_legend = TRUE)
valid_plot <- validation_plot_eqrn(fit_eqrn, uncond_losses_interm, show_legend = FALSE)
plot(train_plot)
plot(valid_plot)
if (!is.null(save_path)) {
  save_myplot(
    plt = train_plot, plt_nm = paste0(save_path, params_string, "_training.pdf"),
    width = 75, height = 75, cairo = FALSE
  )
  save_myplot(
    plt = valid_plot, plt_nm = paste0(save_path, params_string, "_validation.pdf"),
    width = 75, height = 75, cairo = FALSE
  )
}

# Final predictions on X_test
pred_eqrn <- EQRN_predict(fit_eqrn, X_test, prob_lvls_predict, pred_interm, interm_lvl)
pred_eqrn_val <- EQRN_predict(fit_eqrn, X_valid, prob_lvls_predict, valid_quantiles, interm_lvl)
pred_eqrn_train <- EQRN_predict(fit_eqrn, X_train, prob_lvls_predict, intermediate_quantiles, interm_lvl)

nb_prob_lvls_predict <- length(prob_lvls_predict)

pred_eqrn_params <- EQRN_predict_params(fit_eqrn, X_test, pred_interm, return_parametrization = "classical", interm_lvl)
pred_eqrn_params_val <- EQRN_predict_params(fit_eqrn, X_valid, valid_quantiles, return_parametrization = "classical", interm_lvl)
pred_eqrn_params_train <- EQRN_predict_params(fit_eqrn, X_train, intermediate_quantiles, return_parametrization = "classical", interm_lvl)

pred_gbex_params <- gbex_predict(fit_gbex, X_test, to_predict = "par", intermediate_quantiles = pred_interm, interm_lvl = interm_lvl)
pred_gbex_params_val <- gbex_predict(fit_gbex, X_valid, to_predict = "par", intermediate_quantiles = valid_quantiles, interm_lvl = interm_lvl)
pred_gbex_params_train <- gbex_predict(fit_gbex, X_train, to_predict = "par", intermediate_quantiles = intermediate_quantiles, interm_lvl = interm_lvl)

RESULTS <- list(
  pred_eqrn = pred_eqrn, pred_grf_test = pred_grf_test, prob_lvls_predict = prob_lvls_predict,
  train_loss = fit_eqrn$train_loss, valid_loss = fit_eqrn$valid_loss,
  uncond_losses_interm = uncond_losses_interm, params_string = params_string
)
# if(!is.null(save_path)){
#   safe_save_rds(RESULTS, paste0(save_path, "Results_",params_string,".rds"))
# }

cat("\nElapsed time (predictions):\n")
print(Sys.time() - start_time)



# ========== Semi-parametric CI creation ===========

prob_lvl_predict <- prob_lvls_predict
boot_samples <- 1000
boot_size <- n_train_all
CI_type <- "normal"
parallel_strat <- "sequential"
n_workers <- (availableCores() - 1)
err_handling <- "stop"
seed <- 1

pars_trainall <- EQRN_predict_params(fit_eqrn, X_train_all,
  intermediate_quantiles = interm_quantiles_all, return_parametrization = "classical",
  interm_lvl = fit_eqrn$interm_lvl
)
sigmas <- pars_trainall$scales
xis <- pars_trainall$shapes

n_pred <- length(sigmas)
nb_prob_lvls_predict <- length(prob_lvl_predict)
if (nb_prob_lvls_predict != 1) {
  stop("Not implemented error: 'nb_prob_lvls_predict!=1'")
}

predicted_quantiles <- EQRN_predict(fit_eqrn, X_test, prob_lvl_predict, pred_interm, fit_eqrn$interm_lvl)

# random sampling of training observation indices for the bootstrap procedure
if (!file.exists(paste0(save_path, "pred_boot/", "samples_ind.rds")) | force_refit) {
  warning("bootstrap sample re-generated.")
  samples_ind <- semi_param_GPD_bootstrap_inds(
    y = y_train_all, intermediate_quantiles = interm_quantiles_all, scales = sigmas, shapes = xis,
    boot_samples = boot_samples, boot_size = boot_size, seed = seed
  )
  # safe_save_rds(samples_ind, paste0(save_path, "pred_boot/", "samples_ind.rds"))
} else {
  samples_ind <- readRDS(paste0(save_path, "pred_boot/", "samples_ind.rds"))
}

# Semi-parametric bootstrap procedure
if (!file.exists(paste0(save_path, "pred_boot/", params_string, "_pred_matrix.rds")) | force_refit) {
  n_workers <- min(n_workers, n_pred)
  `%fun%` <- set_doFuture_strategy(parallel_strat, n_workers = n_workers)

  # boot_preds <- foreach(b=seq_along(samples_ind), .errorhandling=err_handling, .combine=cbind) %fun% {
  boot_preds <- foreach(b = 1:500, .errorhandling = err_handling, .combine = cbind) %fun% {
    bindsb <- samples_ind[[b]]
    indsb <- bindsb$inds
    yb <- bindsb$yb
    Xb <- X_train_all[indsb, , drop = F]

    # train val split
    n_train_all <- boot_size
    n_test <- nrow(X_test)
    n_valid <- floor(prop_valid * n_train_all)
    n_train <- n_train_all - n_valid

    Xb_train <- Xb[1:n_train, , drop = F]
    yb_train <- yb[1:n_train]
    Xb_val <- Xb[n_train + (1:n_valid), , drop = F]
    yb_val <- yb[n_train + (1:n_valid)]

    Xb_train_all <- Xb[1:(n_train + n_valid), , drop = F]
    yb_train_all <- yb[1:(n_train + n_valid)]

    # Intermediate quantiles fit
    # GRF fit
    fit_grf <- quantile_forest(Xb_train, yb_train,
      num.trees = num.trees, quantiles = quantiles_fit, sample.fraction = sample.fraction, mtry = mtry,
      min.node.size = min.node.size, honesty = honesty, honesty.fraction = honesty.fraction, honesty.prune.leaves = honesty.prune.leaves,
      alpha = alpha, imbalance.penalty = imbalance.penalty, seed = seedGRF
    )
    intermqb_train <- predict(fit_grf, newdata = NULL, quantiles = c(interm_lvl))$predictions
    intermqb_val <- predict(fit_grf, newdata = Xb_val, quantiles = c(interm_lvl))$predictions # TODO fix Xb_val cl everywhere
    intermqb_all <- rbind(intermqb_train, intermqb_val)
    # Predict intermediate quantiles on X_test with GRF
    intermqb_test <- predict(fit_grf, newdata = X_test, quantiles = c(interm_lvl))$predictions

    # EQRN re-fit
    torch_manual_seed(seedT)
    fit_eqrn <- EQRN_load(paste0(path_eqrn, "networks/"), params_string)
    fit_eqrn_b <- EQRN_continue_training(fit_eqrn, Xb_train, yb_train, intermqb_train,
      learning_rate = 1e-5, L2_pen = eqrn_params$L2_pen,
      shape_penalty = eqrn_params$shape_penalty, n_epochs = 500, batch_size = eqrn_params$batch_size,
      X_valid = Xb_val, y_valid = yb_val, quant_valid = intermqb_val, lr_decay = 0.4, patience_decay = 3, min_lr = 1e-5,
      patience_stop = 5
    )
    # fit_eqrn_b <- EQRN_fit_restart(Xb_train, yb_train, intermqb_train, interm_lvl=interm_lvl, number_fits=eqrn_params$nb_restarts, shape_fixed=eqrn_params$shape_fixed,
    #                                net_structure=eqrn_params$net_structure, hidden_fct=eqrn_params$hidden_fct, p_drop=eqrn_params$p_drop, intermediate_q_feature=eqrn_params$intermediate_q_feature,
    #                                learning_rate=eqrn_params$learning_rate, L2_pen=eqrn_params$L2_pen, shape_penalty=eqrn_params$shape_penalty, scale_features=eqrn_params$scale_features,
    #                                n_epochs=eqrn_params$n_epochs, batch_size=eqrn_params$batch_size, X_valid=Xb_val, y_valid=yb_val, quant_valid=intermqb_val,
    #                                lr_decay=eqrn_params$lr_decay, patience_decay=eqrn_params$patience_decay, min_lr=eqrn_params$min_lr, patience_stop=eqrn_params$patience_stop,
    #                                orthogonal_gpd=eqrn_params$orthogonal_gpd)

    EQRN_save(fit_eqrn_b, paste0(save_path, "networks/boot/"), paste0(params_string, "_b_", b))

    pred_eqrnb <- EQRN_predict(fit_eqrn_b, X_test, prob_lvl_predict, intermqb_test, fit_eqrn_b$interm_lvl)

    safe_save_rds(pred_eqrnb, paste0(save_path, "pred_boot/", params_string, "_b_", b, ".rds"))

    matrix(c(pred_eqrnb), ncol = 1)
  }

  safe_save_rds(boot_preds, paste0(save_path, "pred_boot/", params_string, "_pred_matrix.rds"))

  end_doFuture_strategy()
} else {
  boot_preds <- readRDS(paste0(save_path, "pred_boot/", params_string, "_pred_matrix.rds"))
}

library("car")
car::qqPlot(boot_preds[85, ])

# Checking possible alternatives, but central/symmetric CIs are requested.
percentile_CIs <- make_prediction_CIs(predicted_quantiles, boot_preds, confidence, CI_type = "percentile")
basic_CIs <- make_prediction_CIs(predicted_quantiles, boot_preds, confidence, CI_type = "basic")
bc_CIs <- make_prediction_CIs(predicted_quantiles, boot_preds, confidence, CI_type = "bc")
normal_CIs <- make_prediction_CIs(predicted_quantiles, boot_preds, confidence, CI_type = "normal")
bcall_CIs <- make_prediction_CIs(predicted_quantiles, boot_preds, confidence, CI_type = "bcall")

cbind(pred_eqrn, percentile_CIs)
cbind(pred_eqrn, basic_CIs)
cbind(pred_eqrn, bc_CIs)
cbind(pred_eqrn, bcall_CIs)
cbind(pred_eqrn, normal_CIs)

cbind(pred_eqrn, percentile_CIs, basic_CIs, bc_CIs, bcall_CIs, normal_CIs)

sum(pred_eqrn > percentile_CIs[, 1] & pred_eqrn < percentile_CIs[, 2])
sum(pred_eqrn > basic_CIs[, 1] & pred_eqrn < basic_CIs[, 2])
sum(pred_eqrn > bc_CIs[, 1] & pred_eqrn < bc_CIs[, 2])
sum(pred_eqrn > normal_CIs[, 1] & pred_eqrn < normal_CIs[, 2])
sum(pred_eqrn > bcall_CIs[, 1] & pred_eqrn < bcall_CIs[, 2])

summary(c(percentile_CIs[, 2] - percentile_CIs[, 1]))
boxplot(c(percentile_CIs[, 2] - percentile_CIs[, 1]))
summary(c(basic_CIs[, 2] - basic_CIs[, 1]))
boxplot(c(basic_CIs[, 2] - basic_CIs[, 1]))
summary(c(bc_CIs[, 2] - bc_CIs[, 1]))
boxplot(c(bc_CIs[, 2] - bc_CIs[, 1]))
summary(c(normal_CIs[, 2] - normal_CIs[, 1]))
boxplot(c(normal_CIs[, 2] - normal_CIs[, 1]))
summary(c(bcall_CIs[, 2] - bcall_CIs[, 1]))
boxplot(c(bcall_CIs[, 2] - bcall_CIs[, 1]))


if (any(predicted_quantiles != pred_eqrn)) {
  stop("Debug error.")
}

# Choose normal CIs as central/symmetric CIs are requested.
lower_CI <- normal_CIs[, 1]
upper_CI <- normal_CIs[, 2]


# Save final (extreme quantile and 50% quantile CIs) test predictions for Task C1 (in multiple formats)
AnswerC1 <- cbind(pred_eqrn, lower_CI, upper_CI)

colnames(AnswerC1) <- c("quantile", "lower", "upper")

as_tibble(AnswerC1) %>% write_csv(paste0(save_path, "AnswerC1.csv"))
safe_save_rds(AnswerC1, paste0(save_path, "AnswerC1.rds"))

save(AnswerC1, file = paste0(save_path, "AnswerC1.Rdata"))


cat("\nElapsed time (boot CI):\n")
print(Sys.time() - start_time)

# Some basic analysis of the predictions
summary(y_train_all)
boxplot(y_train_all)

summary(c(intermediate_quantiles, valid_quantiles))
boxplot(c(intermediate_quantiles, valid_quantiles))
summary(c(pred_interm))
boxplot(c(pred_interm))

summary(c(pred_eqrn_train, pred_eqrn_val))
boxplot(c(pred_eqrn_train, pred_eqrn_val))
summary(c(pred_eqrn))
boxplot(c(pred_eqrn))


cat("\nRun time:\n")
print((Sys.time() - start_time))
