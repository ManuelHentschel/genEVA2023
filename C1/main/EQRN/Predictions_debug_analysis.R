# This script is a debug and plotting file analysing the C1 predictions from 'EQRN_semiparam_boot_prediction_transf.R'

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

save_path <- "Results/Predictions_debug/"

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



## =============================== END PARAMETERS ===============================


# source(param_file)
check_directory(save_path, recursive = TRUE)

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



start_time <- Sys.time()



# ========== CI 100 analysis ===========

paraboot <- readRDS("Results/EQRN_predictions_param/AnswerC1_alt.rds")
sparboott <- readRDS("Results/EQRN_predictions_boot_transf/AnswerC1.rds")

cbind(paraboot, sparboott)

sum(paraboot[, 1] > paraboot[, 2] & paraboot[, 1] < paraboot[, 3])

sum(sparboott[, 1] > sparboott[, 2] & sparboott[, 1] < sparboott[, 3])

summary(c(paraboot[, 3] - paraboot[, 2]))
summary(c(sparboott[, 3] - sparboott[, 2])) # CI width
boxplot(c(paraboot[, 3] - paraboot[, 2]))
boxplot(c(sparboott[, 3] - sparboott[, 2])) # CI width


summary(c(paraboot[, 1]))
summary(c(sparboott[, 1]))
boxplot(c(paraboot[, 1]))
boxplot(c(sparboott[, 1]))

summary(c(paraboot[, 2]))
summary(c(sparboott[, 2]))
boxplot(c(paraboot[, 2]))
boxplot(c(sparboott[, 2]))


summary(c(paraboot[, 3]))
summary(c(sparboott[, 3]))
boxplot(c(paraboot[, 3]))
boxplot(c(sparboott[, 3]))

my_boxplot <- function(vals, x_tick, y_lab = NULL, expand_pts = NULL) {
  boxplot_col <- my_palette_methods["EQRN"]
  mean_col <- my_palette$red

  bp_plot <- data.frame(Y = vals, X = rep(x_tick, length(vals))) %>%
    ggplot(aes(x = X, y = Y, group = X)) +
    geom_violin(width = 0.75, color = alpha(boxplot_col, 0.3), fill = boxplot_col, alpha = 0.3) +
    stat_boxplot(geom = "errorbar", width = 0.5, color = boxplot_col) +
    geom_boxplot(width = 0.75, color = boxplot_col, alpha = 0.2, outlier.alpha = 1.0) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = mean_col, fill = mean_col) +
    scale_y_continuous(expand = c(0.02, 0)) +
    labs(title = NULL, x = NULL, y = y_lab)

  if (!is.null(expand_pts)) {
    bp_plot <- bp_plot + expand_limits(y = expand_pts)
  }
  return(bp_plot)
}

bp_yobs <- my_boxplot(y_train_all, x_tick = "Train set", y_lab = "Y observations")
bp_preds <- my_boxplot(sparboott[, 1], x_tick = "Test set", y_lab = "Quantile predictions (level 0.9999)", expand_pts = c(80, 250))
bp_width <- my_boxplot(sparboott[, 3] - sparboott[, 2], x_tick = "Test set", y_lab = "50% CI widths", expand_pts = c(4, 12))
plot(bp_yobs)
plot(bp_preds)
plot(bp_width)
save_myplot(bp_yobs, paste0(save_path, "Yobs_boxplot.pdf"), width = 25, height = 80)
save_myplot(bp_preds, paste0(save_path, "Qpred_boxplot.pdf"), width = 25, height = 80)
save_myplot(bp_width, paste0(save_path, "CIwidth_boxplot.pdf"), width = 25, height = 80)




cat("\nRun time:\n")
print((Sys.time() - start_time))
