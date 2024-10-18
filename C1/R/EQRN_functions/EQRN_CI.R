EQRN_predict_CI_paraboot <- function(fit_eqrn, X, prob_lvl_predict, confidence, intermediate_quantiles, interm_lvl = fit_eqrn$interm_lvl,
                                     boot_samples = 1000, boot_size = fit_eqrn$n_excesses, CI_type = c("percentile", "basic"),
                                     parallel_strat = c("sequential", "multisession", "multicore"), n_workers = (availableCores() - 1),
                                     err_handling = c("stop", "remove", "pass"), seed = NULL) {
  CI_type <- match.arg(CI_type)
  parallel_strat <- match.arg(parallel_strat)
  err_handling <- match.arg(err_handling)

  GPD_params_pred <- EQRN_predict_params(fit_eqrn, X, intermediate_quantiles,
    return_parametrization = "classical", interm_lvl
  )
  sigmas <- GPD_params_pred$scales
  xis <- GPD_params_pred$shapes

  n_pred <- length(sigmas)
  nb_prob_lvls_predict <- length(prob_lvl_predict)
  if (nb_prob_lvls_predict != 1) {
    stop("Not implemented error: 'nb_prob_lvls_predict!=1'")
  }

  predicted_quantiles <- GPD_quantiles(prob_lvl_predict, interm_lvl, intermediate_quantiles, sigmas, xis)

  n_workers <- min(n_workers, n_pred)
  `%fun%` <- set_doFuture_strategy(parallel_strat, n_workers = n_workers)

  # if(!is.null(seed)){set.seed(seed)}

  boot_CIs <- foreach(i = 1:n_pred, .errorhandling = err_handling, .combine = rbind) %fun% {
    # Parametric bootstrap resampling
    gpd_samples_i <- matrix(evd::rgpd(boot_samples * boot_size, loc = 0, scale = sigmas[i], shape = xis[i]),
      nrow = boot_samples, ncol = boot_size
    )
    # Quantile estimates
    quants_i <- apply(
      gpd_samples_i, 1,
      function(x) {
        predict_quantiles_from_excesses(x,
          quantiles = prob_lvl_predict,
          interm_lvl = interm_lvl, interm_quantile = intermediate_quantiles[i]
        )
      }
    )
    if (CI_type == "percentile") {
      CI <- matrix(quantile(c(quants_i), probs = c((1 - confidence) / 2, (1 + confidence) / 2), names = FALSE),
        nrow = 1, ncol = 2
      )
    } else {
      stop("Not implemented error 'CI_type'")
      CI <- matrix(c(CId, CIu), nrow = 1, ncol = 2)
    }
    CI
  }

  end_doFuture_strategy()


  return(list(
    predicted_quantiles = predicted_quantiles, scales = sigmas, shapes = xis,
    CIs = boot_CIs
  ))
}



make_prediction_CIs <- function(predicted_quantiles, boot_preds, confidence, CI_type = c("basic", "percentile", "bc", "normal", "bcall"),
                                debug_verb = FALSE) {
  CI_type <- match.arg(CI_type)

  ex_bool <- foreach(i = 1:nrow(boot_preds), .errorhandling = "stop", .combine = c) %do% {
    quants_i <- c(boot_preds[i, ])
    pred_i <- c(predicted_quantiles)[i]
    quants_i > pred_i
  }
  bias_all <- mean(ex_bool)
  if (debug_verb) {
    print(c(bias_all, qnorm(bias_all)))
  }

  boot_CIs <- foreach(i = 1:nrow(boot_preds), .errorhandling = "stop", .combine = rbind) %do% {
    # Quantile estimates
    quants_i <- c(boot_preds[i, ])
    pred_i <- c(predicted_quantiles)[i]
    alpha2 <- (1 - confidence) / 2

    boot_qs_i <- quantile(c(quants_i), probs = c(alpha2, (1 - alpha2)), names = FALSE)
    # boot_qs <- quantile(c(quants_i), probs=c((1-confidence)/2, (1+confidence)/2), names=FALSE)

    if (CI_type == "percentile") {
      CI <- matrix(boot_qs_i, nrow = 1, ncol = 2)
    } else if (CI_type == "basic") {
      CI <- matrix(c(2 * pred_i - boot_qs_i[2], 2 * pred_i - boot_qs_i[1]), nrow = 1, ncol = 2)
    } else if (CI_type == "bc") {
      b_bias <- qnorm(mean(quants_i > pred_i))
      q1 <- pnorm(qnorm(alpha2) - 2 * b_bias)
      q2 <- pnorm(qnorm(1 - alpha2) - 2 * b_bias)
      bc_boot_qs_i <- quantile(c(quants_i), probs = c(q1, q2), names = FALSE)
      if (debug_verb) {
        print(c(mean(quants_i > pred_i), b_bias, q1, q2))
      }
      CI <- matrix(bc_boot_qs_i, nrow = 1, ncol = 2)
    } else if (CI_type == "bcall") {
      q1 <- pnorm(qnorm(alpha2) - 2 * qnorm(bias_all))
      q2 <- pnorm(qnorm(1 - alpha2) - 2 * qnorm(bias_all))
      bc_boot_qs_i <- quantile(c(quants_i), probs = c(q1, q2), names = FALSE)
      if (debug_verb) {
        print(c(q1, q2))
      }
      CI <- matrix(bc_boot_qs_i, nrow = 1, ncol = 2)
    } else if (CI_type == "normal") {
      swidCI <- qnorm(1 - alpha2) * sd(quants_i)
      CI <- matrix(c(pred_i - swidCI, pred_i + swidCI), nrow = 1, ncol = 2)
    } else {
      stop("Not implemented error 'CI_type'")
      CI <- matrix(c(CId, CIu), nrow = 1, ncol = 2)
    }
    CI
  }
  return(boot_CIs)
}


bootstrap_inds <- function(n, boot_samples = 1000, boot_size = n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  samples_ind <- list()
  for (i in seq(boot_samples)) {
    samples_ind[[i]] <- sample.int(n = n, size = boot_size, replace = TRUE)
  }
  return(samples_ind)
}


semi_param_GPD_bootstrap_inds <- function(y, intermediate_quantiles, scales, shapes,
                                          boot_samples = 1000, boot_size = length(c(y)), seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  samples_ind <- bootstrap_inds(n = length(c(y)), boot_samples = boot_samples, boot_size = boot_size, seed = NULL)

  for (inds_i in seq_along(samples_ind)) {
    inds <- samples_ind[[inds_i]]
    yb <- y[inds]
    intqb <- intermediate_quantiles[inds]
    scb <- scales[inds]
    shb <- shapes[inds]

    for (i in seq_along(inds)) {
      if (yb[i] > intqb[i]) {
        yb[i] <- evd::rgpd(1, loc = intqb[i], scale = scb[i], shape = shb[i])
      }
    }
    samples_ind[[inds_i]] <- list(inds = inds, yb = yb)
  }
  return(samples_ind)
}







#' Predict unconditional extreme quantiles using peaks over threshold
#'
#' @param interm_lvl Probability level at which the empirical quantile should be used as the intermediate threshold.
#' @param quantiles Probability levels at which to predict the extreme quantiles.
#' @param interm_quantile
#' @param Y Vector of ("training") observations.
#'
#' @return
#' @export
#'
#' @examples
predict_quantiles_from_excesses <- function(Y, quantiles = c(0.99), interm_lvl, interm_quantile) {
  p0 <- interm_lvl
  t0 <- interm_quantile
  pars <- ismev::gpd.fit(Y, 0, show = FALSE, maxit = 1e6)$mle
  sigma <- pars[1]
  xi <- pars[2]

  q_hat <- GPD_quantiles(quantiles, p0, t0, sigma, xi)

  return(c(q_hat))
}



#' @title EQRN fit function for independent data
#'
#' @description Use the \code{\link{EQRN_fit_restart}} wrapper instead, with \code{data_type="iid"}, for better stability using fitting restart.
#'
#' @param fit_eqrn
#' @param X Matrix of covariates, for training.
#' @param y Response variable vector to model the extreme conditional quantile of, for training.
#' @param intermediate_quantiles Vector of intermediate conditional quantiles at level \code{interm_lvl}.
#' @param learning_rate Initial learning rate for the optimizer during training of the neural network.
#' @param L2_pen L2 weight penalty parameter for regularization during training.
#' @param shape_penalty Penalty parameter for the shape estimate, to potentially regularize its variation from the fixed prior estimate.
#' @param n_epochs Number of training epochs.
#' @param batch_size Batch size used during training.
#' @param X_valid Covariates in a validation set, or \code{NULL}.
#' Used for monitoring validation loss during training, enabling learning-rate decay and early stopping.
#' @param y_valid Response variable in a validation set, or \code{NULL}.
#' Used for monitoring validation loss during training, enabling learning-rate decay and early stopping.
#' @param quant_valid Intermediate conditional quantiles at level \code{interm_lvl} in a validation set, or \code{NULL}.
#' Used for monitoring validation loss during training, enabling learning-rate decay and early stopping.
#' @param lr_decay Learning rate decay factor.
#' @param patience_decay Number of epochs of non-improving validation loss before a learning-rate decay is performed.
#' @param min_lr Minimum learning rate, under which no more decay is performed.
#' @param patience_stop Number of epochs of non-improving validation loss before early stopping is performed.
#' @param tol Tolerance for stopping training, in case of no significant training loss improvements.
#' @param patience_lag The validation loss is considered to be non-improving if it is larger than on any of the previous \code{patience_lag} epochs.
#' @param optim_met DEPRECATED. Optimization algorithm to use during training. \code{"adam"} is the default.
#'
#' @return An EQRN object of classes \code{c("EQRN_iid", "EQRN")}, containing the fitted network,
#' as well as all the relevant information for its usage in other functions.
#' @export
#'
#' @examples
EQRN_continue_training <- function(fit_eqrn, X, y, intermediate_quantiles, learning_rate = 1e-5, L2_pen = 0, shape_penalty = 0, n_epochs = 100, batch_size = 256,
                                   X_valid = NULL, y_valid = NULL, quant_valid = NULL, lr_decay = 1, patience_decay = n_epochs, min_lr = 0,
                                   patience_stop = n_epochs, tol = 1e-6, patience_lag = 1, optim_met = "adam") {
  network <- fit_eqrn$fit_nn # $clone(deep=TRUE)
  interm_lvl <- fit_eqrn$interm_lvl
  intermediate_q_feature <- fit_eqrn$intermediate_q_feature
  X_scaling <- fit_eqrn$X_scaling
  orthogonal_gpd <- fit_eqrn$orthogonal_gpd
  scale_features <- X_scaling$scaling

  fit_eqrn$train_loss_prev <- c(fit_eqrn$train_loss_prev, fit_eqrn$train_loss)
  fit_eqrn$valid_loss_prev <- c(fit_eqrn$valid_loss_prev, fit_eqrn$valid_loss)
  fit_eqrn$n_obs_prev <- c(fit_eqrn$n_obs_prev, fit_eqrn$n_obs)
  fit_eqrn$n_excesses_prev <- c(fit_eqrn$n_excesses_prev, fit_eqrn$n_excesses)
  fit_eqrn$excesses_ratio_prev <- c(fit_eqrn$excesses_ratio_prev, fit_eqrn$excesses_ratio)

  data_excesses <- get_excesses(
    X = X, y = y, quantiles = intermediate_quantiles,
    intermediate_q_feature = intermediate_q_feature, scale_features = scale_features, X_scaling = X_scaling
  )
  Y_excesses <- torch_tensor(data_excesses$Y_excesses, device = device)
  X_feats_excesses <- torch_tensor(data_excesses$X_excesses, device = device)
  # X_scaling <- data_excesses$X_scaling

  # Data Loader
  n_train <- Y_excesses$size()[1]
  trainset <- tensor_dataset(X_feats_excesses, Y_excesses)
  trainloader <- dataloader(trainset, batch_size = batch_size, shuffle = TRUE)

  # Validation dataset (if everything needed is given)
  do_validation <- (!is.null(y_valid) & !is.null(X_valid) & !is.null(quant_valid))
  if (do_validation) {
    data_valid_excesses <- get_excesses(
      X = X_valid, y = y_valid, quantiles = quant_valid,
      intermediate_q_feature = intermediate_q_feature, scale_features = scale_features, X_scaling = X_scaling
    )
    y_valid_ex <- torch_tensor(data_valid_excesses$Y_excesses, device = device)
    X_valid_ex <- torch_tensor(data_valid_excesses$X_excesses, device = device)

    n_valid <- y_valid_ex$size()[1]
    validset <- tensor_dataset(X_valid_ex, y_valid_ex)
    validloader <- dataloader(validset, batch_size = batch_size, shuffle = FALSE)
  }

  # Semi-conditional GPD fit (on y rescaled excesses wrt intermediate_quantiles)
  if (shape_penalty > 0) {
    semicond_gpd_fit <- fit_GPD_unconditional(c(data_excesses$Y_excesses), interm_lvl = NULL, thresh_quantiles = NULL)
  } else {
    semicond_gpd_fit <- NULL
  }

  # Instantiate network
  Dim_in <- ncol(X_feats_excesses)

  # Optimizer
  # optimizer <- setup_optimizer(network, learning_rate, L2_pen, hidden_fct, optim_met=optim_met)
  optimizer <- torch::optim_adam(network$parameters, lr = learning_rate, weight_decay = L2_pen)

  # Train the network
  network$train()

  loss_log_train <- rep(as.double(NA), n_epochs)
  nb_stable <- 0
  if (do_validation) {
    loss_log_valid <- rep(as.double(NA), n_epochs)
    nb_not_improving_val <- 0
    nb_not_improving_lr <- 0
  }
  curent_lr <- learning_rate
  for (e in 1:n_epochs) {
    loss_train <- 0
    coro::loop(for (b in trainloader) {
      b <- fix_dimsimplif(b, Dim_in) # TODO: remove when fixed in torch
      # Forward pass
      net_out <- network(b[[1]])
      # Loss
      loss <- loss_GPD_tensor(net_out, b[[2]],
        orthogonal_gpd = orthogonal_gpd,
        shape_penalty = shape_penalty, prior_shape = semicond_gpd_fit$shape, return_agg = "mean"
      )

      # zero the gradients accumulated in buffers (not overwritten in pyTorch)
      optimizer$zero_grad()
      # Backward pass: compute gradient of the loss with respect to model parameters
      loss$backward()
      # make one optimizer step
      optimizer$step()
      # store loss
      loss_train <- loss_train + (b[[2]]$size()[1] * loss / n_train)$item()
    })
    # Log loss
    loss_log_train[e] <- loss_train
    if (do_validation) {
      network$eval()
      loss_valid <- 0
      coro::loop(for (b in validloader) {
        b <- fix_dimsimplif(b, Dim_in) # TODO: remove when fixed in torch
        valid_out <- network(b[[1]])
        loss <- loss_GPD_tensor(valid_out, b[[2]],
          orthogonal_gpd = orthogonal_gpd,
          shape_penalty = 0, prior_shape = NULL, return_agg = "mean"
        )
        loss_valid <- loss_valid + (b[[2]]$size()[1] * loss / n_valid)$item()
      })
      loss_log_valid[e] <- loss_valid
      # Is the validation loss improving ?
      if (e > patience_lag) {
        if (any(is.na(loss_log_valid[(e - patience_lag):e]))) {
          if (is.nan(loss_log_valid[e])) {
            nb_not_improving_val <- nb_not_improving_val + 1
            nb_not_improving_lr <- nb_not_improving_lr + 1
            cat("NaN validation loss at epoch:", e, "\n")
          }
        } else {
          if (loss_log_valid[e] > (min(loss_log_valid[(e - patience_lag):(e - 1)]) - tol)) {
            nb_not_improving_val <- nb_not_improving_val + 1
            nb_not_improving_lr <- nb_not_improving_lr + 1
          } else {
            nb_not_improving_val <- 0
            nb_not_improving_lr <- 0
          }
          if (abs(loss_log_train[e] - loss_log_train[e - 1]) < tol) {
            nb_stable <- nb_stable + 1
          } else {
            nb_stable <- 0
          }
        }
      }
      # Learning rate decay
      if (curent_lr > min_lr & nb_not_improving_lr >= patience_decay) {
        optimizer <- decay_learning_rate(optimizer, lr_decay)
        curent_lr <- curent_lr * lr_decay
        nb_not_improving_lr <- 0
      }
      if (nb_not_improving_val >= patience_stop) {
        cat(
          "Early stopping at epoch:", e, ", average train loss:", loss_log_train[e],
          ", validation loss:", loss_log_valid[e], ", lr=", curent_lr, "\n"
        )
        break
      }
      network$train()
    } else {
      # If no validation
      if (abs(loss_log_train[e] - loss_log_train[e - 1]) < tol) {
        nb_stable <- nb_stable + 1
      } else {
        nb_stable <- 0
      }
    }
    if (nb_stable >= patience_stop) {
      cat("Early tolerence stopping at epoch:", e, ", average train loss:", loss_log_train[e])
      if (do_validation) {
        cat(", validation loss:", loss_log_valid[e])
      }
      cat(", lr=", curent_lr, "\n")
      break
    }
    # Print progess
    if (e %% 100 == 0 || e == 1) {
      cat("Epoch:", e, "out of", n_epochs, ", average train loss:", loss_log_train[e])
      if (do_validation) {
        cat(", validation loss:", loss_log_valid[e], ", lr=", curent_lr)
      }
      cat("\n")
    }
  }

  network$eval()

  fit_eqrn$fit_nn <- network
  fit_eqrn$train_loss <- loss_log_train[1:e]
  fit_eqrn$n_obs <- length(y)
  fit_eqrn$n_excesses <- n_train
  fit_eqrn$excesses_ratio <- data_excesses$excesses_ratio

  if (do_validation) {
    fit_eqrn$valid_loss <- loss_log_valid[1:e]
  }
  # class(fit_eqrn) <- c("EQRN_iid", "EQRN")

  return(fit_eqrn)
}
