# Credits: Some of the imputation procedures are adapted or inspired from a tutorial by Julie Josse and "R-miss-tastic"
# https://rmisstastic.netlify.app/
# http://juliejosse.com/


library(tidyverse)
library(Amelia)
library(mice)
library(missForest)
library(missMDA)
library(MASS)
library(softImpute)
library(devtools)
# source("../../R/EQRN_functions/utils.R")
library(EQRN)


get_feat_mat <- function(tbl) {
  X <- as.matrix(tbl %>% mutate(Season = as.double(Season)))
  return(X)
}


dat <- read_csv("../../data_wrangled/Amaurot_train.csv") %>% mutate(Season = factor(as.integer(factor(Season, levels = c("S1", "S2"))), levels = c(1, 2)))
test <- read_csv("../../data_wrangled/Amaurot_test.csv") %>% mutate(Season = factor(as.integer(factor(Season, levels = c("S1", "S2"))), levels = c(1, 2)))
dat


colSums(is.na(dat))


print(c(nrow(dat), nrow(drop_na(dat)), nrow(dat) - nrow(drop_na(dat))))
print((nrow(dat) - nrow(drop_na(dat))) / nrow(dat))

alldat <- bind_rows(dat, test) %>%
  rename(Y_feat = Y) %>%
  mutate(Season = as.double(Season))

alldat_sc <- scale(alldat)

X_scaling <- list(
  center_all = attr(alldat_sc, "scaled:center"),
  scale_all = attr(alldat_sc, "scaled:scale"),
  center = attr(alldat_sc, "scaled:center")[-1],
  scale = attr(alldat_sc, "scaled:scale")[-1]
)

safe_save_rds(X_scaling, "../../data_wrangled/other/X_scaling.rds")

X_na <- dat %>%
  rename(Y_feat = Y) %>%
  mutate(Season = as.double(Season)) %>%
  scale(center = X_scaling$center_all, scale = X_scaling$scale_all)
# X_na

train_na_sd <- bind_cols(X_na, dat["Y"])
X_test_sd <- test %>%
  mutate(Season = as.double(Season)) %>%
  scale(center = X_scaling$center, scale = X_scaling$scale) %>%
  as_tibble()
train_na_sd %>% write_csv("../../data_wrangled/train_na_sd.csv")
X_test_sd %>% write_csv("../../data_wrangled/X_test_sd.csv")

## softImpute


# perform softImpute
sft <- softImpute(x = X_na, rank.max = 2, lambda = 0, type = c("als", "svd"))


# compute the factorization
X.sft <- sft$u %*% diag(sft$d) %*% t(sft$v)
# replace missing values by computed values
X.sft[which(!is.na(X_na))] <- X_na[which(!is.na(X_na))]


source("https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/impute/CrossValidation_softImpute.R")
lambda_sft <- cv_sft(X_na)
sft <- softImpute(x = X_na, lambda = lambda_sft)
X.sft <- sft$u %*% diag(sft$d) %*% t(sft$v)
X.sft[which(!is.na(X_na))] <- X_na[which(!is.na(X_na))]
colnames(X.sft) <- colnames(X_na)
bind_cols(X.sft, dat["Y"]) %>% write_csv("../../data_wrangled/imputations/X_sft.csv")
head(X.sft)


## mice



mice_mice <- mice(data = X_na, m = 5, method = "pmm") # contains m=5 completed datasets.
# mice::complete(mice_mice, 1) #get back the first completed dataset of the five available in mice_res


methods(mice)


IMP <- 0
for (i in 1:5) {
  IMP <- IMP + mice::complete(mice_mice, i)
}
X.mice <- IMP / 5 # 5 is the default number of multiple imputations
bind_cols(X.mice, dat["Y"]) %>% write_csv("../../data_wrangled/imputations/X_mice.csv")
head(X.mice)



## missForest


forest <- missForest(xmis = X_na, maxiter = 20, ntree = 200) # , parallelize="forests")


X.forest <- forest$ximp
bind_cols(X.forest, dat["Y"]) %>% write_csv("../../data_wrangled/imputations/X_forest.csv")
head(X.forest)


## missMDA


pca <- imputePCA(X = X_na, ncp = 2, scale = TRUE, method = c("Regularized", "EM"))


ncp.pca <- estim_ncpPCA(X_na, method.cv = "gcv")$ncp
pca <- imputePCA(X_na, ncp = ncp.pca)
X.pca <- pca$comp
bind_cols(X.pca, dat["Y"]) %>% write_csv("../../data_wrangled/imputations/X_pca.csv")
head(X.pca)



# Numerical experiments to compare the different methods

## Synthetic data



MSE_miss <- function(X, Xtrue, mask) {
  return(sqrt(sum((as.matrix(X) * mask - as.matrix(Xtrue) * mask)^2) / sum(mask)))
}



source_url("https://raw.githubusercontent.com/R-miss-tastic/website/master/static/how-to/generate/amputation.R")
HowToImpute <- function(X, perc.list, mecha.list, nbsim) {
  perc_mecha.matrix <- matrix(perc.list, nrow = length(mecha.list) * length(perc.list), ncol = 2)
  perc_mecha.matrix[, 2] <- as.vector(sapply(mecha.list, rep, length(perc.list)))

  results.all <- apply(perc_mecha.matrix, 1, function(perc_mecha) {
    perc <- as.numeric(perc_mecha[1])
    mecha <- perc_mecha[2]

    results.couple <- lapply(1:nbsim, function(iter) {
      XproduceNA <- produce_NA(as.matrix(X), mechanism = mecha, perc.missing = perc)

      XNA <- as.matrix(as.data.frame(XproduceNA$data.incomp))

      ## Mean
      X.mean <- imputeMean(XNA)

      ## MICE
      temp <- mice(XNA, printFlag = FALSE, method = "pmm", remove.collinear = FALSE) # for the predictive mean matching method
      IMP <- 0
      for (i in 1:5) {
        IMP <- IMP + mice::complete(temp, i)
      }
      X.mice <- IMP / 5 # 5 is the default number of multiple imputations

      ## PCA
      ncp.pca <- estim_ncpPCA(XNA)$ncp
      pca <- imputePCA(XNA, ncp = ncp.pca)
      X.pca <- pca$comp

      ## SoftImpute
      lambda_sft <- cv_sft(XNA)
      sft <- softImpute(x = XNA, lambda = lambda_sft, rank.max = min(10, ncol(XNA) - 1))
      X.sft <- sft$u %*% diag(sft$d) %*% t(sft$v)
      X.sft[which(!is.na(XNA))] <- XNA[which(!is.na(XNA))]

      ## RandomForest
      forest <- missForest(XNA, verbose = FALSE)
      X.forest <- forest$ximp


      mse <- sapply(list(X.pca, X.forest, X.mice, X.sft, X.mean), MSE_miss, Xtrue = as.data.frame(X), mask = is.na(XNA))

      cbind.data.frame(mse)
    })

    results <- Reduce("+", results.couple) / length(results.couple)
    rownames(results) <- c("X.pca", "X.forest", "X.mice", "X.soft", "X.mean")
    return(results)
  })

  names(results.all) <- paste0(perc_mecha.matrix[, 1], " ", perc_mecha.matrix[, 2])

  resdf <- as.data.frame(results.all)
  colnames(resdf) <- paste0(perc_mecha.matrix[, 1], " ", perc_mecha.matrix[, 2])

  return(resdf)
}



perc.list <- c(0.1, 0.3, 0.5)
mecha.list <- c("MCAR", "MAR", "MNAR")
res <- HowToImpute(as.matrix(drop_na(as_tibble(X_na))), perc.list = perc.list, mecha.list = mecha.list, nbsim = 2)


res
