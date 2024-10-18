library(tidyverse)
library(evd)
library(ismev)
library(torch)

library(future)
library(doFuture)

library(lubridate)
library(tsibble) # Tidy Temporal Data Frames and Tools
library(feasts) # Feature Extraction and Statistics for Time Series
library(EnvStats)
library(VGAM)

library(latex2exp)
library(egg)
# library(coro)
library(ggpubr)
library(scales)

# devtools::install_github("opasche/EQRN")
library(EQRN)

# coro, mvtnorm, randtoolbox, # stats, tools

# if (cuda_is_available()) {
#   device <- torch_device("cuda")
# } else {
#   device <- torch_device("cpu")
# }

source("R/EQRN_functions/EQRN_CI.R")

source("R/plotting_functions/plot_utils.R")
source("R/plotting_functions/plot_helpers.R")
source("R/plotting_functions/plot_helpers_ts.R")
