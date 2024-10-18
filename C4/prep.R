
### Prep script for C4
# Called from main.R

# Packages
library(graphicalExtremes)
library(memoise)
library(Kmedians)
library(ggplot2)
library(igraph)

# Output directory for figures
PDF_DIR <- './figures'
if(!dir.exists(PDF_DIR)){
  dir.create(PDF_DIR, recursive = TRUE)
}

# Custom plotting wrapper that saves the plot to a PDF file
myPlot <- function(filename, expr, width=5, height=5){
  envir <- parent.frame()
  expr <- substitute(expr)

  filenamePdf <- paste0(filename, '.pdf')
  pdf(
    file = file.path(PDF_DIR, filenamePdf),
    width = width,
    height = height
  )
  tryCatch(
    eval(expr, envir=parent.frame()),
    finally = dev.off()
  )
  invisible(NULL)
}

# Load data
DATA_DIR <- file.path('..', 'data')

U1 <- read.csv(file.path(DATA_DIR, 'UtopulaU1.csv'))
U2 <- read.csv(file.path(DATA_DIR, 'UtopulaU2.csv'))

U <- as.matrix(cbind(U1, U2))

# Set up constants
PHI1 <- 1/300
PHI2 <- 12 * PHI1
S1 <- 5.702133
S2 <- 3.198534
P1 <- 1-PHI1
P2 <- 1-PHI2

MIXED_THRESHOLDS <- c(
  rep(S1, ncol(U1)),
  rep(S2, ncol(U2))
)
S1_THRESHOLDS <- rep(S1, ncol(U))

n <- ncol(U)

diagMask <- diag(n) == 1
nonDiagMask <- !diagMask
upperTriMask <- upper.tri(diagMask)


# Empirical Coefficient of AD
emp_chi_k <- function(data, p=NULL, k=2, m=1000){
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }
  data <- data > 1
  n <- nrow(data)
  d <- ncol(data)
  combs <- combn(seq_len(d), k)
  m <- min(ncol(combs), m)
  combs <- combs[,seq_len(m), drop=FALSE]
  chis <- numeric(m)
  for(i in seq_len(m)){
    cols <- data[,combs[,i]]
    allLarge <- apply(cols, 1, prod)
    ret <- sum(allLarge) / (1/k * sum(cols))
    chis[i] <- ret
  }
  return(chis)
}

# GPD estimation
sgp <- function(q, loc = 0, scale = 1, shape = 0){
  unname((1 + shape * (q - loc)/scale)^(-1/shape))
}
estProbUsingGPD <- function(x, pGpd, xEst){
  # Fit GPD
  ret <- extRemes::fevd(
    x = x,
    type = 'GP',
    threshold = pGpd,
    initial = NULL
  )
  par <- ret$results$par
  
  # Compute loc of GPD
  loc <- unname(quantile(x, pGpd))
  
  # Compute survival probability under GPD
  p0 <- sgp(xEst, loc=loc, scale=par['scale'], shape=par['shape'])
  
  # Combine with threshold exceedance probability
  p1 <- (1-pGpd) * p0
  
  return(p1)
}


# Memoise functions (makes repeated runs of the script faster)
myCache <- cachem::cache_disk(
  'memoise_cache',
  evict = 'fifo'
)
myMemoise <- function(f){
  if('memoised' %in% class(f)){
    return(f)
  }
  return(memoise::memoise(
    f,
    cache = myCache
  ))
}

# Memoise expensive functions:
emp_chi_k <- myMemoise(emp_chi_k)
emp_chi <- myMemoise(emp_chi)
emp_vario <- myMemoise(emp_vario)
estProbUsingGPD <- myMemoise(estProbUsingGPD)
