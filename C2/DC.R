################################################################
################################################################
################################################################
### Data challenge EVA 2023
## Challenge 2 implementation
################################################################
data <- read.csv('../data/Amaurot.csv')
source('quantile_estimate.R')
library(plotly) 
library(latex2exp)
library(mev)
library(fExtremes)
library(gmp)
################################################################

################################################################
## load data
margin <- data$Y
n      <- length(margin)
## Plot tail-index
th1 <- 0.95
n2 <- length(data$Y)
################################################################
## first estimate
qu  <- quantile_estimate(th1,margin,plot=T,return_fit=T)
fit <- qu$fit
u   <- qu$fit$threshold
q0  <- qu$q0 ## ~188 ## First estimate 
##
##
################################################################
## Histogram
x2 <- seq(0.01, max(margin), length = 100)
fun <- evd::dgpd(x2, scale = fit$estimate[1], shape =  fit$estimate[2]) # gpd curve
hist( margin[margin > u] - u , prob = TRUE, col = "white",
     ylim = c(0, max(fun)), breaks = 50,
     main = "Histogram with GPD fit") # Histogram
lines(x2, fun, col = 2, lwd = 2)

################################################################
## As a function of th
################################################################
th     <- c( seq( n*0.9 , n*0.99,by=100)/n)
estimates <- sapply(th, function(k) quantile_estimate(k,margin) )
plot(th,estimates, type = "l") ## for different th.
abline(h=q0, col = 'red')
plot.ts(margin)
abline(h = estimates )
abline(h=q0, col = 'red')
####################################
loss_function <- function(q,qhat){
  val <- NA
  if(0.99*q > qhat) val <- 0.9*(0.99*q - qhat)
  if( abs(q-qhat) <= 0.01*q ) val <- 0
  if( 1.01*q < qhat) val <- 0.1*(qhat - 1.01*q)
  return(val)
}
#################################################################
#################################################################
#################################################################
#################################################################
#### To compute profile like ci.
## Extracts data
year_period <- n/105
## As a function of th
estimates <- sapply(th, function(k) quantile_estimate(k,margin, period=year_period ) )
plot.ts(th,estimates, ylim = c(120,220)) ## plot of estimates
abline(h=q0, col = 'blue', lty = 1)
quantileci        <- sapply(th, function(k) quantile_ci(k,margin,period=year_period ) )
#plot.ts(quantile_ci)
q_eval <- seq( min(quantileci[1,]), max(quantileci[3,]),length=50)
loss   <- sapply(q_eval, function(i)  sapply(q_eval, function(k) loss_function(i,k) ))
##################################################################
#################################################################
## Ploting the loss function

x <- q_eval
y <- q_eval## estimates
data <- expand.grid(x = x, y = y)
data$loss <- apply(data,1, function(k) loss_function(k[1],k[2]) )# as.vector(loss)

gg1 <- ggplot(data, aes(x = x, y = y, fill = loss)) +
  geom_tile() +
  #geom_contour(aes(z = loss ), color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "C") +
  labs(title = "Loss function",
       x = TeX('${q}_T$'),
       y = TeX('$\\widehat{q}_T$')) +
  theme_minimal()
gg1 

################################################################
################################################################
################################################################
################################################################
################################################################
## Validation
## Let T be the ear length
## I'm looking for the event that happens once every 300 hundred years. 
## I have 105 years of training
## For doing 10-fold cross validation, I have sample of sizes 2100
## T = 20 
## first estimate
k_fold     <-  30 
th     <- c( seq( n*0.9 , n*0.99,by=100)/n)
index      <- sample(1:n, n, replace=F )
loss_value <- NULL
loss_curve <- NULL 
loss_validation <- NULL

for(R in k_fold:1){
  ## Extracts data
  kindex_2 <- index[(R-1)*(n/k_fold) +  (1:(n/k_fold))]    ## sort for stationarity
  margin_2 <- margin[kindex_2] 
  index_out<- !1:n%in%kindex_2 
  sorted     <- sort(margin[index_out], decreasing = T)
  n2         <- length(margin_2)
  year_period <- n2/105
  true_approx <- sorted[length(index_out)/(year_period*300)] 
  ## As a function of th
  if(is.null(th))   th <- data.frame(th1 =c(0.95) ) 
  estimates <- sapply(th, function(k) quantile_estimate(k,margin_2, period=year_period ) )
  #plot.ts(th,estimates, ylim = c(120,220))
  #abline(h=true_approx, col = 'blue', lty = 1)
  
  ### profile_like confidence intervals
  quantileci        <- sapply(th, function(k) quantile_ci(k,margin,period=year_period ) )
  #plot.ts(quantile_ci)
  q_eval <- seq( min(quantileci[1,]), max(quantileci[3,]),length=1000)
  #q_eval <- seq( quantile_ci$normal[2], quantile_ci$normal[3],length=1000)
  loss   <- sapply(q_eval, function(i) sapply(q_eval, function(k) loss_function(i,k) )) ## loss(yhat,y)
  #heatmap(loss, Rowv = NA, Colv = NA) #
  par(mfrow=c(3,1))
  ## plot data-based loss vs. real loss
  plot.ts( q_eval, apply( loss, 1, sum )/length(q_eval) , type = "l", 
           main = "estimated loss vs. true loss minimization",
           ylim = c(0,10), ylab = 'loss') ## sum_y loss(yhat,y)
  loss_2 <- sapply(q_eval, function(k) loss_function( true_approx,k ))
  lines(q_eval, loss_2 , lty = 2, col = "blue")
  min <- which.min( apply( loss, 1, sum ) )
  points(q_eval[min],apply( loss, 1, sum )[min]/length(q_eval)  )
  print(min)
  ###
  ## plot estimation vs. true
  plot.ts(margin, ylab = 'estimate', main  = 'quantile estimation')
  abline(h=q_eval[min], col = 'red', lty = 1)
  abline(h=true_approx, col = 'blue', lty = 2)
  print(loss_function( true_approx,q_eval[min] ))
  
  #fine_tune   <- sort( apply( loss[ (min:(min+50)) , ], 1, sum ), index.return = T)
  loss_value  <- rbind(loss_value, sapply( ((min-700):(min+100)) , function(k) loss_function(true_approx,q_eval[k] ) ))
 
  ## validation_line 
  loss_curve      <- rbind(loss_curve, apply( loss, 1, sum )/length(q_eval) )
  loss_validation <- rbind(loss_validation, loss_2)
  
  ## validation line
  plot.ts( sapply(1:800, function(k) (mean(loss_value[,k]))), 
           main = 'Parameter tuning curve',
           ylab = 'loss' , 
           xlab = 'index')
  
}


lambda <- which.min( sapply(1:800, function(k) (mean(loss_value[,k])))) 
## answer : q_eval[min - 700 + lambda] , letting R <- 1 in the above code

