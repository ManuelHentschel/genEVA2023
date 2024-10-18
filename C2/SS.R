################################################################
################################################################
################################################################
### Data challenge EVA 2023
## Simulation Study
################################################################
data <- read.csv('../data/Amaurot.csv')
source('quantile_estimate.R')
library(plotly) 
library(boot)
library(ExtDist)
################################################################

################################################################
## load data
margin <- data$Y
n      <- length(margin)


estimate_loss <- function(k_fold, margin){
  index      <- sample(1:n, n, replace=F )
  loss_value <- NULL
  loss_curve <- NULL 
  loss_validation <- NULL
  for(R in k_fold:1 ){
    ## Extracts data
    kindex_2 <- index[(R-1)*(n/k_fold) +  (1:(n/k_fold))]    ## sort for stationarity
    margin_2 <- margin[kindex_2] 
    index_out<- !1:n%in%kindex_2 
    sorted     <- sort(margin[index_out], decreasing = T)
    n2         <- length(margin_2)
    year_period <- n2/105
    true_approx <- sorted[length(index_out)/(year_period*300)] 
    th <-  data.frame(th1 =c(0.94) ) # data.frame(th1 = c( seq( n2*0.9 , n2*0.99,by=40)/n2)  )
    ## As a function of th
    if(is.null(th$th1))   th <- data.frame(th1 =c(0.94) ) 
    estimates <- sapply(th$th1, function(k) quantile_estimate(k,margin_2, period=year_period ) )
    quantileci        <- sapply(th$th1, function(k) quantile_ci(k,margin,period=year_period,estimates= c(estimates^(.5),estimates^(1.5) )) )
    q_eval <- seq( min(quantileci[1,]), max(quantileci[3,]),length=1000)
    loss   <- sapply(q_eval, function(i) sapply(q_eval, function(k) loss_function(i,k) )) ## loss(yhat,y)
    #heatmap(loss, Rowv = NA, Colv = NA) #
    loss_2 <- sapply(q_eval, function(k) loss_function( true_approx,k ))
    min <- which.min( apply( loss, 1, sum ) )
    print(min)
    ###
    loss_value  <- rbind(loss_value, sapply( ((min-500):(min+100)) , function(k) loss_function(true_approx,q_eval[k] ) ))
    
    ## validation_line 
    loss_curve      <- rbind(loss_curve, apply( loss, 1, sum )/length(q_eval) )
    loss_validation <- rbind(loss_validation, loss_2)
    
    #plot.ts( sapply(1:600, function(k) (mean(loss_value[,k]))))
    
  }
  min_delta <- which.min( sapply(1:600, function(k) (mean(loss_value[,k])))) 
  ##
  year_period <- 200
  estimates <- sapply(th$th1, function(k) quantile_estimate(k,margin_2, period=year_period ) )
  quantileci        <- sapply(th$th1, function(k) quantile_ci(k,margin,period=year_period ) )
  q_eval <- seq( min(quantileci[1,]), max(quantileci[3,]),length=1000)
  loss   <- sapply(q_eval, function(i) sapply(q_eval, function(k) loss_function(i,k) )) ## loss(yhat,y)
  #heatmap(loss, Rowv = NA, Colv = NA) #
  min  <- which.min( apply( loss, 1, sum ) )
  min2 <- (min-500+ min_delta-1)
  
  res <- c('gpd' = estimates, 'fine_tuned' = q_eval[min2])
  
} 
loss_function <- function(q,qhat){
  val <- NA
  if(0.99*q > qhat) val <- 0.9*(0.99*q - qhat)
  if( abs(q-qhat) <= 0.01*q ) val <- 0
  if( 1.01*q < qhat) val <- 0.1*(qhat - 1.01*q)
  return(val)
}


#### Automatic fine-tuning
N <- 3 ## take N larger to reproduce experiments
n <- 7000#21000
rsample   <- function(n) rBurr(n, 1,0.5,2)#rnorm(n))##rt(n,df=4
true      <- qBurr(1-1/(200*300), 1,0.5,2)#qt(1-1/(200*300), df=4)#qnorm(1-1/(200*300))#
k_fold    <-  10#30 
res_SS    <- NULL
loss_SS   <- NULL
for(i in 1:N){
  margin   <- rsample(n)
  res <-estimate_loss(k_fold,margin)
  res_SS   <- rbind(res_SS, res)
  loss_SS  <- rbind(loss_SS, sapply(res, function(k) loss_function(k,true)))
}

boxplot(res_SS/true, ylim = c(0,10))
abline(h=1, lty=2, col ='red')
#boxplot(loss_SS, ylim = c(0,5))
sapply(1:2, function(k) mean(loss_SS[,k]))






