library(data.table)
library(copula)
library(mev)

Rcpp::sourceCpp("C3_functions.cpp")

#define Gumbel CDF
Gumbel_cdf <- function(x){
  exp(-exp(-x))
}
#define GPD CDF
GPD_cdf <- function(u, sigma, xi, x){
  if (1+xi*(x-u)/sigma < 0){
    stop('Evaluation point outside of support')
  }
  if (xi == 0){
    return(1-exp(-(x-u)/sigma))
  } else {
    return(1-(1+xi*(x-u)/sigma)^(-1/xi))
  }
}


############# Task C3
#load data
C3_rawdata <- fread("../data/Coputopia.csv")
dim(C3_rawdata)
sum(complete.cases(C3_rawdata)) # check whether there are missing values

### Exploratory data analysis
# Extremal dependence
u_seq = seq(0.9, 0.998, 0.0005)
chi_mat = matrix(NA, nrow=length(u_seq), ncol=3)
eta_mat = matrix(NA, nrow=length(u_seq), ncol=3)
#computed empirical chi and eta for each pair of variables
for (i in 1:length(u_seq)){
  temp = chieta_matrix(as.matrix(C3_rawdata[,-c(1:2)]), u=u_seq[i])
  chi_mat[i,] = c(temp$chi[1,2], temp$chi[1,3], temp$chi[2,3])
  eta_mat[i,] = c(temp$eta[1,2], temp$eta[1,3], temp$eta[2,3])
}
#obtain chi_bar from eta
chibar_mat = eta_mat*2 - 1

pdf("C3_EDA.pdf", height=4, width=8)
par(mfrow=c(1,2), pty="s")
par(mar = c(5, 3, 0.3, 0.2))

plot(u_seq, chi_mat[,1], xlab=expression(u), ylab=expression(chi(u)),
     ylim=c(0.05,0.3), type="l", col=1, lwd=1.5)
lines(u_seq, chi_mat[,2], type="l", col=2, lwd=1.5)
lines(u_seq, chi_mat[,3], type="l", col=3, lwd=1.5)
legend("topright", c(expression(chi[12](u)),expression(chi[13](u)),
                                           expression(chi[23](u))),
       lty=c(1,1,1), col = c(1,2,3), cex = 1)

plot(u_seq, chibar_mat[,1], xlab=expression(u), ylab=expression(bar(chi)(u)),
     ylim=c(0.1,0.6), type="l", col=1, lwd=1.5)
lines(u_seq, chibar_mat[,2], type="l", col=2, lwd=1.5)
lines(u_seq, chibar_mat[,3], type="l", col=3, lwd=1.5)
legend("topleft", c(expression(bar(chi)[12](u)),expression(bar(chi)[13](u)),
                     expression(bar(chi)[23](u))),
       lty=c(1,1,1), col = c(1,2,3), cex = 1)
dev.off()

# Covariates effect
u <- 0.975 
#check empirical extremal dependence coefficients for subsets of data
C3_rawdata_S1 <- subset(C3_rawdata, Season == "S1")
exdep_S1 <- chieta_matrix(as.matrix(C3_rawdata_S1[,-c(1:2)]), u=u)
C3_rawdata_S2 <- subset(C3_rawdata, Season == "S2")
exdep_S2 <-chieta_matrix(as.matrix(C3_rawdata_S2[,-c(1:2)]), u=u)

exdep_S1$chi
exdep_S2$chi
2*exdep_S1$eta-1
2*exdep_S2$eta-1

quant_atm <- quantile(C3_rawdata$Atmosphere, probs = c(0, 0.3, 0.7, 1))
C3_rawdata_Atm_low <- subset(C3_rawdata, Atmosphere <= quant_atm[2])
exdep_low_atm <- chieta_matrix(as.matrix(C3_rawdata_Atm_low[,-c(1:2)]), u=u)
C3_rawdata_Atm_middle <- subset(C3_rawdata, 
                          Atmosphere > quant_atm[2] & Atmosphere <= quant_atm[3])
exdep_middle_atm <- chieta_matrix(as.matrix(C3_rawdata_Atm_middle[,-c(1:2)]), u=u)
C3_rawdata_Atm_high <- subset(C3_rawdata, Atmosphere > quant_atm[3])
exdep_high_atm <- chieta_matrix(as.matrix(C3_rawdata_Atm_high[,-c(1:2)]), u=u)

exdep_low_atm$chi
exdep_middle_atm$chi
exdep_high_atm$chi
2*exdep_low_atm$eta-1
2*exdep_middle_atm$eta-1
2*exdep_high_atm$eta-1


##### Estimating P1 and P2
# P1: fit a Gaussian copula by eta matching
th <- 0.995
const <- Gumbel_cdf(6)
#compute empirical eta
C3.chieta <- chieta_matrix(as.matrix(C3_rawdata[,-c(1:2)]), th)
#specify the Gaussian copula with the same eta
rho.normal <- 2*C3.chieta$eta[upper.tri(C3.chieta$eta)] - 1
ext.normal.fit.1 <- normalCopula(rho.normal, dim=3, dispstr="un")
l1 <- rep(const, 3)
p1.est <- 1-3*const+pCopula(c(const,const), normalCopula(rho.normal[1]))+
  pCopula(c(const,const), normalCopula(rho.normal[2]))+
  pCopula(c(const,const), normalCopula(rho.normal[3]))-
  pCopula(l1, ext.normal.fit.1)

# P2: fit GPD to min(Y1, Y2)|Y3<m
m <- -log(log(2))
#extract data of min(Y1, Y2)|Y3<m
Y3lessm.ind <- C3_rawdata$Y3<m
minY12.Y3lessm <- pmin(C3_rawdata$Y1[Y3lessm.ind], C3_rawdata$Y2[Y3lessm.ind])

th <- 0.95
q.th <- quantile(minY12.Y3lessm, th)
gpdfit <- fit.gpd(minY12.Y3lessm, threshold = q.th)
p2_gpdfit <- (1-GPD_cdf(q.th, as.numeric(gpdfit$param[1]), 
                      as.numeric(gpdfit$param[2]), 7))*(1-th)*0.5
gpdfit$param
#estimated probability 
p2_gpdfit
