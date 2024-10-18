
### R code for C4

stopifnot(basename(getwd()) == 'C4')

## Source helper functions
source('prep.R')

## Clustering of pairwise interactions
P_CLUSTER <- 0.90
nClusters <- 4 # Chosen after inspection of relevant diagrams

V <- emp_vario(U, p=P_CLUSTER)
chi <- emp_chi(U, p=P_CLUSTER)
coords <- cbind(V[nonDiagMask], chi[nonDiagMask])

d <- dist(coords)
clusters <- hclust(d)
clusterVec <- cutree(clusters, k=nClusters)
cMat <- matrix(0, n, n)
cMat[nonDiagMask] <- clusterVec

centers <- matrix(NA, nClusters, ncol(coords))
for(i in seq_len(nClusters)){
  centers[i,] = colMeans(coords[clusterVec==i,,drop=FALSE])
}


## Figure 5 (left)
plotVvsChi <- function(U, useVario=TRUE){
  if(useVario){
    V <- emp_vario(U, p=P_CLUSTER)
    xlabel <- 'Empirical variogram'
  } else{
    V <- var(U)
    xlabel <- 'Covariance'
  }
  chi <- emp_chi(U, p=P_CLUSTER)
  plot(
    V[nonDiagMask],
    chi[nonDiagMask],
    col=clusterVec,
    xlab = xlabel,
    ylab = 'Coefficient of AD',
    pch = 4
  )
}
myPlot('varioVsChi', plotVvsChi(U))


## Figure 5 (right)
pCount <- 30
pMin <- 0.7
pMax <- 0.999
P_LIST_LINES <- 1-exp(seq(log(1-pMin), log(1-pMax), length.out = pCount))
plotChisPreCluster <- function(U){
  chiList <- lapply(P_LIST_LINES, function(p){
    emp_chi(U, p=p)
  })
  
  chiMats <- list()
  for(i in seq_len(nClusters)){
    cols <- lapply(chiList, function(chi){
      chi[cMat == i & upperTriMask]
    })
    chiMats[[i]] <- do.call(cbind, cols)
  }
  
  plot(NULL, xlim=range(P_LIST_LINES), ylim=c(0,1), ylab='Coefficient of AD', xlab='Threshold')
  for(i in seq_along(chiMats)){
    y <- chiMats[[i]]
    y <- cbind(y, NA)
    x <- rep(c(P_LIST_LINES, NA), nrow(y))
    y <- t(y)
    lines(x, y, 'l', col=i)
  }
}
myPlot('chiVsP_PreCluster', plotChisPreCluster(U))


## Clustering of sites (connected components)
clusterIdDisc <- which.max(centers[,1]) # max variogram => disconnect
adjMat <- 1 * (cMat != clusterIdDisc)

g <- igraph::graph_from_adjacency_matrix(adjMat, mode='undirected')
comps <- igraph::components(g)

vertexClusterId <- comps$membership
nConnComps <- comps$no
connIds <- seq_len(comps$no)
connComps <- lapply(connIds, function(i) which(comps$membership == i))
print('\nconnected components:\n')
print(connComps)

if(nConnComps == 1){
  stop('only one connected component!')
}


## Figure 6 (left)
myPlot(
  'heatmapClustering',
  heatmap(V, symm=TRUE)
)


## Figure 6 (right)
plotIndClusterEntries <- function(){
  chiList <- lapply(P_LIST_LINES, function(p){
    emp_chi(U, p=p)
  })
  cols <- lapply(chiList, function(chi){
    chi[cMat == clusterIdDisc & upperTriMask]
  })
  chiMat <- do.call(cbind, cols)
  y <- cbind(chiMat, NA)
  # y <- matrix(colMeans(y), nrow=1)
  x <- rep(c(P_LIST_LINES, NA), nrow(y))
  y <- t(y)
  plot(
    x,
    y,
    'l',
    col = clusterIdDisc,
    xlab = 'u',
    ylab = 'Coefficient of AD'
  )
  
  lines(P_LIST_LINES, 1-P_LIST_LINES, col=1, lwd = 3)
}
myPlot(
  'disconnectedVsIndependent',
  plotIndClusterEntries()
)


## Figure 7
getSimulExceedanceTable <- function(thresholds, clusterId){
  siteInds <- connComps[[clusterId]]
  isExceedance <- U[,siteInds] > thresholds
  simEx <- rowSums(isExceedance)
  tab <- sapply(seq_along(siteInds), function(i) sum(simEx == i))
  return(tab)
}
plotHistograms <- function(thresholds){
  par0 <- par(
    mfrow=c(length(thresholds), length(connIds)),
    mar = c(2,2,2,2)
  )
  for(sInd in seq_along(thresholds)){
    for(clusterInd in connIds){
      tbl <- getSimulExceedanceTable(thresholds[[sInd]], clusterInd)
      names(tbl) <- seq_along(tbl)
      barplot(
        tbl,
        col=clusterInd
      )
    }
  }
  par(par0)
}
myPlot(
  'jointExceedancesS1AndMixed',
  plotHistograms(list(MIXED_THRESHOLDS, S1_THRESHOLDS)),
  height=3.5,
  width = 7
)


## Compute empirical probabilities per cluster
thresholds <- list(
  rep(S1, ncol(U)),
  MIXED_THRESHOLDS
)
empFreqs <- matrix(NA, nrow=length(thresholds), ncol=length(connIds))

for(sInd in seq_along(thresholds)){
  for(clusterInd in connIds){
    ind <- U[,connComps[[clusterInd]]] > thresholds[[sInd]][clusterInd]
    empFreqs[sInd, clusterInd] <- sum(apply(ind, 1, prod))
  }
}
empProbs <- empFreqs / nrow(U)


## Compute minima per cluster
connComps1 <- lapply(connComps, function(ids) ids[ids <= ncol(U1)])
connComps2 <- lapply(connComps, function(ids) ids[ids > ncol(U1)])
splitComps <- c(connComps1, connComps2)

splitMinima <- matrix(0, nrow=nrow(U), ncol=0)

for(i in seq_along(splitComps)){
  cols <- splitComps[[i]]
  splitMinima <- cbind(splitMinima, apply(U[,cols], 1, min, na.rm=TRUE))
}

minima <- matrix(0, nrow=nrow(U), ncol=0)
for(i in seq_along(connComps)){
  cols <- connComps[[i]]
  minima <- cbind(minima, apply(U[,cols], 1, min, na.rm=TRUE))
}


## Fit GPDs to minima
# For threshold S1
lp0 <- log(1-0.97)
lp1 <- log(1-0.9999)
pVec <- 1 - exp(seq(lp0, lp1, length.out=10))

pEstsS1 <- matrix(0, length(pVec), nConnComps)
for(i in seq_along(pVec)){
  for(j in seq_len(nConnComps)){
    pEstsS1[i, j] <- estProbUsingGPD(minima[,j], pVec[i], S1)
  }
}
estProbsS1 <- colMeans(pEstsS1, na.rm=TRUE)

# For mixed thresholds (only the more extreme part)
lp0 <- log(1-0.9)
lp1 <- log(1-0.9999)
pVec <- 1 - exp(seq(lp0, lp1, length.out=10))

pEstsS1only <- matrix(0, length(pVec), nConnComps)
for(i in seq_along(pVec)){
  for(j in seq_len(nConnComps)){
    pEstsS1only[i, j] <- estProbUsingGPD(splitMinima[,j], pVec[i], S1)
  }
}

estProbsS1only <- colMeans(pEstsS1only, na.rm=TRUE)

# Compute chis between halfs of each cluster
clusterHalfChis <- matrix(0, length(pVec), nConnComps)
for(i in seq_along(pVec)){
  p <- pVec[i]
  clusterChis <- emp_chi(splitMinima, p)
  clusterHalfChis[i,] <- diag(clusterChis[1:nConnComps, (nConnComps+1):(2*nConnComps)])
}
medClusterHalfChis <- apply(clusterHalfChis, 2, median)


## Final computations:
# p2:
useEmp <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
empProbsS1 <- empProbs[1,]
p2 <- prod(empProbsS1[useEmp]) * prod(estProbsS1[!useEmp])

# p1:
useEmp <- c(FALSE, TRUE, TRUE, FALSE, FALSE)
mixedProbs <- numeric(5)

empProbsSmixed <- empProbs[1,]
mixedProbs[useEmp] <- empProbsSmixed[useEmp]

# Fall back to all>S1 for NAs:
isNan <- is.na(estProbsS1only)
mixedProbs[isNan & !useEmp] <- estProbsS1[isNan & !useEmp]

# Use chi as substitute for P(...>S2 | ...>S1)
ind <- !useEmp & !isNan
mixedProbs[ind] <- medClusterHalfChis[ind] * estProbsS1only[ind]

p1 <- prod(mixedProbs)

AnswerC4 <- c(p1, p2)
print(AnswerC4)

cat('Done.\n')

