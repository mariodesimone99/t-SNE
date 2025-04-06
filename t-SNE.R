rm(list=ls())
library(scatterplot3d)
library(beans)

PointPointDist = function(X){
  term = as.matrix(dist(X,method = "euclidean", diag=TRUE, upper=TRUE))
  return(term)
}

Perplexity = function(p_i){
  term = 0
  for (i in 1:length(p_i)) {
    if(p_i[i] != 0){
      term = term + p_i[i]*log2(p_i[i])
    }
  }
  term = 2^(-term)
  return(term)
}

PerplexitySigmas = function(dist_mtx, ppx, eps = 10^-10, iter=1000, lb=10^-20, ub=10000){
  N = nrow(dist_mtx)
  sigmas = rep(0, N)
  #for each point
  for(i in 1:N){
    upper = ub
    lower = lb
    #research
    for (j in 1:iter) {
      sigma = (upper+lower)/2
      p_i = P_ij(dist_mtx[i,], sigma)
      ppx_found = Perplexity(p_i)
      
      if(abs(ppx_found-ppx)<=eps){
        sigmas[i] = sigma
        break
      }
      if(ppx_found>ppx){
        upper=sigma
      }else{
        lower=sigma
      }
      if(j==iter){
        sigmas[i] = sigma
      }
    }
  }
  return(sigmas)
}

P_ij = function(distances, sigmas){
  p = exp(-distances/(2*(sigmas^2)))
  if(is.vector(distances)){
    p[1] = 0
    p=p/sum(p)
  }else{
    diag(p) = 0
    p=p/rowSums(p)
  }
  return(p)
}

SimilarityQ = function(Y){
  dist_mtx = PointPointDist(Y)
  dist_mtx = dist_mtx^2
  q=(1+dist_mtx)^-1
  diag(q) = 0
  q=q/sum(q)
  return(q)
}

SimilarityP = function(X, ppx){
  N = nrow(X)
  distances = PointPointDist(X)
  distances = distances^2
  sigmas = PerplexitySigmas(distances, ppx)
  p = P_ij(distances, sigmas)
  p = (p +t(p))/(2*N)
  return(p)
}

Gradient = function(P, Q, Y){
  N = nrow(Y)
  grad = matrix(0, nrow=N, ncol=ncol(Y))
  dist_mtx = PointPointDist(Y)
  dist_mtx = dist_mtx^2
  for (i in 1:N) {
    term = 0
    for (j in 1:N) {
      if(i != j){
        pdiff = P[i,j]-Q[i,j]
        ydiff = Y[i,]-Y[j,]
        ysum = Y[i,]+Y[j,]
        term=term+(pdiff)*(ydiff)*((1+sum(ysum*ysum))^-1)
      }
    }
    grad[i,] = 4*term
  }
  return(grad)
}

Momentum = function(iter){
  m=0
  if(iter < 250){
    m = 0.5
  }else{
    m=0.8
  }
  return(m)
}

tSNE = function(X, dim, iter, ppx=40, lr=500){
  N = nrow(X)
  Y = array(rep(0, 2*N*dim), c(2, N, dim))
  set.seed(1)
  init_data = rnorm(N*dim, 0, 10^-4)
  Y[1,,] = init_data
  Y[2,,] = init_data
  P = SimilarityP(X, ppx)
  for(i in 1:iter){
    Q = SimilarityQ(Y[2,,])
    grad = Gradient(P, Q, Y[2,,])
    tmp = Y[2,,]
    Y[2,,] = Y[2,,] -lr*grad + Momentum(i)*(Y[2,,]-Y[1,,])
    Y[1,,] = tmp
    if(!(i %% 50)){
      print(paste(i, "iteration"))
    }
  }
  return(Y[2,,])
}

Quality = function(Y, labels){
  dists = PointPointDist(Y)
  levels = as.matrix(labels)
  l = unique(labels)
  for (i in 1:length(l)) {
    levels[which(labels==l[i])]=i
  }
  
  all_levels = combn(unique(levels), 2)
  min_dists = array(rep(0, 3))
  max_dists = min_dists
  for (i in 1:ncol(all_levels)) {
    idx1 = which(levels==all_levels[1, i])
    idx2 = which(levels==all_levels[2, i])
    min_dists[i] = min(dists[idx1,idx2])
    max_dists[i] = max(dists[idx1,idx2])
  }
  return(list("min_dists"=min_dists, "max_dists"=max_dists))
}

set.seed(1)

#First Datset
data(beans)
k=1
nsample =  100
sample_rows = NA
classes = NA
lev = NA
for(i in unique(beans$class)){
  sample_rows = c(sample_rows, sample(which(beans$class==i), nsample))
  print(paste(i))
  classes = c(classes, rep(i, nsample))
  lev = c(lev, rep(k, nsample))
  if(k==3){
    break
  }
  k=k+1
}
lev = lev[-1]
sample_rows=sample_rows[-1]
classes=classes[-1]

Y = tSNE(beans[sample_rows, -17],2,1000, ppx=30)
colors = rainbow(length(unique(classes)))
names(colors) = unique(classes)
par(mgp=c(2.5,1,0))
plot(Y, t='n', main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2, "cex.lab"=1.5)
text(Y, labels=lev, col=colors[classes])
q = Quality(Y, classes)

Y_3d = tSNE(beans[sample_rows, -17],3,1000, ppx=30)
s3d=scatterplot3d(Y_3d, color=colors[classes])
legend(s3d$xyz.convert(0,-60,-100), legend=unique(classes),pch=19, yjust=0, 
       col=seq_along(unique(classes)))
q = Quality(Y_3d, classes)

#Second Dataset
data(iris)
lev = NA
k=0
for (i in unique(iris$Species)) {
  lev = c(lev, rep(k, sum(iris$Species==i)))
  k=k+1
}
lev = lev[-1]
Y1 = tSNE(iris[,-5],2,1000, ppx=30)
colors = rainbow(length(unique(iris$Species)))
names(colors) = unique(iris$Species)
par(mgp=c(2.5,1,0))
plot(Y1, t='n', main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2, "cex.lab"=1.5)
text(Y1, labels=lev, col=colors[iris$Species])
q = Quality(Y1, iris$Species)

Y1_3d = tSNE(iris[,-5],3,1000, ppx=30)
s3d=scatterplot3d(Y1_3d, color=colors[iris$Species])
legend(s3d$xyz.convert(100,-30,-150), legend=unique(iris$Species),pch=19, yjust=0, 
       col=seq_along(unique(iris$Species)))
q = Quality(Y1_3d, iris$Species)