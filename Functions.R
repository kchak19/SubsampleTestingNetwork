# doParallel is required for computation using parallel cores
# Use more parallel cores and large graphs for better results
require(doParallel)
require(data.table)
require(matrixStats)
require(dplyr)

# graph generation function.
graph = function(P)
{
  n = nrow(P)
  A = matrix(0, nrow = n, ncol = n)
  A[col(A) > row(A)] = runif(n*(n-1)/2)
  A = (A + t(A))
  A = (A < P) + 0 ;
  diag(A) = 0
  #A = Matrix(A,sparse=TRUE)
  return(A)
}

# function for estimation : no subsampling, basic ASE.
ASE = function(A, dim)
{
  A.eig = eigen(A)
  A.val = A.eig$values[1:dim]
  A.vec = A.eig$vectors[,1:dim]
  if(dim == 1)
    A.coords = sqrt(A.val) * A.vec
  else
    A.coords = A.vec %*% diag(sqrt(A.val))

  return(A.coords)
}


# procrustes transformation.
procr= function(X,Y)
{
  tmp = t(X) %*% Y
  tmp.svd = svd(tmp)
  W = tmp.svd$u %*% t(tmp.svd$v)
  return(list(d = norm(X%*%W - Y, type = "F"), W = W))
}

ProbM = function(X) {X %*% t(X)}

# function for estimation with subsampling
estSS = function(A,m,k,d, overlap = "random"){
  n = nrow(A)
  b = (n-m)/k
  # Choice of overlap
  if(overlap == "dense"){
    S = apply(A, 1, sum)
    common = tail(order(S),m)
    x = head(order(S),n-m)
    samp = c()
    for(i in 1:k)
    {
      samp = rbind(samp, sample(x, size = b, replace = FALSE))
      x = setdiff(x,samp)
    }
  }
  if(overlap == "random"){
    x = 1:n
    common = sample(x, size = m, replace = FALSE)
    x = setdiff(x,common)
    samp = c()
    for(i in 1:k)
    {
      samp = rbind(samp, sample(x, size = b, replace = FALSE))
      x = setdiff(x,samp)
    }
  }
  
  # reference subsample
  ind1 = c(common,samp[1,])
  A.ref = A[ind1,ind1]
  X.hat.ref = ASE(A.ref,d)
  X.ref.0 = X.hat.ref[1:m,]
  X.ref.1 = X.hat.ref[(m+1):(m+b),]
  
  # other subsamples
  #IMPORTANT : Need to create a parallel cluster with sufficeint cores first
  X.part = foreach(i = 2:k, .combine = rbind) %dopar%
    {
      A.sub = A[c(common,samp[i,]),c(common,samp[i,])]
      X.hat.sub = ASE(A.sub,d)
      X.sub.0 = X.hat.sub[1:m,]
      X.sub.i = X.hat.sub[(m+1):(m+b),]
      
      H = procr(X.sub.0, X.ref.0)$W
      X.sub.i.trans = X.sub.i %*% H
      X.sub.i.trans
    }
  
  X.fin = matrix(0, n, d)
  
  # combine
  for(j in 1:length(common)) {X.fin[common[j],] = X.ref.0[j,]}
  for(j in 1:length(samp[1,])) {X.fin[samp[1,j],] = X.ref.1[j,]}
  for(i in 2:k) for(j in 1:length(samp[i,])) {X.fin[samp[i,j],] = X.part[b*(i-2)+j,]}
  
  return(X.fin)
}

EST = function(A,m,k,d, method = "SS"){
  if(method == "ASE") return(ASE(A,d))
  if(method == "SS") return(estSS(A,m,k,d))
}


# bootstrap testing
testNetwork = function(A, B, m, k, d, bs, method = "SS"){
  Xhat = EST(A,m,k,d, method)
  Yhat = EST(B,m,k,d, method)
  P1hat = ProbM(Xhat)
  P2hat = ProbM(Yhat)
  TS = norm(P1hat - P2hat, type = "F")
  S = rep(NA,bs)
  t1 = proc.time()[3]
  for(boot in 1:bs)
  {
    cat(boot)
    P = (P1hat + P2hat)/2
    A1 = graph(P)
    B1 = graph(P)
    Xb = EST(A1,m,k,d, method)
    Yb = EST(B1,m,k,d, method)
    Q1 = ProbM(Xb)
    Q2 = ProbM(Yb)
    f = norm(Q1-Q2, type = "F")
    S[boot] = f
  }
  t2 = proc.time()[3]
  T1 = (t2-t1)
  pval = (sum(S > TS))/bs
  return(list("pval" = pval, "time" = T1))
}


# multiple testing
multTestNetwork = function(A, B, m, k, d, bs){
  n = nrow(A)
  b = (n-m)/k
  x = 1:n
  common = sample(x, size = m, replace = FALSE)
  x = setdiff(x,common)
  
  samp = c()
  for(i in 1:k)
  {
    samp = rbind(samp, sample(x, size = b, replace = FALSE))
    x = setdiff(x,samp)
  }
  T1 = proc.time()[3]
  pval2 = foreach(j = 1:k, .combine=c) %dopar%
    {
      A10 = A[c(common,samp[j,]),c(common,samp[j,])]
      B10 = B[c(common,samp[j,]),c(common,samp[j,])]
      n0 = nrow(A10)
      Bp = testNetwork(A10, B10, m, k, d, bs,"ASE")
      Bp$pval
    }
  pval = max(pval2)
  T2 = proc.time()[3]
  ps = list("pval" = pval, "time" = T2-T1)
  return(ps)
}

