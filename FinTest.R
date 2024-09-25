#setwd("~/Desktop/JCGS Codes")
source("Functions.R")

n = 500              #no of vertices
m =  20            #size of common region
k = 5              #number of subgraphs
bs = 50
d = 2                          #rank of X
eps = 0             #alternative level is epsilon

#minimum 100 iterations recommended for level
iter1 = 100
iter = 100
#b = (n-m)/k
alpha = 0.05/k

X = matrix(runif(n*d),n,d)
P = ProbM(X)
P = P/max(P)
diag(P) = 0

Pe = ProbM(X+eps)
Pe = Pe/max(Pe)
diag(Pe) = 0

registerDoParallel(6)

## Without subsampling

power = rep(0,iter)
times = rep(0,iter)
tm1 = proc.time()[3]
for(r in 1:iter)
{
  A1 = graph(P)
  A2 = graph(Pe)
  Bp = testNetwork(A1, A2, m, k, d, bs,"ASE")
  times[r] = Bp$time
  power[r] = Bp$pval
}
ls1 = (iter-length(which(power>0.05)))/iter
tm2 = proc.time()[3]
Tm1 = tm2-tm1


## With subsampling

power = rep(0,iter1)
times = rep(0,iter1)
tm1 = proc.time()[3]
for(r in 1:iter1)
{
  A1 = graph(P)
  A2 = graph(Pe)
  Bp = testNetwork(A1, A2, m, k, d, bs,"SS")
  times[r] = Bp$time
  power[r] = Bp$pval
}
ls2 = (iter1-length(which(power>0.05)))/iter1
tm2 = proc.time()[3]
Tm2 = tm2-tm1


# Multiple Testing

x = 1:n
b = (n-m)/k
common = sample(x, size = m, replace = FALSE)
x = setdiff(x,common)
samp = c()
for(i in 1:k)
{
  samp = rbind(samp, sample(x, size = b, replace = FALSE))
  x = setdiff(x,samp)
}

lvlS = c()
tp1 = proc.time()[3]
for(r in 1:iter1)
{
  A1 = graph(P)
  A2 = graph(Pe)
  tk3 = 0
  pval2 = foreach(j = 1:k, .combine = c) %dopar%
    {
      A10 = A1[c(common,samp[j,]),c(common,samp[j,])]
      A20 = A2[c(common,samp[j,]),c(common,samp[j,])]
      Bp = testNetwork(A10, A20, m, k, d, bs,"ASE")
      Bp$pval
    }
  lvlS = rbind(lvlS, pval2)
}
L1 = c()
for(j in 1:k)
{
  L = (iter1-length(which(lvlS[,j]>alpha)))/iter1
  L1 = c(L1,L)
}
l1 = max(L1)
tp2 = proc.time()[3]
Tp1 = tp2-tp1

stopImplicitCluster()

Powers = c(ls1, ls2, l1)   #power values
Times = c(Tm1, Tm2, Tp1)   #times

result.test = data.frame(Powers = Powers, Times = Times)
rownames(result.test) = c("without ss", "with ss", "multtest")
write.csv(result.test, "FinTest.csv", row.names = F)
