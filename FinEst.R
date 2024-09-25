#setwd("~/Desktop/JCGS Codes")
source("Functions.R")

n = 5000              #no of vertices
ms =  c(40,100,200,500)            #size of common region
ks = c(5,10,20)              #number of subgraphs

d = 2                #rank of X
repl = 1            #number of replicates


registerDoParallel(6)

N1 = T1 = c()

time2 = var2 = matrix(0,length(ms),length(ks))
norm2 = c()


X = matrix(runif(n*d),n,d)
P = X %*% t(X)
P = P/max(P)
diag(P) = 0


for(r1 in 1:repl)
{
  A = graph(P)
  m = ms[1]; k = ks[1]
  # Without sub-sampling
  t1 = proc.time()[3]
  X.hat = EST(A,m,k,d,"ASE")
  t2 = proc.time()[3]
  P.hat = ProbM(X.hat)
  
  N1 = c(N1, norm(P-P.hat,type = "F")/norm(P, type="F"))
  T1 = c(T1, t2-t1)
  
  # With subsampling
  N2 = T2 = V2 = matrix(0,length(ms),length(ks))
  
  for(i in 1:length(ms))
    for(j in 1:length(ks)){
      
      m = ms[i]
      k = ks[j]
      b = (n-m)/k
      
      n2 = Time = c()
      for(r2 in 1:repl)
      {
        t3 = proc.time()[3]
        X.fin = EST(A,m,k,d,"SS")
        t4 = proc.time()[3]
        P.fin = ProbM(X.fin)
        Time = c(Time, t4-t3)
        n2 = c(n2, norm(P-P.fin,type = "F")/norm(P, type= "F"))
      }
      
      Time = mean(Time)
      n2.mean = mean(n2)
      n2.var = repl*var(n2)
      
      
      N2[i,j] = n2.mean
      T2[i,j] = Time
      V2[i,j] = n2.var
    }
  
  var2 = var2+V2
  norm2 = append(norm2, list(N2))
  time2 = time2+T2
}

V1 = var(N1)
N1 = mean(N1)
T1 = mean(T1)

T2 = time2/repl
N2 = Reduce("+",norm2)/repl
vb = sapply(1:(length(ms)*length(ks)), function(i) {
  vec = sapply(norm2, '[[', i)
  var(vec)*repl
})

dim(vb) = c(length(ms),length(ks))
var2 = (var2+vb)/repl

stopImplicitCluster()


rownames(N2) = paste0("overlap","=",ms)
colnames(N2) = paste0("NumSubsamp","=",ks)

rownames(T2) = paste0("overlap","=",ms)
colnames(T2) = paste0("NumSubsamp","=",ks)

rownames(var2) = paste0("overlap","=",ms)
colnames(var2) = paste0("NumSubsamp","=",ks)

result = list(WithoutSubsampErr=N1, WithoutSubsampTime=T1, WithoutSubsampVar=V1,
              WithSubsampErr=N2, WithSubsampTime=T2, WithSubsampVar=var2)

sink("FinEst.txt")
print(result)
sink()

