#setwd("~/Desktop/JCGS Codes")
source("Functions.R")

dat = fread("friendfeed_ita.mpx", header = F)
friend = as.data.frame(dat)

likes = filter(friend, V3=="like")
comments = filter(friend, V3=="comment")
follows = filter(friend, V3=="follow")

n=21006

A.lk = A.cm = A.fl = matrix(0,n,n)

N = paste0("A",1:n)

for(i in 1:nrow(likes)){
  v1 = likes[i,1]
  v2 = likes[i,2]
  A.lk[match(v1,N),match(v2,N)] = 1
  A.lk[match(v2,N),match(v1,N)] = 1
}
diag(A.lk)=0
sum2(A.lk)/2

for(i in 1:nrow(comments)){
  v1 = comments[i,1]
  v2 = comments[i,2]
  A.cm[match(v1,N),match(v2,N)] = 1
  A.cm[match(v2,N),match(v1,N)] = 1
}
diag(A.cm)=0
sum2(A.cm)/2

for(i in 1:nrow(follows)){
  v1 = follows[i,1]
  v2 = follows[i,2]
  A.fl[match(v1,N),match(v2,N)] = 1
  A.fl[match(v2,N),match(v1,N)] = 1
}
diag(A.fl)=0
sum2(A.fl)/2

registerDoParallel(6)

m = 1006            #size of common region
k = 25              #number of subgraphs
b = (n-m)/k
d = 2                #rank of X
bs = 20           #no of booststrap resamples


# with subsampling

Lik.Com = testNetwork(A.lk, A.cm, m, k, d, bs, method = "SS")
#Com.Fol = testNetwork(A.cm, A.fl, m, k, d, bs, method = "SS")
#Lik.Fol = testNetwork(A.lk, A.fl, m, k, d, bs, method = "SS")

# multiple testing

Lik.Com.mult = multTestNetwork(A.lk, A.cm, m, k, d, bs)
#Com.Fol.mult = multTestNetwork(A.cm, A.fl, m, k, d, bs)
#Lik.Fol.mult = multTestNetwork(A.lk, A.fl, m, k, d, bs)

pvals = c(Lik.Com$pval, Lik.Com.mult$pval)
times = c(Lik.Com$time, Lik.Com.mult$time)

stopImplicitCluster()

result = data.frame(pvals = pvals, times = times)
rownames(result) = c("with ss", "multtest")
write.csv(result, "FriendFeedSS.csv", row.names = F)


