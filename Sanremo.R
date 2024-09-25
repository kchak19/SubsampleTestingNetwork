#setwd("~/Desktop/JCGS Codes")
source("Functions.R")

Sanremo = read.table("Sanremo2016_final_multiplex.edges", quote="\"", comment.char="")
Interactions = filter(Sanremo, V4>=2)

Retweets = filter(Interactions, V1==1)
Mentions = filter(Interactions, V1==2)
Replies = filter(Interactions, V1==3)

lowWInd = fread("lowWeightIndex.csv", header=T)
highWInd = setdiff(1:56562, data.matrix(lowWInd))
n = length(highWInd)

node.number = function(node) {which(highWInd == node)}

registerDoParallel(6)

A = matrix(0,n,n)
for(i in 1:nrow(Retweets)){
  v1 = node.number(Retweets[i,3])
  v2 = node.number(Retweets[i,2])
  A[v1,v2] = 1
  A[v2,v1] = 1
}
diag(A) = 0
sum2(A)/2

B = matrix(0,n,n)
for(i in 1:nrow(Mentions)){
  v1 = node.number(Mentions[i,3])
  v2 = node.number(Mentions[i,2])
  B[v1,v2] = 1
  B[v2,v1] = 1
}
diag(B) = 0
sum2(B)/2

C = matrix(0,n,n)
for(i in 1:nrow(Replies)){
  v1 = node.number(Replies[i,3])
  v2 = node.number(Replies[i,2])
  C[v1,v2] = 1
  C[v2,v1] = 1
}
diag(C) = 0
sum2(C)/2

m = 1080
k = 25
d = 2
bs = 20

# with subsampling

Ret.Men = testNetwork(A, B, m, k, d, bs, method = "SS")
#Men.Rep = testNetwork(B, C, m, k, d, bs, method = "SS")
#Ret.Rep = testNetwork(A, C, m, k, d, bs, method = "SS")

# multiple testing

Ret.Men.mult = multTestNetwork(A, B, m, k, d, bs)
#Men.Rep = multTestNetwork(B, C, m, k, d, bs)
#Ret.Rep = multTestNetwork(A, C, m, k, d, bs)

pvals = c(Ret.Men$pval, Ret.Men.mult$pval)
times = c(Ret.Men$time, Ret.Men.mult$time)

stopImplicitCluster()

result = data.frame(pvals = pvals, times = times)
rownames(result) = c("with ss", "multtest")
write.csv(result, "SanremoSS.csv", row.names = F)


