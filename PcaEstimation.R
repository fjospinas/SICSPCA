rm(list = ls())
library(mirt)
library(FactoMineR)
data = as.matrix(read.table(file = '/home/mirt/Comparaciones_Montenegro/3PL/Datasets/Test_10_1_1000.csv',header=T))

#Z = (data -  apply(X = data,MARGIN = 2,mean)) / apply(X = data,MARGIN = 2,sd)
#Z2 = (t(Z) %*% Z) / 1000
#eig = eigen(Z2)
#eig
ngroups = 20
pca.math <- PCA(data,ncp=3)
pca.math$eig[1,1]
Ycoor = pca.math$ind$coord[,1]
theta = (Ycoor - mean(Ycoor))/sd(Ycoor)
cut = seq(0,ngroups) / ngroups
limits = quantile(theta,cut)
limits[1] = limits[1] - .1
limits[ngroups + 1] = limits[ngroups +1] + .1
length(limits)
index = list()
f = numeric(ngroups)
for(i in 1:ngroups){
  indG = list(unique(which((theta > limits[i]) & (theta <= limits[i+1]))))
  index = append(index,indG)
  f[i] = length(indG[[1]])
}

r = matrix(0,nrow = ngroups,ncol = ncol(data))
for(i in 1:ngroups){
  grupoAct = index[[i]]
  dataGrupo = data[grupoAct,]
  r[i, ] = colSums(dataGrupo)
}
r
data[index[[1]],]
