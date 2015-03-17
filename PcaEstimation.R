#rm(list = ls())
library(mirt)
library(FactoMineR)
data = as.matrix(read.table(file = '/home/mirt/Comparaciones_Montenegro/3PL/Datasets/Test_10_1_1000.csv',header=T))

setwd("/home/mirt/Git/GrupoSICS/dev/SICSPCA/")
dyn.load("pasoe3.so")
dyn.load("pasom3.so")

start.andrade.d = function(datos,thetas){
  I = ncol(datos)
  P = nrow(datos)
  m = 5
  
  #scores
  T = apply(datos,1,sum)
  
  #correlacion biserial
  corr.bis = rep(NA,I) 
  for(i in 1:I){
    corr.bis[i]  = cor(datos[,i],T,method="pearson")    
  } 
  #a inicial
  a.ini = sqrt(corr.bis^2/(1-(corr.bis^2)))
  #proporcion de aciertos
  pi = as.vector(apply(datos,2,sum) / P)
  #b inicial
  b.ini = -(qnorm(pi) / corr.bis)
  d.ini = -b.ini * a.ini
  
  #c inicial
  matB = matrix(b.ini,nrow = P,ncol = I,byrow = T)
  matTheta = matrix(theta,nrow = P,ncol= I,byrow = F)
  c.ini = colMeans(matTheta < matB & data == 1)
  
  ini.andrade = matrix(c(a.ini,d.ini,c.ini),ncol=3)
  colnames(ini.andrade) = c("a","d","c")
  ini.andrade
}

patrones = function(datos){
  items = ncol(datos)
  comprim = apply(datos,MARGIN=1,FUN=paste,collapse="/")
  freq = table(comprim)
  pats = names(freq)
  freq = as.vector(freq)
  pats = as.numeric(unlist(lapply(pats,FUN=strsplit,split="/")))
  pats = matrix(pats,ncol=items,byrow=T)
  pats = cbind(pats,freq)
  pats
}

indexPat = function(data,pats){
  comprimData = apply(data,MARGIN=1,FUN=paste,collapse="/")
  comprimPats = apply(pats[,1:ncol(data)],MARGIN=1,FUN=paste,collapse="/")
  index = lapply(comprimPats,FUN = function(x) {which(x == comprimData)})
  index
}

LL2 = function(zita.vec,R,fvec,pt.cuad,nitems){
  .Call("Loglik",zita.vec,R,fvec,pt.cuad)
}

gradLoglik = function(zita.vec,R,fvec,pt.cuad,nitems){
  .Call("grad",zita.vec,R,fvec,pt.cuad)
}

calculoThetaG = function(theta,index){
  sapply(X = 1:ngroups,FUN = function(x){median(theta[index[[x]]])})
}

#Z = (data -  apply(X = data,MARGIN = 2,mean)) / apply(X = data,MARGIN = 2,sd)
#Z2 = (t(Z) %*% Z) / 1000
#eig = eigen(Z2)
#eig

pats = patrones(data)
nitems = ncol(data)
indexPat = indexPat(data = data,pats = pats)


ngroups = 40
pca.math <- PCA(data,ncp=3)
pca.math$eig[1,1]
Ycoor = pca.math$ind$coord[,1]
theta = (Ycoor - mean(Ycoor))/sd(Ycoor)
ini.vals = start.andrade.d(datos = data,thetas = theta)
ini.vals[,3] = qlogis(ini.vals[,3]) 
zita = ini.vals
cut = seq(0,ngroups) / ngroups

for(pp in 1:100){
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
  if(length(grupoAct) == 1){
    print(grupoAct)
    r[i, ] = dataGrupo
  }else{
    r[i, ] = colSums(dataGrupo)
  }  
}



thetaG = calculoThetaG(theta,index)


  zita.vec = as.vector(zita)

  opt = optim(par=zita.vec,fn=LL2,gr=gradLoglik,method= "BFGS",R=r,fvec=f,pt.cuad=thetaG,nitems=nitems,control=list(maxit=20),hessian = T)
print(opt$value)

zita = matrix(opt$par,nrow=nitems,byrow=F)
zita[,1] = ifelse(abs(zita[,1]) > 10, ini.vals[,1], zita[,1])
zita[,2] = ifelse(abs(zita[,2]) > 40, ini.vals[,2], zita[,2])
zita[,3] = ifelse(abs(zita[,3]) > 600, ini.vals[,3], zita[,3])
print(zita)
w.cuad = sapply(1:length(index),function(X) length(index[[X]]))
w.cuad = w.cuad / nrow(data)

patsTheta = scoresEAP(pats = pats,zita = zita,w.cuad = w.cuad,thetaG = thetaG)

for(mm in 1:nrow(patsTheta)){
  theta[indexPat[[mm]]] = patsTheta[mm,ncol(data) +1]        
}
print(summary(theta))
}