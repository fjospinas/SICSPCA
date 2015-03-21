#rm(list = ls())
library(mirt)
library(FactoMineR)
data = as.matrix(read.table(file = '/home/mirt/Comparaciones_Montenegro/3PL/Datasets/Test_10_1_1000.csv',header=T))

setwd("/home/mirt/Git/GrupoSICS/dev/SICSPCA/")
#dyn.load("pasoe3.so")
dyn.load("pasom3.so")

#Función de probabilidad
gg = function(a,d, cp,  theta){
  exp(cp)/(1+exp(cp))+ (1-(exp(cp)/(1+exp(cp))))*(1 + exp(-(a*theta+ d)))^(-1)
}

#Probabilidad del patrón
ggpat = function(pat,zita,nitems,cuad){
  p = rep(0,nitems)
  for(i in 1:nitems){
    p[i] = gg(a = zita[i,1],d = zita[i,2],cp = zita[i,3],theta = cuad)
    if(pat[i] == 0){
      p[i] = 1 - p[i] 
    }
  }
  prod(p)
}

#Funcion que calcula scores por EAP
scoresEAP = function(pats = pats,zita = zita,w.cuad = w.cuad,thetaG = thetaG){
  
  npats = nrow(pats)
  nitems = ncol(pats) -1
  
  ncuads = length(thetaG)
  
  thetaEst = rep(0,npats)
  for(j in 1:npats){
    sumNum = 0
    sumDen = 0
    for(k in 1:ncuads){
      pat = pats[j,][1:nitems]
      sumNum = sumNum + (thetaG[k] * w.cuad[k] * ggpat(pat = pat,zita = zita,nitems = nitems,cuad = thetaG[k]))
      sumDen = sumDen + (w.cuad[k] * ggpat(pat = pat,zita = zita,nitems = nitems,cuad = thetaG[k]))
    }
    thetaEst[j] = sumNum / sumDen
  }
  
  ret = cbind(pats[,-ncol(pats)],thetaEst)
  colnames(ret) = c(sapply(X = 1:nitems,FUN = function(x){paste("Item",x,sep = " ")}),"Score")
  ret
}

#valoes iniciales de Andrade modificados para que retornen d (en vez de b)
#calculan el c basado en los thetas inferiores al b que responden el item
#el c lo restorna en la escala [0,1]
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

#Función que calcula Patrones
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

#Función que indexa patrones para devolverl los scores una vez calculados 
#sobre los patrones
indexPat = function(data,pats){
  comprimData = apply(data,MARGIN=1,FUN=paste,collapse="/")
  comprimPats = apply(pats[,1:ncol(data)],MARGIN=1,FUN=paste,collapse="/")
  index = lapply(comprimPats,FUN = function(x) {which(x == comprimData)})
  index
}

#Logverosimilitu escrita en Cpp
LL2 = function(zita.vec,R,fvec,pt.cuad,nitems){
  .Call("Loglik",zita.vec,R,fvec,pt.cuad)
}

#Gradiente de logverosimilutud escrito en Cpp
gradLoglik = function(zita.vec,R,fvec,pt.cuad,nitems){
  .Call("grad",zita.vec,R,fvec,pt.cuad)
}

#Calcula el candidato de cada grupo de la partición
#en este caso es la mediana
calculoThetaG = function(theta,index){
  sapply(X = 1:ngroups,FUN = function(x){median(theta[index[[x]]])})
}

#Calculo directo de eigen (No funciona)
#Z = (data -  apply(X = data,MARGIN = 2,mean)) / apply(X = data,MARGIN = 2,sd)
#Z2 = (t(Z) %*% Z) / 1000
#eig = eigen(Z2)
#eig

#variables que no cambian durante la estimación
pats = patrones(data)
nitems = ncol(data)
indexPat = indexPat(data = data,pats = pats)

#número de grupos 
ngroups = 40

#Thetas iniciales
pca.math <- PCA(data,ncp=3)
pca.math$eig[1,1]
Ycoor = pca.math$ind$coord[,1]
theta = (Ycoor - mean(Ycoor))/sd(Ycoor)

#Valores iniciales de andrade
ini.vals = start.andrade.d(datos = data,thetas = theta)
ini.vals[,3] = qlogis(ini.vals[,3]) 
zita = ini.vals

#cuantiles para particionar el theta
cut = seq(0,ngroups) / ngroups

###################
# Bucle principal #
###################

#inicialmente se ponen ciclos fijos para determinar paso a paso la convergencia
for(ciclos in 1:100){
  
  #Se definen los limites de la partición
  limits = quantile(theta,cut)
  limits[1] = limits[1] - .1
  limits[ngroups + 1] = limits[ngroups +1] + .1
  
  #La variable Index es una lista que conserva los indices de la partición
  #es decir si la partición es de 40 grupos la lista tiene 40 vectores,
  #el vector 1 son los indices de los individuos que pertenecen a el grupo 1 y así sucesivamente.
  #se usa una lista debido a que los grupos no son necesariamente del mismo tamaño y por ende los 
  #vectores de indices tampoco
  index = list()
  
  #Se consyruye f y se llena Index
  f = numeric(ngroups)
  for(i in 1:ngroups){
    indG = list(unique(which((theta > limits[i]) & (theta <= limits[i+1]))))
    index = append(index,indG)
    f[i] = length(indG[[1]])
  }

  #Se construye matriz R
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

  #Se calculan los candidatos de cada grupo (la mediana)
  thetaG = calculoThetaG(theta,index)

  #Se transforma zita en un vector para poder ser usado por el optimizador
  zita.vec = as.vector(zita)
  opt = optim(par=zita.vec,fn=LL2,gr=gradLoglik,method= "BFGS",R=r,fvec=f,pt.cuad=thetaG,nitems=nitems,control=list(maxit=20),hessian = T)

  #Se imprime la logverosimilitud evaluada en el optim para ver cómo va la convergencia
  print(opt$value)

  #Se transfocma zita nuevamente en una matriz y se acota
  zita = matrix(opt$par,nrow=nitems,byrow=F)
  zita[,1] = ifelse(abs(zita[,1]) > 10, ini.vals[,1], zita[,1])
  zita[,2] = ifelse(abs(zita[,2]) > 40, ini.vals[,2], zita[,2])
  zita[,3] = ifelse(abs(zita[,3]) > 600, ini.vals[,3], zita[,3])
  print(zita)
  
  #Se calculan los pesos de cada grupo basado en los indices
  w.cuad = sapply(1:length(index),function(X) length(index[[X]]))
  w.cuad = w.cuad / nrow(data)

  #Se calculan los thetas para los patrones
  patsTheta = scoresEAP(pats = pats,zita = zita,w.cuad = w.cuad,thetaG = thetaG)

  #se expanden los thetas a la matriz original
  for(mm in 1:nrow(patsTheta)){
    theta[indexPat[[mm]]] = patsTheta[mm,ncol(data) +1]        
  }
  
  #se imprime el summary de los thetas para ver como va
  print(summary(theta))
}
