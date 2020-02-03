library(astrolibR) # gcirc
source('Init_Fct_Set_detect_sub_v3.R')

fct_Distsphe2 <- function(RaDec_deg){
  n=length(RaDec_deg[,1])
  Distsphe <- matrix(nrow = n, ncol =n)
  for (ind in seq(from=1, to=n-1)){
    for (jind in seq(from=ind+1, to=n)){
      Distsphe[ind, jind] <- astrolibR::gcirc(2,RaDec_deg[ind,1],RaDec_deg[ind,2],RaDec_deg[jind,1],RaDec_deg[jind,2])
      Distsphe[jind,ind]<-Distsphe[ind, jind] 
    }
    Distsphe[ind, ind]<-NA
  }
  Distsphe[n, n]<-NA
  return(Distsphe/3600) # routine gcirc gives a value back in arcsec --> convert here in degree}
  
}

MstDist<- function(distMatrix){
  nedge<-length(distMatrix[1,])-1
  v1<-1
  matAux=distMatrix[v1,]
  v2<-which.min(matAux)
  d<-matAux[v2]
  tree<-cbind.data.frame(v1,v2,d)  
  intree<-c(v1,v2)
  while (length(intree) <= nedge) {
    w1<-which(distMatrix==min(distMatrix[intree,-intree],na.rm = TRUE),arr.ind = TRUE)[1,]
    if(is.element(w1[1], intree)){
      v1<-w1[[1]]
    }else{
      v1<-w1[[2]]
    }
    matAux=distMatrix[v1,]
    matAux[intree]<-NA
    v2<-which.min(matAux)

    d<-matAux[v2]
    temp<-cbind.data.frame(v1,v2,d)
    tree<-rbind.data.frame(tree,temp)
    intree<-c(intree,v2)
  }
  return(tree)
}

paramQ<- function(distMatrix){
 
  s <- mean(distMatrix, na.rm = TRUE)
  
  area<-1
  req<-sqrt(area/pi)
  
  
  sbar=s/req
  np=length(distMatrix[,1])
  ms<-MstDist(distMatrix)
  mdMST<-mean(ms$d)
  print(mdMST)
  mbar=mdMST/(sqrt(area*np)/(np-1))
  print(mbar)
  
  q=mbar/sbar
  
  
  return(cbind.data.frame(q, mdMST, mbar, s, area, req, sbar))
}

regionDensity2D <- function(distMatrix){ 
  nneigh <- fct_nns(distMatrix,6)$knd
  meanDist <- mean(nneigh)
  densFromMean <- 1.151082*5/(pi*meanDist^2)

  return(densFromMean)
}
calculateEps2D<-function(nn1, rho){

funOPC<- function(r){
  densEmp<-density(nn1$knd, from = 0)
  fdens<-approxfun(densEmp)
  return(fdens(r)/ProbNN_RandTheo_D2_k1_r(r,rho))
}


r=seq(from=min(nn1$knd), to=max(nn1$knd), length.out=5000)
l10<-log10(funOPC(r))
l10[is.na(l10)] <- 0
indi=which(diff(sign(l10))<0)
eps<-mean(r[c(indi[1], indi[1]+1)])
return(eps)
}

calculateNmin2D<-function(nn1,rho, eps){
  alpha=99.85
  tol=1.e-5

    P_val <- 1-alpha/100
    Test_P_val <-1
    rseq1 <- NN_RandTheo_D2_k1_seqAut(rho, N=35000,nsig=1)$rseq
    if(max(rseq1)<eps){
      print(eps)
      print(max(rseq1))
      rseq1<-seq(from=min(nn1$knd), to=max(nn1$knd), length.out=35000)
    }
    k <-0 
    # while(Test_P_val >= P_val){  
    while((Test_P_val -P_val)>tol){
      k<-k+1 
      print('npoints')
      print(k)
      Test_P_val <- CumuNNEmpir_RandTheo_seq_val(eps,rseq1,k,rho,2)
      print(Test_P_val)
    }
    Nmin <-k + 1
  return(Nmin)
}
