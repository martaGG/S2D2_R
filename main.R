source('auxFuncts.R')
source('Init_Fct_Set_detect_sub_v3.R')
library(fpc)
input=read.csv('inputExample.csv', as.is = TRUE)
rd=read.csv(input$filename)

if(input$dim==2){
  #calc q parameter
  Distsphe=fct_Distsphe2(rd)
  parQ<-paramQ(Distsphe)$q
  print(paste('Q parameter:', as.character(parQ)))
  if(is.na(input$eps)){
    if(parQ<input$Qlim){
      #calculate density
      rho=regionDensity2D(Distsphe)
      #first nearest neighbour distribution
      nn1=fct_nns(Distsphe,1)
      #calculate epsilon
      eps=calculateEps2D(nn1,rho)
      print(paste('epsilon:', as.character(eps)))
      
      #calculate nmin
      Nmin=calculateNmin2D(nn1,rho,eps)
      print(paste('Nmin:', as.character(Nmin)))
      
      #dbscan: structure detection
      diag(Distsphe)<-0
      dbs <- fpc::dbscan(Distsphe,eps, MinPts=Nmin,method="dist")
      #output
      print(dbs)
      spl<-strsplit(input$filename, '\\.')
      ilin<-max(length(spl[[1]])-1,1)
      fout<-paste(spl[[1]][ilin],'out.dat', sep='_')
      write.table(cbind.data.frame(rd,dbs$cluster), fout)
    }else{
      print('To search for small significant structures the regions must be structured')
      print('The Q parameter criteria does not hold:')
      print(parQ)
    }
  }else if(input$eps>0){
    if(input$Nmin>0){
      #calculate density
      rho=regionDensity2D(Distsphe)
      #first nearest neighbour distribution
      nn1=fct_nns(Distsphe,1)
      eps=input$eps
      Nmin=input$Nmin
      #dbscan
      diag(Distsphe)<-0
      dbs <- fpc::dbscan(Distsphe,eps, MinPts=Nmin,method="dist")
      #output
      print(dbs)
      spl<-strsplit(input$filename, '.')
      ilin<-max(length(spl[[1]])-1,1)
      fout<-paste(spl[[1]][ilin],'out.dat', sep='_')
      write.table(cbind.data.frame(rd,dbs$cluster), fout)
    }else{
      print("please input also an minimum number of points for dbscan")
    }
  }else{
    print("not a valid epsilon value. Should be NA for automatic calculation")
    print("or a number >0 for input in dbscan")
  }
    
}else{
  print("only 2D for now")
}    
