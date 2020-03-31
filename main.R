source('auxFuncts.R')
source('Init_Fct_Set_detect_sub_v3.R')
library(fpc)
library(stats)
input=read.csv('inputExample.csv', as.is = TRUE)
rd=read.csv(input$filename)
print(paste('read file:', input$filename))
spl<-strsplit(input$filename, '\\.')
ilin<-max(length(spl[[1]])-1,1)
nam=spl[[1]][ilin]
if(input$dim==2){
  #calc q parameter
  print(input$coord)
  if(input$coord =='Ra Dec'){
    print('Great circle distance')
    Distsphe=fct_Distsphe2(rd)
    }
  else{
    print('Euclidean distance')
    Distsphe=as.matrix(dist(rd))
    diag(Distsphe)<-NA
  }
  parQ<-paramQ(Distsphe)$q
  print(paste('Q parameter:', as.character(parQ)))
  if(is.na(input$eps)){
    print('Calculating eps and Nmin for sifnificance')
    if(parQ<input$Qlim){
      #calculate density
      rho=regionDensity2D(Distsphe)
      #first nearest neighbour distribution
      nn1=fct_nns(Distsphe,1)
      #calculate epsilon
      eps=calculateEps2D(nn1,rho)
      print(paste('epsilon:', as.character(eps)))
      
      #calculate nmin
      if(is.na(input$signif)){
        signif=99.85 
      }else if ((input$signif>0)&(input$signif<100)){
        signif=input$signif
      }else{
        print('Not a valid significance percentage value')
        print('using default')
        signif=99.85
      }
      print(paste('Significance level:', as.character(signif)))
      Nmin=calculateNmin2D(nn1,rho,eps,signif)
      print(paste('Nmin:', as.character(Nmin)))
      
      #dbscan: structure detection
      diag(Distsphe)<-0
      dbs <- fpc::dbscan(Distsphe,eps, MinPts=Nmin,method="dist")
      #output
      print(dbs)

      fout<-paste(nam,'out', sep='.')
      write.table(cbind.data.frame(rd,dbs$cluster), fout)
      pdf(paste(nam,'pdf', sep='.'))
      if(input$coord =='Ra Dec'){
        plot(rd, xlim=c(max(rd[,1]),min(rd[,1])), pch=19, col='gray', cex=0.8, main=nam, xlab='RA (deg)', ylab='Dec (deg)')
      }else{
        plot(rd,  pch=19, col='gray', cex=0.8, main=nam, xlab='X', ylab='Y')
      }
      if(max(dbs$cluster)>0){
        inNests<-dbs$cluster > 0
        library(viridis)
        vd<-viridis(max(dbs$cluster))
        points(rd[inNests,],col=vd[dbs$cluster[inNests]],pch=19)}
      dev.off()
    }else{
      print('To search for small significant structures the regions must be structured')
      print(paste('The Q parameter value is:', as.character(parQ)))
      print(paste('With  user suplied limit of: ', as.character(input$Qlim)))
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

      fout<-paste(nam,'out', sep='.')
      write.table(cbind.data.frame(rd,dbs$cluster), fout)
      pdf(paste(nam,'pdf', sep='.'))
      if(input$coord =='Ra Dec'){
        plot(rd, xlim=c(max(rd[,1]),min(rd[,1])), pch=19, col='gray', cex=0.8, main=nam, xlab='RA (deg)', ylab='Dec (deg)')
      }else{
        plot(rd,  pch=19, col='gray', cex=0.8, main=nam, xlab='X', ylab='Y')
      }
      if(max(dbs$cluster)>0){
        inNests<-dbs$cluster > 0
        library(viridis)
        vd<-viridis(max(dbs$cluster))
        points(rd[inNests,],col=vd[dbs$cluster[inNests]],pch=19)}
      dev.off()
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
