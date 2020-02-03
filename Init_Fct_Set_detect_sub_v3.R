
################################################################################################
# Function fct_rho (test_rho,RaDec,WinPoly=NULL) - Compute mean density
	# Input for a polygonal window to estimate the mean density 4 points with Ra and Dec Winpoly[4,2]
	# what kind of mean density input
	#
	# Test_rho <- 1: given estimate of the mean density from the user
	# Test_rho <- 2: Estimate from the whole region (OK if no big empty space)
	# Test_rho <- 3:  estimate from a given polygonal window
	# Test_rho <- 4: estimate from the chull of the polyonal window (default)
# 

# DistEucli2D (xyz,NATrue=TRUE): Function distance euclidienne 2D

#fct_nns (sd.dist,k): kth Nearest neigbor  distance, from distance matrix sd.dist
		# BEWARE: the diagonal of the distance matrix should be NA

#fct_Distsphe (RaDec_deg,method=1): Spherical distance between two points on the sky
	# input coord in degree RaDec_deg[1:n,1:2] 1->RA; 2->DEC;
	# return: arc distance in degree
	# method=1 use gcirc function of astrolibR package 
	# method=2 uses in cartesian coordinate compute DX,DY,DZ then the associated distance D = sqrt(DX^2+DY^2+DZ^2); Arc=2 asin(D/2);	
	# method=3 uses the atan function
	# BEWARE: methods 2 and 3 need to be tested (bugs)

# fct_RA_deg (coord): Convert coordinates RA from coord[h, min, sec] to ra_deg[degree] 

# fct_DEC_deg (coord): Convert coordinates DEC from coord[deg, arcmin, arcsec] to dec_deg[degree] 

# fct_ReadData (FileDir,FileName,nrow,ncol): Read data 

###########################################################################################################
#The following functions are a set of functions related to the computation of theoretical random nearest neighbo statistics

# Note: the distribution of the 1-NNS does not depend strongly on the assumption of the infiniteness, so it does not depend strongly on the shape not the surface. It depends mainly on the value of the mean density rho.

# NN_RandTheo_D2_k1_seqAut (rho, rseq=NULL,N=1000,nsig=1): distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial sequence
# NN_RandTheo_D2_k1_seq (rho, rseq): : distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, based on the given rseq the radial sequence
# NN_RandTheo_D2_k1_seqlogAut (rho, logrseq=NULL,N=1000,nsig=3):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial logarithm sequence
# NN_RandTheo_D2_k1_seqlog (rho, logrseq):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# LogNN_RandTheo_D2_k1_seqlog (rho,logrseq): LOGARITHM of the distribution (log PDF) of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# ProbNN_RandTheo (r,k,rho,Dim): PDF of the k-NNS in space of dimension Dim for a given r 
# ProbNN_RandTheo_seq (rseq,k,rho,Dim):  PDF of the k-NNS in space of dimension Dim for a given radial sequence rseq 
# CumuNNEmpir_RandTheo_seq_val (rval,rseq,k,rho,Dim): Empirical cumulative function of the PDF of the k-NNS in space of dimension Dim for a given rval
# ProbNN_RandTheo_D2 (r,k,rho): PDF of the k-NNS in space of dimension 2 for a given r 
# ProbNN_RandTheo_D2_k1_r (r,rho): PDF of the first-NNS in space of dimension 2 for a given r 
# ProbNN_RandTheo_D2_k_r (r,rho,k):  PDF of the k-NNS in space of dimension 2 for a given r  
# CumuNN_RandTheo_D2_k1_r (r,rho): Theoretical cumulative function of the PDF of the first NNS in space of dimension 2 for a given r
# CumuNN_RandTheo_D2_k2_r (r,rho): Theoretical cumulative function of the PDF of the second NNS in space of dimension 2 for a given r
# CumuNN_RandTheo_D2_k3_r (r,rho): Theoretical cumulative function of the PDF of the third NNS in space of dimension 2 for a given r
# ProbNN_RandTheo_D2_k1_rlog (rlog, rho): PDF of the first-NNS in space of dimension 2 for a given radius logarithm rlog 
# LogProbNN_RandTheo_D2_k1_rlog (rlog,rho): Log PDF of the first-NNS in space of dimension 2 for a given radius logarithm rlog 




##########################################################################################

#
# Function fct_rho (test_rho,RaDec,WinPoly=NULL) - Compute mean density


##########################################################################################

# Function fct_rho (test_rho,RaDec,WinPoly=NULL)
# Compute mean density
#
# 
# Input for a polygonal window to estimate the mean density 4 points with Ra and Dec Winpoly[4,2]
# what kind of mean density input
#
# Test_rho <- 1: given estimate of the mean density from the user
# Test_rho <- 2: Estimate from the whole region (OK if no big empty space) - default
# Test_rho <- 3:  estimate from a given polygonal window
# Test_rho <- 4: estimate from the chull of the polyonal window 
# 

fct_rho<- function (test_rho,RaDec,WinPoly=NULL){
 	if(test_rho==2){
 		Win <- RaDec[chull(RaDec),]
 		WinPoly <- as(RaDec[chull(RaDec),], "gpc.poly")
 		#WinPoly_area <- area.poly(WinPoly)
 		test1 <- WinPoly@pts[[1]]$x
 		test2 <-WinPoly@pts[[1]]$y
 		test2l <-test2[length(test2):1]
 		test1l <-test1[length(test1):1]
 		OwinPoly<- owin(poly=list(cbind(test1l,test2l)))
 		nb <- length(which(inside.owin(RaDec[,1],RaDec[,2],as.owin(OwinPoly))==TRUE))
		area <- area.owin(OwinPoly)
 		rho <- nb/area
 		df <-data.frame(rho,nb,area)
 		return(df)}
 		
 		if(test_rho==3){
 		nb <- length(which(inside.owin(RaDec[,1],RaDec[,2],as.owin(WinPoly))==TRUE))
		area <- area.owin(WinPoly)
 		rho <- nb/area
 		df <-data.frame(rho,nb,area)
 		return(df)}

 		if(test_rho==4){
 			Pts <- which(inside.owin(RaDec[,1],RaDec[,2],as.owin(WinPoly))==TRUE)
 			Win <- RaDec[Pts,][chull(RaDec[Pts,]),]
 			WinPoly2 <- as(Win, "gpc.poly")
 			test1 <- WinPoly2@pts[[1]]$x
 			test2 <-WinPoly2@pts[[1]]$y
 			test2l <-test2[length(test2):1]
 			test1l <-test1[length(test1):1]
 			OwinPoly<- owin(poly=list(cbind(test1l,test2l)))	
 			nb <- length(which(inside.owin(RaDec[,1],RaDec[,2],as.owin(OwinPoly))==TRUE))
			area <- area.owin(OwinPoly)
 			rho <- nb/area
 			df <-data.frame(rho,nb,area)
 		return(df)}		 		
 }		
################################################################################################
#DistEucli2D (xyz,NATrue=TRUE): Function distance euclidienne 2D
#

# Function distance euclidienne 2D
#
DistEucli2D<-function(xyz,NATrue=TRUE){n <- length(xyz)/2
	dist2D <- matrix(ncol=n, nrow=n)
	for (i in 1:n){
		for (j in 1:n){
			dist2D[i,j] <-sqrt((xyz[j,1]-xyz[i,1])**2+(xyz[j,2]-xyz[i,2])**2)
		}
	}
			if (NATrue==TRUE){for (i in 1:n){dist2D[i,i] <-NA}}
	return(dist2D)
	}

################################################################################################
#fct_nns (sd.dist,k): kth Nearest neigbor  distance, from distance matrix sd.dist
# BEWARE: the diagonal of the distance matrix should be NA

# kth Nearest neigbor distance 
# BEWARE the diagonal of the distance matrix should be NA
#


fct_nns <- function(sd.dist,k){
k.for.nn <- k
n = nrow(sd.dist)
knn.mat = matrix(0, ncol = k.for.nn, nrow = n)
knd.mat = knn.mat
knn = vector("numeric", length= n)
knd = knn

for(i in 1:n){
  knn.mat[i,] = order(sd.dist[i,])[1:k.for.nn]
  knd.mat[i,] = sd.dist[i,knn.mat[i,]]
    knn[i] = knn.mat[i,k.for.nn]
  knd[i] = knd.mat[i,k.for.nn]
}
df <-data.frame(knd,knn)
return(df)}

################################################################################################
#fct_Distsphe <- function(RaDec_deg,method=1): Spherical distance between two points on the sky
	# input coord in degree RaDec_deg[1:n,1:2] 1->RA; 2->DEC;
	# return: arc distance in degree
	# method=1 use gcirc function of astrolibR package 
	# method=2 uses in cartesian coordinate compute DX,DY,DZ then the associated distance D = sqrt(DX^2+DY^2+DZ^2); Arc=2 asin(D/2);	
	# method=3 uses the atan function
	# BEWARE: methods 2 and 3 need to be tested (bugs)

#fct_Distsphe <- function(RaDec_deg,method=1): Spherical distance between two points on the sky#
# interpoint distance compute on the sky seen as a sphere 
# input coord in degree RaDec_deg[1:n,1:2] 1->RA; 2->DEC;
# return: arc distance in degree
# method=1 use gcirc function of astrolibR package
#
# ATTENTION: do not use yet methods 2 and 3 (bugs to be fixed) 
# method=2 uses in cartesian coordinate compute DX,DY,DZ then the associated distance D = sqrt(DX^2+DY^2+DZ^2); Arc=2 asin(D/2);
# method=3 uses the atan function
 

fct_Distsphe <- function(RaDec_deg,method=1){
	n=length(RaDec_deg[,1])
	Distsphe <- matrix(nrow = n, ncol =n)
	if (method==1){
	for (ind in seq(from=1, to=n)){
		for (jind in seq(from=1, to=n)){
			 Distsphe[ind, jind] <- astrolibR::gcirc(2,RaDec_deg[ind,1],RaDec_deg[ind,2],RaDec_deg[jind,1],RaDec_deg[jind,2])}
		Distsphe[ind, ind]<-NA}
		return(Distsphe/3600) # routine gcirc gives a value back in arcsec --> convert here in degree}
	if (method==2){
		RaDec_rad=RaDec_deg/180*pi	
		for (id in seq(from=1, to=n)){
			for (jd in seq(from=1, to=n)){
				DX <- cos(RaDec_rad[jd,2])*cos(RaDec_rad[jd,1])-cos(RaDec_rad[id,2])*cos(RaDec_rad[id,1])
				DY <- cos(RaDec_rad[jd,2])*sin(RaDec_rad[jd,1])-cos(RaDec_rad[id,2])*sin(RaDec_rad[id,1])
				DZ <- sin(RaDec_rad[jd,2])-sin(RaDec_rad[id,2])
				C <- sqrt(DX**2+DY**2+DZ**2)
				Distsphe[id, jd] <- 2*asin(C/2)
				Distsphe[id, jd]/pi*180}
				Distsphe[id, id]<-NA}
				return(Distsphe)}
	if (method==3){
		RaDec_rad=RaDec_deg/180*pi	
		for (id in seq(from=1, to=n)){
			for (jd in seq(from=1, to=n)){			
				X<- sqrt((sin(RaDec_rad[id,2])*sin(RaDec_rad[jd,2])) + (cos(RaDec_rad[id,2])*cos(RaDec_rad[jd,2])*cos(RaDec_rad[jd,1]-RaDec_rad[id,1])) )
				Y<- sqrt( (cos(RaDec_rad[jd,2])*sin(RaDec_rad[jd,1]-RaDec_rad[id,1]))**2
		 			+(cos(RaDec_rad[id,2])*sin(RaDec_rad[jd,2])-sin(RaDec_rad[id,2])*cos(RaDec_rad[jd,2])*cos(RaDec_rad[jd,1]-RaDec_rad[id,1]))**2)
		 		Distsphe[id, jd] <- atan(Y/X)	
		 	Distsphe[id, jd]/pi*180}
		Distsphe[id, id]<-NA}
	return(Distsphe)}
	}}
		




################################################################################################

# fct_RA_deg (coord): Convert coordinates RA from coord[h, min, sec] to ra_deg[degree] 

# fct_DEC_deg (coord): Convert coordinates DEC from coord[deg, arcmin, arcsec] to dec_deg[degree] 


#
# Convert coordinates from RA h, min, sec to degree DEC deg,arcmin, arcsec to degree
# RA_deg, DEC_deg, RADEC_deg (coord())
#return

	fct_RA_deg <- function(coord){
	ra_deg <- (coord[1]+coord[2]/60+coord[3]/3600)/24*360
return(ra_deg)}

	fct_DEC_deg <- function(coord){
	dec_deg <-  (coord[4]+coord[5]/60+coord[6]/3600)
return(dec_deg)}

	fct_RADEC_deg <- function(coord){
	ra_deg <- (coord[1]+coord[2]/60+coord[3]/3600)/24*360
	dec_deg <-  (coord[4]+coord[5]/60+coord[6]/3600)
	df <-data.frame(ra_deg,dec_deg)
return(df)}


################################################################################################

# fct_ReadData (FileDir,FileName,nrow,ncol): Read data 

#
fct_ReadData <- function(FileDir,FileName,nrow,ncol){
FilDirDef <- getwd()
setwd(FileDir)
data <- matrix(scan(FileName, n = nrow*ncol, what = character()), nrow, ncol, byrow = TRUE)
setwd(FilDirDef)
return(data)}
#



###########################################################################################################
# COLLECTION OF FUNCTIONS RELATED TO NEAREST NEIGHBOR IN RANDOM SPATIAL DISTRIBUTION
#	depend on the dimension of the space (dim) and the kth nearest neighbor
# Functions that give the Moments of the first nearest neighbor distribution
#
# Function NNLog_RandTheo_D2_k1_seqlog (logrseq,rho)
#
#  give as a return the Probability distribution of the first nearest neighbor, the distance sequence given as an input in Log10, in the case of spatial random distribution in an infinite medium.
#
# Note: the distribution of the 1-NNS does not depend strongly on the assumption of the infiniteness, so it does not depend strongly on the shape not the surface. It depends mainly on the value of the density
#
########################################################################################################### This function return automatically the distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process
# Values NN_RandTheo_D2_k1_seqAut

# NN_RandTheo_D2_k1_seqAut (rho, rseq=NULL,N=1000,nsig=1): distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process

 NN_RandTheo_D2_k1_seqAut <- function(rho, rseq=NULL,N=1000,nsig=1)
{	
	if (is.null(rseq)){
		MeanNN <-1/2/sqrt(rho)
		SdNN <- sqrt((4-pi)/(4*pi*rho))
		MaxNN <- MeanNN+nsig*SdNN
		dNN <-MaxNN/(N-1)
		rseq<- seq(0,MaxNN,dNN)}
	Prob_seq <-vector(length=length(rseq))
	for (ind in 1:length(rseq)){
	r <-rseq[ind]
	Prob_seq[ind] <-ProbNN_RandTheo_D2_k1_r(r,rho)
	}
	Prob <- Prob_seq
	df <- data.frame(rseq,Prob)
	return(df)}


########################################################################################################### This function return automatically the 
# NN_RandTheo_D2_k1_seq (rho, rseq): : distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, based on the given rseq the radial sequence
	
NN_RandTheo_D2_k1_seq <- function(rho, rseq)
{	Prob_seq <-vector(length=length(rseq))
	for (ind in 1:length(rseq)){
	r <-rseq[ind]
	Prob_seq[ind] <-ProbNN_RandTheo_D2_k1_r(r,rho)
	}
	return(Prob_seq)}	
	

########################################################################################################### This function return automatically the 
# NN_RandTheo_D2_k1_seqlogAut (rho, logrseq=NULL,N=1000,nsig=3):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial logarithm sequence

NN_RandTheo_D2_k1_seqlogAut <- function(rho, logrseq=NULL,N=1000,nsig=3)
{	
	if (is.null(logrseq)){
		MeanNN <-log10(1/2/sqrt(rho))
		SdNN <- log10(sqrt((4-pi)/(4*pi*rho)))
		MaxNN <- MeanNN+nsig*abs(SdNN)
		MinNN <- MeanNN-nsig*abs(SdNN)
		dNN <-(MaxNN-MinNN)/(N-1)
		logrseq<- seq(MinNN,MaxNN,dNN)}

	Prob_seq <-vector(length=length(logrseq))
	for (ind in 1:length(logrseq)){
	logr <-logrseq[ind]
	Prob_seq[ind] <-ProbNN_RandTheo_D2_k1_rlog(logr,rho)
	}
	logr <-logrseq
	Prob <- Prob_seq
	df <- data.frame(logr,Prob)
	return(df)}	
	

###########################################################################################################

# NN_RandTheo_D2_k1_seqlog (rho, logrseq):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence

	
NN_RandTheo_D2_k1_seqlog <- function(rho, logrseq)
{	Prob_seq <-vector(length=length(logrseq))
	for (ind in 1:length(logrseq)){
	logr <-logrseq[ind]
	Prob_seq[ind] <-ProbNN_RandTheo_D2_k1_rlog(logr,rho)
	}
	return(Prob_seq)}	
	


###########################################################################################################

# LogNN_RandTheo_D2_k1_seqlog (rho,logrseq): LOGARITHM of the distribution (log PDF) of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence


#  give as a return the Probability in Log10 that the first nearest neighbor is located at a value rlog, the distance given as an input in Log10  in the case of spatial random distribution in an infinite medium.
#
LogNN_RandTheo_D2_k1_seqlog <- function(rho,logrseq)
{	logProb_seq <-vector(length=length(logrseq))
	for (ind in 1:length(logrseq)){
	rlog <-logrseq[ind]
	logProb_seq[ind] <-LogProbNN_RandTheo_D2_k1_rlog(rlog,rho)
	}
	return(logProb_seq)}	
	



#

###########################################################################################################

# ProbNN_RandTheo (r,k,rho,Dim): PDF of the k-NNS in space of dimension Dim

ProbNN_RandTheo <- function(r,k,rho,Dim)
{cd <- pi**(Dim/2)/gamma(Dim/2+1)
Vol_r <- cd*r**Dim
	if (k <100) {Prob <- Dim/r*(rho*Vol_r)**k/gamma(k)*exp(-rho*Vol_r)}
	if (k >=100) {Prob <-0}
return(Prob)}


###########################################################################################################

# ProbNN_RandTheo_seq (rseq,k,rho,Dim):  PDF of the k-NNS in space of dimension Dim for a given radial sequence rseq 

ProbNN_RandTheo_seq <- function(rseq,k,rho,Dim)
{	
	Prob_seq <-vector(length=length(rseq))
	for (ind in 1:length(rseq)){
	r <- rseq[ind]
	Prob_seq[ind] <- ProbNN_RandTheo (r,k,rho,Dim)}
return(Prob_seq)}

###########################################################################################################
 
# CumuNNEmpir_RandTheo_seq_val (rval,rseq,k,rho,Dim): Empirical cumulative function of the PDF of the k-NNS in space of dimension Dim for a given rval

CumuNNEmpir_RandTheo_seq_val <- function(rval,rseq,k,rho,Dim)
{
	indmax <- which(rseq >=rval)[1]
	ProbCumu <-0
	for (ind in 2:indmax){
		ProbCumu <- ProbCumu + (rseq[ind]-rseq[ind-1])*ProbNN_RandTheo (rseq[ind],k,rho,Dim)
	}	
return(ProbCumu)
}


###########################################################################################################

# ProbNN_RandTheo_D2 (r,k,rho): PDF of the k-NNS in space of dimension 2 for a given r 

ProbNN_RandTheo_D2 <- function(r,k,rho)
{if (k <100) {Prob <- 2*(pi*rho)**k/gamma(k)*r**(2*k-1)*exp(-rho*pi*r**2)}
	if (k >100) {Prob <-0}
return(Prob)}

###########################################################################################################
#This function return automatically the distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process
# Values NN_RandTheo_D2_k1_seqAut
# Note: the distribution of the 1-NNS does not depend strongly on the assumption of the infiniteness, so it does not depend strongly on the shape not the surface. It depends mainly on the value of the mean density

# NN_RandTheo_D2_k1_seqAut (rho, rseq=NULL,N=1000,nsig=1): distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial sequence
# NN_RandTheo_D2_k1_seq (rho, rseq): : distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, based on the given rseq the radial sequence
# NN_RandTheo_D2_k1_seqlogAut (rho, logrseq=NULL,N=1000,nsig=3):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial logarithm sequence
# NN_RandTheo_D2_k1_seqlog (rho, logrseq):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# LogNN_RandTheo_D2_k1_seqlog (rho,logrseq): LOGARITHM of the distribution (log PDF) of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# ProbNN_RandTheo (r,k,rho,Dim): PDF of the k-NNS in space of dimension Dim for a given r 
# ProbNN_RandTheo_seq (rseq,k,rho,Dim):  PDF of the k-NNS in space of dimension Dim for a given radial sequence rseq 
# CumuNNEmpir_RandTheo_seq_val (rval,rseq,k,rho,Dim): Empirical cumulative function of the PDF of the k-NNS in space of dimension Dim for a given rval
# ProbNN_RandTheo_D2 (r,k,rho): PDF of the k-NNS in space of dimension 2 for a given r 
# ProbNN_RandTheo_D2_k1_r (r,rho): PDF of the first-NNS in space of dimension 2 for a given r 

ProbNN_RandTheo_D2_k1_r <- function(r,rho)
{	
	Prob <-2*pi*rho*r*exp(-pi*rho*(r**2))
	return(Prob)}

###########################################################################################################
#This function return automatically the distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process
# Values NN_RandTheo_D2_k1_seqAut
# Note: the distribution of the 1-NNS does not depend strongly on the assumption of the infiniteness, so it does not depend strongly on the shape not the surface. It depends mainly on the value of the mean density

# NN_RandTheo_D2_k1_seqAut (rho, rseq=NULL,N=1000,nsig=1): distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial sequence
# NN_RandTheo_D2_k1_seq (rho, rseq): : distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, based on the given rseq the radial sequence
# NN_RandTheo_D2_k1_seqlogAut (rho, logrseq=NULL,N=1000,nsig=3):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial logarithm sequence
# NN_RandTheo_D2_k1_seqlog (rho, logrseq):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# LogNN_RandTheo_D2_k1_seqlog (rho,logrseq): LOGARITHM of the distribution (log PDF) of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# ProbNN_RandTheo (r,k,rho,Dim): PDF of the k-NNS in space of dimension Dim for a given r 
# ProbNN_RandTheo_seq (rseq,k,rho,Dim):  PDF of the k-NNS in space of dimension Dim for a given radial sequence rseq 
# CumuNNEmpir_RandTheo_seq_val (rval,rseq,k,rho,Dim): Empirical cumulative function of the PDF of the k-NNS in space of dimension Dim for a given rval
# ProbNN_RandTheo_D2 (r,k,rho): PDF of the k-NNS in space of dimension 2 for a given r 
# ProbNN_RandTheo_D2_k1_r (r,rho): PDF of the first-NNS in space of dimension 2 for a given r 
# ProbNN_RandTheo_D2_k_r (r,rho,k):  PDF of the k-NNS in space of dimension 2 for a given r  

ProbNN_RandTheo_D2_k_r <- function(r,rho,k)
{	
	Prob <-2*pi*rho*r*exp(-pi*rho*(r**2))
	return(Prob)}


###########################################################################################################
#This function return automatically the distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process
# Values NN_RandTheo_D2_k1_seqAut
# Note: the distribution of the 1-NNS does not depend strongly on the assumption of the infiniteness, so it does not depend strongly on the shape not the surface. It depends mainly on the value of the mean density

# NN_RandTheo_D2_k1_seqAut (rho, rseq=NULL,N=1000,nsig=1): distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial sequence
# NN_RandTheo_D2_k1_seq (rho, rseq): : distribution of the first nearest neighbor in 2D based on the mean density of the spatial random (Poisson) process, based on the given rseq the radial sequence
# NN_RandTheo_D2_k1_seqlogAut (rho, logrseq=NULL,N=1000,nsig=3):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, compute automatically the radial logarithm sequence
# NN_RandTheo_D2_k1_seqlog (rho, logrseq):  distribution of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# LogNN_RandTheo_D2_k1_seqlog (rho,logrseq): LOGARITHM of the distribution (log PDF) of the first nearest neighbor LOGARITHM  in 2D based on the mean density of the spatial random (Poisson) process, based on a given radial logarithm sequence
# ProbNN_RandTheo (r,k,rho,Dim): PDF of the k-NNS in space of dimension Dim for a given r 
# ProbNN_RandTheo_seq (rseq,k,rho,Dim):  PDF of the k-NNS in space of dimension Dim for a given radial sequence rseq 
# CumuNNEmpir_RandTheo_seq_val (rval,rseq,k,rho,Dim): Empirical cumulative function of the PDF of the k-NNS in space of dimension Dim for a given rval
# ProbNN_RandTheo_D2 (r,k,rho): PDF of the k-NNS in space of dimension 2 for a given r 
# ProbNN_RandTheo_D2_k1_r (r,rho): PDF of the first-NNS in space of dimension 2 for a given r 
# ProbNN_RandTheo_D2_k_r (r,rho,k):  PDF of the k-NNS in space of dimension 2 for a given r  
# CumuNN_RandTheo_D2_k1_r (r,rho): Theoretical cumulative function of the PDF of the first NNS in space of dimension 2 for a given r

CumuNN_RandTheo_D2_k1_r <- function(r,rho)
{	
	Cumu <-1-exp(-pi*rho*r*r)
	return(Cumu)}

###########################################################################################################

# CumuNN_RandTheo_D2_k2_r (r,rho): Theoretical cumulative function of the PDF of the second NNS in space of dimension 2 for a given r

CumuNN_RandTheo_D2_k2_r <- function(r,rho)
{	
	Cumu <- 1-(pi*rho*r*r)*exp(-pi*r**2)
	return(Cumu)}

###########################################################################################################

# CumuNN_RandTheo_D2_k3_r (r,rho):: Theoretical cumulative function of the PDF of the third NNS in space of dimension 2 for a given r
	
CumuNN_RandTheo_D2_k3_r <- function(r,rho)
{	
	Cumu <- 1/2*( exp(-pi*rho*r*r) *( -pi*rho*r*r*(pi*rho*r*r+2)-2 )  +2) 
	return(Cumu)}
		

###########################################################################################################

# ProbNN_RandTheo_D2_k1_rlog (rlog, rho): PDF of the first-NNS in space of dimension 2 for a given radius logarithm rlog 

ProbNN_RandTheo_D2_k1_rlog <- function(rlog, rho)
{	r <-10**rlog
	Prob <-2*log(10)*pi*rho*r**2*exp(-pi*rho*r**2)
	return(Prob)}

###########################################################################################################

# LogProbNN_RandTheo_D2_k1_rlog (rlog,rho): Log PDF of the first-NNS in space of dimension 2 for a given radius logarithm rlog 

LogProbNN_RandTheo_D2_k1_rlog <- function(rlog,rho)
{	r <-10**rlog
	Log_Prob <- log10(2*log(10)*pi*rho)+2*r-(pi*rho/log(10))*r**2
	return(Log_Prob)}
	
