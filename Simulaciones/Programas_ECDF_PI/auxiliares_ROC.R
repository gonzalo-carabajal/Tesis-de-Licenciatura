

####################################################
# ROC CONDITIONAL
# ESTIMATORS BASED ON ECDF
####################################################
## INPUT
## yD      	: biomarker in the diseased population (length nD)
## yH      	: biomarker in the healthy population (length nH)
## XCOMP_D	: covariates of the diseased population
##		  must be compositional and a matrix of 
##		  dimension nD*k
## XCOMP_H	: covariates of the healthy population
##		  must be compositional and a matrix of 
##		  dimension nH*k
## grilla_x	: matrix of compositional data where the ROC is to be evaluated
##		: dimension nx*k
##
## grilla_p	: grid of points where the ROC is evaluated
## ------------------------------------------------------
## A LINEAR REGRESSION MODEL IS FITTED
## ------------------------------------------------------
## OUTPUT
##
## ROC_CL_x_grilla	: Conditional ROC curve 
##				   conditional to grilla_x 
##				   evaluated at grilla_p
## AUC_CL_x			: Conditional AUC curve 
##				   conditional to grilla_x 
## YOUDEN_CL_x		: Conditional Youden Index
## ROC_PI_x_grilla	: PLUG--IN conditional ROC curve 
##				   assuming a binormal model
##				   conditional to grilla_x 
##				   evaluated at grilla_p
## AUC_PI_x			: PLUG--IN conditional AUC curve 
##				   assuming a binormal model
##				   conditional to grilla_x 
## YOUDEN_PI_x		: PLUG--IN conditional Youden Index
##				    assuming a binormal model

## beta_D			: estimated slope in the ilr space
##				   for the diseased population
## beta_H			: estimated slope in the ilr space
##				   for the healthy population
## betaCOMP_D		: estimated slope in the simplex
##				   for the diseased population
## betaCOMP_H		: estimated slope in the simplex
##				   for the healthy population
## alfa_D			: estimated intercept for the diseased population
## alfa_H			: estimated intercept for the healhy population
## sigma_D			: estimated errors SD for the diseased population
## sigma_H			: estimated errors SD for the diseased population
## ------------------------------------------------------


CURVA_ROC_ECDF <- function(Y_D,XCOMP_D,Y_H, XCOMP_H, grilla_x, grilla_p){
  	ene_D=length(Y_D)
  	ene_H=length(Y_H)
  
  	x_D <- ilr(XCOMP_D)
  	x_H <- ilr(XCOMP_H)
  	regresion_D  <- lm(Y_D~x_D)
  	regresion_H <- lm(Y_H~x_H)
  
  	residuos_D  <-  regresion_D$residuals      #Estimacion de los residuos
  	residuos_H <- regresion_H$residuals  #Estimacion de los residuos
  
  
  	beta_D <- unname(regresion_D$coef)
  	beta_H<- unname(regresion_H$coef) 
  
  
	alfa_D <-beta_D[1]
	alfa_H <-beta_H[1]

	betaCOMP_D <- ilrInv(beta_D[-1])
	betaCOMP_H <- ilrInv(beta_H[-1])
  
  	SD_D <- sqrt(sum(residuos_D^2)/(ene_D-2))  
  	SD_H<- sqrt(sum(residuos_H^2)/(ene_H-2))
  
  	residuos_D_std<- residuos_D/SD_D
  	residuos_H_std<- residuos_H/SD_H
  
  
  
  	largox <- dim(grilla_x)[1]
  	ROC_CL_x_grilla<-matrix(rep(NA,largox*length(grilla_p)),nrow=largox,ncol=length(grilla_p))
  	ce_CL_x <- YOUDEN_CL_x <- AUC_CL_x<-  rep(NA,largox)

  

  	ROC_PI_x_grilla<-matrix(rep(NA,largox*length(grilla_p)),nrow=largox,ncol=length(grilla_p))
  	YOUDEN_PI_x <- AUC_PI_x<-  rep(NA,largox)

  
  	for(k in 1:largox){
    		x=ilr(grilla_x[k,])
    		Z=c(1,x)
    		#################################################################
    		# definiendo Z=c(1,X) puedo poner  (beta_H-beta_D)%*% Z
    		#################################################################
    
    		dif_predict <-as.numeric((beta_H-beta_D)%*% Z)  ##### OJO ESTO ES -(mu_D-mu_H)=mu_H-mu_D
    
    		for(i in 1:length(grilla_p)){
      		p<-grilla_p[i]
      		cuantil_H <- unname(quantile(residuos_H_std,probs=1-p))
      
      		punto<- cuantil_H*(SD_H/SD_D) + dif_predict/SD_D 
      		ROC_CL_x_grilla[k,i]<- 1-mean(1*(residuos_D_std<=punto))
      	
			if(p!=0 & p!=1){
				cuantil <- qnorm(1-p)
				punto<- cuantil*(SD_H/SD_D) + dif_predict/SD_D 
				ROC_PI_x_grilla[k,i] <- 1- pnorm(punto)
     			 }
  

			if(p==1){
				ROC_PI_x_grilla[k,i] <-   1
				ROC_CL_x_grilla[k,i] <- 1
			}
			if(p==0){
				ROC_PI_x_grilla[k,i] <-  ROC_CL_x_grilla[k,i] <- 0
			}
    		}
    
		AUC_CL_x[k]<-mean(ROC_CL_x_grilla[k,])
    		YOUDEN_CL_x[k] <- max(ROC_CL_x_grilla[k,]- grilla_p)
    		cualindice <- which.max(ROC_CL_x_grilla[k,]- grilla_p)
    		cualp <- grilla_p[cualindice]
    		cual_cuantil_H <- unname(quantile(residuos_H_std,probs=1-cualp))
    		ce_CL_x[k] <- beta_H%*% Z + SD_H*cual_cuantil_H


		AUC_PI_x[k]<-mean(ROC_PI_x_grilla[k,])
    		YOUDEN_PI_x[k] <- max(ROC_PI_x_grilla[k,]- grilla_p)
    	
 
  	} 
  	return(list(ROC_CL_x_grilla=ROC_CL_x_grilla, AUC_CL_x=AUC_CL_x, YOUDEN_CL_x=YOUDEN_CL_x, ce_CL_x=ce_CL_x , 
		ROC_PI_x_grilla=ROC_PI_x_grilla, AUC_PI_x=AUC_PI_x, YOUDEN_PI_x=YOUDEN_PI_x,
    		beta_D=beta_D, beta_H=beta_H, betaCOMP_D=betaCOMP_D, betaCOMP_H=betaCOMP_H, alfa_D=alfa_D, alfa_H=alfa_H, sigma_D=SD_D, sigma_H=SD_H))
}
 
 



####################################################
# ROC CONDITIONAL
# ESTIMATORS BASED ON SMOOTHED ONES
####################################################
## INPUT
## yD      	: biomarker in the diseased population (length nD)
## yH      	: biomarker in the healthy population (length nH)
## XCOMP_D	: covariates of the diseased population
##		  must be compositional and a matrix of 
##		  dimension nD*k
## XCOMP_H	: covariates of the healthy population
##		  must be compositional and a matrix of 
##		  dimension nH*k
## grilla_x	: matrix of compositional data where the ROC is to be evaluated
##		: dimension nx*k
##
## grilla_p	: grid of points where the ROC is evaluated
## ------------------------------------------------------
## A LINEAR REGRESSION MODEL IS FITTED
## ------------------------------------------------------
## OUTPUT
##
## ROC_SMOOTH_x_grilla	: SMOOTHED conditional ROC curve 
##				   conditional to grilla_x 
##				   evaluated at grilla_p
## AUC_SMOOTH_x			: SMOOTHED conditional AUC curve 
##				   conditional to grilla_x 
## YOUDEN_SMOOTH_x		: SMOOTHED conditional Youden Index
##				    assuming a binormal model

## beta_D			: estimated slope in the ilr space
##				   for the diseased population
## beta_H			: estimated slope in the ilr space
##				   for the healthy population
## betaCOMP_D		: estimated slope in the simplex
##				   for the diseased population
## betaCOMP_H		: estimated slope in the simplex
##				   for the healthy population
## alfa_D			: estimated intercept for the diseased population
## alfa_H			: estimated intercept for the healhy population
## sigma_D			: estimated errors SD for the diseased population
## sigma_H			: estimated errors SD for the diseased population
## ------------------------------------------------------


CURVA_ROC_SMOOTH <- function(Y_D,XCOMP_D,Y_H, XCOMP_H, grilla_x, grilla_p, ache=NULL, nucleo="epan"){
  	ene_D=length(Y_D)
  	ene_H=length(Y_H)
  
  	x_D <- ilr(XCOMP_D)
  	x_H <- ilr(XCOMP_H)
  	regresion_D  <- lm(Y_D~x_D)
  	regresion_H <- lm(Y_H~x_H)
  
  	residuos_D  <-  regresion_D$residuals      #Estimacion de los residuos
  	residuos_H <- regresion_H$residuals  #Estimacion de los residuos
  
  
  	beta_D <- unname(regresion_D$coef)
  	beta_H<- unname(regresion_H$coef) 
  
  
	alfa_D <-beta_D[1]
	alfa_H <-beta_H[1]

	betaCOMP_D <- ilrInv(beta_D[-1])
	betaCOMP_H <- ilrInv(beta_H[-1])
  
  	SD_D <- sqrt(sum(residuos_D^2)/(ene_D-2))  
  	SD_H<- sqrt(sum(residuos_H^2)/(ene_H-2))
  
  	residuos_D_std<- residuos_D/SD_D
  	residuos_H_std<- residuos_H/SD_H
  
  
  
  	largox <- dim(grilla_x)[1]
 

  	ROC_SMOOTH_x_grilla<-matrix(rep(NA,largox*length(grilla_p)),nrow=largox,ncol=length(grilla_p))
  	YOUDEN_SMOOTH_x <- AUC_SMOOTH_x<-  rep(NA,largox)

 
  
  	for(k in 1:largox){
    		x=ilr(grilla_x[k,])
    		Z=c(1,x)
    		#################################################################
    		# definiendo Z=c(1,X) puedo poner  (beta_H-beta_D)%*% Z
    		#################################################################
    
    		dif_predict <-as.numeric((beta_H-beta_D)%*% Z)  ##### OJO ESTO ES -(mu_D-mu_H)=mu_H-mu_D
    
		    			 
  		#################################################################
    		#PARA DEFINIR EL SMOOTH
		# DEFINO LAS PSEUDO--VARIABLES W
		# W= residuos_D_std * SD_D/SD_H + (mu_D-mu_H)/SD_H
		#  =  residuos_D/SD_H - dif_predict/SD_H
    		#################################################################
 
		punto_SMOOTH<- residuos_D/SD_H  - dif_predict /SD_H 
			
		W<-rep(NA, ene_D) #Las observaciones de la pseudovariable W
      

		for(jw in 1:ene_D){
        		W[jw]<-1-mean(1*(residuos_H_std<=punto_SMOOTH[jw]))
      
		}

		

    		for(i in 1:length(grilla_p)){
      		p<-grilla_p[i]
 
			if(p!=0 & p!=1){
     			 
  			#################################################################
    			#DEFINO EL SMOOTH
    			#################################################################
   

				if(is.null(ache)){
					cte.cn=1+1.8*ene_D^(-1/5)
					ache= cte.cn*sqrt(5*p*(1-p))/sqrt(2*ene_H)
				}

				argumento <- (p-W)/ache
				ROC_SMOOTH_x_grilla[k,i] <- mean(nucleo_int(argumento, kernel=nucleo))
 
			}

			if(p==1){
				ROC_SMOOTH_x_grilla[k,i] <- 1
				
			}
			if(p==0){
				ROC_SMOOTH_x_grilla[k,i] <- 0
			}
    		}
    
	
		AUC_SMOOTH_x[k]<-mean(ROC_SMOOTH_x_grilla[k,])
    		YOUDEN_SMOOTH_x[k] <- max(ROC_SMOOTH_x_grilla[k,]- grilla_p)

  	} 
  	return(list(ROC_SMOOTH_x_grilla=ROC_SMOOTH_x_grilla, AUC_SMOOTH_x=AUC_SMOOTH_x, YOUDEN_SMOOTH_x=YOUDEN_SMOOTH_x,
    		beta_D=beta_D, beta_H=beta_H, betaCOMP_D=betaCOMP_D, betaCOMP_H=betaCOMP_H, alfa_D=alfa_D, alfa_H=alfa_H, sigma_D=SD_D, sigma_H=SD_H))
}
 
