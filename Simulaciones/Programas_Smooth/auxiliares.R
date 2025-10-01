
######################################
# ROC NO CONDICIONAL
######################################
## INPUT
## yD      	: biomarker in the diseased population
## yH      	: biomarker in the healthy population
## p		: grid of points where the ROC is evaluated
## ------------------------------------------------------
## OUTPUT
##
##  the ROC evaluated at  p
## ------------------------------------------------------


ROC.clasica.dim1<- function(yD,yH,p=seq(0,1,length=101)) {
  	xD.emp <- ecdf(yD)
	lp <- length(p)
	roc.cl = rep(NA, length=lp)
	for(i in 1:lp){
      	punto_p<-p[i]
		roc.cl[i] <- 1- xD.emp(quantile(yH,1-punto_p))
		if(punto_p==1){roc.cl[i]  <- 1}
		if(punto_p==0){roc.cl[i]  <- 0}
	}

	return(roc.cl)
}


####################################################
# ROC CONDITIONAL
# ESTIMATORS BASED ON ECDF and SMOOTHED ONES
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
## ROC_SMOOTH_x_grilla	: SMOOTHED conditional ROC curve 
##				   conditional to grilla_x 
##				   evaluated at grilla_p
## AUC_SMOOTH_x			: SMOOTHED conditional AUC curve 
##				   conditional to grilla_x 
## YOUDEN_SMOOTH_x		: SMOOTHED conditional Youden Index
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


CURVA_ROC_BOTH <- function(Y_D,XCOMP_D,Y_H, XCOMP_H, grilla_x, grilla_p, ache=NULL, nucleo="epan"){
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

  	ROC_SMOOTH_x_grilla<-matrix(rep(NA,largox*length(grilla_p)),nrow=largox,ncol=length(grilla_p))
  	YOUDEN_SMOOTH_x <- AUC_SMOOTH_x<-  rep(NA,largox)

  	ROC_PI_x_grilla<-matrix(rep(NA,largox*length(grilla_p)),nrow=largox,ncol=length(grilla_p))
  	YOUDEN_PI_x <- AUC_PI_x<-  rep(NA,largox)

  
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
      		cuantil_H <- unname(quantile(residuos_H_std,probs=1-p))
      
      		punto<- cuantil_H*(SD_H/SD_D) + dif_predict/SD_D 
      		ROC_CL_x_grilla[k,i]<- 1-mean(1*(residuos_D_std<=punto))
      	
			if(p!=0 & p!=1){
				cuantil <- qnorm(1-p)
				punto<- cuantil*(SD_H/SD_D) + dif_predict/SD_D 
				ROC_PI_x_grilla[k,i] <- 1- pnorm(punto)
     			 
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
				ROC_PI_x_grilla[k,i] <- ROC_SMOOTH_x_grilla[k,i] <- 1
				ROC_CL_x_grilla[k,i] <- 1
			}
			if(p==0){
				ROC_PI_x_grilla[k,i] <- ROC_SMOOTH_x_grilla[k,i] <- ROC_CL_x_grilla[k,i] <- 0
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
    	
		AUC_SMOOTH_x[k]<-mean(ROC_SMOOTH_x_grilla[k,])
    		YOUDEN_SMOOTH_x[k] <- max(ROC_SMOOTH_x_grilla[k,]- grilla_p)

  	} 
  	return(list(ROC_CL_x_grilla=ROC_CL_x_grilla, AUC_CL_x=AUC_CL_x, YOUDEN_CL_x=YOUDEN_CL_x, ce_CL_x=ce_CL_x , 
		ROC_PI_x_grilla=ROC_PI_x_grilla, AUC_PI_x=AUC_PI_x, YOUDEN_PI_x=YOUDEN_PI_x,
    		ROC_SMOOTH_x_grilla=ROC_SMOOTH_x_grilla, AUC_SMOOTH_x=AUC_SMOOTH_x, YOUDEN_SMOOTH_x=YOUDEN_SMOOTH_x,
    		beta_D=beta_D, beta_H=beta_H, betaCOMP_D=betaCOMP_D, betaCOMP_H=betaCOMP_H, alfa_D=alfa_D, alfa_H=alfa_H, sigma_D=SD_D, sigma_H=SD_H))
}
 

############################
# EPANECHNIKOV KERNEL INTEGRATED 
#########################

nucleo_EPAN=function(x){
	a=0.75*(x-x^3/3+2/3)
	b=a*(abs(x)<1)+1*(1<=x)
	return(b)
}


nucleo_int<- function(x, kernel="gaussian"){

	 if(kernel=="epan"){
   		 nucleo_int<-nucleo_EPAN(x)
  	}

	if(kernel=="gaussian"){
    		nucleo_int <- pnorm(x) 
    	}

	return(nucleo_int)
}

	
#########################################
#Kernel Distribution Function Estimation
#########################################

  ##INPUTS
  ## datos: observaciones sobre las cuales se realiza la estimación
  ## kernel: núcleo elegido
  ## h: ventana 
  ##OUTPUT
  ## dist: kernel cumulative dist estimated function



dist_acum_ker<-function(datos, kernel="gaussian", h=NULL){
  
  	n<-length(datos)
  	if(is.null(h)){ #Si no se pasa como parámetro, calculo ventana de Silverman
    		S<-sd(datos)
    		IQR<-IQR(datos)
    		h<-1.06*min(c(S, IQR/1.349))*n^{-1/5}
  	}
	dist<-function(t){
    		res<-rep(NA, length(t))
    		for (i in 1:length(t)){
      		arg<-(t[i]-datos)/h
      		res[i]<- mean(nucleo_int(arg))
    		}
    		return(res)
  	}
  
  	return(dist)
}



#####################################################
# TRUE CONDITIONAL ROC ASSUMING NORMAL ERRORS
#####################################################
## INPUT
## XCOMP	:vector of compositional data where the ROC is to be evaluated
##		: length k
## beta_d	: true slope in the simplex (length k)
##				   for the diseased population (length k)
## beta_h	: true slope in the simplex
##				   for the healthy population
## alfa_d	: true intercept for the diseased population
## alfa_h	: true intercept for the healhy population
## sigma_d	: true SD of the errors  for the diseased population
## sigma_h	: true SD of the errors for the diseased population

## p		: real valued point where the ROC is evaluated
## ------------------------------------------------------
## A LINEAR REGRESSION MODEL IS ASSUMED
## ------------------------------------------------------
## OUTPUT
##
## roc	: Conditional ROC curve conditional to XCOMP 
##		  evaluated at p
## ------------------------------------------------------


roc_x_verdadera<-function(p,XCOMP, alfa_d, alfa_h, beta_d, beta_h, sigma_d, sigma_h,dist_eps=list(type_H = 'gaussian', type_D = 'gaussian',gl_H=NA,gl_D=NA)){
	###########################################################################
	# XCOMP, beta_d and beta_h : vectors of compositions of the same dimension
	# alfa_d and alfa_h correspond to the intercept
	# p is a real number where the ROC_x is computed
	###########################################################################

	mu_d <- alfa_d + prod.interno.atch(beta_d,XCOMP)
	mu_h <- alfa_h + prod.interno.atch(beta_h,XCOMP)

	diferencia <- (mu_h-mu_d)/sigma_d
	cociente <- sigma_h/sigma_d
      	if(p!=0 & p!=1){

		if(dist_eps$type_H == 'gaussian'){
			cuantil <- qnorm(1-p)
		}

		if(dist_eps$type_H == 'Student'){
			gradosH = dist_eps$gl_H
			cuantil <- qt(1-p, df=gradosH)*sqrt((gradosH-2)/gradosH)
		}

		if(dist_eps$type_H == 'logistic'){
			cuantil <- qlogis(1-p, 0,scale=sqrt(3)/pi)
		}

		if(dist_eps$type_H == 'uniform'){
			cuantil <- qunif(1-p, -sqrt(3) ,sqrt(3))
		}

		##########################################
		# ARGUMENTO DEL G_D
		##########################################

		arg <- diferencia+cociente*cuantil

		##########################################
		# DEFINO G_D SEGUN DISTRIBUCION
		##########################################


		if(dist_eps$type_D == 'gaussian'){
			roc <- 1- pnorm(arg)
		}

		if(dist_eps$type_D == 'Student'){
			gradosD = dist_eps$gl_D
			corrijo = sqrt(gradosD/(gradosD-2))
			roc <- 1- pt(arg *corrijo , df=gradosD)
		}

		if(dist_eps$type_D == 'logistic'){
			roc <- 1- plogis(arg , 0 , scale=sqrt(3)/pi)
		}

		if(dist_eps$type_D == 'uniform'){
			roc <- 1- punif(arg , -sqrt(3) ,sqrt(3))
		}

	}
	if(p==1){roc <- 1}
	if(p==0){roc <- 0}

	 
  	return(roc)
}


auc_x_verdadera<-function(XCOMP, alfa_d, alfa_h, beta_d, beta_h, sigma_d, sigma_h){
	###########################################################################
	# XCOMP, beta_d and beta_h : vectors of compositions of the same dimension
	# alfa_d and alfa_h correspond to the intercept
	# p is a real number where the ROC_x is computed
	###########################################################################

	mu_d <- alfa_d + prod.interno.atch(beta_d,XCOMP)
	mu_h <- alfa_h + prod.interno.atch(beta_h,XCOMP)

	diferencia <- (mu_h-mu_d)
	denominador <- sqrt(sigma_d^2+sigma_h^2)
       
	argumento <- diferencia/denominador
	auc <- 1-pnorm(argumento)
  	return(auc)
}

#####################################################
# ATCHINSON INNER PRODUCT
#####################################################
## INPUT
## x	:vector of compositional data length d
## y	:vector of compositional data length d
## ------------------------------------------------------
## OUTPUT
## ------------------------------------------------------
## Returns the atchinson inner product between x and y
## ------------------------------------------------------

########################################
# ONLY ADMITS VECTORS, NOT MATRICES
########################################
prod.interno.atch <- function(x,y){
	d <- length(x)
	producto=NA
	if(length(y)!=d){
		print("ERROR")
	}
	if(length(y)==d ){
		if(min(x)!=0 & min(y)!=0){
			clrx <- as.vector(unname(clr(x)))
			clry <- as.vector(unname(clr(y)))
			producto <- sum(clrx*clry)

			##############################################
			# ES IGUAL A sum(ilr(x)*ilr(y))
			##############################################
		}
	}
		return(producto)

}

####################################
# PARA HACER PELICULA
###################################

######################################################################
# ESTA FUNCION ES COMO movie3d PERO FACILITA PEGAR EN LATEX
######################################################################
 peli<-   function (f, duration, dev = cur3d(), ..., fps = 10, movie = "movie", 
    frames = movie, dir =  getwd(), convert = NULL, clean = TRUE, 
    verbose = TRUE, top = !rgl.useNULL(), type = "gif", startTime = 0, 
    webshot = TRUE) 
{
    olddir <- setwd(dir)
    on.exit(setwd(olddir))
    for (i in round(startTime * fps):(duration * fps)) {
        time <- i/fps
        if (cur3d() != dev) 
            set3d(dev)
        stopifnot(cur3d() != 0)
        args <- f(time, ...)
        subs <- args$subscene
        if (is.null(subs)) 
            subs <- currentSubscene3d(dev)
        else args$subscene <- NULL
        for (s in subs) par3d(args, subscene = s)
        filename <- paste(frames,i, '.png', sep='')
        if (verbose) {
            cat(gettextf("Writing '%s'\r", filename))
            flush.console()
        }
        if (top) 
            rgl.bringtotop()
        snapshot3d(filename = filename, webshot = webshot)
    }
    cat("\n")
    if (.Platform$OS.type == "windows") 
             system <- shell
    if (is.null(convert) && requireNamespace("magick", quietly = TRUE)) {
        m <- NULL
        for (i in round(startTime * fps):(duration * fps)) {
            filename <- paste(frames,i, '.png', sep='')
            frame <- magick::image_read(filename)
            if (is.null(m)) 
                m <- frame
            else m <- c(m, frame)
            if (clean) 
                unlink(filename)
        }
        m <- magick::image_animate(m, fps = fps, loop = 1, dispose = "previous")
        magick::image_write(m, paste0(movie, ".", type))
        return(invisible(m))
    }
    else if (is.null(convert)) {
        warning("R package 'magick' is not installed; trying external package.")
        convert <- TRUE
    }
    if (is.logical(convert) && convert) {
        progname <- "magick"
        version <- try(system2(progname, "--version", stdout = TRUE, 
            stderr = TRUE), silent = TRUE)
        if (inherits(version, "try-error") || !length(grep("ImageMagick",
version))) {
            progname <- "convert"
            version <- try(system2(progname, "--version", stdout = TRUE, 
                stderr = TRUE), silent = TRUE)
        }
        if (inherits(version, "try-error") || !length(grep("ImageMagick", 
            version))) 
            stop("'ImageMagick' not found")
        filename <- paste0(movie, ".", type)
        if (verbose) 
            cat(gettextf("Will create: %s\n", file.path(dir, 
                filename)))
        convert <- paste(progname, "-delay 1x%d %s*.png %s.%s")
    }
    if (is.character(convert)) {
        convert <- sprintf(convert, fps, frames, movie, type, 
            duration, dir)
        if (verbose) {
            cat(gettextf("Executing: '%s'\n", convert))
            flush.console()
        }
        system(convert)
        if (clean) {
            if (verbose) 
                cat(gettext("Deleting frames\n"))
            for (i in 0:(duration * fps)) {
                filename <- sprintf("%s%03d.png", frames, i)
                unlink(filename)
            }
        }
    }
    invisible(convert)
} 
