
library(compositions)
library(Compositional)

library("Ternary")
library("ggtern")

source('auxiliares.R')
source('auxiliares_ROC.R')

source('genero_unif.R')

simulacion_ECDF_ROC_UNIF<- function(Nrep=1000, grilla_p, ejex,ejey, nD, nH, 
                                       beta_D,beta_H, sigma_D, sigma_H, 
                                       alfa_D, alfa_H, 
                                       a1D=c(-2,2), a2D=c(-2,2), 
					a1H=c(-2,2),a2H=c(-2,2), 
				       dist_eps=list(type_H = 'gaussian', type_D = 'gaussian',gl_H=NA,gl_D=NA),		
				       carpeta='ROC_UNIF',
					archivo ='simu-ROC_UNIF.RData') {
  
  	largox=length(ejex)
	largoy=length(ejey)

	 
	lpes=length(grilla_p)


  	####################################
  	# CAMBIO A LA CARPETA DE GUARDADO
  	####################################
   

  	if (!dir.exists(carpeta)) dir.create(carpeta)
  
  	setwd(carpeta)
  
  	###########################################
  	###########################################
	grilla_region_x <- matrix(0,nrow=largox*largoy,ncol=3) 

	coord.tern.regionx <- matrix(0,nrow=largox ,ncol=largoy) 
	coord.tern.regiony <- matrix(0,nrow=largox ,ncol=largoy) 


	coordx_region  <- coordy_region<- rep(NA, length=largox*largoy)

	nx <- 0;
	j <- 1;
	while(j<=largox){
   		for (k in 1:largoy){
       		punto <- c(ejex[j], ejey[k]);
        		X.pred <- as.vector( unname(ilrInv(punto)))
	 		pepe <- as.vector(CoordinatesToXY(X.pred))
	 		grilla_region_x[k+nx,1:3] <-  X.pred # COORDENADAS COMPOSICIONALES
	 			
			####################################
			# GUARDO LAS COORDENADAS TERNARIAS
			# DE DOS FORMAS PARA DESPUES GRAFICAR
			###################################
			coordx_region[k+nx]<- pepe[1] 
			coordy_region[k+nx]<-pepe[2]

			coord.tern.regionx[j,k] <- pepe[1]
			coord.tern.regiony[j,k] <- pepe[2]
   		}
   	nx <- nx+largoy;
   	j <- j+1;
	}

  	#########################################
  	# GRAFICO GRILLA
  	#########################################
  	pdf('ternary-plot-grid.pdf',bg='transparent')
  
  	par(mar = rep(0.2, 4))
  
  	TernaryPlot(axis.labels = seq(0, 1, by = 0.1))
  
  	AddToTernary(graphics::points,  grilla_region_x[,1:3], col="gold", pch = 1, cex = 0.5)


  	dev.off()
  
  	############################################
  	#ROC VERDADERA
  	############################################

	dimension_ROC_region <- c(largox,largoy,lpes)

	TRUE_ROC_x_region<-array(rep(NA,largox*largoy*lpes),dim= dimension_ROC_region)

	TRUE_AUC_x_region<-  matrix(NA, ncol=largox,nrow=largoy) 

 
	for(j in 1:largox){
		for(k in 1:largoy){
    	 		punto <- c(ejex[j], ejey[k]);
       			punto_x <- as.vector( unname(ilrInv(punto)))
 			for(jp in 1:lpes){
      				p<-grilla_p[jp]
		
				TRUE_ROC_x_region[j,k,jp] <- roc_x_verdadera(p,XCOMP=punto_x, 
				alfa_d=alfa_D, alfa_h=alfa_H, beta_d=beta_D, beta_h=beta_H, sigma_d=sigma_D, 
				sigma_h=sigma_H, dist_eps=dist_eps)
			}
		TRUE_AUC_x_region[j,k] <- mean(TRUE_ROC_x_region[j,k,])
		}
	} 

	############################################
  	# GUARDO AUC REAL y EJES
  	############################################

	archivox<- 'ejex-modelo1.txt'

	archivoy<- 'ejey-modelo1.txt'
 

	write(t(ejex),file=archivox,ncolumns=largox,append=T)

	write(t(ejey),file=archivoy,ncolumns=largoy,append=T)


	ARCHIVO_AUC_REAL <- 'AUC_REAL_modelo1.txt'

	vector_AUC <- as.vector(TRUE_AUC_x_region)

	lvec_AUC <- length(vector_AUC)
	write(t(vector_AUC),file=ARCHIVO_AUC_REAL,ncolumns=lvec_AUC)



  	############################################
  	# ARCHIVOS DE GUARDADO
  	############################################

	ARCHIVO_EST <- 'Est_beta_alfa_modelo1_unif.txt'

	ARCHIVO_AUC_EST <- 'AUC_est_modelo1_unif.txt'

	ARCHIVO_MSE_ROC <- 'MSE_ROC_x_modelo1_unif.txt'
 	ARCHIVO_MAXDIF_ROC <- 'MAXDIF_ROC_x_modelo1_unif.txt'


 	ARCHIVO_AUC_EST_PI <- 'AUC_PI_est_modelo1_unif.txt'

	ARCHIVO_MSE_ROC_PI <- 'MSE_PI_ROC_x_modelo1_unif.txt'
 	ARCHIVO_MAXDIF_ROC_PI <- 'MAXDIF_PI_ROC_x_modelo1_unif.txt'


  	############################################
  	#Contenedores de los estimadores (en forma de matrices)
  	############################################
	
	dimension_Arreglo_est<-c(Nrep, largox,largoy, lpes)
  	ROC_CL_est <- array(NA,dim= dimension_Arreglo_est)
  	AUC_CL_est <- array(NA,dim= dimension_Arreglo_est[1:3])
  
	MSE_x_simu <- MAXDIF_x_simu <- array(NA,dim= dimension_Arreglo_est[1:3])
  
	MSE_simu <- MAXDIF_simu <- rep(NA, length=Nrep)

	MSE_AUC <- MAXDIF_AUC<- rep(NA, length=Nrep)


	ROC_PI_est <- array(NA,dim= dimension_Arreglo_est)
  	AUC_PI_est <- array(NA,dim= dimension_Arreglo_est[1:3])
  
	MSE_x_simu_PI <- MAXDIF_x_simu_PI <- array(NA,dim= dimension_Arreglo_est[1:3])
  
	MSE_simu_PI <- MAXDIF_simu_PI <- rep(NA, length=Nrep)

	MSE_AUC_PI <- MAXDIF_AUC_PI<- rep(NA, length=Nrep)

  	############################################
  	#EMPIEZA LA SIMULACIÓN
  	############################################
  	t1=Sys.time()
  
  	 
  
  	for(iter in 1:Nrep){
    
		print(iter)
		#####################
		#GENERO LOS DATOS  
                ####################

    		datos<-genero_unif(iter,nD, nH, alfa_D, alfa_H, beta_D, beta_H, a1D, a2D, a1H,a2H,  sigma_D, sigma_H, dist_eps=dist_eps)
    
		y_d<-datos$yD
    		xcomp_d<-datos$xcomp_D
    		y_h<-datos$yH
    		xcomp_h<-datos$xcomp_H
    
    		#################################
		#ESTIMADOR ROC CONDICIONAL
		#################################
    
		ROC_cond<-CURVA_ROC_ECDF(Y_D=y_d, XCOMP_D= xcomp_d,Y_H=y_h, XCOMP_H=xcomp_h, 
                           grilla_x=grilla_region_x[,1:3], grilla_p=grilla_p)
  
    		est_ROC_cond<-ROC_cond$ROC_CL_x_grilla
    		est_AUC<-ROC_cond$AUC_CL_x
 
		est_PI_ROC_cond<-ROC_cond$ROC_PI_x_grilla
    		est_PI_AUC<-ROC_cond$AUC_PI_x


		beta_est_COMP_D<- ROC_cond$betaCOMP_D 
		beta_est_COMP_H<- ROC_cond$betaCOMP_H 

		alfa_est_D<- ROC_cond$alfa_D 
		alfa_est_H<- ROC_cond$alfa_H
 
   		#######################################################
		#Guardo la info de la estimación de esta iteración
   		#######################################################
		############################################
		# ESTIMADORES DE ALFA y BETA
		###########################################
		vector_est <- c(iter,alfa_est_D,alfa_est_H,beta_est_COMP_D,beta_est_COMP_H)

		lvec <- length(vector_est)
		write(t(vector_est),file=ARCHIVO_EST,ncolumns=lvec,append=T)


		############################################
		# ESTIMADORES DE LA ROC y AUC
		############################################

		nx <- 0;
		j <- 1;
		while(j<=largox){
   			for (k in 1:largoy){
				ROC_CL_est[iter,j,k,] <- est_ROC_cond[k+nx,]  
				AUC_CL_est[iter,j,k] <- est_AUC[k+nx] 
			 
 				ROC_PI_est[iter,j,k,] <- est_PI_ROC_cond[k+nx,]  
				AUC_PI_est[iter,j,k] <- est_PI_AUC[k+nx] 

   			}
   			nx <- nx+largoy;
   			j <- j+1;
		}


   		#######################################################
		#Guardo errores de estimación para cada x de la grilla
		#######################################################

		#######################################
		# BASADO EN EMPIRICA
		######################################

    		diferencia.ROC <- TRUE_ROC_x_region-ROC_CL_est[iter, , ,]
    		MSE_x<- MAXDIF_x <-matrix(NA, ncol=largox,nrow=largoy) 
    		for(j in 1:largox){
			for(k in 1:largoy){
    	      			MSE_x[j,k] <- mean(diferencia.ROC[j,k,]^2)
      				MAXDIF_x[j,k] <- max(abs(diferencia.ROC[j,k,])) #Dist de Kolmogorov
    			}
		}

		MSE_x_simu[iter,, ] <- MSE_x
		MAXDIF_x_simu[iter,, ] <- MAXDIF_x

   
		MSE_simu[iter] <- mean(diferencia.ROC^2)
		MAXDIF_simu[iter] <-  max(abs(diferencia.ROC))
		
	
		diferencia.AUC <- TRUE_AUC_x_region-AUC_CL_est[iter, , ]
		MSE_AUC[iter] <- mean(diferencia.AUC^2)
		MAXDIF_AUC[iter] <- max(abs(diferencia.AUC))
 

		#######################################
		# GUARDO EN ARCHIVOS
		#######################################



		vector_AUC_est <- as.vector(AUC_CL_est[iter, , ] )

		lvec_AUC <- length(vector_AUC_est)
		write(t(vector_AUC_est),file=ARCHIVO_AUC_EST,ncolumns=lvec_AUC,append=T)

		vector_MSE_x <- as.vector(MSE_x_simu[iter, , ] )

		lvec_MSE <- length(vector_MSE_x)
		write(t(vector_MSE_x),file=ARCHIVO_MSE_ROC,ncolumns=lvec_MSE,append=T)

		vector_MAXDIF_x <- as.vector(MAXDIF_x_simu[iter, , ] )

		lvec_MAXDIF <- length(vector_MAXDIF_x)
		write(t(vector_MAXDIF_x),file=ARCHIVO_MAXDIF_ROC,ncolumns=lvec_MAXDIF,append=T)


		#######################################
		# BASADO EN PLUG-IN BINORMAL
		######################################

    		diferencia.ROC_PI <- TRUE_ROC_x_region-ROC_PI_est[iter, , ,]
    		MSE_x_PI<- MAXDIF_x_PI <-matrix(NA, ncol=largox,nrow=largoy) 
    		for(j in 1:largox){
			for(k in 1:largoy){
    	      			MSE_x_PI[j,k] <- mean(diferencia.ROC_PI[j,k,]^2)
      				MAXDIF_x_PI[j,k] <- max(abs(diferencia.ROC_PI[j,k,])) #Dist de Kolmogorov
    			}
		}

		MSE_x_simu_PI[iter,, ] <- MSE_x_PI
		MAXDIF_x_simu_PI[iter,, ] <- MAXDIF_x_PI

   
		MSE_simu_PI[iter] <- mean(diferencia.ROC_PI^2)
		MAXDIF_simu_PI[iter] <-  max(abs(diferencia.ROC_PI))
		
	
		diferencia.AUC_PI <- TRUE_AUC_x_region-AUC_PI_est[iter, , ]
		MSE_AUC_PI[iter] <- mean(diferencia.AUC_PI^2)
		MAXDIF_AUC_PI[iter] <- max(abs(diferencia.AUC_PI))
 

		#######################################
		# GUARDO EN ARCHIVOS
		#######################################



		vector_AUC_PI <- as.vector(AUC_PI_est[iter, , ] )

		lvec_AUC_PI <- length(vector_AUC_PI)
		write(t(vector_AUC_PI),file=ARCHIVO_AUC_EST_PI,ncolumns=lvec_AUC_PI,append=T)

		vector_MSE_x_PI <- as.vector(MSE_x_simu_PI[iter, , ] )

		lvec_MSE_PI <- length(vector_MSE_x_PI)
		write(t(vector_MSE_x_PI),file=ARCHIVO_MSE_ROC_PI,ncolumns=lvec_MSE_PI,append=T)

		vector_MAXDIF_x_PI <- as.vector(MAXDIF_x_simu_PI[iter, , ] )

		lvec_MAXDIF_PI <- length(vector_MAXDIF_x_PI)
		write(t(vector_MAXDIF_x_PI),file=ARCHIVO_MAXDIF_ROC_PI,ncolumns=lvec_MAXDIF_PI,append=T)

 


      	}

  	t2=Sys.time()
  
  	tiempo=t2-t1
  	print(c('tardo= ', tiempo))

 
	
	save(nD, nH, alfa_D, alfa_H, beta_D, beta_H, 
            	a1D, a2D, a1H,a2H, 
		sigma_D, sigma_H, grilla_p, ejex, ejey,
		grilla_region_x, coordx_region, coordy_region,
		coord.tern.regionx, coord.tern.regiony,
		TRUE_ROC_x_region, ROC_CL_est, ROC_PI_est,
		TRUE_AUC_x_region,AUC_CL_est,  AUC_PI_est,
		MSE_x_simu, MAXDIF_x_simu, MSE_simu, MAXDIF_simu, MSE_AUC, MAXDIF_AUC,
		MSE_x_simu_PI, MAXDIF_x_simu_PI, MSE_simu_PI, MAXDIF_simu_PI, MSE_AUC_PI, MAXDIF_AUC_PI,
		 file=archivo)  

   }
