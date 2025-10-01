
library(compositions)
library(Compositional)

library("Ternary")
library("ggtern")


source('auxiliares.R')
source('auxiliares_ROC.R')

source('genero.R')

simulacion_SMOOTH_ROC_Dirichlet<- function(Nrep=1000, grilla_p, ejex,ejey, nD, nH, 
                                       beta_D,beta_H, sigma_D, sigma_H, 
                                       alfa_D, alfa_H, 
                                       alphadiriD, alphadiriH, 
				       dist_eps=list(type_H = 'gaussian', type_D = 'gaussian',gl_H=NA,gl_D=NA),
				       carpeta='ROC_prueba',archivo ='simu-ROC_prueba.RData', ache=NULL, nucleo="epan") {
  
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
				alfa_d=alfa_D, alfa_h=alfa_H, beta_d=beta_D, beta_h=beta_H, 
				sigma_d=sigma_D, sigma_h=sigma_H, dist_eps=dist_eps)
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


	ARCHIVO_AUC_REAL <- 'AUC_REAL_modelo1_unif.txt'

	vector_AUC <- as.vector(TRUE_AUC_x_region)

	lvec_AUC <- length(vector_AUC)
	write(t(vector_AUC),file=ARCHIVO_AUC_REAL,ncolumns=lvec_AUC)



  	############################################
  	# ARCHIVOS DE GUARDADO
  	############################################


	ARCHIVO_EST <- 'Est_beta_alfa_modelo1.txt'
 
	ARCHIVO_AUC_EST_SMOOTH <- 'AUC_SMOOTH_est_modelo1.txt'

	ARCHIVO_MSE_ROC_SMOOTH <- 'MSE_SMOOTH_ROC_x_modelo1.txt'
 	ARCHIVO_MAXDIF_ROC_SMOOTH <- 'MAXDIF_SMOOTH_ROC_x_modelo1.txt'

 

  	############################################
  	#Contenedores de los estimadores (en forma de matrices)
  	############################################
	
	dimension_Arreglo_est<-c(Nrep, largox,largoy, lpes)
 
	ROC_SMOOTH_est <- array(NA,dim= dimension_Arreglo_est)
  	AUC_SMOOTH_est <- array(NA,dim= dimension_Arreglo_est[1:3])
  
	MSE_x_simu_SMOOTH <- MAXDIF_x_simu_SMOOTH <- array(NA,dim= dimension_Arreglo_est[1:3])
  
	MSE_simu_SMOOTH <- MAXDIF_simu_SMOOTH <- rep(NA, length=Nrep)

	MSE_AUC_SMOOTH <- MAXDIF_AUC_SMOOTH<- rep(NA, length=Nrep)


  	############################################
  	#EMPIEZA LA SIMULACIÓN
  	############################################
  	t1=Sys.time()
  
  	 
  
  	for(iter in 1:Nrep){
    
		print(iter)
		#####################
		#GENERO LOS DATOS  
                ####################

      		datos<-genero(iter,nD, nH, alfa_D, alfa_H, beta_D, beta_H, 
                  alphadiriD, alphadiriH, sigma_D, sigma_H, dist_eps=dist_eps)
 
		y_d<-datos$yD
    		xcomp_d<-datos$xcomp_D
    		y_h<-datos$yH
    		xcomp_h<-datos$xcomp_H
    
    		#################################
		#ESTIMADOR ROC CONDICIONAL
		#################################
    
		ROC_cond<-CURVA_ROC_SMOOTH(Y_D=y_d, XCOMP_D= xcomp_d,Y_H=y_h, XCOMP_H=xcomp_h, 
                           grilla_x=grilla_region_x[,1:3], grilla_p=grilla_p, ache=ache, nucleo=nucleo)
 
		est_SMOOTH_ROC_cond<-ROC_cond$ROC_SMOOTH_x_grilla
    		est_SMOOTH_AUC<-ROC_cond$AUC_SMOOTH_x
 
   		
		
   		
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
			 
				ROC_SMOOTH_est[iter,j,k,] <- est_SMOOTH_ROC_cond[k+nx,]  
				AUC_SMOOTH_est[iter,j,k] <- est_SMOOTH_AUC[k+nx] 

   			}
   			nx <- nx+largoy;
   			j <- j+1;
		}


   		#######################################################
		#Guardo errores de estimación para cada x de la grilla
		#######################################################

		

		#######################################
		# BASADO EN Smoothing 
		######################################

    		diferencia.ROC_SMOOTH <- TRUE_ROC_x_region-ROC_SMOOTH_est[iter, , ,]
    		MSE_x_SMOOTH<- MAXDIF_x_SMOOTH <-matrix(NA, ncol=largox,nrow=largoy) 
    		for(j in 1:largox){
			for(k in 1:largoy){
    	      			MSE_x_SMOOTH[j,k] <- mean(diferencia.ROC_SMOOTH[j,k,]^2)
      				MAXDIF_x_SMOOTH[j,k] <- max(abs(diferencia.ROC_SMOOTH[j,k,])) #Dist de Kolmogorov
    			}
		}

		MSE_x_simu_SMOOTH[iter,, ] <- MSE_x_SMOOTH
		MAXDIF_x_simu_SMOOTH[iter,, ] <- MAXDIF_x_SMOOTH

   
		MSE_simu_SMOOTH[iter] <- mean(diferencia.ROC_SMOOTH^2)
		MAXDIF_simu_SMOOTH[iter] <-  max(abs(diferencia.ROC_SMOOTH))
		
	
		diferencia.AUC_SMOOTH <- TRUE_AUC_x_region-AUC_SMOOTH_est[iter, , ]
		MSE_AUC_SMOOTH[iter] <- mean(diferencia.AUC_SMOOTH^2)
		MAXDIF_AUC_SMOOTH[iter] <- max(abs(diferencia.AUC_SMOOTH))
 

		#######################################
		# GUARDO EN ARCHIVOS
		#######################################



		vector_AUC_SMOOTH <- as.vector(AUC_SMOOTH_est[iter, , ] )

		lvec_AUC_SMOOTH <- length(vector_AUC_SMOOTH)
		write(t(vector_AUC_SMOOTH),file=ARCHIVO_AUC_EST_SMOOTH,ncolumns=lvec_AUC_SMOOTH,append=T)

		vector_MSE_x_SMOOTH <- as.vector(MSE_x_simu_SMOOTH[iter, , ] )

		lvec_MSE_SMOOTH <- length(vector_MSE_x_SMOOTH)
		write(t(vector_MSE_x_SMOOTH),file=ARCHIVO_MSE_ROC_SMOOTH,ncolumns=lvec_MSE_SMOOTH,append=T)

		vector_MAXDIF_x_SMOOTH <- as.vector(MAXDIF_x_simu_SMOOTH[iter, , ] )

		lvec_MAXDIF_SMOOTH <- length(vector_MAXDIF_x_SMOOTH)
		write(t(vector_MAXDIF_x_SMOOTH),file=ARCHIVO_MAXDIF_ROC_SMOOTH,ncolumns=lvec_MAXDIF_SMOOTH,append=T)



      	}

  	t2=Sys.time()
  
  	tiempo=t2-t1
  	print(c('tardo= ', tiempo))

 

	
	save(nD, nH, alfa_D, alfa_H, beta_D, beta_H, 
            	alphadiriD, alphadiriH, sigma_D, sigma_H, grilla_p, ejex, ejey,
		grilla_region_x, coordx_region, coordy_region,
		coord.tern.regionx, coord.tern.regiony,
		TRUE_ROC_x_region, ROC_SMOOTH_est,
		TRUE_AUC_x_region, AUC_SMOOTH_est, 
		MSE_x_simu_SMOOTH, MAXDIF_x_simu_SMOOTH, MSE_simu_SMOOTH, MAXDIF_simu_SMOOTH, MSE_AUC_SMOOTH, MAXDIF_AUC_SMOOTH,
		 file=archivo)  

   }
