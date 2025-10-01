rm(list=ls())
 
source('simula_BOTH_unif_GRACIELA.R')

#Armo grilla_p
grilla_p=seq(0, 1, length=101)
 
#####################################
#Defino los parámetros
# ESTE MODELO se llama Modelo1
#####################################
nD<-50
nH<-50

beta_D<-c(0.1, 0.3, 0.6)
beta_H<-c(0.05, 0.55, 0.4)
 

sigma_D<- 1
sigma_H <- 1


beta_D<-  c(0.1, 0.3, 0.6)  
beta_H<- c(0.05, 0.55, 0.4)  

alfa_D <- 4
alfa_H <- 1
 
a1D=c(-2,2)

a2D=c(-2,2)

a1H=c(-2,2)

a2H=c(-2,2)
  

largo=50

ejex.inf <- max(a1D[1],a1H[1])
ejex.sup <- min(a1D[2],a1H[2])


ejey.inf <- max(a2D[1],a2H[1])
ejey.sup <- min(a2D[2],a2H[2])

ejex= seq(ejex.inf, ejex.sup,length=largo) 
ejey= seq(ejey.inf, ejey.sup,length=largo)  


Nrep=1000

 

########################################
# STUDENT-STUDENT
########################################


dist_eps =list(type_H ='Student', type_D='Student', gl_H=4, gl_D=4)

#############################
# DONDE GUARDO
# NOMBRE ARCHIVO
############################



carpeta=paste('ROC_modelo1_UNIF_nD_', nD, '_n_H_', nH, '_Nrep_',Nrep,'_',
	dist_eps$type_H,'_', dist_eps$type_D, '_',
	dist_eps$gl_H,'_', dist_eps$gl_D,sep='')

archivo =paste('simu-ROC_modelo1_UNIF_',
		dist_eps$type_H,'_', dist_eps$type_D, '_',
		dist_eps$gl_H,'_', dist_eps$gl_D,
		'.RData',sep='')

 



t1 <- Sys.time()

simulacion_BOTH_ROC_UNIF(Nrep=Nrep, grilla_p=grilla_p, ejex=ejex,ejey=ejey, 
			nD=nD, nH=nH, 
                        beta_D=beta_D,beta_H=beta_H, 
			sigma_D=sigma_D, sigma_H=sigma_H, 
			alfa_D=alfa_D, alfa_H=alfa_H, 
			a1D=a1D, a2D=a2D, a1H=a1H,a2H=a2H,
			dist_eps=dist_eps,
			carpeta=carpeta, archivo=archivo, ache=NULL, nucleo="epan")  
 

 
t2 <- Sys.time()
t2-t1



 

########################################
# Logistic-Logistic
########################################


dist_eps =list(type_H ='logistic', type_D='logistic', gl_H=NA, gl_D=NA)

#############################
# DONDE GUARDO
# NOMBRE ARCHIVO
############################



carpeta=paste('ROC_modelo1_UNIF_nD_', nD, '_n_H_', nH, '_Nrep_',Nrep,'_',
	dist_eps$type_H,'_', dist_eps$type_D, '_',
	dist_eps$gl_H,'_', dist_eps$gl_D,sep='')

archivo =paste('simu-ROC_modelo1_UNIF_',
		dist_eps$type_H,'_', dist_eps$type_D, '_',
		dist_eps$gl_H,'_', dist_eps$gl_D,
		'.RData',sep='')

 



t1 <- Sys.time()

simulacion_BOTH_ROC_UNIF(Nrep=Nrep, grilla_p=grilla_p, ejex=ejex,ejey=ejey, 
			nD=nD, nH=nH, 
                        beta_D=beta_D,beta_H=beta_H, 
			sigma_D=sigma_D, sigma_H=sigma_H, 
			alfa_D=alfa_D, alfa_H=alfa_H, 
			a1D=a1D, a2D=a2D, a1H=a1H,a2H=a2H,
			dist_eps=dist_eps,
			carpeta=carpeta, archivo=archivo, ache=NULL, nucleo="epan")  
 

 
t2 <- Sys.time()
t2-t1

 
