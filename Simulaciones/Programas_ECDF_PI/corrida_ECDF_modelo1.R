rm(list=ls())
 
source('simula_ECDF_GRACIELA.R')

#Armo grilla_p
grilla_p=seq(0, 1, length=101)
 
#####################################
#Defino los parámetros
# ESTE MODELO se llama Modelo1
#####################################
nD<-300
nH<-300

beta_D<-c(0.1, 0.3, 0.6)
beta_H<-c(0.05, 0.55, 0.4)
 

sigma_D<- 1
sigma_H <- 1


beta_D<-  c(0.1, 0.3, 0.6)  
beta_H<- c(0.05, 0.55, 0.4)  

alfa_D <- 4
alfa_H <- 1
 

alphadiriD <- c(3,7,10)
alphadiriH <- c(3,7,10)
 
  

largo=50

ejex= seq(-1.4,3.2,length=largo) 
ejey= seq(-0.8,2.6,length=largo)  

Nrep=1000


dist_eps =list(type_H ='gaussian', type_D='gaussian', gl_H=NA, gl_D=NA)

#############################
# DONDE GUARDO
# NOMBRE ARCHIVO
############################


carpeta=paste('ROC_modelo1_ECDF_nD_', nD, '_n_H_', nH, '_',
	dist_eps$type_H,'_', dist_eps$type_D, '_',
	dist_eps$gl_H,'_', dist_eps$gl_D,sep='')

archivo =paste('simu-ROC_modelo1_ECDF_',
		dist_eps$type_H,'_', dist_eps$type_D, '_',
		dist_eps$gl_H,'_', dist_eps$gl_D,
		'.RData',sep='')



t1 <- Sys.time()

simulacion_ECDF_ROC_Dirichlet(Nrep=Nrep, grilla_p=grilla_p, ejex=ejex,ejey=ejey, 
				nD=nD, nH=nH, 
                        beta_D=beta_D,beta_H=beta_H, 
				sigma_D=sigma_D, sigma_H=sigma_H, 
				alfa_D=alfa_D, alfa_H=alfa_H, 
				alphadiriD=alphadiriD, alphadiriH=alphadiriH, 
				dist_eps=dist_eps,
				carpeta=carpeta, archivo=archivo)  
 
 
t2 <- Sys.time()
t2-t1

 




 



