rm(list=ls())
 
source('simula_ECDF_unif_GRACIELA.R')

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


dist_eps =list(type_H ='gaussian', type_D='gaussian', gl_H=NA, gl_D=NA)



Nrep=1000


#############################
# DONDE GUARDO
# NOMBRE ARCHIVO
############################


carpeta=paste('ROC_modelo1_UNIF_ECDF_nD_', nD, '_n_H_', nH, 
	dist_eps$type_H,'_', dist_eps$type_D, '_',
	dist_eps$gl_H,'_', dist_eps$gl_D,sep='')

archivo =paste('simu-ROC_modelo1_UNI_ECDF_',
		dist_eps$type_H,'_', dist_eps$type_D, '_',
		dist_eps$gl_H,'_', dist_eps$gl_D,
		'.RData',sep='')

 


t1 <- Sys.time()

simulacion_ECDF_ROC_UNIF(Nrep=Nrep, grilla_p=grilla_p, ejex=ejex,ejey=ejey, 
			nD=nD, nH=nH, 
                        beta_D=beta_D,beta_H=beta_H, 
			sigma_D=sigma_D, sigma_H=sigma_H, 
			alfa_D=alfa_D, alfa_H=alfa_H, 
			a1D=a1D, a2D=a2D, a2H=a1H,a2H=a2H,
			dist_eps=dist_eps,
			carpeta=carpeta, archivo=archivo)
 

 

t2 <- Sys.time()
t2-t1

 




 



