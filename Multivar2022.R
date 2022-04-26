##Librerias necesarias (instalar con calma)
library(MASS)
library(rgl)
library(plotrix)
library(Rtsne)
#library(umap)
library(uwot) #otra implementacion
library(ica)
library(vegan)
library(ggplot2)
library(factoextra)
library(rattle)
library(e1071)
library(vegan)


#setwd("~/Documentos/Estadistica")

##El problema de la multidimensionalidad
covarianza1=matrix(byrow = T, c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
media1=c(1,1,1)
covarianza2=matrix(byrow = T, c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
media2=c(10,1,1)
covarianza3=matrix(byrow = T, c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
media3=c(1,10,1)

aleatorio1<-mvrnorm(n=100,mu=media1,Sigma=covarianza1)
aleatorio2<-mvrnorm(n=100,mu=media2,Sigma=covarianza2)
aleatorio3<-mvrnorm(n=100,mu=media3,Sigma=covarianza3)

tratamiento=as.factor(c(rep("tratamiento1",100),rep("tratamiento2",100),rep("tratamiento3",100)))

poblacion=rbind(aleatorio1,aleatorio2,aleatorio3)
colnames(poblacion)<-c("Variable1","Variable2","Variable3")
tabla=cbind(tratamiento,as.data.frame(poblacion))

dim(poblacion)

plot3d(x=poblacion[,1],xlab="x",
       y=poblacion[,2],ylab="y",
       z=poblacion[,3],zlab="z",
       type="p",col=as.numeric(tratamiento))

pairs(poblacion,col=as.numeric(tratamiento))

ggplot(tabla,aes(x=Variable1,y=Variable2))+
  geom_point(aes(col=tratamiento))+
  scale_color_manual(values=c("blue","red","yellow"))+
  geom_density2d(aes(col=tratamiento),bins=10)

ggplot(tabla,aes(x=Variable2,y=Variable3))+
  geom_point(aes(col=tratamiento))+
  scale_color_manual(values=c("blue","red","yellow"))+
  geom_density2d(aes(col=tratamiento),bins=10)

ggplot(tabla,aes(x=Variable2,y=Variable3))+
  geom_point(aes())+
  geom_density2d(aes(col=..level..),bins=10)

############################################################
###########Conclusiones ####################################
observaciones<-iris[,c(1:3)] 
tratamiento<-as.factor(iris$Species)
plot3d(x=observaciones[,1],
         y=observaciones[,2],
         z=observaciones[,3],
         xlim=c(-10,20),ylim=c(-15,20),zlim=c(-5,10),
       type="p",col=as.numeric(tratamiento))

covarianza=cov(observaciones)
media=apply(observaciones,2,mean)

aleatorio<-mvrnorm(n=1000,mu=media,Sigma=covarianza)

plot3d(x=aleatorio[,1],
       y=aleatorio[,2],
       z=aleatorio[,3],col="red",
       xlim=c(-10,20),ylim=c(-15,20),zlim=c(-5,10),add=T)

lines3d(x=c(-10,20),y=c(0,0),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(-15,20),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(0,0),z=c(-5,10),add=T)

observaciones_centrado<-scale(observaciones,scale = F)
#observaciones_centrado<-as.matrix(observaciones)

plot3d(x=observaciones_centrado[,1],
       y=observaciones_centrado[,2],
       z=observaciones_centrado[,3],col=as.numeric(tratamiento),
       xlim=c(-10,20),ylim=c(-15,20),zlim=c(-5,10),
       type="p")
lines3d(x=c(-10,20),y=c(0,0),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(-15,20),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(0,0),z=c(-5,10),add=T)

#https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca

covarianza2=cov(observaciones_centrado)
propios=eigen(covarianza2)

L_matriz=diag(propios$values,3,3)
V=as.matrix(propios$vectors)

round(covarianza2,3)==round(V%*%L_matriz%*%t(V),3)

arrow3d(p0=c(0,0,0),p1=5*V[,1],add=T,col="red")
arrow3d(p0=c(0,0,0),p1=5*V[,2],add=T,col="green")
arrow3d(p0=c(0,0,0),p1=5*V[,3],add=T,col="yellow")

acos(sum(V[,1]*V[,3])/(sqrt(sum(V[,1] * V[,1])) * sqrt(sum(V[,3]* V[,3]))))

proyeccion<-observaciones_centrado%*%V
plot(proyeccion[,1],proyeccion[,2],col=as.numeric(tratamiento),pch=16)


## Singular Value decomposition ##
#https://www.youtube.com/watch?v=mBcLRGuAFUk
plot3d(x=observaciones_centrado[,1],
       y=observaciones_centrado[,2],
       z=observaciones_centrado[,3],col="blue",
       xlim=c(-10,20),ylim=c(-15,20),zlim=c(-5,10),
       type="p")
lines3d(x=c(-10,20),y=c(0,0),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(-15,20),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(0,0),z=c(-5,10),add=T)

descomposicion<-svd(observaciones_centrado)

reconst<-descomposicion$u%*%diag(descomposicion$d)%*%t(descomposicion$v)

plot(descomposicion$u%*%diag(descomposicion$d),col=as.numeric(tratamiento),pch=16)
plot(descomposicion$u,col=as.numeric(tratamiento),pch=16)


dim(observaciones_centrado)
dim(descomposicion$u)
dim(diag(descomposicion$d))
dim(descomposicion$v)

round(observaciones_centrado,3)==round(reconst,3)

round(descomposicion$v,3)==round(V,3)

observaciones_centrado%*%descomposicion$v
descomposicion$u%*%diag(descomposicion$d)

plot(descomposicion$v,col=as.numeric(tratamiento),pch=16)

## GENERAR DATOS CON CURVATURA
theta=runif(500,0,1.4*pi)#500 números de 0 a pi radianes (theta)
r=rep(1,500)#500 números 1 (r)
ys=r*sin(theta)#y coordenadas rectangulares
xs=r*cos(theta)#x coordenadas rectangulares
plot(ys,xs)#graficar
zs=rnorm(500,2.5,1)#500 valores z
trata2=vector(length = 500)
trata2[zs>=2.5]<-"grupo1"
trata2[!zs>=2.5]<-"grupo2"
trata2<-as.factor(trata2)
trata3=vector(length = 500)
trata3[ys>=0.6]<-"grupo1"
trata3[!ys>=0.6]<-"grupo2"
trata3<-as.factor(trata3)
trata4=vector(length = 500)
trata4[xs>=0]<-"grupo1"
trata4[!xs>=0]<-"grupo2"
trata4<-as.factor(trata4)

coordenadas<-matrix(byrow = T,nrow =3,ncol=500, c(xs,ys,zs))
dim(coordenadas)
plot3d(x=coordenadas[1,],
       y=coordenadas[2,],
       z=coordenadas[3,],
       type="p",col=as.numeric(trata4))

Mz=function(angulo){
  matrix(byrow = T, ncol = 3, nrow = 3, c(cos(angulo),-sin(angulo),0,
                                          sin(angulo),cos(angulo),0,
                                          0,0,1))
}
My=function(angulo){
  matrix(byrow = T, ncol = 3, nrow = 3, c(cos(angulo),0,sin(angulo),
                                          0,1,0,
                                          -sin(angulo),0,cos(angulo)))
}
Mx=function(angulo){
  matrix(byrow = T, ncol = 3, nrow = 3, c(1,0,0,
                                          0,cos(angulo),-sin(angulo),
                                          0,sin(angulo),cos(angulo)))
}

coordenadas2=Mz(0.5)%*%My(0.8)%*%Mx(0)%*%coordenadas
plot3d(x=coordenadas2[1,],xlab="x",xlim = c(-4,4),
       y=coordenadas2[2,],ylab="y",ylim = c(-4,4),
       z=coordenadas2[3,],zlab="z",zlim = c(-4,4),
       type="p",col=as.numeric(trata4))

lines3d(x=c(-4,4),y=c(0,0),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(-4,4),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(0,0),z=c(-4,4),add=T)

coordenadas3<-scale(t(coordenadas2),scale = F)

descomposicion<-svd(coordenadas3)
reconst<-descomposicion$u%*%diag(descomposicion$d)%*%t(descomposicion$v)
round(coordenadas3,3)==round(reconst,3)

descomposicion$v

plot3d(x=coordenadas3[,1],xlab="x",xlim = c(-4,4),
       y=coordenadas3[,2],ylab="y",ylim = c(-4,4),
       z=coordenadas3[,3],zlab="z",zlim = c(-4,4),
       type="p",col=as.numeric(trata4))

lines3d(x=c(-4,4),y=c(0,0),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(-4,4),z=c(0,0),add=T)
lines3d(x=c(0,0),y=c(0,0),z=c(-4,4),add=T)

arrow3d(p0=c(0,0,0),p1=5*descomposicion$v[,1],add=T,col="red")
arrow3d(p0=c(0,0,0),p1=5*descomposicion$v[,2],add=T,col="green")
arrow3d(p0=c(0,0,0),p1=5*descomposicion$v[,3],add=T,col="yellow")

proyeccion2<-descomposicion$u%*%diag(descomposicion$d)
proyeccion3<-coordenadas3%*%descomposicion$v
plot(proyeccion2[,1],proyeccion2[,2],col=as.numeric(trata4),pch=16)
plot(proyeccion3[,1],proyeccion3[,2],col=as.numeric(trata4),pch=16)

########### Conclusiones ###############
########################################
###Reducción dimensional no lineal
#Multidimensional scaling (MDS) Principal Coordinates Analysis
matriz_distancia<-dist(t(coordenadas2),method = "euclidean",diag=T,upper = T)
dim(as.matrix(matriz_distancia))
MDS_result <- cmdscale(matriz_distancia, eig = TRUE, k = 2)
plot(MDS_result$points[,1],MDS_result$points[,2],col=as.numeric(trata4),pch=16)

##Independent component analysis
library(ica)
ica_resultado <- icafast(t(coordenadas2), 2,
                           center = TRUE, maxit = 100,
                           tol = 1e-6)
plot(ica_resultado$Y[,1],ica_resultado$Y[,2],col=as.numeric(trata4),pch=16)

##t-Distributed Stochastic Neighbor Embedding (t-SNE)
#https://www.analyticsvidhya.com/blog/2017/01/t-sne-implementation-r-python/
#https://towardsdatascience.com/t-distributed-stochastic-neighbor-embedding-t-sne-bb60ff109561
tsne_resultado<- Rtsne(t(coordenadas2),pca = FALSE, perplexity = 12,theta = 0.0,
                       max_iter=2000)
plot(tsne_resultado$Y[,1],tsne_resultado$Y[,2],col=as.numeric(trata4),pch=16)


#Uniform manifold approximation and projection (UMAP)
#https://umap-learn.readthedocs.io/en/latest/
#https://cran.r-project.org/web/packages/umap/vignettes/umap.html
umap_resultado = umap(t(coordenadas2))
plot(umap_resultado$layout[,1],umap_resultado$layout[,2],col=as.numeric(trata4),pch=16)

library(uwot) #otra implementacion
umap_resultado2<- umap(t(coordenadas2),n_neighbors = 15,min_dist = 0.2, spread = 1)
plot(umap_resultado2[,1],umap_resultado2[,2],col=as.numeric(trata4),pch=16)

########################################################
############Conclusiones################################
##########EXPLORACIÓN BASE DE DATOS#####################
matriz_metabo<-read.csv("metabo_hojas.csv",header = T)
dim(matriz_metabo)
head(names(matriz_metabo))
##[1] "Hoja"         "Inoculo"      "Combinado"    "X50.83771282"
hoja=as.factor(matriz_metabo$Hoja)
matriz1=scale(as.matrix(matriz_metabo[,-c(1)]),scale=F)
dim(matriz1)
##Singular value decomposition (para PCA)
pca=svd(matriz1)
dim(pca$u)
dim(diag(pca$d))
dim(pca$v)
proyeccion=pca$u%*%diag(pca$d)
plot(pca$u,col=as.numeric(hoja),pch=16)
plot(proyeccion,col=as.numeric(hoja),pch=16)##estos son las unidades experimentales (llamados scores)
plot(pca$v)##estos son los metabolitos en el plano ordenado (llamados loadings)
head(pca$u)

plot(pca$u,col=as.numeric(hoja),pch=16)
plot(proyeccion[,1],proyeccion[,2],col=as.numeric(hoja),pch=16)

contribucion=pca$d/sum(pca$d)
plot(contribucion,type = "b")
plot(cumsum(contribucion), type="b")

plot3d(x=proyeccion[,1],xlab="x",
       y=proyeccion[,2],ylab="y",
       z=proyeccion[,3],zlab="z",
       type="p",col=as.numeric(hoja))

elipse_H0 = ellipse3d(cov(proyeccion[hoja=="H0",]), 
                    centre=c(mean(proyeccion[hoja=="H0",1]), 
                             mean(proyeccion[hoja=="H0",2]), 
                             mean(proyeccion[hoja=="H0",3])), level = 0.75)

elipse_H1 = ellipse3d(cov(proyeccion[hoja=="H1",]), 
                      centre=c(mean(proyeccion[hoja=="H1",1]), 
                               mean(proyeccion[hoja=="H1",2]), 
                               mean(proyeccion[hoja=="H1",3])), level = 0.75)

elipse_H2 = ellipse3d(cov(proyeccion[hoja=="H2",]), 
                      centre=c(mean(proyeccion[hoja=="H2",1]), 
                               mean(proyeccion[hoja=="H2",2]), 
                               mean(proyeccion[hoja=="H2",3])), level = 0.75)

shade3d(elipse_H0, col = "#D95F02", alpha = 0.1, lit = FALSE)
wire3d(elipse_H0, col = "#D95F02",  lit = FALSE)

shade3d(elipse_H1, col = "#00FF00", alpha = 0.1, lit = FALSE)
wire3d(elipse_H1, col = "#00FF00",  lit = FALSE)

shade3d(elipse_H2, col = "#EE6AA7", alpha = 0.1, lit = FALSE)
wire3d(elipse_H2, col = "#EE6AA7",  lit = FALSE)

###con la traspuesta
pca2=svd(t(matriz1))
dim(t(matriz1)) #tratamientos=columnas, filas= metabolitos
dim(pca2$u)
dim(diag(pca2$d))
dim(pca2$v)
proyeccion2=pca2$u%*%diag(pca2$d)
plot(pca2$u)##ahora son los scores (metabolitos) PARA LA MATRIZ TRANSPUESTA
plot(pca2$v,col=as.numeric(hoja),pch=16)#ahora son los loadings(unidades experim)
head(pca2$v)#filas=tratamientos,columas:variables latentes (eigenmetabolitos)
#plot(proyeccion2[,1],proyeccion2[,2],col=as.numeric(hoja),pch=16)
#ahora lo que se debe de interpretar es la matriz v final (pca2$v)


library(vegan)
pca2<-rda(matriz1,scale=F)
plot(pca2$CA$u[,1],pca2$CA$u[,2],col=as.numeric(hoja),pch=16)##Es lo mismo

umap_resultado2<- umap(n_components = 3,matriz1,n_neighbors = 20,min_dist = 0.2, spread = 1)
plot(umap_resultado2[,1],umap_resultado2[,2],col=as.numeric(hoja),pch=16,cex=1)

plot3d(x=umap_resultado2[,1],xlab="x",
       y=umap_resultado2[,2],ylab="y",
       z=umap_resultado2[,3],zlab="z",
       type="p",col=as.numeric(hoja))

####################Conclusiones##################################
##################################################################
#### Un problema de clasificación######
plot(proyeccion[,1],proyeccion[,2],col=as.numeric(hoja),pch=16)
reducido3d<-data.frame(hoja,proyeccion[,1:3])
plot3d(x=reducido3d$X1,xlab="x",
       y=reducido3d$X2,ylab="y",
       z=reducido3d$X3,zlab="z",
       type="p",col=as.numeric(hoja))

##Support vector machine 
###https://www.datacamp.com/community/tutorials/support-vector-machines-r
library(e1071)
svm_tune <- tune(svm, hoja~., data=reducido3d ,kernel ="linear", ##Probar radial
                 ranges = list(cost=c(0.001, 0.01,0.1, 1, 10, 100)))
svm_tune
svmfit = svm(hoja ~ ., data = reducido3d, kernel = "linear", 
             cost = 1, scale = FALSE)

plot(svmfit,data=reducido3d,X2~X3)

prueba<-reducido3d[,-1]

make.grid = function(x, n = 20) {
  grange = apply(x, 2, range)
  x1 = seq(from = grange[1,1], to = grange[2,1], length = n)
  x2 = seq(from = grange[1,2], to = grange[2,2], length = n)
  x3 = seq(from = grange[1,3], to = grange[2,3], length = n)
  expand.grid(X1 = x1, X2 = x2, X3= x3)
}

xgrid=make.grid(prueba)
ygrid = predict(svmfit, xgrid)
plot(xgrid[,1],xgrid[,2], col = c("red","blue","yellow")[as.numeric(ygrid)], pch = 20, cex = .2)
points(proyeccion[,1],proyeccion[,2],col=as.numeric(hoja),pch=16)

plot3d(x=xgrid[,1],xlab="x",
       y=xgrid[,2],ylab="y",
       z=xgrid[,3],zlab="z",
       type="p",col = c("red","blue","yellow")[as.numeric(ygrid)],add=T)

###### k- medias ###########################
#https://uc-r.github.io/kmeans_clustering
reducido3d2<-reducido3d
ggplot(reducido3d2,aes(x=X1,y=X2,color=hoja))+
  geom_point()+
  stat_ellipse(geom="polygon",aes(fill=hoja),
               alpha=0.25,linetype=2)

row.names(reducido3d2)<-paste(reducido3d2$hoja,1:length(reducido3d2$hoja))
kmedias<-kmeans(reducido3d2[,-c(1,3)],centers = 3,nstart = 25)##puede cambiarse a reducido3d[,-c(1,3)]
fviz_cluster(kmedias,data=reducido3d2[,-1])
cluster=as.factor(kmedias$cluster)
reducido3d2<-cbind(reducido3d2,cluster)

ggplot(reducido3d2,aes(x=X1,y=X2,color=hoja))+
  geom_point()+
  stat_ellipse(geom="polygon",aes(fill=hoja),
               alpha=0.25,linetype=2)

ggplot(reducido3d2,aes(x=X1,y=X2,group=cluster))+
  geom_point(aes(col=cluster))+
  stat_ellipse(aes(col=cluster))


###################################################
#Linear Discriminant Analysis
#https://www.displayr.com/linear-discriminant-analysis-in-r-an-introduction/
#http://www.sthda.com/english/articles/36-classification-methods-essentials/146-discriminant-analysis-essentials-in-r/
#https://www.geeksforgeeks.org/linear-discriminant-analysis-in-r-programming/
#https://www.rpubs.com/dvallslanaquera/lda
#https://rpubs.com/ranvirkumarsah/LDA
#Both Linear Discriminant Analysis (LDA) and Principal Component Analysis (PCA) 
#are linear transformation techniques that are commonly used for 
#dimensionality reduction. PCA can be described as an “unsupervised” 
#algorithm, since it “ignores” class labels and its goal is to find the 
#directions (the so-called principal components) that maximize the variance 
#in a dataset. In contrast to PCA, LDA is “supervised” and computes 
#the directions (“linear discriminants”) that will represent the axes 
#that that maximize the separation between multiple classes.
matriz_metabo<-read.csv("metabo_hojas.csv",header = T)
dim(matriz_metabo)
head(names(matriz_metabo))
##[1] "Hoja"         "Inoculo"      "Combinado"    "X50.83771282"
hoja=as.factor(matriz_metabo$Hoja)
matriz1=scale(as.matrix(matriz_metabo[,-c(1)]),scale=F)
dim(matriz1)
library(MASS)
lda(Hoja~.,data=matriz_metabo)
##El error es porque tenemos factores con varianza cero, por lo que son 
##constantes. NO poseen valor informativo
posicion<-round(apply(matriz_metabo,2,var),7)!=0 ##redondeado a 7 decimales
posicion[1]<-TRUE
length(posicion)
matriz_sin<-matriz_metabo[,posicion]

lda_result<-lda(Hoja~.,data=matriz_sin)
plot(lda_result)
lda_componentes<-predict(lda_result)
ldahist(lda_componentes$x[,1],g=matriz_sin$Hoja)

proyectados_lda=data.frame(Hoja=matriz_sin$Hoja,lda=lda_componentes$x)

ggplot(proyectados_lda,aes(x=lda.LD1,y=lda.LD2))+
  geom_point(aes(col=Hoja))+
  stat_ellipse(aes(col=Hoja))

##En realidad debería ser aplicado para bases más pequeñas.
#NO es muy conveniente para bases donde el numero de variables
#supera al de observaciones, como en este caso
###################################################
######Conclusiones#################################
###################################################
library(rattle)
rattle()
#bosques aleatorios
#redes neuronales

