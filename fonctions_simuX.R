simuX=function(n,p=100,J=50,lambda=(1:J)^(-2),loi='norm',...)
{
  if (loi=='norm')  ksi=rnorm(n*J) else ksi=runif(n*J,-sqrt(3),sqrt(3))
  ksi=matrix(ksi,n,J)*matrix(rep(sqrt(lambda),n),n,J,byrow=TRUE)
  x=matrix(rep(seq(0,1,length.out=p),J),J,p,byrow=TRUE)
  Jm=matrix(rep(seq(0.5,J-0.5,by=1),p),J,p) 
  base=sqrt(2)*sin(pi*Jm*x)
  pertinit <- matrix(rep(rnorm(n),p),n,p,byrow=FALSE)
  ksi%*%base + pertinit
}

simuX_FV2000 = function(n,p=100,...){
  # prend en entree la taille n de l'echantillon, le nombre p de points de discretisation.
  # retourne un echantillon X[i,j]=X_i(t_j) de courbes aleatoires simulees selon la loi 
  #  X(t) = a + b*t + c*exp(t) + sin(d*t), a~U(0,100), b~U(-30,30), c~U(0,10), d~U(1,3)
  # (Ferraty et Vieu 2000)(t_j)_{1<=j<=p} 
  
  a <- runif(n,0,100)
  b <- runif(n,-30,30)
  c <- runif(n,0,10)
  d <- runif(n,1,3)
  
  x <- seq(0,1,length.out=p)
  Mx <- matrix(rep(x,n),n,p,byrow=TRUE)
  matrix(rep(a,p),n,p) + matrix(rep(b,p),n,p)*Mx + matrix(rep(c,p),n,p)*exp(Mx) + sin(matrix(rep(d,p),n,p)*Mx) 
}
### Classe X_AFKV2008
#Constructeur
simuX_AFKV2008 <- function(coeffs = runif(3)){
  obj <- coeffs
  class(obj) <- "X_AFKV2008"
  obj
}
as.function.X_AFKV2008 <- function(X){
  function (x){X[1]*cos(2*pi*x)+X[2]*sin(4*pi*x)+2*X[3]*(x-0.25)*(x-0.5)}
}

plot.X_AFKV2008 <- function(X,x = seq(0,1,length.out=100),col=1){
  plot(x,as.function(X)(x),type='l',col=col)
}

points.X_AFKV2008 <- function(X,x = seq(0,1,length.out=100),col=1){
  points(x,as.function(X)(x),type='l',col=col)
}

simuXpol=function(n,p,J,varj){
  ksi<-rnorm(n*J)
  ksi=matrix(ksi,n,J)*matrix(rep(sqrt(varj),n),n,J,byrow=TRUE)
  x=matrix(rep(seq(0,1,length.out=p),J),J,p,byrow=TRUE)
  Jm=matrix(rep(0:(J-1),p),J,p)
  base=x^Jm
  ksi%*%base
}

MvB=function(n,p)
{
 X=rnorm(n*p)/sqrt(p)
 X=matrix(X,ncol=p)
 MvB=matrix(0,ncol=p,nrow=n)
 for (i in 1:n)
 {
  MvB[i,]=cumsum(X[i,])
 }
 MvB
}

SimulGauss=function(n,p,theta=0.2,sigma=1,alpha=1.5)
{
### Simulation d'un processus de covariance Gaussienne
### de parametre theta
### n trajectoires sur [0,1]
### p points de discretisation equidistants
X=rnorm(n*p)
X=matrix(X,nrow=n,ncol=p)
Gamma=matrix(0,ncol=p,nrow=p)
temps=seq(0,1,length=p)
Temps.mat=Gamma
for (j in 1:p)
{
 Temps.mat[,j]=abs(temps-temps[j])^alpha
}
Gamma=sigma^2*exp(-Temps.mat/theta)+(1e-08)*diag(p)
Gamma.sqrt=chol(Gamma)
mat=X%*%Gamma.sqrt
list(S=mat,Gamma=Gamma)
}

simuX_per=function(n,p,J,lambda=(1:J)^(-2),loi='norm')
{
  I<-(J-1)/2
  if (loi=='norm')  ksi=rnorm(n*J) else ksi=runif(n*J,-sqrt(3),sqrt(3))
  ksi=matrix(ksi,n,J)*matrix(rep(sqrt(lambda),n),n,J,byrow=TRUE)
  x=seq(0,1,length.out=p)
  base=matrix(rep(0,p*J),J,p)
  base[1,]=1
  base[seq(2,2*I,by=2),]=sqrt(2)*cos(2*pi*(1:I)%*%t(x))
  base[seq(3,J,by=2),]=sqrt(2)*sin(2*pi*(1:I)%*%t(x))
  list(coeffs=ksi , X=ksi%*%base)
}
