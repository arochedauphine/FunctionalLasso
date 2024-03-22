tikho <- function(X,Y,prods,rho,beta_init=NULL,niter=100,calpha=1,expalpha=1,plots=FALSE){
  
  
  alpha = calpha*(1:niter)^(-expalpha)
  
  n = length(Y)
  
  res = vector("list",niter)
  
  if (is.null(beta_init)){
    betacourant = init(X,n)
    betacourantAVG = init(X,n)
  } else {
    betacourant = beta_init
    betacourantAVG = beta_init
  }
  
  p = length(betacourant)
  
  x = seq(0,1,length.out=length(beta_init[[1]]))
  
  predcourant = ps_tot(betacourant,X,prods)
  rescourant = as.vector(Y-predcourant)
  erreur_moy = rep(mean(rescourant^2),niter)
  erreur_moyAVG = rep(mean(rescourant^2),niter)
  #print(c("erreur moyenne =",mean(rescourant^2)))
  
  res[[1]] = list(beta=beta_init,preds=predcourant)
  
  
  for (j in 2:niter){
    for (k in 1:p){
      gradk = as.vector(-2*crossprod(rescourant,X[[k]])/n+2*rho*betacourant[[k]])
      betacourant[[k]] = betacourant[[k]] - alpha[k]*gradk
      betacourantAVG[[k]] = ((niter-1)/niter)*betacourantAVG[[k]] + betacourant[[k]]/niter
    }
    predcourant = ps_tot(betacourant,X,prods)
    rescourant = as.vector(Y-predcourant)
    predcourantAVG = ps_tot(betacourantAVG,X,prods)
    rescourantAVG = as.vector(Y-predcourantAVG)
    res[[j]] = list(beta=betacourant,preds=predcourant,betaAVG=betacourantAVG)
    erreur_moy[j] = mean(rescourant^2)
    erreur_moyAVG[j] = mean(rescourantAVG^2)
  }
  if (plots){
    plot(1:niter,erreur_moy,type='l')
  }
  res=res[[niter]]
}

tikho.CV <- function(X,Y,prods,rho_grid,beta_init=NULL,niter=100,V=5,calpha=1,expalpha=1){
  
  n <- length(Y)
  size <- n/V
  p <- length(X)
  errormean <- rep(NA,length(rho_grid))
  for (j in 1:length(rho_grid)){
    rho = rho_grid[j]
    print(c("rho=",rho))
    predict <- rep(NA,n)
    for (f in 1:V){
      validationIndex <- seq((f - 1) * size + 1, f * size) 
      Xvalid = vector("list",p)
      for (k in 1:p){
        if (is.vector(X[[k]])) {
          Xvalid[[k]] = X[[k]][-validationIndex]
        } else {
          Xvalid[[k]] = X[[k]][-validationIndex,]
        }
      }
      for (i in validationIndex){
        Xpred = vector("list",p)
        for (k in 1:p){
          if (is.vector(X[[k]])) {
            Xpred[[k]] = X[[k]][i]
          } else {
            Xpred[[k]] = X[[k]][i,]
          }
        }
        res = tikho(Xvalid,Y[-validationIndex],prods,rho,beta_init,niter,calpha,expalpha)
        predict[i] <- ps_tot(res$beta,Xpred,prods)
      }
      errormean[j] <- mean((predict-Y)^2,na.rm=TRUE)^(1/2)
    }
    hatrho <- rho_grid[which.min(errormean)] 
  }
  list(res=tikho(X,Y,prods,hatrho,beta_init,niter,calpha,expalpha),rho=hatrho)
}
