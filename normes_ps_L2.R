ps <- function(f,g){
  
  # Calculate the scalar products of two functions f et g (in discretized or function form)
  if(is.function(f)){
    f <- f(seq(0,1,length.out=100))
  }
  if(is.function(g)){
    g <- g(seq(0,1,length.out=100))
  }
  if (is.vector(f)){
    nbpt<-length(f)
    if (is.vector(g)){
      rep<-as.numeric((1/(2*(nbpt-1)))*(f[1]*g[1]+f[nbpt]*g[nbpt]+2*sum(f[2:(nbpt-1)]*g[2:(nbpt-1)])))
    } else { 
      nbg<-NROW(g)
      rep<-as.vector((1/(2*(nbpt-1)))*(2*tcrossprod(f,g)-f[1]*(g[,1])-f[nbpt]*(g[,nbpt])))}
  }else{  if (is.vector(g)){
    nbpt<-length(g)
    rep<-as.vector((1/(2*(nbpt-1)))*(2*tcrossprod(g,f)-g[1]*(f[,1])-g[nbpt]*(f[,nbpt])))
    
  } else {
    nbpt<-NCOL(g)
    rep<-(1/(2*(nbpt-1)))*(2*tcrossprod(f,g)-tcrossprod(f[,1],g[,1])-tcrossprod(f[,nbpt],g[,nbpt]))}
    
  }
  
  return(rep)
}





norm <- function(X,...){
  # Input : a matrix X=(X_i(t_j))_{1<= i <= n, 1 <= j<= p})) or a function
  # Output : vector containing the approximations of the norms of X_i
  if(is.vector(X)) {
    N <- sqrt(ps(X,X))
  } else {
    N <- sqrt(diag(ps(X,X)))
  }
  N
}

