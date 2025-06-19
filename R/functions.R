# some functions used in the algorithm
#### MCP 
sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))
  thval <- 1 - th/xn
  
  if(xn ==0)
  {
    value <- xn
  }else{
    value <- thval*((thval) >0)*x
  }
}
mcp <- function(x,lam,nu,gam)
{
  temp <- gam*lam
  xn <- sqrt(sum(x^2))
  if(xn <= temp)
  {
    z <- sfun(x,lam/nu) / (1 - 1/(gam*nu))
  }else{
    z <- x
  }
  z <- as.matrix(z)
  return(z)
}



#### get group ####
getgroup = function(deltam, n, tol = 1e-2)
{
  p = nrow(deltam)
  b2value =sqrt(colMeans(deltam^2))
  b2value[b2value <= tol] = 0
  
  d2 = matrix(0, n, n)
  for(j in 1:(n-1))
  {
    indexj1 = (2*n -j)*(j-1)/2 + 1
    indexj2 = indexj1 + n - j - 1
    d2[(n - indexj2 + indexj1):n,j] = b2value[indexj1:indexj2]
  }
  d2 = t(d2) + d2
  
  
  ngj = 1:n
  groupest = rep(0,n)
  j = 0
  
  while (length(ngj) >0)
  {
    j = j + 1
    gj = (1:n)[d2[ngj[1],] ==0]
    indexj = ngj %in% gj
    gj = ngj[indexj]
    ngj = ngj[!indexj]
    groupest[gj] = j * rep(1,length(gj))
  }
  
  
  return(groupest)
  
}



#### quick inverse ###
inverseR <- function(indexy, x,nu = 1)
{
  uindex <- unique(indexy)
  n <- length(uindex)
  p <- ncol(x)
  Ip <- diag(1,p,p)
  
  DB <- matrix(0,p,p)
  
  A0B <- matrix(0, n*p, p)
  matinv <- matrix(0, n*p, n*p)
  
  for(i in 1:n)
  {
    xmati <- x[indexy==uindex[i],,drop=FALSE]
    matj <- solve(t(xmati)%*%xmati + nu*n*Ip)
    DB <- DB  + matj
    indexi <- ((i-1)*p+1):(i*p)
    A0B[indexi,] <- matj
    matinv[indexi,indexi] <- matj
  }
  
  temp1 <- solve(1/nu*Ip - DB)
  
  A1 <- A0B %*% temp1 %*% t(A0B)
  
  
  matinv <- matinv + A1 
  return(matinv)
}


inverseR2 <- function(indexy,z, x,nu = 1)
{
  uindex <- unique(indexy)
  n <- length(uindex)
  p <- ncol(x)
  q <- ncol(z)
  Ip <- diag(1,p,p)
  
  # D^A(-1)B ###
  DB <- matrix(0,p,p)
  zA0z <- matrix(0, q, q)
  A0z <- matrix(0,p,q)
  
  A0B <- matrix(0, n*p, p)
  A0Z <- matrix(0, n*p, q)
  matinv <- matrix(0, n*p, n*p)
  
  for(i in 1:n)
  {
    xmati <- x[indexy==uindex[i],,drop=FALSE]
    zmati <- z[indexy==uindex[i],,drop=FALSE]
    matj <- solve(t(xmati)%*%(xmati) + nu*n*Ip)
    DB <- DB  + matj
    zA0z <- zA0z + t(zmati)%*%(zmati) -  t(zmati)%*%(xmati) %*% matj%*% t(xmati)%*%(zmati)
    A0z <- A0z + matj%*% t(xmati)%*%(zmati)
    indexi <- ((i-1)*p+1):(i*p)
    A0B[indexi,] <- matj
    A0Z[indexi,] <- matj%*%t(xmati)%*%(zmati)
    matinv[indexi,indexi] <- matj
    
  }
  
  temp1 <- solve(1/nu*Ip - DB - A0z %*% solve(zA0z) %*%t(A0z))
  
  A1 <- A0B %*% temp1 %*% t(A0B)
  A2 <- A0Z %*% solve(zA0z)%*% t(A0z)
  A3 <- A2 %*%temp1
  A22 <- A3 %*% t(A2)
  A12 <- A3 %*% t(A0B)
  
  matinv <- matinv+ A0Z %*% solve(zA0z) %*% t(A0Z) + A1 + A22 + A12 + t(A12)
  return(matinv)
}



