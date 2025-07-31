#' initial values functions
#'@param dom vector for domain
#'@param y response
#'@param x if needs intercept, need to add before fitting the model
#'@param index_d the column index for different coef
#'@export

initial_nested <- function(dom, y, x, index_d, wts = NULL,
                           lam0 = 0.001)
{

  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ncx <- ncol(x)
  len_d <- length(index_d)


  if(is.null(wts))
  {
    wts <- rep(1,n0)
  }else{
    wts <- nr * wts/sum(wts)
  }

  x <- as.matrix(x)*sqrt(wts)
  y <- y*sqrt(wts)

  x_d <- x[,index_d, drop = FALSE] ## subgroup part

  # for X expand
  Xm <- xexpand(dom, uindex, n0, x_d, len_d = len_d, cluster = 1:nr)

  if(len_d < ncx)
  {
    z <-  x[,-(index_d), drop = FALSE] ## common part

    ztz <- solve(t(z) %*% z)%*%t(z)
    Qz <- diag(1, n0, n0) - z%*% ztz
    Xty <- t(Xm)%*% Qz %*% y
    Xinv <- inverseR2(dom, z, x_d, lam0)

    beta00_vec <- Xinv %*% Xty
    beta00 <-  matrix(beta00_vec,ncol= len_d, byrow = TRUE)
    eta00 <- ztz %*% (y - Xm %*% beta00_vec)

    cluster0 <- kmeans(beta00, centers = sqrt(nr), nstart = 100)$cluster

    ### redo above again
    ng <- length(unique(cluster0))
    Ip <- diag(1,len_d,len_d)
    D <- matrix(0,ng*(ng-1)/2,ng)
    for(j in 1:(ng-1))
    {
      indexj <- (2*ng-j)*(j-1)/2
      indexvj <- indexj + (1:(ng-j))
      D[indexvj,j] <- 1
      D[cbind(indexvj,(j+1):ng)] <- -1
    }
    Am <- D %x% Ip
    AtA <- t(D)%*%D %x% Ip


    Xm <- xexpand(dom, uindex, n0, x_d, len_d = len_d, cluster = cluster0)
    Xty <- t(Xm)%*% Qz %*% y
    Xinv <- solve(t(Xm) %*% Qz %*% Xm + lam0 *AtA)
    beta00_vec <- Xinv %*% Xty
    beta00 <-  matrix(beta00_vec,ncol= len_d, byrow = TRUE)
    eta00 <- ztz %*% (y - Xm %*% beta00_vec)

    beta00 <- beta00[cluster0, ,drop = FALSE]

    out <- list(eta = eta00, betam = beta00)
  }

  if(len_d == ncx)
  {
    Xty <- t(Xm) %*% y
    Xinv <- inverseR(dom, x_d, lam0)

    cluster0 <- kmeans(beta00, centers = sqrt(nr), nstart = 100)$cluster

    ### redo above again
    ng <- length(unique(cluster0))
    Ip <- diag(1,len_d,len_d)
    D <- matrix(0,ng*(ng-1)/2,ng)
    for(j in 1:(ng-1))
    {
      indexj <- (2*ng-j)*(j-1)/2
      indexvj <- indexj + (1:(ng-j))
      D[indexvj,j] <- 1
      D[cbind(indexvj,(j+1):ng)] <- -1
    }
    Am <- D %x% Ip
    AtA <- t(D)%*%D %x% Ip


    Xm <- xexpand(dom, uindex, n0,  x_d, len_d = len_d, cluster = cluster0)
    Xty <- t(Xm) %*% y
    Xinv <- solve(t(Xm) %*% Xm + lam0 *AtA)

    beta00_vec <- Xinv %*% Xty
    beta00 <-  matrix(beta00_vec,ncol= len_d, byrow = TRUE)

    beta00 <- beta00[cluster0, ,drop = FALSE]
    out <- list(eta = NULL, betam = beta00)
  }

  return(out)
}

xexpand <- function(dom, uindex, n0, x_d, len_d, cluster)
{
  ng <- length(unique(cluster))
  Ux <- matrix(0, n0, ng*len_d)
  for(j in 1:ng)
  {
    uindexj <- uindex[cluster==j]
    indexj <- dom %in% uindexj
    Ux[indexj,(len_d*(j-1) + 1) : (len_d*j)] <- x_d[indexj,]
  }
  return(Ux)
}
