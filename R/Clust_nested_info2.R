#' Clust_nested linear model  with sampling weights, warm start

#'@param wts is the sampling weights
#'@param dom vector for domain
#'@param y response
#'@param x if needs intercept, need to add before fitting the model
#'@param index_d the column index for different coef
#'@export

Clust_nested_info2 <- function(dom, y, x, wts, index_d, lam,
                              init0 = NULL, vv0 = NULL,
                              nu = 1, gam = 3, lam0 = 0.001,
                              maxiter= 500,
                              tol = 1e-3)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  nbar <- mean(as.numeric(table(dom)))
  ncx <- ncol(x)

  wts_tilde <- wts/sum(wts)* nr

  x <- as.matrix(x)

  ## sparse function
  D <- matrix(0,nr*(nr-1)/2,nr)
  for(j in 1:(nr-1))
  {
    indexj <- (2*nr-j)*(j-1)/2
    indexvj <- indexj + (1:(nr-j))
    D[indexvj,j] <- 1
    D[cbind(indexvj,(j+1):nr)] <- -1
  }

  Hmat <- as.matrix(model.matrix(~0+as.factor(dom))) ## for random effect
  colnames(Hmat) <- NULL

  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)

  x_d <- x[,index_d, drop = FALSE] ## subgroup part

  # for X expand
  Xm <- matrix(0, n0, nr*len_d)
  for(i in 1:nr)
  {
    Xm[dom == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[dom == uindex[i],]
  }

  if(is.null(init0))
  {
    init0 <- initial_nested(dom, y, x, index_d, lam0 = lam0)
  }
  if(is.null(vv0))
  {
    vv0 <- rep(0, nr)
  }

  if(len_d < ncx)
  {
    z <-  x[,-(index_d), drop = FALSE] ## common part

    #### initial values
    eta00 <- init0$eta
    beta00 <- init0$betam
    beta00_vec <- c(t(beta00))

    deltam <- D %*% beta00
    um <-  matrix(0, nr*(nr-1)/2 , len_d)
    vv_cur <- vv0 ### random effects
    tau00 <- mean((y - z%*% eta00- Xm %*% beta00_vec - vv_cur[dom])^2)

    eta_cur <- eta00
    beta_cur <- c(t(beta00))
    tau0_cur <- tau00/2
    tau1_cur <- tau00/2


    for(j in 1:maxiter)
    {
      wm <- 1/tau0_cur* wts_tilde

      ## update eta and beta
      eta_new <- solve(t(wm * z) %*% z) %*% t(wm * z) %*% (y - Xm %*% beta_cur - vv_cur[dom])

      Xty <- t(wm * Xm) %*% (y - z %*%eta_new -  vv_cur[dom])
      Xinv <- inverseR(dom, x_d*sqrt(wm), nu = nu)
      beta_new <-  Xinv %*% (Xty + nu*c(t(deltam -  um/nu) %*% D))
      betam_new <- matrix(beta_new, nr, len_d, byrow = TRUE)

      ## update vv ##
      Vinv <- 1/(as.numeric(by(as.numeric(wm*nbar), dom, sum)) + 1/tau1_cur)
      vv_new <- Vinv * (t(wm *nbar* Hmat) %*% (y - z %*% eta_new - Xm %*% beta_new))


      #### update tau0 and tau1 ##
      tau0_new <- sum(wts_tilde*(y - z %*% eta_new - Xm %*% beta_new - vv_new[dom])^2)/nr
      tau1_new <- mean(Vinv + vv_new^2)

      ###
      psim <- D %*% betam_new + um/nu
      deltam_new <- t(matrix(sapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam)),
                             ncol=nrow(psim)))
      um <- um + nu * (D %*% betam_new - deltam_new)

      rm <- sqrt(sum((D %*% betam_new  - deltam_new)^2))

      deltam <- deltam_new
      vv_cur <- vv_new
      beta_cur <- beta_new
      eta_cur <- eta_new
      tau0_cur <- tau0_new
      tau1_cur <- tau1_new

      if(rm <= tol*sqrt(nr*len_d))
      {break}
    }

    cluster_est <- getgroup(t(deltam), n = nr)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    l_pql <- n0*log(tau0_new) + sum(vv_new^2)/tau1_new

    out <- list(eta = eta_new, beta = betam_new,
                tau0 = tau0_new, tau1 = tau1_new,
                rand = vv_new,
                betac = beta_c, cluster = cluster_est, loss = l_pql,
                deltam = deltam, niters = j, init0 = init0)

  }

  if(len_d == ncx)
  {
    #### initial values

    beta00 <- init0$betam
    beta00_vec <- c(t(beta00))


    deltam <- D %*% beta00
    um <-  matrix(0, nr*(nr-1)/2 , len_d)
    vv_cur <- vv0 ### random effects
    tau00 <- mean((y - Xm %*% beta00_vec - vv_cur[dom])^2)

    beta_cur <- c(t(beta00))
    tau0_cur <- tau00/2
    tau1_cur <- tau00/2


    for(j in 1:maxiter)
    {
      wm <- 1/tau0_cur* wts_tilde

      ## update  beta

      Xty <- t(wm * Xm) %*% (y -  vv_cur[dom])
      Xinv <- inverseR(dom, x_d*sqrt(wm), nu = nu)
      beta_new <-  Xinv %*% (Xty + nu*c(t(deltam -  um/nu) %*% D))
      betam_new <- matrix(beta_new, nr, len_d, byrow = TRUE)

      ## update vv ##
      Vinv <- 1/(as.numeric(by(as.numeric(wm*nbar), dom, sum)) + 1/tau1_cur)
      vv_new <- Vinv * (t(wm *nbar* Hmat) %*% (y - Xm %*% beta_new))

      #### update tau0 and tau1 ##
      tau0_new <- sum(wts_tilde*(y - Xm %*% beta_new - vv_new[dom])^2)/nr
      tau1_new <- mean(Vinv + vv_new^2)

      ###
      psim <- D %*% betam_new + um/nu
      deltam_new <- t(matrix(sapply(1:nrow(psim),function(xx) mcp(psim[xx,],lam,nu,gam)),
                             ncol=nrow(psim)))
      um <- um + nu * (D %*% betam_new - deltam_new)

      rm <- sqrt(sum((D %*% betam_new  - deltam_new)^2))

      deltam <- deltam_new
      vv_cur <- vv_new
      beta_cur <- beta_new
      tau0_cur <- tau0_new
      tau1_cur <- tau1_new

      if(rm <= tol*sqrt(nr*len_d))
      {break}
    }

    cluster_est <- getgroup(t(deltam), n = nr)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    l_pql <- n0*log(tau0_new) + sum(vv_new^2)/tau1_new

    out <- list(eta = NULL, beta = betam_new,
                tau0 = tau0_new, tau1 = tau1_new,
                rand = vv_new,
                betac = beta_c, cluster = cluster_est, loss = l_pql,
                deltam = deltam, niters = j, init0 = init0)
  }


  return(out)


}

#'@export
Clust_nested_info_bic2 <- function(dom, y, x, wts, index_d, lamvec,
                                  init0 = NULL, vv0 = NULL,
                                 nu = 1, gam = 3, lam0 = 0.001,
                                 maxiter= 500,
                                 tol = 1e-3)
{
  nlam <- length(lamvec)
  K_est <- rep(0, nlam)
  loss_value <- rep(0, nlam)
  out <- list()
  m <- length(unique(dom))
  n0 <- length(dom)
  len_d <- length(index_d)

  if(is.null(init0))
  {
    init0 <- initial_nested(dom, y, x, index_d, lam0 = lam0)
  }
  if(is.null(vv0))
  {
    vv0 <- rep(0, m)
  }


  for(jj in 1:nlam)
  {
    resj <- Clust_nested_info2(dom, y, x, wts, index_d, lam = lamvec[jj],
                            init0 = init0, vv0 = vv0,
                             nu = nu, gam = gam, lam0 = lam0,
                             maxiter = maxiter, tol = tol)
    init0$eta <- resj$eta
    init0$betam <- resj$beta
    vv0 <- resj$rand
    out[[jj]] <- resj
    K_est[jj] <- length(unique(resj$cluster))
    loss_value[jj] <- resj$loss

  }
  bic_df <- data.frame(lambda= lamvec, Khat = K_est,
                       loss = loss_value)

  bic_df <- bic_df %>%
    mutate(mbic = loss + log(m)*log(n0)*(Khat*len_d))

  ### fitted

  res_fitted <- out[[which.min(bic_df$mbic)]]

  outlist <- list(out = out, bic = bic_df, fitted = res_fitted)
  return(outlist)
}
