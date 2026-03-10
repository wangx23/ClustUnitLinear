#' bootstrap_infor is used to estimate MSE.

#'@param dom vector for domain.
#'@param y response.
#'@param x if needs intercept, need to add before fitting the model.
#'@param wts sampling weights.
#'@param index_d the column index for different coefficients.
#'@param Xbar Xbar is a m * (p+1) matrix, p is the number of column of x; the first column is the dom index.
#'@param fit object from est_area_clut.
#'@param useparallel if using parallel.
#'@param ncores number of cores.
#'@param B bootstrap replicates.
#'@param conf confidence level.
#'@import lme4
#'@import doParallel
#'@import parallel
#'@import foreach
#'@export

bootstrap_info <- function(dom, y, x, wts = NULL, index_d, Xbar,
                           fit, useparallel = FALSE, ncores = 4,
                            B = 200, conf = 0.95,
                           seed = 123456)
{
  m <- length(unique(dom))
  n0 <- length(dom)
  len_d <- length(index_d)
  ncx <- ncol(x)
  ns <- as.vector(table(dom))


  cluster_est <- fit$fit$fitted$cluster ### estimated clusters
  est_area <- fit$est  ### estimated

  betam <- fit$fit$fitted$betac[cluster_est,drop=FALSE,-1] ### estimated betam
  eta_est <- fit$fit$fitted$eta  #### estimated eta
  tau1_est <- fit$fit$fitted$tau1
  tau0_est <- fit$fit$fitted$tau0


  Xbar0 <- Xbar[,-1]
  Xbar_d <- Xbar0[,index_d, drop = FALSE]
  x <- as.matrix(x)
  x_d <- x[,index_d, drop = FALSE] ## subgroup part

  if(len_d == ncx)
  {
    refitd <- fit$refit
    if(!is.null(refitd))
    {
      betam <- matrix(refitd$est[,1], ncol = len_d, byrow= TRUE)[cluster_est,,drop=FALSE]
      tauvec_est <-  as.data.frame(VarCorr(refitd$fit))$vcov
      tau1_est <- tauvec_est[1]
      tau0_est <- tauvec_est[2]
    }

    muhat <- rowSums(Xbar * betam)
    muhat_obs <- rowSums(x * betam[dom,, drop = FALSE])
  }

  if(len_d < ncx)
  {
    refitd <- fit$refit
    if(!is.null(refitd))
    {
      estvec <- refitd$est[,1]
      tauvec_est <-  as.data.frame(VarCorr(refitd$fit))$vcov

      eta_est <- estvec[1:(ncx-len_d)]
      betam <- matrix(estvec[-(1:(ncx-len_d))],
                      ncol = len_d, byrow= TRUE)[cluster_est,,drop=FALSE]
      tau1_est <- tauvec_est[1]
      tau0_est <- tauvec_est[2]
    }

    Zbar_d <- Xbar0[,-index_d, drop = FALSE]
    z <-  x[,-(index_d), drop = FALSE] ## common part
    muhat <- Zbar_d %*% eta_est + rowSums(Xbar_d * betam)
    muhat_obs <- z %*% eta_est + rowSums(x_d * betam[dom,,drop = FALSE]) ## for simulation

  }

  ##### bootstrap start here ####
  cutoff <- qnorm(1-(1-conf)/2)
  lamvec <- fit$fit$bic$lambda

  if(!useparallel)
  {
    mse_estB1 <- matrix(NA, B, m)
    for(b in 1:B)
    {
      set.seed(seed+b)

      tryCatch({
        vib <- rnorm(m)*sqrt(tau1_est) ### random effects
        thetab <- muhat + vib ## target
        yb <- muhat_obs + vib[dom] + rnorm(n0) * sqrt(tau0_est) ### sampled


        ### refit the model

        resb <- est_area_clust(dom, y = yb, x = x,  wts = wts,
                               index_d = index_d, Xbar = Xbar,
                               lamvec = lamvec, algorithm = 2)
        estb <- resb$est

        mse_estB <- (estb - thetab)^2

      }, error = function(e)
      {
        mse_estB[b,] <- NA
      })

      print(b)
    }

    Bt <- sum(rowSums(is.na(mse_estB)) ==0)
    mse_estB1 <- mse_estB[rowSums(is.na(mse_estB))==0,]
  }

  if(useparallel)
  {


    subfun <- function(bb)
    {
      set.seed(seed+bb)
      vib <- rnorm(m)*sqrt(tau1_est) ### random effects
      thetab <- muhat + vib ## target
      yb <- muhat_obs + vib[dom] + rnorm(n0) * sqrt(tau0_est) ### sampled


      ### refit the model

      resb <- est_area_clust(dom, y = yb, x = x,  wts = wts,
                             index_d = index_d, Xbar = Xbar,
                             lamvec = lamvec, algorithm = 2)
      estb <- resb$est

      mse_bb <- (estb - thetab)^2
      return(mse_bb)
    }

    packages <- c("dplyr","magrittr","lme4")

    cl = makeCluster(ncores)
    registerDoParallel(cl)
    out_list <- foreach(mm = 1:B,
                        .packages = packages,
                        .errorhandling = "remove") %dopar% subfun(mm)
    stopCluster(cl)

    Bt <- length(out_list)
    mse_estB <- matrix(0, Bt, m)
    for(b in 1:Bt)
    {
      outb <- out_list[[b]]
      mse_estB[b,] <- outb
    }
  }


  mse_est <- colMeans(mse_estB)
  ci_nomial <- cbind(est_area[,1] - cutoff*sqrt(mse_est),
                      est_area[,1] + cutoff*sqrt(mse_est))
  colnames(ci_nomial) <- c("lower","upper")


  outlist <- list(mse = mse_est,
                  ci_nomial = ci_nomial,
                  Bt = Bt)
  return(outlist)
}

