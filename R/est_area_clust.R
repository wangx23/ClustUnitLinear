#' est_area_clust is used to estimate domain means based on unit level models

#'@param dom vector for domain.
#'@param y response.
#'@param x if needs intercept, need to add before fitting the model.
#'@param wts sampling weights.
#'@param index_d the column index for different coefficients.
#'@param Xbar Xbar is a m * (p+1) matrix, p is the number of column of x; the first column is the dom index.
#'@param lamvec a vector of lambda values
#'@import lme4
#'@export

est_area_clust <- function(dom, y, x, wts = NULL, index_d, Xbar,
                           lamvec, algorithm = 2,
                           init0 = NULL, vv0 = NULL)
{

  m <- length(unique(dom))
  n0 <- length(dom)
  len_d <- length(index_d)
  ncx <- ncol(x)
  ns <- as.vector(table(dom))

  ybar <- aggregate(y, by = list(dom = dom), mean)[,2]
  xbar <- as.matrix(rowsum(x,group = dom)/ns)

  ### without wts
  if(is.null(wts))
  {
    if(algorithm==1)
    {
      res <- Clust_nested_error_bic(dom, y, x, index_d, lamvec)
    }
    if(algorithm==2)
    {
      res <- Clust_nested_error_bic2(dom, y, x, index_d, lamvec,
                                     init0 = init0, vv0 =vv0)
    }
    cluster_est <- res$fitted$cluster
    fitted <- res$fitted
    eta_est <- fitted$eta
    betam_est <- fitted$beta
    rand_est <- fitted$rand
  }

  if(!is.null(wts))
  {
    if(algorithm==1)
    {
      res <- Clust_nested_info_bic(dom, y, x, wts, index_d, lamvec)
    }
    if(algorithm==2)
    {
      res <- Clust_nested_info_bic2(dom, y, x, wts, index_d, lamvec,
                                     init0 = init0, vv0 =vv0)
    }
    cluster_est <- res$fitted$cluster
    fitted <- res$fitted
    eta_est <- fitted$eta
    betam_est <- fitted$beta
    rand_est <- fitted$rand
  }

  Xbar0 <- Xbar[,-1]
  Xbar_d <- Xbar0[,index_d, drop = FALSE]
  xbar_d <- xbar[,index_d, drop = FALSE]

  if(len_d < ncx)
  {
    # refit if can
    refitted <- NULL
    tryCatch({
      refitted <- fit_clust_nested(dom, y = y, x = x, wts = wts, index_d = index_d,
                                   cluster = cluster_est)
      rand_est <- refitted$rand
      coef_est <- refitted$est[,1]
      eta_est <- coef_est[1:(ncx-len_d)]
      betam_est <- matrix(coef_est[-(1:(ncx-len_d))], ncol = len_d, byrow= TRUE)[cluster_est,,drop=FALSE]
    }, error = function(e) {
    })


    Zbar_d <- Xbar0[,-index_d, drop = FALSE]
    Ybar_est <- Zbar_d %*% eta_est + rowSums(Xbar_d * betam_est) + rand_est
  }

  if(len_d == ncx)
  {
    refitted <- NULL
    tryCatch({
      refitted <- fit_clust_nested(dom, y = y, x = x, wts = wts, index_d = index_d,
                                   cluster = cluster_est)
      rand_est <- refitted$rand
      coef_est <- refitted$est[,1]
      betam_est <- matrix(coef_est, ncol = len_d, byrow= TRUE)[cluster_est,,drop=FALSE]
    }, error = function(e) {
    })

    Ybar_est <- rowSums(Xbar_d * betam_est) + rand_est
  }


  out <- list(est = Ybar_est, fit = res, refit = refitted)
  return(out)

}
