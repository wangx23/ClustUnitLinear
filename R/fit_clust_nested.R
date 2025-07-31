#' fit the LMM model for a given cluster structure

#'@param dom vector for domain
#'@param y response
#'@param x if needs intercept, need to add before fitting the model
#'@param index_d the column index for different coef
#'@param wts is the sampling weights
#'@import lme4
#'@export

fit_clust_nested <- function(dom, y, x, wts = NULL, index_d,
                             cluster)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  len_d <- length(index_d)

  x <- as.matrix(x)

  x_d <- x[,index_d, drop = FALSE]

  ng <- length(unique(cluster))
  Ux <- matrix(0, n0, ng*len_d)
  for(j in 1:ng)
  {
    uindexj <- uindex[cluster==j]
    indexj <- dom %in% uindexj
    Ux[indexj,(len_d*(j-1) + 1) : (len_d*j)] <- x_d[indexj,]
  }

  names_d <- paste(rep(paste("x",1:len_d,sep=""),ng), rep(paste("c",1:ng,sep=""), each = len_d),sep="_")

  if(len_d == ncol(x))
  {
    Ux <- Ux
    colnames(Ux) <- names_d
  }else{
    Ux <- cbind(x[,-index_d, drop = FALSE], Ux)
    colnames(Ux) <- c(paste("z",1:(ncol(x)-len_d),sep=""), names_d)
  }


  df_f <- as.data.frame(cbind(y, Ux))

  for1 <- as.formula(paste("y", "~0+", paste(c(colnames(Ux),"(1|dom)"), collapse = "+")))

  if(!is.null(wts))
  {
    res1 <- lmer(for1, data = df_f, weights = n0*wts/sum(wts))
  }else{
    res1 <- lmer(for1, data = df_f)
  }


  randv <- ranef(res1)$dom[,1]
  est_mat <- summary(res1)$coefficients

  out <- list(est = est_mat, rand = randv, fit = res1)

  return(out)
}
