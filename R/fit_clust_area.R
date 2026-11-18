######## fit the model based on a give cluster structure and given vardir ####
## cluster should have two columns

#' Fit a Fay-Herriot area-level model for a given cluster structure
#'
#' Fits a Fay-Herriot model in which the domain-specific covariate effects are
#' allowed to differ across a supplied clustering of the areas, and returns the
#' EBLUP area estimates together with their MSE.
#'
#' @param dom Vector identifying the area (domain) of each observation.
#' @param y Vector of (direct) responses.
#' @param x Covariate matrix.
#' @param index_d Column indices of `x` whose effects vary by cluster.
#' @param vardir Vector of known sampling variances.
#' @param cluster Integer vector giving the cluster label of each unique area.
#' @return A list with the clustered area estimates (`area_est_ada`), the
#'   fitted model object (`fit`), the MSE (`mse`) and residuals (`resids`).
#' @export
fit_clust_area <- function(dom, y, x, index_d, vardir, cluster)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  len_d <- length(index_d)
  
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
  
  df_f <- as.data.frame(cbind(y, vardir, Ux))
  
  forFH <-as.formula(paste("y~0+",paste(colnames(Ux), collapse = "+")))
  res_mse <- mseFH(forFH, vardir = vardir, data = df_f, method ="ML", MAXITER = 500)
  res_fh <- res_mse$est
  mse <- res_mse$mse
  
  
  out <- list(area_est_ada =res_fh$eblup,
              fit =res_fh$fit,
              mse = mse,
              resids = y - Ux %*% res_fh$fit$estcoef[,1])
  return(out)
}


#' Fit a spatial (SAR) Fay-Herriot model for a given cluster structure
#'
#' Spatial Fay-Herriot counterpart of [fit_clust_area()] in which area
#' effects follow a simultaneous autoregressive (SAR) process defined by a
#' proximity matrix.
#'
#' @inheritParams fit_clust_area
#' @param adjW Spatial proximity (adjacency) matrix passed to [sae::mseSFH()].
#' @return A list with the clustered area estimates (`area_est_ada`), the
#'   fitted model object (`fit`), the MSE (`mse`) and residuals (`resids`).
#' @export
fit_clustSAR_area <- function(dom, y, x, index_d, vardir, adjW, cluster)
{
  n0 <- length(y)
  uindex <- unique(dom)
  nr <- length(uindex)
  ns <- as.numeric(table(dom)) # number of replicates
  len_d <- length(index_d)
  
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
    #Ux <- Xm%*%W
    Ux <- Ux
    colnames(Ux) <- names_d
  }else{
    #Ux <- cbind(x[,-index_d, drop = FALSE], Xm%*%W)
    Ux <- cbind(x[,-index_d, drop = FALSE], Ux)
    colnames(Ux) <- c(paste("z",1:(ncol(x)-len_d),sep=""), names_d)
  }
  
  df_f <- as.data.frame(cbind(y, vardir, Ux))
  
  forFH <-as.formula(paste("y~0+",paste(colnames(Ux), collapse = "+")))
  res_mse <- mseSFH(forFH, vardir = vardir,  proxmat = adjW,
                    data = df_f, method ="REML", MAXITER = 500)
  res_fh <- res_mse$est
  mse <- res_mse$mse
  
  
  out <- list(area_est_ada =res_fh$eblup,
              fit =res_fh$fit,
              mse = mse, 
              resids = y - Ux %*% res_fh$fit$estcoef[,1])
  return(out)
}