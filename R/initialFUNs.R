#### initial values ####

#' Ridge-type initial values for the clustered area-level algorithm
#'
#' Computes starting values for the area-specific coefficient matrix used to
#' initialise the ADMM fusion algorithms, via a graph-penalised ridge fit.
#'
#' @param dom Vector identifying the area (domain) of each observation.
#' @param y Vector of (direct) responses.
#' @param x Covariate matrix.
#' @param index_d Column indices of `x` whose effects vary by area.
#' @param graph0 An `igraph` object describing the area adjacency graph.
#' @param lam0 Ridge penalty parameter.
#' @param k0 Unused; retained for backward compatibility.
#' @return A matrix of initial area-specific coefficients (one row per area).
#' @export
initialFUNs <- function(dom, y, x, index_d, graph0, lam0 = 0.001, k0 = 10)
{
  n0 <- length(y)
  uindex <- unique(dom)
  m <- length(uindex)
  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)
  
  Hm <- incidentGen(graph0, m)
  AtA <- (t(Hm) %*% Hm) %x% Ip
  
  x_d <- x[,index_d, drop = FALSE]
  Xm <- matrix(0, n0, m*len_d)
  for(i in 1:m)
  {
    Xm[dom == uindex[i],(len_d*(i-1) + 1) : (len_d*i)] <- x_d[dom == uindex[i],]
  }
  
  
  if(len_d == ncol(x))
  {
    beta00 <- matrix(solve(t(Xm) %*% Xm + lam0*AtA) %*% t(Xm) %*% y,
                     ncol = len_d, byrow = TRUE) 
  }
  
  if(len_d < ncol(x))
  {
    z <- x[,-index_d, drop = FALSE]
    ztz <- solve(t(z) %*% z)%*%t(z)
    Qz <- diag(1, n0, n0) - z%*% ztz
    
    beta00 <- matrix(solve(t(Xm) %*% Qz%*% Xm + lam0*AtA) %*% t(Xm)%*% Qz %*% y, 
                     ncol = len_d, byrow = TRUE)
  }
  
  
  
  return(beta00)
}

