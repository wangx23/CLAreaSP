# some functions used in the algorithm
#### MCP and SCAD ####
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

### penalty function 
mcp_pen <- function(x, lam, gam)
{
  xn <- sqrt(sum(x^2))
  if (xn <=  gam * lam){
    lam * xn - xn ^ 2 / (2 * gam) 
  }else{
    gam * lam^2 / 2
  }
}


mcp_vec <- function(psim, lam, nu, gam)
{
  norms <- sqrt(rowSums(psim^2))
  temp  <- gam * lam
  # soft-threshold scale for rows where MCP applies
  sf    <- pmax(1 - (lam / nu) / pmax(norms, .Machine$double.eps), 0) / (1 - 1 / (gam * nu))
  apply_mcp <- norms <= temp
  out <- psim
  out[apply_mcp, ] <- sf[apply_mcp] * psim[apply_mcp, , drop = FALSE]
  out
}




### linear regression model without vardir ###
lm_clust_reg <- function(dom, y, x, index_d, cluster)
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
  
  # df_f <- as.data.frame(cbind(y, Ux))
  
  coefu <- as.numeric(solve(t(Ux) %*% Ux + 0.001*diag(ncol(Ux))) %*% t(Ux) %*% y)
  # forlm <-as.formula(paste("y~0+",paste(colnames(Ux), collapse = "+")))
  # res_lm <- lm(forlm, data = df_f)
  # 
  # return(coef(res_lm))
  return(coefu)
}


### linear regression model with vardir ###
## use iterative to estimate everything
lm_clust_reg_ridge <- function(dom, y, x, vardir, index_d, cluster)
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
  
  # df_f <- as.data.frame(cbind(y, Ux))
  # forlm <-as.formula(paste("y~0+",paste(colnames(Ux), collapse = "+")))
  # res_lm <- lm(forlm, data = df_f)
  # return(coef(res_lm))
  
  coefu <- as.numeric(solve(t(Ux*wt) %*% Ux + 0.001*diag(ncol(Ux))) %*% t(Ux*wt) %*% y)
  
  return(coefu)
}




### get clusters 

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



get_sp_cluster_graph <- function(graph0,m,deltam, tol = 1e-4) {
  E(graph0)$weight <- apply(deltam, 1, norm, "2")
  link_matrix <- igraph::as_data_frame(graph0) %>%
    dplyr::mutate(edge = row_number()) %>%
    filter(weight < tol) %>%
    dplyr::select(from, to)
  graph_fused <- graph_from_data_frame(link_matrix, vertices = 1:m, directed = FALSE)
  components(graph_fused)$membership
}


incidentGen <- function(graph, m) {
  # return the incident matrix of a graph
  linkMatrix <- igraph::as_data_frame(graph) %>%
    dplyr::mutate(edge = row_number()) %>% dplyr::select(edge, from, to)
  linkMatrix$from <- linkMatrix$from %>% as.integer()
  linkMatrix$to <- linkMatrix$to %>% as.integer()
  H <- as(matrix(0, dim(linkMatrix)[1], m), "dgCMatrix")
  H[as.matrix(linkMatrix[, c(1, 2)])] = 1
  H[as.matrix(linkMatrix[, c(1, 3)])] = -1
  H
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



