#### clustered area level model with SAR(1) model
### graph algorithm based on admm

## for a give lambda

#' ADMM graph-fusion fit for a spatial (SAR) clustered area-level model
#'
#' Spatial counterpart of [Clust_area_graph()]: fuses area-specific
#' coefficients along a graph with an MCP penalty while modelling the area
#' effects with a simultaneous autoregressive (SAR) structure, for a single
#' fusion penalty value.
#'
#' @param dom Vector identifying the area (domain) of each observation.
#' @param y Vector of (direct) responses.
#' @param x Covariate matrix.
#' @param index_d Column indices of `x` whose effects vary by area.
#' @param vardir Vector of known sampling variances.
#' @param graph0 An `igraph` object describing the area adjacency graph.
#' @param betam0 Matrix of initial area-specific coefficients.
#' @param tau0,rho0 Optional starting values for the SAR variance and spatial
#'   parameters; estimated internally when `NULL`.
#' @param adjW Spatial proximity (adjacency) matrix for the SAR process.
#' @param lam Fusion penalty value.
#' @param nu MCP scale parameter.
#' @param gam MCP concavity parameter.
#' @param lam0 Ridge parameter for the within-step solves.
#' @param maxiter Maximum number of ADMM iterations.
#' @param tol Convergence / fusion tolerance.
#' @return A list with the fused coefficients, estimated cluster labels, the
#'   estimated spatial parameter, area estimates and loss.
#' @export
ClustSAR_area_graph <- function(dom, y, x, index_d, vardir, graph0,
                              betam0, tau0 = NULL, rho0 = NULL, 
                              adjW, nu = 1, gam = 3, lam = 0.5,
                              lam0 = 0.001, maxiter= 500, tol = 1e-4)
{
  m <- length(dom)

  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)
  Im <- diag(m)
  
  ### H matrix for the graph
  Hm <- incidentGen(graph0, m)
  Am <- Hm %x% Ip
  AtA <- (t(Hm) %*% Hm) %x% Ip
  
  nL <- nrow(Hm) ### number of edges

  #### fixed design information ########
  x_d <- x[,index_d, drop = FALSE]
  Xm <- bdiag(lapply(seq_len(m), function(i) x_d[i, , drop = FALSE]))

  beta_cur <- c(t(betam0))
  
  nu <- 1
  ##### two cases ###
  ## case 1
  if(len_d == ncol(x))
  {
    ##### initial of others 
    deltam <- Hm %*% betam0
    vm <-  matrix(0, nL , len_d)
    
    adjWt <- t(adjW)  

    if(is.null(tau0) | is.null(rho0))
    {
      resids0 <- as.numeric(y - Xm %*% beta_cur)
      resj <- optim(par = c(mean(resids0^2), 0), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids0,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999))
      tau_cur <- resj$par[1]
      rho_cur <- resj$par[2]
    }else{
      tau_cur <- tau0
      rho_cur <- rho0
    }

    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta (fix 2: Cholesky)
      Cm <- (Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW)
      Cm_chol <- chol(Cm)
      Cm_inv   <- chol2inv(Cm_chol)
      Sm <- diag(vardir) + tau_cur * Cm_inv
      Sm_chol <- chol(Sm)
      Sm_inv   <- chol2inv(Sm_chol)
      Wm <- Sm_inv / (m)/mean(1/vardir)
      Xmt <- t(Xm) %*% Wm

      # nu <- mean(1/vardir)

      beta_new <- solve(Xmt %*% Xm + nu * AtA, Xmt %*% y + nu * c(as.matrix(t(deltam - vm/nu) %*% Hm)))
      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)

      resids_new <- as.numeric(y - Xm %*% beta_new)
      resj <- optim(par = c(tau_cur, rho_cur), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids_new,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999),
                    control = list(maxit = 20))
      tau_new <- resj$par[1]
      rho_new <- resj$par[2]

      Dbeta <- Hm %*% betam_new
      psim  <- Dbeta + vm/nu
      deltam_new <- mcp_vec(psim, lam, nu, gam)

      vm <- vm + nu * (Dbeta - deltam_new)

      diff2 <- sqrt(sum((Dbeta - deltam_new)^2)) / sqrt(nL * len_d)

      deltam   <- deltam_new
      tau_cur  <- tau_new
      rho_cur  <- rho_new
      beta_cur <- beta_new

      if(diff2 < tol) {break}
    }


    cluster_est <- get_sp_cluster_graph(graph0, m, deltam)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    Cm_chol_f <- chol((Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW))
    Cm_inv_f  <- chol2inv(Cm_chol_f)
    Sm_f      <- diag(vardir) + tau_cur * Cm_inv_f
    Sm_inv_f  <- chol2inv(chol(Sm_f))

    muhat     <- rowSums(x * beta_c[cluster_est, -1])
    ratio_est <- tau_cur * Cm_inv_f %*% Sm_inv_f
    area_est  <- muhat + ratio_est %*% (y - muhat)

    out <- list(beta = betam_new, tau = tau_cur, rho = rho_cur,
                betac = beta_c, cluster = cluster_est,
                muhat = muhat, area_est = area_est,
                vardir = vardir, deltam = deltam, niters = j)

  }
  if(len_d < ncol(x))
  {
    # having same part 
    # initial of beta 
    z <- x[,-index_d, drop = FALSE]
    ztz <- solve(t(z) %*% z)%*%t(z)
    Qz <- diag(1, m, m) - z%*% ztz

    eta_cur <- ztz %*% (y - Xm %*% beta_cur)

    ##### initial of others
    deltam <- Hm %*% betam0   # nL x len_d (consistent with vm)
    vm <-  matrix(0, nL, len_d)

    adjWt <- t(adjW)

    if(is.null(tau0) | is.null(rho0))
    {
      resids0 <- as.numeric(y - z %*% eta_cur - Xm %*% beta_cur)
      resj <- optim(par = c(mean(resids0^2), 0), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids0,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999))
      tau_cur <- resj$par[1]
      rho_cur <- resj$par[2]
    }else{
      tau_cur <- tau0
      rho_cur <- rho0
    }

    ### iterations
    for(j in 1:maxiter)
    {
      Cm_chol <- chol((Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW))
      Cm_inv   <- chol2inv(Cm_chol)
      Sm       <- diag(vardir) + tau_cur * Cm_inv
      Sm_inv   <- chol2inv(chol(Sm))
      Wm       <- Sm_inv / (m )/mean(1/vardir)
      ztz      <- solve(t(z) %*% Wm %*% z) %*% t(z) %*% Wm
      Qz       <- Wm - Wm %*% z %*% ztz

      #nu <- mean(1/vardir)

      beta_new <- solve(t(Xm) %*% Qz %*% Xm + nu * AtA,
                        t(Xm) %*% Qz %*% y + nu * c(as.matrix(t(deltam - vm/nu) %*% Hm)))

      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)
      eta_new   <- as.numeric(ztz %*% (y - Xm %*% beta_new))

      resids_new <- as.numeric(y - z %*% eta_new - Xm %*% beta_new)
      resj <- optim(par = c(tau_cur, rho_cur), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids_new,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999),
                    control = list(maxit = 20))
      tau_new <- resj$par[1]
      rho_new <- resj$par[2]

      Dbeta <- as.matrix(Hm %*% betam_new)
      psim  <- Dbeta + vm/nu
      deltam_new <- mcp_vec(psim, lam, nu, gam)

      vm <- vm + nu * (Dbeta - deltam_new)

      diff2 <- sqrt(sum((Dbeta - deltam_new)^2)) / sqrt(nL * len_d)

      deltam  <- deltam_new
      beta_cur <- beta_new
      tau_cur  <- tau_new
      rho_cur  <- rho_new

      if(diff2 < tol) {break}
    }

    cluster_est <- get_sp_cluster_graph(graph0, m, deltam)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    Cm_chol_f <- chol((Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW))
    Cm_inv_f  <- chol2inv(Cm_chol_f)
    Sm_f      <- diag(vardir) + tau_cur * Cm_inv_f
    Sm_inv_f  <- chol2inv(chol(Sm_f))

    muhat     <- z %*% eta_new + rowSums(x_d * beta_c[cluster_est, -1])
    ratio_est <- tau_cur * Cm_inv_f %*% Sm_inv_f
    area_est  <- muhat + ratio_est %*% (y - muhat)

    out <- list(beta = betam_new, eta = eta_new,
                betac = beta_c,
                tau = tau_cur,
                rho = rho_cur,
                cluster = cluster_est,
                muhat = muhat, area_est = area_est, vardir = vardir,
                deltam = deltam, niters = j)

  }
  
  out$beta00 <- betam0
  out$loss <- lossgradSARfun(c(tau_cur, rho_cur), vardir, resids_new, adjW, adjWt)
  out$diff_max <- diff2
  
  
  return(out)
  
}


### this is the function for bic, 
ClustSAR_area_graph_fit <- function(dom, y, x, index_d, vardir, graph0,
                                    Hm, AtA, x_d, Xm, Am,
                                betam0, tau0 = NULL, rho0 = NULL,
                                adjW, nu = 1, gam = 3, lam = 0.5,
                                lam0 = 0.001, maxiter= 500, tol = 1e-4)
{
  m <- length(dom)

  len_d <- length(index_d)
  Ip <- diag(1,len_d,len_d)
  Im <- diag(m)

  nL <- nrow(Hm) ### number of edges
  
  
  beta_cur <- c(t(betam0))
  ##### two cases ###
  ## case 1
  if(len_d == ncol(x))
  {
    ##### initial of others 
    deltam <- Hm %*% betam0
    vm <-  matrix(0, nL , len_d)
    
    adjWt <- t(adjW)  
    
    if(is.null(tau0) | is.null(rho0))
    {
      resids0 <- as.numeric(y - Xm %*% beta_cur)
      resj <- optim(par = c(mean(resids0^2), 0), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids0,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999))
      tau_cur <- resj$par[1]
      rho_cur <- resj$par[2]
    }else{
      tau_cur <- tau0
      rho_cur <- rho0
    }
    
    ### iterations
    for(j in 1:maxiter)
    {
      # update beta and eta (fix 2: Cholesky)
      Cm <- (Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW)
      Cm_chol <- chol(Cm)
      Cm_inv   <- chol2inv(Cm_chol)
      Sm <- diag(vardir) + tau_cur * Cm_inv
      Sm_chol <- chol(Sm)
      Sm_inv   <- chol2inv(Sm_chol)
      Wm <- Sm_inv / (m)/mean(1/vardir)
      Xmt <- t(Xm) %*% Wm
      
      #nu <- mean(1/vardir)
      
      beta_new <- solve(Xmt %*% Xm + nu * AtA, Xmt %*% y + nu * c(as.matrix(t(deltam - vm/nu) %*% Hm)))
      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)
      
      resids_new <- as.numeric(y - Xm %*% beta_new)
      resj <- optim(par = c(tau_cur, rho_cur), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids_new,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999),
                    control = list(maxit = 20))
      tau_new <- resj$par[1]
      rho_new <- resj$par[2]
      
      Dbeta <- Hm %*% betam_new
      psim  <- Dbeta + vm/nu
      deltam_new <- mcp_vec(psim, lam, nu, gam)
      
      vm <- vm + nu * (Dbeta - deltam_new)
      
      diff2 <- sqrt(sum((Dbeta - deltam_new)^2)) / sqrt((nL) * len_d)
      
      deltam   <- deltam_new
      tau_cur  <- tau_new
      rho_cur  <- rho_new
      beta_cur <- beta_new
      
      if(diff2 < tol) {break}
    }
    
    
    cluster_est <- get_sp_cluster_graph(graph0, m, deltam)
    
    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)
    
    Cm_chol_f <- chol((Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW))
    Cm_inv_f  <- chol2inv(Cm_chol_f)
    Sm_f      <- diag(vardir) + tau_cur * Cm_inv_f
    Sm_inv_f  <- chol2inv(chol(Sm_f))
    
    muhat     <- rowSums(x * beta_c[cluster_est, -1])
    ratio_est <- tau_cur * Cm_inv_f %*% Sm_inv_f
    area_est  <- muhat + ratio_est %*% (y - muhat)
    
    out <- list(beta = betam_new, tau = tau_cur, rho = rho_cur,
                betac = beta_c, cluster = cluster_est,
                muhat = muhat, area_est = area_est,
                vardir = vardir, deltam = deltam, niters = j)
    
  }
  if(len_d < ncol(x))
  {
    # having same part 
    # initial of beta 
    z <- x[,-index_d, drop = FALSE]
    ztz <- solve(t(z) %*% z)%*%t(z)
    Qz <- diag(1, m, m) - z%*% ztz
    
    eta_cur <- ztz %*% (y - Xm %*% beta_cur)
    
    ##### initial of others 
    deltam <- Hm %*% betam0
    vm <-  matrix(0, nL , len_d)
    
    adjWt <- t(adjW)  
    
    if(is.null(tau0) | is.null(rho0))
    {
      resids0 <- as.numeric(y - z %*% eta_cur - Xm %*% beta_cur)
      resj <- optim(par = c(mean(resids0^2), 0), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids0,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999))
      tau_cur <- resj$par[1]
      rho_cur <- resj$par[2]
    }else{
      tau_cur <- tau0
      rho_cur <- rho0
    }
    
    ### iterations
    for(j in 1:maxiter)
    {
      Cm_chol <- chol((Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW))
      Cm_inv   <- chol2inv(Cm_chol)
      Sm       <- diag(vardir) + tau_cur * Cm_inv
      Sm_inv   <- chol2inv(chol(Sm))
      Wm       <- Sm_inv / (m)/mean(1/vardir)
      ztz      <- solve(t(z) %*% Wm %*% z) %*% t(z) %*% Wm
      Qz       <- Wm - Wm %*% z %*% ztz
      
      #nu <- mean(1/vardir)
      
      beta_new <- solve(t(Xm) %*% Qz %*% Xm + nu * AtA,
                        t(Xm) %*% Qz %*% y + nu * c(as.matrix(t(deltam - vm/nu) %*% Hm)))
      
      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)
      eta_new   <- as.numeric(ztz %*% (y - Xm %*% beta_new))
      
      resids_new <- as.numeric(y - z %*% eta_new - Xm %*% beta_new)
      resj <- optim(par = c(tau_cur, rho_cur), fn = lossgradSARfun,
                    gr = gradSARfun_new, vardir = vardir, resids = resids_new,
                    adjW = adjW, adjWt = adjWt,
                    method = "L-BFGS-B",
                    lower = c(0, -0.999), upper = c(mean(y^2), 0.999),
                    control = list(maxit = 20))
      tau_new <- resj$par[1]
      rho_new <- resj$par[2]
      
      Dbeta <- as.matrix(Hm %*% betam_new)
      psim  <- Dbeta + vm/nu
      deltam_new <- mcp_vec(psim, lam, nu, gam)
      
      vm <- vm + nu * (Dbeta - deltam_new)
      
      diff2 <- sqrt(sum((Dbeta - deltam_new)^2)) / sqrt((nL) * len_d)
      
      deltam  <- deltam_new
      beta_cur <- beta_new
      tau_cur  <- tau_new
      rho_cur  <- rho_new
      
      if(diff2 < tol) {break}
    }
    
    cluster_est <- get_sp_cluster_graph(graph0, m, deltam)
    
    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)
    
    Cm_chol_f <- chol((Im - rho_cur * adjWt) %*% (Im - rho_cur * adjW))
    Cm_inv_f  <- chol2inv(Cm_chol_f)
    Sm_f      <- diag(vardir) + tau_cur * Cm_inv_f
    Sm_inv_f  <- chol2inv(chol(Sm_f))
    
    muhat     <- z %*% eta_new + rowSums(x_d * beta_c[cluster_est, -1])
    ratio_est <- tau_cur * Cm_inv_f %*% Sm_inv_f
    area_est  <- muhat + ratio_est %*% (y - muhat)
    
    out <- list(beta = betam_new, eta = eta_new,
                betac = beta_c,
                tau = tau_cur,
                rho = rho_cur,
                cluster = cluster_est,
                muhat = muhat, area_est = area_est, vardir = vardir,
                deltam = deltam, niters = j)
    
  }
  
  out$beta00 <- betam0
  out$loss <- lossgradSARfun(c(tau_cur, rho_cur), vardir, resids_new, adjW, adjWt)
  out$diff_max <- diff2
  
  
  return(out)
  
}


###bic #####
#' SAR graph-fusion path with BIC selection grid
#'
#' Runs the spatial (SAR) graph-fusion fit over a grid of fusion penalties and
#' records the estimated number of clusters and loss at each value.
#'
#' @inheritParams ClustSAR_area_graph
#' @param lamvec Vector of fusion penalty values to evaluate.
#' @return A list with `res_df` (grid of `lambda`, `Khat`, `loss`) and `out`
#'   (the fitted object for each penalty value).
#' @export
ClustSAR_area_graph_bic <- function(dom, y, x, index_d, vardir, graph0,
                                    adjW, lamvec,
                                    nu = 1, gam = 3, lam0 = 0.001,
                                    maxiter= 500, tol = 1e-4)
{
  m <- length(unique(dom))
  nlam <- length(lamvec)
  len_d <- length(index_d)
  
  ### H matrix for the graph
  Ip <- diag(1,len_d,len_d)
  Hm <- incidentGen(graph0, m)
  Am <- Hm %x% Ip
  AtA <- (t(Hm) %*% Hm) %x% Ip
  
  adjWt <- t(adjW)  
  
  #### fixed design information ########
  x_d <- x[,index_d, drop = FALSE]
  Xm <- bdiag(lapply(seq_len(m), function(i) x_d[i, , drop = FALSE]))
  

  K_est <- rep(0, nlam)
  loss_value <- rep(0, nlam)
  out <- list()
  
  betam0 <- initialFUNs(dom = dom, y = y, x =x,
                       index_d = index_d, graph0, 
                       lam0 = lam0)
  res00 <- fit_clustSAR_area(dom, y, x, index_d, vardir, adjW, cluster = rep(1,m))
  
  tau0 <- res00$fit$refvar
  rho0 <- res00$fit$spatialcorr

  
  for(j in 1:length(lamvec))
  {
    resj <- ClustSAR_area_graph_fit(dom, y, x, index_d, vardir, graph0 = graph0,
                                    Hm = Hm, AtA = AtA, x_d = x_d, Xm = Xm,
                                    betam0 = betam0,tau0 = tau0, rho0 = rho0,
                                    adjW = adjW,
                                    nu = nu, gam = gam, lam = lamvec[j],
                                    lam0 = lam0, maxiter= maxiter, tol = tol)
    loss_value[j] <- resj$loss
    K_est[j] <- length(unique(resj$cluster))
    out[[j]] <- resj
  }
  
  
  
  res_df <- data.frame(lambda = lamvec,Khat = K_est, 
                       loss = loss_value)
  
  return(list(res_df = res_df, out = out))
}


#### estimator ####


#' Area-level estimator via SAR graph fusion with BIC model selection
#'
#' Top-level spatial estimator: fits the SAR graph-fusion path with
#' [ClustSAR_area_graph_bic()], selects the penalty minimising a modified BIC,
#' and returns the selected clustering together with the (re-fitted) area
#' estimates.
#'
#' @inheritParams ClustSAR_area_graph_bic
#' @return A list containing the BIC table, the selected fit, the estimated
#'   cluster labels and the area estimates.
#' @export
Est_clustSAR_area_graph <- function(dom, y, x, index_d, vardir, graph0,
                                    adjW, lamvec, nu = 1, gam = 3,
                                    lam0 = 0.001, maxiter= 500, tol = 1e-4)
{
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  m <- length(unique(dom))
  
  res_bic <- ClustSAR_area_graph_bic(dom=dom, y, x, index_d, vardir, 
                                       graph0= graph0, adjW = adjW,
                                       lamvec = lamvec, 
                                       nu = nu, gam = gam, 
                                       lam0 = lam0, maxiter= maxiter,  
                                       tol = tol)

  
  bic_df <- res_bic$res_df %>% 
    mutate(mbic = loss + log(m)*log(m)*(Khat*len_d + nz +1))
  
  res_fm <- res_bic$out[[which.min(bic_df$mbic)]]
  area_estm <- res_fm$area_est
  
  
  fit_fm <- NULL
  tryCatch({
    fit_fm <- fit_clustSAR_area(dom = dom, y = y, x = x, index_d = index_d,
                                vardir = vardir, adjW = adjW, cluster = res_fm$cluster)
    area_estm <- fit_fm$area_est_ada
  }, error = function(e) {
    area_estm <- res_fm$area_est
  })
  
  

  out <- list(bic_df = bic_df,
              resfm = res_fm,
              cluster_estm = res_fm$cluster,
              area_estm = area_estm,
              refit = fit_fm,
              out_bic = res_bic)
  
  return(out)
  
}


.sar_cache <- new.env(parent = emptyenv())

lossgradSARfun <- function(xx, vardir, resids, adjW, adjWt)
{
  tau <- xx[1]; rho <- xx[2]
  Im  <- diag(1, length(vardir))
  adjW  <- as.matrix(adjW)
  adjWt <- as.matrix(adjWt)

  Cm      <- (Im - rho * adjWt) %*% (Im - rho * adjW)
  Cm_chol <- chol(Cm)
  Cm_inv  <- chol2inv(Cm_chol)
  Sm      <- diag(vardir) + tau * Cm_inv
  Sm_chol <- chol(Sm)
  Sm_inv  <- chol2inv(Sm_chol)

  # cache intermediates for gradient reuse
  .sar_cache$xx     <- xx
  .sar_cache$Cm_inv <- Cm_inv
  .sar_cache$Sm_inv <- Sm_inv
  .sar_cache$adjWt  <- adjWt
  .sar_cache$adjW   <- adjW

  log_det <- 2 * sum(log(diag(Sm_chol)))
  Sm_inv_r <- Sm_inv %*% resids
  loss <- as.numeric(log_det + t(resids) %*% Sm_inv_r)
  return(loss)
}

gradSARfun_new <- function(xx, vardir, resids, adjW, adjWt)
{
  # reuse cache if same (tau,rho) was just evaluated by lossgradSARfun
  if (!is.null(.sar_cache$xx) && identical(xx, .sar_cache$xx)) {
    Cm_inv <- .sar_cache$Cm_inv
    Sm_inv <- .sar_cache$Sm_inv
    adjWt  <- .sar_cache$adjWt
    adjW   <- .sar_cache$adjW
  } else {
    tau <- xx[1]; rho <- xx[2]
    Im  <- diag(1, length(vardir))
    adjW  <- as.matrix(adjW)
    adjWt <- as.matrix(adjWt)
    Cm     <- (Im - rho * adjWt) %*% (Im - rho * adjW)
    Cm_inv <- chol2inv(chol(Cm))
    Sm     <- diag(vardir) + xx[1] * Cm_inv
    Sm_inv <- chol2inv(chol(Sm))
  }
  tau <- xx[1]; rho <- xx[2]

  Sm_inv_r <- Sm_inv %*% resids

  # Sm_inv, Cm_inv (and Cm_inv_dC_Cm_inv) are symmetric, so traces reduce to
  # elementwise products (O(m^2) instead of forming the full O(m^3) matrix product),
  # and quadratic forms reduce to vector chains via Sm_inv_r = Sm_inv %*% resids.
  value_tau1 <- sum(Sm_inv * Cm_inv)
  value_tau2 <- as.numeric(crossprod(Sm_inv_r, Cm_inv %*% Sm_inv_r))
  value_tau  <- value_tau1 - value_tau2

  dC <- adjWt - 2 * rho * adjWt %*% adjW + adjW   # dC/d(rho)
  Cm_inv_dC_Cm_inv <- Cm_inv %*% dC %*% Cm_inv
  value_rho1 <- tau * sum(Sm_inv * Cm_inv_dC_Cm_inv)
  value_rho2 <- tau * as.numeric(crossprod(Sm_inv_r, Cm_inv_dC_Cm_inv %*% Sm_inv_r))
  value_rho  <- value_rho1 - value_rho2

  return(c(value_tau, value_rho))
}










