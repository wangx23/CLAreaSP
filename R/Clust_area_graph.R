#### clustered area level model without SAR model
### delta = A beta

#' ADMM graph-fusion fit for a clustered area-level model
#'
#' Core ADMM solver that fuses area-specific coefficients along the edges of a
#' supplied graph using an MCP penalty, for a single fusion penalty value.
#'
#' @param dom Vector identifying the area (domain) of each observation.
#' @param y Vector of (direct) responses.
#' @param x Covariate matrix.
#' @param index_d Column indices of `x` whose effects vary by area.
#' @param vardir Vector of known sampling variances.
#' @param graph0 An `igraph` object describing the area adjacency graph.
#' @param betam0 Matrix of initial area-specific coefficients.
#' @param rho ADMM augmented-Lagrangian parameter.
#' @param gam MCP concavity parameter.
#' @param lam Fusion penalty value.
#' @param lam0 Ridge parameter for the within-step solves.
#' @param maxiter Maximum number of ADMM iterations.
#' @param tol Convergence / fusion tolerance.
#' @return A list with the fused coefficients, estimated cluster labels,
#'   area estimates and loss for the given penalty value.
#' @export
Clust_area_graph <- function(dom, y, x, index_d, vardir, graph0, betam0,
                           rho = 1, gam = 3, lam = 0.5,
                           lam0 = 0.001, maxiter= 500, tol = 1e-4)
{

  m <- length(dom)
  len_d <- length(index_d)
  Ip <- Matrix::Diagonal(len_d)          # sparse identity

  rho <- 1

  ### graph and H matrix (kept sparse) ###
  Hm <- incidentGen(graph0, m)
  AtA <- Matrix::kronecker(Matrix::crossprod(Hm), Ip)  # (t(Hm) %*% Hm) %x% Ip

  nL <- nrow(Hm) ### number of edges

  #### fixed design information: Xm as a sparse block matrix ########
  x_d <- x[,index_d, drop = FALSE]
  Xm <- Matrix::sparseMatrix(
    i = rep(seq_len(m), each = len_d),
    j = as.integer(unlist(lapply(seq_len(m), function(i) (i-1)*len_d + seq_len(len_d)))),
    x = as.vector(t(x_d)),
    dims = c(m, m*len_d))

  beta_cur <- c(t(betam0))

  ##### two cases ###
  ## case 1
  if(len_d == ncol(x))
  {
    ##### initial of others
    deltam <- as.matrix(Hm %*% betam0)
    vm <-  matrix(0, nL, len_d)

    resids0 <- as.numeric(y - Xm %*% beta_cur)

    lossj <- function(xx)
    {
      lossfun2(xx, vardir, resids0)
    }

    sig2_cur <- optimize(lossj, interval = c(0,mean(resids0^2)))$minimum

    ### iterations
    chol_cache <- NULL   # reuse the symbolic factorisation across iterations
    for(j in 1:maxiter)
    {
      # update beta
      wm <- 1/(m*(vardir + sig2_cur))/mean(1/vardir)

      # XtWX is block-diagonal, AtA is the graph Laplacian %x% Ip: both sparse
      Amat <- Matrix::forceSymmetric(Matrix::crossprod(Xm, wm * Xm) + rho * AtA)
      rhs  <- as.numeric(Matrix::crossprod(Xm, wm * y)) +
              rho * as.numeric(Matrix::crossprod(deltam - vm/rho, Hm))  # len_d x m, area-major

      beta_new <- tryCatch({
        if(is.null(chol_cache)) chol_cache <- Matrix::Cholesky(Amat, LDL = FALSE)
        else                    chol_cache <- Matrix::update(chol_cache, Amat)
        as.numeric(Matrix::solve(chol_cache, rhs))
      }, error = function(e) as.numeric(Matrix::solve(Amat, rhs)))

      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)

      # update sigma_v^2
      resids_new <- as.numeric(y - Xm %*% beta_new)

      lossj <- function(xx)
      {
        lossfun2(xx, vardir, resids_new)
      }

      sig2_new <- optimize(lossj, interval = c(0,mean(resids_new^2)))$minimum

      Dbeta <- as.matrix(Hm %*% betam_new)
      psim <- Dbeta + vm/rho
      deltam_new <- mcp_vec(psim, lam, rho, gam)   # vectorised MCP over edges

      vm <- vm + rho * (Dbeta - deltam_new)

      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      sig2_cur <- sig2_new

      if(diff2 < tol)
      {break}


    }


    cluster_est <- get_sp_cluster_graph(graph0, m, deltam)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)


    muhat <- rowSums(x * beta_c[cluster_est, -1])
    ratio_est <- sig2_cur/(vardir + sig2_cur)
    area_est <- ratio_est * y + (1-ratio_est)*muhat

    out <- list(beta = betam_new, sig2 = sig2_cur,
                betac = beta_c, cluster = cluster_est,
                muhat = muhat, area_est = area_est,
                vardir = vardir,deltam = deltam, niters = j)

  }
  if(len_d < ncol(x))
  {
    # having same part
    # initial of eta
    z <- x[,-index_d, drop = FALSE]
    ztz <- solve(t(z) %*% z)%*%t(z)

    eta_cur <- ztz %*% (y - Xm%*%beta_cur)

    ##### initial of others
    deltam <- as.matrix(Hm %*% betam0)   # nL x len_d (consistent with vm)
    vm <-  matrix(0, nL , len_d)

    resids0 <- as.numeric(y - z %*% eta_cur - Xm %*% beta_cur)

    lossj <- function(xx)
    {
      lossfun2(xx, vardir, resids0)
    }

    sig2_cur <- optimize(lossj, interval = c(0,mean(resids0^2)))$minimum


    ### iterations
    chol_cache <- NULL   # reuse the symbolic factorisation across iterations
    for(j in 1:maxiter)
    {
      # update beta and eta.
      # The system matrix is  M = Xm' Qz Xm + rho*AtA  with the projection
      #   Qz = Wm - Wm z (z'Wz)^-1 z'W.  Written out,  M = S - U C^-1 U'  where
      #   S = Xm'Wm Xm + rho*AtA  (sparse: block-diagonal data term + Laplacian%x%Ip),
      #   U = Xm'Wm z (p x nz),  C = z'Wm z (nz x nz).
      # The z-part is a rank-nz correction, so a Woodbury update lets us keep a
      # cached sparse Cholesky of S and never form the dense p x p solve or n x n Qz.
      wm <- 1/m*(1/(vardir + sig2_cur)) / mean(1/vardir)

      S    <- Matrix::forceSymmetric(Matrix::crossprod(Xm, wm * Xm) + rho * AtA)
      U    <- as.matrix(Matrix::crossprod(Xm, wm * z))   # p x nz
      Cmat <- as.matrix(Matrix::crossprod(z,  wm * z))   # nz x nz
      pen  <- rho * as.numeric(Matrix::crossprod(deltam - vm/rho, Hm))  # len_d x m, area-major
      r    <- as.numeric(Matrix::crossprod(Xm, wm * y)) + pen -
              U %*% solve(Cmat, as.numeric(Matrix::crossprod(z, wm * y)))

      beta_new <- tryCatch({
        if(is.null(chol_cache)) chol_cache <- Matrix::Cholesky(S, LDL = FALSE)
        else                    chol_cache <- Matrix::update(chol_cache, S)
        SinvU <- as.matrix(Matrix::solve(chol_cache, U))   # p x nz
        Sinvr <- as.numeric(Matrix::solve(chol_cache, r))  # p
        G     <- Cmat - crossprod(U, SinvU)                # nz x nz Schur complement
        as.numeric(Sinvr + SinvU %*% solve(G, crossprod(U, Sinvr)))
      }, error = function(e) {
        # dense fallback (original formulation)
        Qz <- Matrix::Diagonal(x = wm) - (wm * z) %*% solve(Cmat, t(wm * z))
        as.numeric(solve(as.matrix(Matrix::crossprod(Xm, Qz %*% Xm) + rho * AtA),
                         as.numeric(Matrix::crossprod(Xm, Qz %*% y)) + pen))
      })

      betam_new <- matrix(beta_new, m, len_d, byrow = TRUE)
      vresid    <- as.numeric(y - Xm %*% beta_new)
      eta_new   <- solve(Cmat, as.numeric(Matrix::crossprod(z, wm * vresid)))  # (z'Wz)^-1 z'W(y-Xb)


      # update sigma_v^2
      resids_new <- as.numeric(vresid - z %*% eta_new)

      lossj <- function(xx)
      {
        lossfun2(xx, vardir, resids_new)
      }

      sig2_new <- optimize(lossj, interval = c(0,mean(resids_new^2)))$minimum

      Dbeta <- as.matrix(Hm %*% betam_new)
      psim <- Dbeta + vm/rho

      deltam_new <- mcp_vec(psim, lam, rho, gam)   # vectorised MCP over edges

      vm <- vm + rho * (Dbeta - deltam_new)

      diff2 <- max(abs(Dbeta- deltam_new))

      deltam <- deltam_new
      sig2_cur <- sig2_new
      beta_cur <- beta_new

      if(diff2 < tol)
      {break}
    }

    cluster_est <- get_sp_cluster_graph(graph0, m, deltam)

    beta_c <- aggregate(betam_new, by = list(cluster_est), FUN = mean)

    muhat <- as.numeric(z %*% eta_new) +  rowSums(x_d * beta_c[cluster_est,-1])
    ratio_est <- sig2_cur/(vardir + sig2_cur)
    area_est <- ratio_est * y + (1-ratio_est)*muhat

    out <- list(beta = betam_new, eta = eta_new,
                betac = beta_c, sig2 = sig2_cur,
                cluster = cluster_est,
                muhat = muhat, area_est = area_est, vardir = vardir,
                deltam = deltam, niters = j)

  }

  out$beta00 <- betam0
  out$loss <- lossfun2(sig2_new, vardir, resids_new)
  out$diff_max <- diff2

  return(out)
}


###bic 
#' Graph-fusion path with BIC selection grid for clustered area models
#'
#' Runs [Clust_area_graph()] over a grid of fusion penalties and records the
#' estimated number of clusters and loss at each value.
#'
#' @inheritParams Clust_area_graph
#' @param lamvec Vector of fusion penalty values to evaluate.
#' @return A list with `res_df` (grid of `lambda`, `Khat`, `loss`) and `out`
#'   (the fitted object for each penalty value).
#' @export
Clust_area_graph_bic <- function(dom, y, x, index_d, vardir,
                                   graph0, lamvec,
                                   rho = 1, gam = 3, lam0 = 0.01,
                                   maxiter= 500, tol = 1e-4)
{
  m <- length(dom)
  nlam <- length(lamvec)
  len_d <- length(index_d)
  
  
  K_est <- rep(0, nlam)
  loss_value <- rep(0, nlam)
  out <- list()
  
  betam0 <- initialFUNs(dom = dom, y = y, x =x,
                       index_d = index_d, graph0, 
                       lam0 = lam0)
  
  
  
  for(j in 1:length(lamvec))
  {
    resj <- Clust_area_graph(dom, y, x, index_d, vardir, 
                               graph0 = graph0, betam0 = betam0,
                               rho = rho, gam = gam, lam = lamvec[j],
                               lam0 = lam0, maxiter= maxiter,
                               tol = tol)

    loss_value[j] <- resj$loss
    K_est[j] <- length(unique(resj$cluster))
    out[[j]] <- resj
  }
  
  
  
  res_df <- data.frame(lambda = lamvec,
                       Khat = K_est, 
                       loss = loss_value)
  
  return(list(res_df = res_df, out = out))
}



##### area estimate #####

#' Area-level estimator via graph fusion with BIC model selection
#'
#' Top-level estimator: fits the graph-fusion path with
#' [Clust_area_graph_bic()], selects the penalty minimising a modified BIC, and
#' returns the selected clustering together with the (re-fitted) area estimates.
#'
#' @inheritParams Clust_area_graph_bic
#' @return A list containing the BIC table, the selected fit, the estimated
#'   cluster labels and the area estimates.
#' @export
Est_clust_area_graph <- function(dom, y, x, index_d, vardir,
                                graph0,lamvec,
                                rho = 1, gam = 3, 
                                lam0 = 0.001, maxiter= 500, tol = 1e-4)
{
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  m <- length(dom)
  
  res_bic <- Clust_area_graph_bic(dom=dom, y, x, index_d, vardir, 
                                  graph0= graph0,lamvec = lamvec, 
                                  rho = rho, gam = gam, 
                                  lam0 = lam0, maxiter= maxiter,  
                                  tol = tol)
  
  bic_df <- res_bic$res_df %>% 
    mutate(mbic = loss + log(m)*log(m)*(Khat*len_d + nz +1))
  
  res_fm <- res_bic$out[[which.min(bic_df$mbic)]]
  area_estm <- res_fm$area_est
  
  
  tryCatch({
    fit_fm <- fit_clust_area(dom = dom, y = y, x = x, index_d = index_d,
                             vardir = vardir, cluster = res_fm$cluster)
    area_estm <- fit_fm$area_est_ada
  }, error = function(e) {
    area_estm <- res_fm$area_est
  })
  

  out <- list(bic_df = bic_df,
              resfm = res_fm,
              cluster_estm = res_fm$cluster,
              area_estm = area_estm,
              out_bic = res_bic)
  
  return(out)
  
}


## grad of function of sigma_v^2 
# resids y - zeta - xbeta 
gradfun2 <- function(sigv2, vardir, resids)
{
  value1 <- sum(1/(vardir + sigv2))
  value2 <- sum(resids^2/(vardir + sigv2)^2)
  value <- value1 - value2
  return(value)
}

lossfun2 <- function(sigv2, vardir, resids)
{
  value1 <-  sum(log(vardir + sigv2))
  value2 <- sum(resids^2/(vardir + sigv2))
  value <- value1 + value2
  return(value)
}



  



