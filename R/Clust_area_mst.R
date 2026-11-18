#### MST adaptive approach ####

#' MST fusion path with BIC selection grid for clustered area models
#'
#' Runs the minimum-spanning-tree (MST) adaptive fusion algorithm over a grid
#' of penalty values, returning the fitted object and the estimated number of
#' clusters and loss at each value.
#'
#' @param dom Vector identifying the area (domain) of each observation.
#' @param y Vector of (direct) responses.
#' @param x Covariate matrix.
#' @param index_d Column indices of `x` whose effects vary by area.
#' @param vardir Vector of known sampling variances.
#' @param graph0 An `igraph` object describing the area adjacency graph.
#' @param lamvec Vector of fusion penalty values to evaluate.
#' @param rho ADMM augmented-Lagrangian parameter.
#' @param gam MCP concavity parameter.
#' @param lam0 Ridge parameter for the within-step solves.
#' @param maxiter Maximum number of ADMM iterations.
#' @param tol Convergence / fusion tolerance.
#' @return A list with `res_df` (grid of `lambda`, `Khat`, `loss`) and `out`
#'   (the fitted object for each penalty value).
#' @export
Clust_area_mst_bic <- function(dom, y, x, index_d, vardir,
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
  
  ### initial mst 
  
  graph_mst0 <- mst(graph0)
  
  
  betam0 <- initialFUNs(dom = dom, y = y, x =x,
                       index_d = index_d, graph0,
                       lam0 = 0.001)
  
  
  
  for(j in 1:length(lamvec))
  {
    resj <- Clust_area_graph(dom, y, x, index_d, vardir, 
                             graph0 = graph_mst0, betam0 = betam0,
                             rho = rho, gam = gam, lam = lamvec[j],
                             lam0 = lam0, maxiter= maxiter,
                             tol = tol)
    
    graph_updated <- graph0
    E(graph_updated)$weight <- apply(incidentGen(graph0, m) %*% as.matrix(resj$betac[resj$cluster,-1]), 1, norm, "2")
    
    betam0_up <- as.matrix(resj$betac[resj$cluster,-1])
    graph_mst_updated <- mst(graph_updated)
    
    
    resj <- Clust_area_graph(dom, y, x, index_d, vardir, 
                             graph0 = graph_mst_updated, betam0 = betam0_up,
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


### area level estimator ###

#' Area-level estimator via MST fusion with BIC model selection
#'
#' Top-level estimator: fits the MST fusion path with [Clust_area_mst_bic()],
#' selects the penalty minimising a modified BIC, and returns the selected
#' clustering together with the (re-fitted) area estimates.
#'
#' @inheritParams Clust_area_mst_bic
#' @return A list containing the BIC table (`bic_df`), the selected fit
#'   (`resfm`), the estimated cluster labels (`cluster_estm`), the area
#'   estimates (`area_estm`) and the full BIC object (`out_bic`).
#' @export
Est_clust_area_mst <- function(dom, y, x, index_d, vardir,
                                 graph0, lamvec, 
                                 rho = 1, gam = 3, 
                                 lam0 = 0.001, maxiter= 500, tol = 1e-4)
{
  len_d <- length(index_d)
  nz <- ncol(x) - len_d
  m <- length(dom)
  
  res_bic <- Clust_area_mst_bic(dom=dom, y, x, index_d, vardir, 
                                  graph0= graph0, lamvec = lamvec, 
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
