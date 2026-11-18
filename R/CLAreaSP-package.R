#' CLAreaSP: Spatially Clustered Effects Area Level Models for Survey Data
#'
#' Spatially clustered area-level small area estimation under survey
#' sampling. The package implements ADMM-based graph and minimum-spanning-tree
#' fusion algorithms with MCP penalties, Fay-Herriot and spatial (SAR)
#' Fay-Herriot fitting, BIC model selection, and area-level predictors for
#' clustered random effects.
#'
#' @importFrom dplyr filter mutate select row_number
#' @importFrom igraph as_data_frame mst E "E<-" components graph_from_data_frame
#' @importFrom magrittr %>%
#' @import Matrix
#' @importClassesFrom Matrix dgCMatrix
#' @import sae
#' @keywords internal
"_PACKAGE"
