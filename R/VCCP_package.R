#' VCCP: Detect change points in the functional connectivity by Vine Copula Change Point Model
#'
#' The vccp package implements the Vine Copula Change Point (VCCP)
#' methodology for the estimation of the number and location of multiple
#' change points in the network structure of multivariate time series.
#' The method uses vine copulas and adapted binary segmentation algorithm to
#' identify multiple change points and a likelihood ratio test for inference.
#' The vine copulas allow for various forms of dependence between time series
#' including tail, symmetric and asymmetric dependence. The functions have been
#' extensively tested on simulated multivariate time series data and fMRI data.
#' For details on the VCCP methodology, please see Xiong & Cribben (2021).
#'
#' @section VCCP functions:
#' \link{mvn.sim.2.cps}, \link{getTestPlot} and \link{vccp.fun}
#' @examples
#' # See examples in the function ``vccp.fun''.
#'
#' @section Author(s):
#'  Xin Xiong, Ivor Cribben (\email{cribben@@ualberta.ca})
#' @section References:
#'  "Beyond linear dynamic functional connectivity: a vine copula change point model", Xiong and Cribben (2021), preprint.
#' @docType package
#' @name VCCP
NULL
#> NULL
#' @importFrom VineCopula RVineStructureSelect
#' @importFrom CDVine CDVineCopSelect
#' @importFrom CDVine CDVineBIC
