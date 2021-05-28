#' vccp: Detect multiple change points in the vine copula structure of multivariate time series by Vine Copula Change Point Model
#'
#' The vccp package implements the Vine Copula Change Point (VCCP)
#' methodology for the estimation of the number and location of multiple
#' change points in the vine copula structure of multivariate time series.
#' The method uses vine copulas, various state-of-the-art segmentation methods
#' to identify multiple change points, and a likelihood ratio test or the
#' stationary bootstrap for inference. The vine copulas allow for various forms
#' of dependence between time series including tail, symmetric and asymmetric
#' dependence. The functions have been extensively tested on simulated multivariate
#' time series data and fMRI data. For details on the VCCP methodology, please see
#' Xiong & Cribben (2021).
#'
#' @section vccp functions:
#' \link{mvn.sim.2.cps}, \link{getTestPlot} and \link{vccp.fun}
#' @examples
#' # See examples in the function vccp.fun.
#'
#' @section Author(s):
#'  Xin Xiong, Ivor Cribben (\email{cribben@@ualberta.ca})
#' @section References:
#'  "Beyond linear dynamic functional connectivity: a vine copula change point model", Xiong and Cribben (2021), bioRxiv 2021.04.25.441254.
#' @docType package
#' @name vccp
NULL
#> NULL
#' @importFrom VineCopula RVineStructureSelect
#' @importFrom VineCopula D2RVine
#' @importFrom VineCopula RVineCopSelect
#' @importFrom VineCopula RVineVuongTest
