#' VCCP: Detect change points in the functional connectivity by Vine Copula Change Point Model
#'
#' The VCCP package uses Vine Copula Change Points (VCCP) model to estimate the number and location of multiple
#' change points in the functional connectivity structure of multivariate
#' time series. Possible binary segmentation methods of the VCCP model
#' include NBS, OBS, MOSUM, and WBS. You can also check the
#' significance of detected candidates by Stationary Bootstrap or
#' Vuong test method. This package provides one main function
#' to implement the VCCP model with two auxiliary functions to
#' generate multivariate normal data with 2 change points
#' and plot the result.
#'
#' @section VCCCP functions:
#' \link{random.mvn.simulate.2.changes}, \link{getTestPlot} and \link{vccp.fun}
#' @examples
#' # See examples for the function ``vccp.fun''.
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
