#' Plot output from the VCCP model
#'
#'  This function plots the change points in the network structure between multivariate time series detected by the VCCP model.
#'
#' @param vccp_result A list generated from \code{\link{vccp.fun}}.
#' @return No return value, called for a plotting purpose.
#' @export
#' @examples
#' \donttest{
#' ## Simulate MVN data with 2 change points
#' data = cbind(1:180, mvn.sim.2.cps(180,8,seed=101))
#' ## Change point detection using VCCP (it may take several minutes to complete...)
#' result = vccp.fun(data, "NBS", test = "V")
#' ## Plot the result
#' getTestPlot(result)
#'
#' result.2 = vccp.fun(data, "NBS", test = "B")
#' ## Plot the result
#' getTestPlot(result.2)
#' }
#' @seealso \code{\link{vccp.fun}}
getTestPlot <- function(vccp_result){
  if(length(vccp_result)==0){
    return(cat("Inappropriate VCCP model specification!"))
  }else{
    if(dim(vccp_result$test_df)[2]==4){
      GetTestPlot.Vuong(vccp_result$test_df, vccp_result$T,
                        vccp_result$sig_alpha)
    }else{
      GetTestPlot.Boot(vccp_result$test_df, vccp_result$T)
    }
  }
}
