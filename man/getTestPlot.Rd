% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTestPlot.R
\name{getTestPlot}
\alias{getTestPlot}
\title{Plot output from the VCCP model}
\usage{
getTestPlot(vccp_result)
}
\arguments{
\item{vccp_result}{A list generated from \code{\link{vccp.fun}}.}
}
\value{
No return value, called for a plotting purpose.
}
\description{
This function plots the change points in the network structure between multivariate time series detected by the VCCP model.
}
\examples{
\donttest{
## Simulate MVN data with 2 change points
data = cbind(1:180, mvn.sim.2.cps(180,8,seed=101))
## Change point detection using VCCP (it may take several minutes to complete...)
result = vccp.fun(data, "NBS", test = "V")
## Plot the result
getTestPlot(result)

result.2 = vccp.fun(data, "NBS", test = "B")
## Plot the result
getTestPlot(result.2)
}
}
\seealso{
\code{\link{vccp.fun}}
}
