#' Multiple change point detection in the vine copula structure of multivariate time series
#'
#' This function detects multiple change points in the vine
#' copula structure of a multivariate time series using
#' vine copulas, various state-of-the-art segmentation methods to identify
#' multiple change points, and a likelihood ratio test or the stationary bootstrap
#' for inference.
#'
#' The time series \code{X_t} is of dimensionality p and we are
#' looking for changes in the vine copula structure between
#' the different time series components \code{X_{t}^{(1)}, X_{t}^{(2)},
#'  ..., X_{t}^{(p)}}. VCCP uses vine copulas, various state-of-the-art
#'  segmentation methods to identify multiple change points,
#'  and a likelihood ratio test or the stationary bootstrap for inference.
#'
#' @param X A numerical matrix representing the multivariate
#' time series, with the columns representing its components.
#' If multiple subjects are included (panel data), vertically
#' stack the subject data and identify timestamps of each subject in the first column.
#'
#' @param method A character string, which defines the
#'  segmentation method. If \code{method} = "NBS", which is the
#'  default method, then the adapted binary segmentation is used.
#'  Similarly, if \code{method} = "OBS", "MOSUM" or "WBS", then binary
#'  segmentation, MOSUM and wild binary segmentation are used, respectively.
#'
#' @param delta A positive integer number with default value equal to 30.
#'  It is used to define the minimum distance acceptable between
#'  change points. In general, \code{delta} >= 5*ncol(X))
#'  is recommended to ensure sufficient data when estimating the
#'  vine copula model.
#'
#' @param G A positive real number between 0 and 1 with default value equal to 0.1.
#'  It is used to define the moving sum bandwidth relative to \code{T} in MOSUM when
#'  \code{method} = "MOSUM" is chosen. Alternatively, a positive integer
#'  less than half of the time series length can be set to define the absolute bandwidth.
#'
#' @param M A positive integer with default value equal to floor(9*log(T)) (T is the length of the time series).
#'  It represents the number of sub-samples in WBS when
#'  \code{method} = "WBS" is chosen.
#'
#' @param test A character string, which defines the inference
#'  method used. If \code{test} = "V", which is the default method,
#'  the Vuong test is performed. If \code{test} = "B", the
#'  stationary bootstrap is performed.
#'
#' @param CDR A character string, which defines the vine structure.
#'  If \code{CDR} = "D", which is the default method,
#'  a D-vine is used. Similarly, if \code{CDR} = "C" or \code{CDR}
#'  = "R", a C-vine or an R-vine is used, respectively.
#'
#' @param trunc_tree A positive integer, which defines the level
#'  of truncation for the vine copula. If \code{trunc_tree} = "NA",
#'  which is the default value, the Vine contains \code{dim(X)[2]-2}
#'  levels of trees.
#'
#' @param family_set A positive integer, which defines the bivariate copula
#'  family. If \code{familyset} = 1, which is the default value, only the
#'  Gauss copula is selected and VCCP detects change points in
#'  the linear correlation graph. Coding of pair-copula
#'  families is the same as in \code{\link[VineCopula]{BiCop}}.
#'
#' @param pre_white A positive integer, which defines whether
#'  the data is pre-whitened. If \code{pre-white} = 0, which is the
#'  default value, no pre-whitening is performed. If
#'  \code{pre_white} = 1, an autoregressive time series model
#'  (method: yule-walker) is used to preprocess the raw data.
#'
#' @param ar_num A positive integer, which defines the maximum
#'  order of model to fit to preprocess the data (see \code{pre_white}).
#'  If \code{ar_num} = 1, which is the default value, then an AR(1)
#'  model is fit to the data.
#'
#' @param p A positive real number between 0 and 1 which is
#'  defined as the block size in the stationary bootstrap
#'  method (\code{rgeom(T,p)}) if \code{test} = "B" is chosen.
#'  If \code{p} = 0.3, which is the default value, each resampled block
#'  has 1/0.3 time points on average.
#'
#' @param N A positive integer, which defines the number
#'  of the stationary bootstrap resamples used. The default value is \code{N} = 100.
#'
#' @param sig_alpha A positive real number between 0 and 1, which
#'  defines the significance level of the inference test.
#'  The default values is 0.05.
#'
#' @return A list with the following components:
#'
#' \tabular{ll}{
#'  \code{loc_of_cpts} \tab The locations of the detected change points. \cr
#'  \code{no_of_cpts} \tab The number of detected change points. \cr
#'  \code{test_df} \tab A dataframe containing the test result.  \cr
#'  \code{compute_time} \tab Time (in minutes) to run \code{vccp.fun}. \cr
#'  \code{T} \tab The length of the time series data. \cr
#'  \code{sig_alpha} \tab The significance level for the inference test. \cr
#' }
#'
#'
#' @export
#' @examples
#' \donttest{
#' ## Simulate MVN data with 2 change points
#' data <- cbind(1:180, mvn.sim.2.cps(180, 8, seed = 101))
#' T <- 180
#' ## Change point detection using VCCP (it may take several minutes to complete...)
#' result.NV <- vccp.fun(data, method = "NBS", delta = 30, test = "V")
#' ## Plot the results
#' getTestPlot(result.NV)
#' #title("VCCP: NBS + Vuong")
#'
#' ## Change point detection using NBS and stationary bootstrap for inference
#' result.NB <- vccp.fun(data, method = "NBS", delta = 30, test = "B")
#' ## Plot the results
#' getTestPlot(result.NB)
#' title("VCCP: NBS + Stationary Bootstrap")
#' }
#' @seealso  \code{\link{getTestPlot}}
#' @section Author(s):
#'  Xin Xiong, Ivor Cribben (\email{cribben@@ualberta.ca})
#' @section References:
#'  "Beyond linear dynamic functional connectivity: a vine copula change point model", Xiong and Cribben (2021), bioRxiv 2021.04.25.441254.
vccp.fun <- function(X, method = 'NBS', delta = 30, G = 0.1, M = NA, test = "V", CDR = "D", trunc_tree = NA,
                   family_set = 1, pre_white = 0, ar_num = 1,
                   p = 0.3, N = 100, sig_alpha = 0.05) {
  if (method != "NBS" & method != "OBS" & method != "MOSUM" & method != "WBS"){
    stop("You can only specify method as 'NBS', 'OBS', 'MOSUM' or 'WBS'!")
  }else{
    if (CDR != "C" & CDR != "D" & CDR != "R") {
      stop("You can only specify CDR as 'C', 'D' or 'R'!")
    } else {
      if (test != "V" & test != 'B') {
       stop("You can only specify test as 'B' or 'V'!")
      } else {
        re.list = list()
        t = proc.time()
        result = switch (method,
          "NBS" = VC_NBS(X, delta, test, CDR, trunc_tree,
                         family_set, pre_white, ar_num, p, N, sig_alpha),
          "OBS" = VC_OBS(X, delta, test, CDR, trunc_tree,
                         family_set, pre_white, ar_num, p, N, sig_alpha),
          "MOSUM" = VC_MOSUM(X, delta, G, test, CDR, trunc_tree,
                         family_set, pre_white, ar_num, p, N, sig_alpha),
          "WBS" = VC_WBS(X, delta, M, test, CDR, trunc_tree,
                         family_set, pre_white, ar_num, p, N, sig_alpha))
        compute = (proc.time() - t)[3]/60
        re.list = list("loc_of_cpts"=result$t,
                       "no_of_cpts" = length(result$t),
                       "test_df" = result,
                       "compute_time" = compute,
                       "T"=length(unique(X[, 1])),
                       "sig_alpha"=sig_alpha)
          return(re.list)
      }
    }
  }
}
