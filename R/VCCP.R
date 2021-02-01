#' Multiple change point detection in the vine copula structure of multivariate time series
#'
#' \code{vccp.fun} detects multiple change points in the vine
#' copula structure of a multivariate time series using
#' vine copulas, various segmentation methods and a
#' likelihood ratio or Stationary Boostrap test for inference.
#'
#' The time series $X_t$ is of dimensionality p and we are
#' looking for changes in the vine copula structure between
#' the different time series components $X_{t}^{(1)}$, $X_{t}^{(2)}$,
#'  ..., $X_{t}^{(p)}$. We use vine copulas, various segmentation
#'  methods and a likelihood ratio test for inference.
#'  The function has been extensively tested on fMRI data.
#'
#' @param X A numerical matrix representing the multivariate
#' time series, with the columns representing its components
#' and rows representing time. If multiple subjects are included,
#' vertically stack the subject data and identify timestamps
#' of each subject in the first column.
#'
#' @param method A character string, which defines the
#'  segmentation method. If \code{method} = "NBS", which is the
#'  default method, then the new binary segmentation is used.
#'  Similarly, "OBS", "MOSUM" and "WBS" represent binary
#'  segmentation, MOSUM and wild binary segmentation,
#'  respectively.
#'
#' @param delta A positive integer number with default value equal to 30.
#'  It is used to define the minimum distance acceptable between
#'  detected change points. Normally \code{delta} >= 5*ncol(X))
#'  is recommended to ensure sufficient data when estimating the
#'  vine copula model.
#'
#' @param G A positive real number between 0 and 1 with default value equal to 0.1.
#'  It is used to define the moving sum bandwidth relative to \code{T} in MOSUM when method = ”MOSUM” is chosen.
#'  Alternatively, an positive integer less than \code{T/2} can be set to define the absolute bandwith.
#'
#' @param M A positive integer with default value floor(9*log(T).
#'  It represents the number of sub-samples in WBS when
#'  \code{method}="WBS" is chosen.
#'
#' @param test A character string, which defines the inference
#'  method. If \code{test} = "V", which is the default method,
#'  the Vuong test is performed. If \code{test} = "B", the
#'  Stationary Bootstrap is performed.
#'
#' @param CDR A character string, which defines the type of Vine
#'  Copula. If \code{CDR} = "D", which is the default method,
#'  a D-vine is used. Similarly, if \code{CDR} = "C" or \code{CDR}
#'  = "R", C-vine or R-vine copula is used.
#'
#' @param trunc_tree A positive integer, which defines the level
#'  of truncation for the vine copula. If \code{trunc_tree} = "NA",
#'  which is the default value, the Vine contains \code{dim(X)[2]-2}
#'  levels of trees.
#'
#' @param family_set A positive integer, which defines the bivariate copula
#'  family. If \code{familyset} = 1, which is the default value, only the
#'  Gauss copula is selected and VCCP detects change points in
#'  the linear correlation network. Coding of pair-copula
#'  families is the same as in \code{\link[VineCopula]{BiCop}}.
#'
#' @param pre_white A positive integer, which defines whether
#'  to pre-whiten the data. If \code{pre-white} = 0, which is the
#'  default value, no pre-whitening is performed. If
#'  \code{pre_white} =1, an autoregressive time series model
#'  (method: yule-walker) is used to preprocess the raw data.
#'
#' @param ar_num A positive integer, which defines the maximum
#'  order of model to fit to preprocess the data (see \code{pre_white}).
#'  If \code{ar_num} = 1, which is the default values, then an AR(1)
#'  model is fit to the data.
#'
#' @param p A positive real number between 0 and 1 which is
#'  defined to control the block size in Stationary Boostrap
#'  method (\code{rgeom(T,p)}) if \code{test} = "B" is chosen.
#'  If \code{p}=0.3, which is the default, each block resample
#'  1/0.3 time points.
#'
#' @param N A positive integer, which defines the number
#'  of Stationary Bootstrap resamples. The default value is \code{N=100}.
#'
#' @param sig_alpha positive real number between 0 and 1, which
#'  defines the significance level of the inference test.
#'  The default values is 0.05.
#'
#' @return A list with the following components:
#'
#' \tabular{ll}{
#'  \code{change.points} \tab The locations of the detected change points. \cr
#'  \code{no.of.cpts} \tab The number of the detected change points. \cr
#'  \code{test.df} \tab A dataframe containing the test result. If \code{test="B"}, the dataframe contains 5 columns.
#'  The first column contains possible change point candidates;
#'  the second one corresponds to the reduced BIC values (left VC + right VC - all VC);
#'  the third and the fourth columns are the lower and upper bound of reduced BIC values calculated by the
#'  Stationary Bootsrtap test; and the fifth one is the inference result. If \code{test="V"}, the dataframe contains 4 columns.
#'  The first column contains possible change point candidates;
#'  the second and the third one correspond to the P-values of the left and right Vuong tests
#'  with Schwarz correction; and the fourth one is the inference result.
#'  \cr
#'  \code{compute.time} \tab Time (in minutes), to run \code{vccp.fun}. \cr
#'  \code{T} \tab The total length of the time series data. \cr
#'  \code{sig_alpha} \tab The significance level of the inference test. \cr
#' }
#'
#'
#' @export
#' @examples
#' data <- cbind(1:180, random.mvn.simulate.2.changes(180, 8, seed = 101))
#' T <- 180
#'
#'
#' result.NV <- vccp.fun(data, method = "NBS", delta = 30, test = "V")
#' getTestPlot(result.NV)
#' title("VCCP: NBS + Vuong")
#'
#'
#' result.NB <- vccp.fun(data, method = "NBS", delta = 30, test = "B")
#' getTestPlot(result.NB)
#' title("VCCP: NBS + Stationary Bootstrap")
#'
#' @seealso  \code{\link{getTestPlot}}
#' @section Author(s):
#'  Xin Xiong, Ivor Cribben (\email{cribben@@ualberta.ca})
#' @section References:
#'  "Beyond linear dynamic functional connectivity: a vine copula change point model", Xiong and Cribben (2021), preprint.
vccp.fun <- function(X, method = 'NBS', delta = 30, G = 0.1, M = NA, test = "V", CDR = "D", trunc_tree = NA,
                   family_set = 1, pre_white = 0, ar_num = 1,
                   p = 0.3, N = 100, sig_alpha = 0.05) {
  if (method != "NBS" & method != "OBS" & method != "MOSUM" & method != "WBS"){
    return(cat("You can only specify method as 'NBS', 'OBS', 'MOSUM' or 'WBS'!"))
  }else{
    if (CDR != "C" & CDR != "D" & CDR != "R") {
      return(cat("You can only specify CDR as 'C', 'D' or 'R'!"))
    } else {
      if (test != "V" & test != 'B') {
        return(cat("You can only specify test as 'B' or 'V'!"))
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
        re.list = list("change.points"=result$t,
                       "no.of.cpts" = length(result$t),
                       "test.df" = result,
                       "compute.time" = compute,
                       "T"=length(unique(X[, 1])),
                       "sig_alpha"=sig_alpha)
          return(re.list)
      }
    }
  }
}
