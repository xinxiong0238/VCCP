#' Find and test change points in connection network by VCCP (MUSUM) model
#'
#' \code{VC_MOSUM} returns a dataframe containing the inference result
#' of significant change points in the connection network among the multi/single-subject
#' time series data. The inference method can be selected from the Stationary
#' Bootstrap and the Vuong test.
#'
#' @param X_raw A matrix with at least three columns. The first one must
#'  be integers indicating time. If multiple subjects are included,
#'  you shoud stack their data vertically with the corresponding
#'  timestamps shown in the first column. The other two or more columns are
#'  the multi-dimensional time series data (n*T times p).
#'
#' @param delta An integer indicating the least distance between two
#'  change points. This parameter is vital for the VC_MOSUM method
#'  and must be specified. Normally \code{delta >= 5*(dim(X_raw)[2]-1)} is recommended to ensure sufficient
#'  data when fitting the VC model.
#'
#' @param G An integer for the moving sum bandwidth in the MOSUM model.
#'  G should be less than \code{T/2}. Alternatively, a number between
#'  0 and 0.5 describing the moving sum bandwidth relative to \code{T}
#'  can be given. Default: \code{G=0.1}. This parameter greatly influences
#'  the performance of the VC_MOSUM model.
#'
#' @param test An upper letter specifying the inference test type of the VCCP model:
#'  "\code{V}" = Vuong test (default),
#'  "\code{B}" = Stationary Bootstrap.
#'
#' @param CDR An upper letter specifying the tree type of the vine model:
#'  "\code{D}" = D-Vine (default),
#'  "\code{C}" = C-Vine,
#'  "\code{R}" = R-Vine.
#'
#' @param trunc_tree An integer or NA; level of truncation for the vine model.
#'  If NA is specified (default), the vine contains \code{dim(X_raw)[2]-2} levels of trees.
#'
#' @param family_set An integer vector of pair-copula families to select from.
#'  The vector has to include at least one pair-copula family that allows for
#'  positive and one that allows for negative dependence.
#'  If \code{familyset} = 1 (default), only Gauss copula is selected and VCCP can only
#'  detect changes in the linear correlation network.
#'  Coding of pair-copula families is the same as in \code{\link[VineCopula]{BiCop}}.

#' @param p A decimal between 0 and 1; controlling the sampling block size in
#'  Stationary Boostrap method (\code{rgeom(T,p)}) if you choose \code{test=="B"}. The larger p (default=0.3) is, the fewer time
#'  points each bootstrap sample contains.
#'
#' @param N An integer; specifying the number of Stationary Bootstrap samples (default=100) if you choose \code{test=="B"}.
#'
#' @param pre_white,ar_num Integers helping fit (method: yule-walker) an
#' autoregressive time series model to preprocess the raw data. If \code{pre_white}=0(default),
#' no ar model is fitted. If \code{pre_white}=1, ar_num is the maximum order of model to fit.
#'
#' @param sig_alpha A decimal between 0 and 1; significance level of the inference test.
#'
#' @return A dataframe. If you choose \code{test=="B"}, the dataframe contains 5 columns.
#'  The first column contains possible change point candiates;
#'  the second one corresponds to the reduced BIC values (left VC + right VC - all_VC);
#'  the third and the fourth columns are the lower and upper bound of reduced BIC values calculated by the
#'  Stationary Bootsrtap test; and the fifth one is the inference result.
#'
#'  If you choose \code{test=="V"}, the dataframe contains 4 columns.
#'  The first column contains possible change point candiates;
#'  the second and the third one correspond to the P-values of the left and right Vuong tests
#'  with Schwarz correction; and the fourth one is the inference result.
#'
#' @export
#' @examples
#' data <- cbind(1:180, random.mvn.simulate.2.changes(180, 8, seed = 101))
#' result <- VC_MOSUM(data, 30, test = "V")
#' @seealso  \code{\link{VCCP}}
VC_MOSUM <- function(X_raw, delta, G = 0.1, test = "V", CDR = "D", trunc_tree = NA,
                     family_set = 1, pre_white = 0, ar_num = 1,
                     p = 0.3, N = 100, sig_alpha = 0.05) {
  if (CDR != "C" & CDR != "D" & CDR != "R") {
    return(cat("You can only specify CDR as 'C', 'D' or 'R'!"))
  } else {
    if (test == "V") {
      infer <- VC_MOSUM_Vuong(
        X_raw, delta, G, CDR, trunc_tree,
        family_set, pre_white, ar_num,
        sig_alpha
      )
      return(infer)
    } else {
      if (test == "B") {
        infer <- VC_MOSUM_Boot(
          X_raw, delta, G, CDR,
          trunc_tree, family_set,
          pre_white, ar_num, p, N, sig_alpha
        )
        return(infer)
      } else {
        return(cat("You can only specify test as 'B' or 'V'!"))
      }
    }
  }
}
