MultiInd <- function(X, start, end) {
  T <- length(unique(X[, 1]))
  subnum <- dim(X)[1] / T
  sam <- c()
  for (i in 1:subnum) {
    sam <- rbind(sam, X[(start + (i - 1) * T):(end - 1 + (i - 1) * T), ])
  }
  return(sam[, -1])
}

#' Find changes in the correlation network by the VC+NBS method
#'
#' \code{VC_NBS_FindPoints} returns possible change points in the correlation network
#'  among the time series data.
#'
#'  The time series data can come from one subject, or from multiple subjects.
#'  The VC model are fitted using \code{\link[VineCopula]{RVineStructureSelect}} or
#'  \code{\link[CDVine]{CDVineCopSelect}}.
#'
#'
#' @param X_raw A matrix with at least three columns. The first one must
#'  be integers indicating time. If multiple subjects are included,
#'  you shoud stack their data vertically with the corresponding
#'  timestamps shown in the first column. The other two or more columns are
#'  the multi-dimensional time series data.
#'
#' @param delta An integer indicating the least distance between two
#'  change points. This parameter is vital for the New Binary Segmentation (NBS) method
#'  and must be specified. Normally \code{delta >= 5*(dim(X_raw)[2]-1)} is recommended to ensure sufficient
#'  data when fitting the VC model.
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
#'
#' @param pre_white,ar_num Integers helping fit (method: yule-walker) an
#' autoregressive time series model to preprocess the raw data. If \code{pre_white}=0(default),
#' no ar model is fitted. If \code{pre_white}=1, ar_num is the maximum order of model to fit.
#'
#' @return A two-element list. The first one contains a integer vector of
#'  detected possible change points. The second one is the preprocessed time
#'  series data (transforming data of each dimension to be univariately uniformly distributed
#'  by \code{\link[VineCopula]{pobs}})
#' @export
#' @examples
#' data = cbind(1:180, random.mvn.simulate.2.changes(180,8,seed=101))
#' result = VC_NBS_FindPoints(data, 30)
#' @seealso  \code{\link{TestPoints_Boot}}, \code{\link{TestPoints_Vuong}}
#'
#' @importFrom VineCopula RVineStructureSelect
#' @importFrom CDVine CDVineCopSelect
#' @importFrom CDVine CDVineBIC
VC_NBS_FindPoints <- function(X_raw, delta, CDR = "D", trunc_tree = NA, family_set = 1, pre_white = 0, ar_num = 1) {
  T <- length(unique(X_raw[, 1]))
  subnum <- dim(X_raw)[1] / T
  cut_point0 <- c(1, T + 1)
  BIC_cut <- 1
  k <- 0
  X <- X_raw
  if (pre_white == 1) {
    for (i in 1:subnum) {
      armodel <- stats::ar(X_raw[((i - 1) * T + 1):(i * T), -1], FALSE, ar_num)
      ar_resid <- armodel$resid
      ar_resid[which(is.na(ar_resid) == TRUE)] <- 0
      X[((i - 1) * T + 1):(i * T), -1] <- ar_resid
    }
  }
  X[, -1] <- VineCopula::pobs(X[, -1])
  while (BIC_cut > 0) {
    k <- k + 1
    cat(paste("Binary search, round",k,"..."),fill = TRUE)
    BIC_de <- rep(0, T)
    for (i in delta:(T - delta)) {
      if (sum(is.element((i - delta + 1):(i + delta), cut_point0)) != 0) {
        next
      } else {
        t_point <- i
        t_start <- cut_point0[which.max(i <= cut_point0) - 1]
        t_end <- cut_point0[which.max(i <= cut_point0)]
        if (CDR == "R") {
          BIC_de[i] <- RVineStructureSelect(MultiInd(X, t_start, t_end),
            familyset = unlist(family_set),
            selectioncrit = "BIC", method = "mle",
            indeptest = TRUE, trunclevel = trunc_tree
          )$BIC -
            (RVineStructureSelect(MultiInd(X, t_start, t_point),
              familyset = unlist(family_set),
              selectioncrit = "BIC", method = "mle",
              indeptest = TRUE, trunclevel = trunc_tree
            )$BIC +
              RVineStructureSelect(MultiInd(X, t_point, t_end),
                familyset = unlist(family_set),
                selectioncrit = "BIC", method = "mle",
                indeptest = TRUE, trunclevel = trunc_tree
              )$BIC)
          #print(c(BIC_de[i], i))
        } else {
          if (CDR == "C") {
            BIC_de[i] <- RVineStructureSelect(MultiInd(X, t_start, t_end),
              familyset = unlist(family_set), type = 1,
              selectioncrit = "BIC", method = "mle",
              indeptest = TRUE, trunclevel = trunc_tree
            )$BIC -
              (RVineStructureSelect(MultiInd(X, t_start, t_point),
                familyset = unlist(family_set), type = 1,
                selectioncrit = "BIC", method = "mle",
                indeptest = TRUE, trunclevel = trunc_tree
              )$BIC +
                RVineStructureSelect(MultiInd(X, t_point, t_end),
                  familyset = unlist(family_set), type = 1,
                  selectioncrit = "BIC", method = "mle",
                  indeptest = TRUE, trunclevel = trunc_tree
                )$BIC)
            #print(c(BIC_de[i], i))
          } else {
            D_1 <- CDVineCopSelect(MultiInd(X, t_start, t_end),
              type = 2, familyset = unlist(family_set),
              selectioncrit = "BIC", indeptest = TRUE
            )
            D_2 <- CDVineCopSelect(MultiInd(X, t_start, t_point),
              type = 2, familyset = unlist(family_set),
              selectioncrit = "BIC", indeptest = TRUE
            )
            D_3 <- CDVineCopSelect(MultiInd(X, t_point, t_end),
              type = 2, familyset = unlist(family_set),
              selectioncrit = "BIC", indeptest = TRUE
            )
            BIC_de[i] <- CDVineBIC(MultiInd(X, t_start, t_end),
              type = 2,
              family = D_1$family, par = D_1$par, par2 = D_1$par2
            )$BIC -
              (
                CDVineBIC(MultiInd(X, t_start, t_point),
                  type = 2,
                  family = D_2$family, par = D_2$par, par2 = D_2$par2
                )$BIC +
                  CDVineBIC(MultiInd(X, t_point, t_end),
                    type = 2,
                    family = D_3$family, par = D_3$par, par2 = D_3$par2
                  )$BIC
              )
            #print(c(BIC_de[i], i))
          }
        }
      }
    }
    BIC_cut <- max(BIC_de)
    cut_new_point <- which.max(BIC_de)
    if(BIC_cut > 0){
      cat(paste("Find candidate",k,": t =",cut_new_point,"\n \n"),fill = TRUE)
    }else{
      cat(paste("No more candidate is found. \n \n"),fill = TRUE)
    }
    cut_point0 <- sort(c(cut_point0, cut_new_point))
  }
  return(list(c(cut_point0[cut_point0 != cut_new_point & cut_point0 != T + 1 & cut_point0 != 1]), X))
}

