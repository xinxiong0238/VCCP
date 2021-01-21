Genind <- function(ind) {
  ind_series <- c()
  for (i in 1:dim(ind)[1]) {
    ind_series <- c(ind_series, ind[i, 1]:ind[i, 2])
  }
  return(ind_series)
}

MultiGenpSeudoSample <- function(p, T, X, t_start, t_end) {
  tSam <- stats::rgeom(T, p) + 1
  subnum <- dim(X)[1] / T
  GeoSam <- tSam[1:which.min(cumsum(tSam) <= (t_end - t_start))]
  I_series <- sample(t_end:t_start, length(GeoSam), replace = T)
  ind <- cbind(I_series, I_series + GeoSam - 1)
  Y <- c()
  Z <- rbind(X, X)
  for (i in 1:subnum) {
    Z1 <- Z[Genind(ind) + (i - 1) * T, ]
    Y <- rbind(Y, Z1[1:(t_end - t_start), ])
  }
  return(Y)
}

MultiSeudoInd <- function(seudo, t, start, end, subnum) {
  sam <- c()
  for (i in 1:subnum) {
    sam <- rbind(sam, seudo[(start + (i - 1) * t):(end - 1 + (i - 1) * t), ])
  }
  return(sam[, -1])
}


#' @importFrom VineCopula RVineStructureSelect
#' @importFrom CDVine CDVineCopSelect
#' @importFrom CDVine CDVineBIC
Multi_CDR_NewTestPoint <- function(t_point, t_start, t_end, X, delta, CDR, trunc_tree, family_set, p, N, sig_alpha) {
  T <- length(unique(X[, 1]))
  subnum <- dim(X)[1] / T
  seudo_BIC <- dis_cut_BIC <- c()
  family_set <- unlist(family_set)
  if (CDR == "R") {
    BIC0 <- RVineStructureSelect(MultiInd(X, t_start, t_end),
      familyset = family_set,
      selectioncrit = "BIC", method = "mle",
      indeptest = TRUE, trunclevel = trunc_tree
    )$BIC
    test_BIC <- BIC0 - (RVineStructureSelect(MultiInd(X, t_start, t_point),
      familyset = family_set,
      selectioncrit = "BIC", method = "mle",
      indeptest = TRUE, trunclevel = trunc_tree
    )$BIC +
      RVineStructureSelect(MultiInd(X, t_point, t_end),
        familyset = family_set,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )$BIC)
    pb <- utils::txtProgressBar(style = 3)
    for (i in 1:N) {
      utils::setTxtProgressBar(pb, i / N)
      seudo <- MultiGenpSeudoSample(p, T, X, t_start, t_end)
      t <- t_end - t_start
      seudo_BIC[i] <- RVineStructureSelect(MultiSeudoInd(seudo, t, 1, t_point - t_start + 1, subnum),
        familyset = family_set,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )$BIC +
        RVineStructureSelect(MultiSeudoInd(seudo, t, (t_point - t_start + 1), (t_end - t_start + 1), subnum),
          familyset = family_set,
          selectioncrit = "BIC", method = "mle",
          indeptest = TRUE, trunclevel = trunc_tree
        )$BIC
      seudo_BIC0 <- RVineStructureSelect(MultiSeudoInd(seudo, t, 1, t_end - t_start + 1, subnum),
        familyset = family_set,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )$BIC
      dis_cut_BIC[i] <- seudo_BIC0 - seudo_BIC[i]
    }
    close(pb)
  } else {
    if (CDR == "C") {
      BIC0 <- RVineStructureSelect(MultiInd(X, t_start, t_end),
        familyset = family_set, type = 1,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )$BIC
      test_BIC <- BIC0 - (RVineStructureSelect(MultiInd(X, t_start, t_point),
        familyset = family_set, type = 1,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )$BIC +
        RVineStructureSelect(MultiInd(X, t_point, t_end),
          familyset = family_set, type = 1,
          selectioncrit = "BIC", method = "mle",
          indeptest = TRUE, trunclevel = trunc_tree
        )$BIC)
      pb <- utils::txtProgressBar(style = 3)
      for (i in 1:N) {
        utils::setTxtProgressBar(pb, i / N)
        seudo <- MultiGenpSeudoSample(p, T, X, t_start, t_end)
        t <- t_end - t_start
        seudo_BIC[i] <- RVineStructureSelect(MultiSeudoInd(seudo, t, 1, t_point - t_start + 1, subnum),
          familyset = family_set, type = 1,
          selectioncrit = "BIC", method = "mle",
          indeptest = TRUE, trunclevel = trunc_tree
        )$BIC +
          RVineStructureSelect(MultiSeudoInd(seudo, t, (t_point - t_start + 1), (t_end - t_start + 1), subnum),
            familyset = family_set, type = 1,
            selectioncrit = "BIC", method = "mle",
            indeptest = TRUE, trunclevel = trunc_tree
          )$BIC
        seudo_BIC0 <- RVineStructureSelect(MultiSeudoInd(seudo, t, 1, t_end - t_start + 1, subnum),
          familyset = family_set, type = 1,
          selectioncrit = "BIC", method = "mle",
          indeptest = TRUE, trunclevel = trunc_tree
        )$BIC
        dis_cut_BIC[i] <- seudo_BIC0 - seudo_BIC[i]
      }
      close(pb)
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
      test_BIC <- CDVineBIC(MultiInd(X, t_start, t_end),
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
      pb <- utils::txtProgressBar(style = 3)
      for (i in 1:N) {
        utils::setTxtProgressBar(pb, i / N)
        seudo <- MultiGenpSeudoSample(p, T, X, t_start, t_end)
        t <- t_end - t_start
        seudo_D1 <- CDVineCopSelect(MultiSeudoInd(seudo, t, 1, t_point - t_start + 1, subnum),
          type = 2, familyset = unlist(family_set),
          selectioncrit = "BIC", indeptest = TRUE
        )
        seudo_D2 <- CDVineCopSelect(MultiSeudoInd(seudo, t, (t_point - t_start + 1), (t_end - t_start + 1), subnum),
          type = 2, familyset = unlist(family_set),
          selectioncrit = "BIC", indeptest = TRUE
        )
        seudo_D0 <- CDVineCopSelect(MultiSeudoInd(seudo, t, 1, t_end - t_start + 1, subnum),
          type = 2, familyset = unlist(family_set),
          selectioncrit = "BIC", indeptest = TRUE
        )
        dis_cut_BIC[i] <- CDVineBIC(MultiSeudoInd(seudo, t, 1, t_end - t_start + 1, subnum),
          type = 2,
          family = seudo_D0$family, par = seudo_D0$par, par2 = seudo_D0$par2
        )$BIC -
          (
            CDVineBIC(MultiSeudoInd(seudo, t, 1, t_point - t_start + 1, subnum),
              type = 2,
              family = seudo_D1$family, par = seudo_D1$par, par2 = seudo_D1$par2
            )$BIC +
              CDVineBIC(MultiSeudoInd(seudo, t, (t_point - t_start + 1), (t_end - t_start + 1), subnum),
                type = 2,
                family = seudo_D2$family, par = seudo_D2$par, par2 = seudo_D2$par2
              )$BIC
          )
      }
      close(pb)
    }
  }

  return(c(test_BIC, stats::quantile(dis_cut_BIC, sig_alpha), stats::quantile(dis_cut_BIC, 1 - sig_alpha)))
}


#' Test change points in connection network by Stationary Bootstrap
#'
#' \code{TestPoints_Boot} returns a dataframe containing the inference result
#' of significant change points in the connection network among the multi/single-subject
#' time series data. It can be used after detecting possible change points
#'  by \code{\link{VC_NBS_FindPoints}}.
#'
#' If you get the change point candidates from \code{\link{VC_NBS_FindPoints}},
#' parameters used in \code{TestPoints_Boot} are strongly recommended to be consistent
#' with the ones specified when detecting possible candidates (i.e., \code{delta},
#' \code{CDR}, \code{trunc_tree}, \code{family_set}, \code{pre_white} and \code{ar_num}).
#' You can also combine the detection (NBS) and inference process (Stationary Bootstrap or Vuong Test)
#' together by \code{\link{VC_NBS}}.
#'
#' @param v_t_point A integer vector of change point candidates, normally
#' from the first element of the list generated by \code{\link{VC_NBS_FindPoints}}.
#'
#' @param X_raw A matrix with at least three columns. The first one must
#'  be integers indicating time. If multiple subjects are included,
#'  you shoud stack their data vertically with the corresponding
#'  timestamps shown in the first column. The other two or more columns are
#'  the multi-dimensional time series data (n*T times p).
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
#' @param p A decimal between 0 and 1; controlling the sampling block size in
#'  Stationary Boostrap method (\code{rgeom(T,p)}). The larger p (default=0.3) is, the fewer time
#'  points each bootstrap sample contains.
#'
#' @param N An integer; specifying the number of Stationary Bootstrap samples (default=100).
#'
#' @param pre_white,ar_num Integers helping fit (method: yule-walker) an
#' autoregressive time series model to preprocess the raw data. If \code{pre_white}=0(default),
#' no ar model is fitted. If \code{pre_white}=1, ar_num is the maximum order of model to fit.
#'
#' @param sig_alpha A decimal between 0 and 1; significance level of the Stationary Bootsrtap test.
#'
#' @return A dataframe. The first column contains possible change point candiates;
#'  the second one corresponds to the reduced BIC values (left VC + right VC - all_VC);
#'  the third and the fourth columns are the lower and upper bound of reduced BIC values calculated by the
#'  Stationary Bootsrtap test; and the fifth one is the inference result.
#' @export
#' @examples
#' data = cbind(1:180, random.mvn.simulate.2.changes(180,8,seed=101))
#' result = VC_NBS_FindPoints(data, 30)
#' inference = TestPoints_Boot(result[[1]], data, delta=30)
#' @seealso  \code{\link{VC_NBS_FindPoints}}, \code{\link{TestPoints_Vuong}}, \code{\link{VC_NBS}}
#'
#' @importFrom VineCopula RVineStructureSelect
#' @importFrom CDVine CDVineCopSelect
#' @importFrom CDVine CDVineBIC
TestPoints_Boot <- function(v_t_point, X_raw, delta, CDR = "D", trunc_tree = NA, family_set = 1,
                            pre_white = 0, ar_num = 1, p = 0.3, N = 100, sig_alpha = 0.05) {
  T <- length(unique(X_raw[, 1]))
  subnum <- dim(X_raw)[1] / T
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

  test_result <- as.data.frame(matrix(0, length(v_t_point), 5))
  if (length(v_t_point) == 0) {
    return(test_result)
  } else {
    test_result[, 1] <- v_t_point
    a <- c(1, v_t_point, T + 1)
    full_point <- a[!duplicated(a)]
    for (i in 1:length(v_t_point)) {
      cat(paste("Test for candidate", i, ": t =", v_t_point[i]), fill = TRUE)
      test_result[i, 2:4] <- Multi_CDR_NewTestPoint(
        full_point[i + 1], full_point[i],
        full_point[i + 2], X, delta, CDR, trunc_tree, family_set, p, N, sig_alpha
      )
      if (test_result[i, 2] > test_result[i, 4]) {
        test_result[i, 5] <- "significant"
        cat("Significant! \n\n")
      } else {
        test_result[i, 5] <- "not significant"
        cat("Insignificant! \n\n")
      }
    }
    names(test_result) <- c("t", "reduced BIC", "Lower_CI", "Upper_CI", "judgement")
    return(test_result)
  }
}


#' @export
Boot_GetTestPlot <- function(test_result, T) {
  if (dim(test_result)[1] == 0) {
    plot(1:T,1:T,type = "n",yaxt="n",ylab = NA)
    graphics::text("No candidate is found.")
  }
  else {
    reduced_BIC <- lwr <- upr <- rep(0, T)
    reduced_BIC[test_result[, 1]] <- test_result[, 2]
    plot(1:T, reduced_BIC,
      type = "h",
      ylim = c(min(test_result[, 2:4]) - 3, max(test_result[, 2:4] + 10, 0)), lwd = 2
    )
    graphics::points(test_result[, 1], test_result[, 4], pch = 20, col = "red")
    graphics::points(test_result[, 1], test_result[, 3], pch = 20, col = "red")
    if (sum(test_result[, 5] == "significant") != 0) {
      graphics::text(test_result[test_result[, 5] == "significant", 1], test_result[test_result[, 5] == "significant", 2] + 5,
        labels = test_result[test_result[, 5] == "significant", 1], cex = 0.8
      )
    }
  }
}


