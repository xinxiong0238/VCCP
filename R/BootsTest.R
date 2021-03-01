Genind <- function(ind) {
  ind_series <- c()
  for (i in 1:dim(ind)[1]) {
    ind_series <- c(ind_series, ind[i, 1]:ind[i, 2])
  }
  return(ind_series)
}
#' @importFrom VineCopula RVineStructureSelect
#' @importFrom VineCopula D2RVine
#' @importFrom VineCopula RVineCopSelect
#' @importFrom VineCopula RVineVuongTest
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


Multi_CDR_NewTestPoint <- function(t_point, t_start, t_end, X, delta, CDR, trunc_tree, family_set, p, N, sig_alpha) {
  T <- length(unique(X[, 1]))
  subnum <- dim(X)[1] / T
  p_X = dim(X)[2] - 1
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
    #pb <- utils::txtProgressBar(style = 3)
    for (i in 1:N) {
      #utils::setTxtProgressBar(pb, i / N)
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
    #close(pb)
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
      #pb <- utils::txtProgressBar(style = 3)
      for (i in 1:N) {
        #utils::setTxtProgressBar(pb, i / N)
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
      #close(pb)
    } else {
      D_X = D2RVine(1:p_X, family = rep(0, p_X*(p_X-1)/2),
                    par = rep(0, p_X*(p_X-1)/2),
                    par2 = rep(0, p_X*(p_X-1)/2))

      test_BIC <- RVineCopSelect(MultiInd(X, t_start, t_end), selectioncrit = "BIC",
                                  method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                                  indeptest = TRUE, trunclevel = trunc_tree)$BIC -
        (
          RVineCopSelect(MultiInd(X, t_start, t_point), selectioncrit = "BIC",
                         method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                         indeptest = TRUE, trunclevel = trunc_tree)$BIC +
            RVineCopSelect(MultiInd(X, t_point, t_end), selectioncrit = "BIC",
                           method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                           indeptest = TRUE, trunclevel = trunc_tree)$BIC
        )

      for (i in 1:N) {
        seudo <- MultiGenpSeudoSample(p, T, X, t_start, t_end)
        t <- t_end - t_start
        dis_cut_BIC[i] <- RVineCopSelect(MultiSeudoInd(seudo, t, 1, t_end - t_start + 1, subnum),
                                         selectioncrit = "BIC",method = "mle",
                                         Matrix = D_X$Matrix,familyset = unlist(family_set),
                                         indeptest = TRUE, trunclevel = trunc_tree)$BIC -
          (
            RVineCopSelect(MultiSeudoInd(seudo, t, (t_point - t_start + 1), (t_end - t_start + 1), subnum),
                           selectioncrit = "BIC",method = "mle",
                           Matrix = D_X$Matrix,familyset = unlist(family_set),
                           indeptest = TRUE, trunclevel = trunc_tree)$BIC +
              RVineCopSelect(MultiSeudoInd(seudo, t, (t_point - t_start + 1), (t_end - t_start + 1), subnum),
                             selectioncrit = "BIC",method = "mle",
                             Matrix = D_X$Matrix,familyset = unlist(family_set),
                             indeptest = TRUE, trunclevel = trunc_tree)$BIC
          )
      }
    }
  }

  return(c(test_BIC, stats::quantile(dis_cut_BIC, sig_alpha), stats::quantile(dis_cut_BIC, 1 - sig_alpha)))
}




TestPoints.Boot <- function(v_t_point, X_raw, delta, CDR = "D", trunc_tree = NA, family_set = 1,
                            pre_white = 0, ar_num = 1, p = 0.3, N = 100, sig_alpha = 0.05) {
  T <- length(unique(X_raw[, 1]))
  subnum <- dim(X_raw)[1] / T
  X <- X_raw
  v_t_point = sort(v_t_point)
  if(sum(c(v_t_point,T)-c(0,v_t_point)<dim(X_raw)[2]-1)>0) stop("Some candidates are too close to each other or to the boundary (distance < p)")
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
    message("Perform stationary bootstrap test on candidates...")
    for (i in 1:length(v_t_point)) {
      test_result[i, 2:4] <- Multi_CDR_NewTestPoint(
        full_point[i + 1], full_point[i],
        full_point[i + 2], X, delta, CDR, trunc_tree, family_set, p, N, sig_alpha
      )
      if (test_result[i, 2] > test_result[i, 4]) {
        test_result[i, 5] <- "significant"
      } else {
        test_result[i, 5] <- "not significant"
      }
    }
    names(test_result) <- c("t", "reduced BIC", "Lower_CI", "Upper_CI", "judgement")
    return(test_result)
  }
}



GetTestPlot.Boot <- function(test_result, T) {
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


