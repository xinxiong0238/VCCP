#' @importFrom VineCopula RVineVuongTest
#' @importFrom CDVine CDVineVuongTest
Vuong_Multi_CDR_NewTestPoint <- function(t_point, t_start, t_end, X, delta, CDR, trunc_tree, family_set) {
  T <- length(unique(X[, 1]))
  subnum <- dim(X)[1] / T
  family_set <- unlist(family_set)
  if (CDR == "R") {
    RV0 <- RVineStructureSelect(MultiInd(X, t_start, t_end),
      familyset = family_set,
      selectioncrit = "BIC", method = "mle",
      indeptest = TRUE, trunclevel = trunc_tree
    )
    RV1 <- RVineStructureSelect(MultiInd(X, t_start, t_point),
      familyset = family_set,
      selectioncrit = "BIC", method = "mle",
      indeptest = TRUE, trunclevel = trunc_tree
    )
    RV2 <- RVineStructureSelect(MultiInd(X, t_point, t_end),
      familyset = family_set,
      selectioncrit = "BIC", method = "mle",
      indeptest = TRUE, trunclevel = trunc_tree
    )
    Vuong_test_left <- RVineVuongTest(MultiInd(X, t_start, t_point), RV0, RV1)
    Vuong_test_right <- RVineVuongTest(MultiInd(X, t_point, t_end), RV0, RV2)
    re <- matrix(c(
      sign(Vuong_test_left$statistic.Schwarz) * Vuong_test_left$p.value.Schwarz,
      sign(Vuong_test_right$statistic.Schwarz) * Vuong_test_right$p.value.Schwarz
    ), 1, 2)
    colnames(re) <- c("left Schwarz p", "right Schwarz p")
  } else {
    if (CDR == "C") {
      RV0 <- RVineStructureSelect(MultiInd(X, t_start, t_end),
        familyset = family_set, type = 1,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )
      RV1 <- RVineStructureSelect(MultiInd(X, t_start, t_point),
        familyset = family_set, type = 1,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )
      RV2 <- RVineStructureSelect(MultiInd(X, t_point, t_end),
        familyset = family_set, type = 1,
        selectioncrit = "BIC", method = "mle",
        indeptest = TRUE, trunclevel = trunc_tree
      )
      Vuong_test_left <- RVineVuongTest(MultiInd(X, t_start, t_point), RV0, RV1)
      Vuong_test_right <- RVineVuongTest(MultiInd(X, t_point, t_end), RV0, RV2)
      re <- matrix(c(
        sign(Vuong_test_left$statistic.Schwarz) * Vuong_test_left$p.value.Schwarz,
        sign(Vuong_test_right$statistic.Schwarz) * Vuong_test_right$p.value.Schwarz
      ), 1, 2)
      colnames(re) <- c("left Schwarz p", "right Schwarz p")
    } else {
      d <- dim(X)[2] - 1
      D_0 <- CDVineCopSelect(MultiInd(X, t_start, t_end),
        type = 2, familyset = unlist(family_set),
        selectioncrit = "BIC", indeptest = TRUE
      )
      D_1 <- CDVineCopSelect(MultiInd(X, t_start, t_point),
        type = 2, familyset = unlist(family_set),
        selectioncrit = "BIC", indeptest = TRUE
      )
      D_2 <- CDVineCopSelect(MultiInd(X, t_point, t_end),
        type = 2, familyset = unlist(family_set),
        selectioncrit = "BIC", indeptest = TRUE
      )
      Vuong_test_left <- CDVineVuongTest(MultiInd(X, t_start, t_point), 1:d, 1:d,
        D_0$family, D_1$family,
        D_0$par, D_1$par, D_0$par2, D_1$par2,
        Model1.type = 2, Model2.type = 2
      )
      Vuong_test_right <- CDVineVuongTest(MultiInd(X, t_point, t_end), 1:d, 1:d,
        D_0$family, D_2$family,
        D_0$par, D_2$par, D_0$par2, D_2$par2,
        Model1.type = 2, Model2.type = 2
      )
      re <- matrix(c(
        sign(Vuong_test_left$statistic.Schwarz) * Vuong_test_left$p.value.Schwarz,
        sign(Vuong_test_right$statistic.Schwarz) * Vuong_test_right$p.value.Schwarz
      ), 1, 2)
      colnames(re) <- c("left Schwarz p", "right Schwarz p")
    }
  }

  return(re)
}



#' @export
TestPoints_Vuong <- function(v_t_point, X_raw, delta, CDR = "D", trunc_tree = NA, family_set = 1,
                             pre_white = 0, ar_num = 1, sig_alpha = 0.05) {
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
  test_result <- as.data.frame(matrix(0, length(v_t_point), 4))
  if (length(v_t_point) == 0) {
    return(test_result)
  } else {
    test_result[, 1] <- v_t_point
    a <- c(1, v_t_point, T + 1)
    full_point <- a[!duplicated(a)]
    for (i in 1:length(v_t_point)) {
      cat(paste("Test for candidate", i, ": t =", v_t_point[i]), fill = TRUE)
      test_result[i, 2:3] <- Vuong_Multi_CDR_NewTestPoint(
        full_point[i + 1], full_point[i],
        full_point[i + 2], X, delta, CDR, trunc_tree, family_set
      )
      if (test_result[i, 2] < 0 & abs(test_result[i, 2]) < sig_alpha &
        test_result[i, 3] < 0 & abs(test_result[i, 3]) < sig_alpha
      ) {
        test_result[i, 4] <- "significant"
        cat("Significant! \n\n")
      } else {
        test_result[i, 4] <- "not significant"
        cat("Insignificant! \n\n")
      }
    }
    names(test_result) <- c("t", "left Schwarz p", "right Schwarz p", "judgement")
    return(test_result)
  }
}

#' @export
GetTestPlot_Vuong <- function(test_result, T, sig_alpha = 0.05) {
  if (dim(test_result)[1] == 0) {
    plot(1:T,1:T,type = "n",yaxt="n",ylab = NA)
    graphics::text("No candidate is found.")
  }
  else {
    LeftSch <- RightSch <- rep(0, T)
    for (ii in 1:dim(test_result)[1]) {
      if (1 / abs(test_result[ii, 2]) > 1e100) test_result[ii, 2] <- -1e-100
      if (1 / abs(test_result[ii, 3]) > 1e100) test_result[ii, 3] <- -1e-100
    }

    LeftSch[test_result[, 1]] <- log(abs(1 / test_result[, 2]))
    RightSch[test_result[, 1]] <- log(abs(1 / test_result[, 3]))
    y_min <- min(LeftSch, RightSch)
    y_max <- max(LeftSch, RightSch)
    plot(1:T, rep(1, T),
      ylab = c("log(1/p_value)"), type = "n", cex = .5, pch = 16,
      ylim = c(max(0.1, y_min - 0.1), y_max + 6)
    )

    n_point <- dim(test_result)[1]
    for (i in 1:n_point) {
      graphics::lines(data.frame(
        x = c(test_result[i, 1], test_result[i, 1]),
        y = c(
          log(1 / abs(min(test_result[i, 2:3]))),
          log(1 / abs(max(test_result[i, 2:3])))
        )
      ))
    }

    graphics::points(test_result[test_result[, 2] < 0, 1], log(-1 / test_result[test_result[, 2] < 0, 2]), pch = 20, col = "blue")
    graphics::points(test_result[test_result[, 3] < 0, 1], log(-1 / test_result[test_result[, 3] < 0, 3]), pch = 20, col = "green")

    graphics::points(test_result[test_result[, 2] > 0, 1], log(1 / test_result[test_result[, 2] > 0, 2]), pch = 4, col = 1)
    graphics::points(test_result[test_result[, 3] > 0, 1], log(1 / test_result[test_result[, 3] > 0, 3]), pch = 4, col = 1)

    graphics::abline(h = c(0, log(1 / sig_alpha)), lwd = 1)
    graphics::text(70, log(1 / sig_alpha) - 0.3, labels = paste(
      "alpha =", sig_alpha,
      "log(1/sig_alpha) =", round(log(1 / sig_alpha), digits = 3)
    ), cex = .6)
    txt <- c("left Schwarz p", "right Schwarz p", "sign_stat>0")
    graphics::legend("topright", legend = txt, col = c("blue", "green", 1), pch = c(16, 16, 4), cex = .6)
    if (sum(test_result[, 4] == "significant") != 0) {
      sig_re <- test_result[test_result[, 4] == "significant", ]
      log_inv_sig_re <- log(-1 / sig_re[, c(-1, -4)])
      p_max <- apply(log_inv_sig_re, 1, max)
      graphics::text(test_result[test_result[, 4] == "significant", 1], p_max + 1,
        labels = test_result[test_result[, 4] == "significant", 1], cex = 0.8
      )
    }
  }
}
