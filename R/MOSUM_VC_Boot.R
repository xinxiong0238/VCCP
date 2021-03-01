#' @importFrom VineCopula RVineStructureSelect
#' @importFrom VineCopula D2RVine
#' @importFrom VineCopula RVineCopSelect
#' @importFrom VineCopula RVineVuongTest
MultiGenXLocal <- function(X_raw, delta, CDR, trunc_tree, family_set, pre_white, ar_num) {
  T <- length(unique(X_raw[, 1]))
  subnum <- dim(X_raw)[1] / T
  cut_point0 <- c(1, T + 1)
  p_X = dim(X_raw)[2] - 1
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
  BIC_de <- rep(0, T)
  family_set <- unlist(family_set)
  for (i in (delta + 1):(T - delta - 1)) {
    if (CDR == "R") {
      BIC_de[i] <- RVineStructureSelect(MultiInd(X, (i - delta), (i + delta)),
                                        familyset = family_set,
                                        selectioncrit = "BIC", method = "mle",
                                        indeptest = TRUE, trunclevel = trunc_tree
      )$BIC -
        (RVineStructureSelect(MultiInd(X, (i - delta), i),
                              familyset = family_set,
                              selectioncrit = "BIC", method = "mle",
                              indeptest = TRUE, trunclevel = trunc_tree
        )$BIC +
          RVineStructureSelect(MultiInd(X, i, (i + delta)),
                               familyset = family_set,
                               selectioncrit = "BIC", method = "mle",
                               indeptest = TRUE, trunclevel = trunc_tree
          )$BIC)
      # print(c(BIC_de[i],i))
    } else {
      if (CDR == "C") {
        BIC_de[i] <- RVineStructureSelect(MultiInd(X, (i - delta), (i + delta)),
                                          familyset = family_set, type = 1,
                                          selectioncrit = "BIC", method = "mle",
                                          indeptest = TRUE, trunclevel = trunc_tree
        )$BIC -
          (RVineStructureSelect(MultiInd(X, (i - delta), i),
                                familyset = family_set, type = 1,
                                selectioncrit = "BIC", method = "mle",
                                indeptest = TRUE, trunclevel = trunc_tree
          )$BIC +
            RVineStructureSelect(MultiInd(X, i, (i + delta)),
                                 familyset = family_set, type = 1,
                                 selectioncrit = "BIC", method = "mle",
                                 indeptest = TRUE, trunclevel = trunc_tree
            )$BIC)
        # print(c(BIC_de[i],i))
      } else {
        D_X = D2RVine(1:p_X, family = rep(0, p_X*(p_X-1)/2),
                      par = rep(0, p_X*(p_X-1)/2),
                      par2 = rep(0, p_X*(p_X-1)/2))
        BIC_de[i] <- RVineCopSelect(MultiInd(X, (i - delta), (i + delta)), selectioncrit = "BIC",
                                    method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                                    indeptest = TRUE, trunclevel = trunc_tree)$BIC -
          (
            RVineCopSelect(MultiInd(X, (i - delta), i), selectioncrit = "BIC",
                           method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                           indeptest = TRUE, trunclevel = trunc_tree)$BIC +
              RVineCopSelect(MultiInd(X, i, (i + delta)), selectioncrit = "BIC",
                             method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                             indeptest = TRUE, trunclevel = trunc_tree)$BIC
          )
        # print(c(BIC_de[i],i))
      }
    }
  }
  return(list(BIC_de, X))
}



FindLocalMax <- function(BIC_re, range_set) {
  local_max <- c()
  k <- 1
  for (i in 1:(length(range_set) - 1)) {
    t_start <- range_set[i]
    t_end <- range_set[i + 1]
    if (max(BIC_re[t_start:(t_end - 1)]) > 0) {
      local_max[k] <- which.max(BIC_re[t_start:t_end]) + t_start - 1
      k <- k + 1
    }
  }
  return(local_max)
}


VC_MOSUM_Boot <- function(X_raw, delta, G = 0.1, CDR = "D", trunc_tree = NA, family_set = 1, pre_white = 0, ar_num = 1,
                          p = 0.3, N = 100, sig_alpha = 0.05) {
  message("MOSUM search ...")
  re <- MultiGenXLocal(X_raw, delta, CDR, trunc_tree, family_set, pre_white, ar_num)
  T <- length(unique(X_raw[, 1]))
  BIC_re <- re[[1]]
  te <- mosum::mosum(BIC_re, G)
  plot(te, display = "data", sub = "MOSUM segmentation")
  a <- sort(c(1, T + 1, te$cpts))
  ini_po <- a[!duplicated(a)]
  points0 <- FindLocalMax(BIC_re, ini_po)
  if (length(points0) == 0) {
    message("No candidate is found.")
    return(matrix(0, 0, 5))
  } else {
    if (length(points0) == 1) {
      tema_NRVBS <- TestPoints.Boot(
        points0, X_raw, delta, CDR, trunc_tree, family_set,
        pre_white, ar_num, p, N, sig_alpha
      )
      return(tema_NRVBS)
    } else {
      points_new <- points <- points0[order(-BIC_re[points0])]
      j <- 1
      while (j <= length(points_new)) {
        d_p <- abs(points_new[j] - points_new)
        points_new <- points_new[d_p >= delta | d_p == 0]
        j <- j + 1
      }
      points_new <- unique(sort(points_new))
      tema_NRVBS <- TestPoints.Boot(
        points_new, X_raw, delta, CDR, trunc_tree, family_set,
        pre_white, ar_num, p, N, sig_alpha
      )

      return(tema_NRVBS)
    }
  }
}
