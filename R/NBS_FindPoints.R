MultiInd <- function(X, start, end) {
  T <- length(unique(X[, 1]))
  subnum <- dim(X)[1] / T
  sam <- c()
  for (i in 1:subnum) {
    sam <- rbind(sam, X[(start + (i - 1) * T):(end - 1 + (i - 1) * T), ])
  }
  return(sam[, -1])
}


#' @export
#' @importFrom VineCopula RVineStructureSelect
#' @importFrom CDVine CDVineCopSelect
#' @importFrom CDVine CDVineBIC

NBS.VC.FindPoints <- function(X_raw, delta, CDR = "D", trunc_tree = NA, family_set0 = 1, pre_white = 0, ar_num = 1) {
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
  } else {
    X <- X_raw
  }
  X[, -1] <- VineCopula::pobs(X[, -1])
  while (BIC_cut > 0) {
    k <- k + 1
    print(k)
    BIC_de <- rep(0, T)
    for (i in delta:(T - delta)) {
      if (sum(is.element((i - delta + 1):(i + delta), cut_point0)) != 0) {
        next
      } else {
        t_point <- i
        t_start <- cut_point0[which.max(i <= cut_point0) - 1]
        t_end <- cut_point0[which.max(i <= cut_point0)]
        if (CDR == "R") {
          family_set0 <- c(0, family_set0)
          family_set <- family_set0
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
          print(c(BIC_de[i], i))
        } else {
          if (CDR == "C") {
            family_set <- family_set0
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
            print(c(BIC_de[i], i))
          } else {
            family_set <- family_set0
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
            print(c(BIC_de[i], i))
          }
        }
      }
    }
    BIC_cut <- max(BIC_de)
    cut_new_point <- which.max(BIC_de)
    cut_point0 <- sort(c(cut_point0, cut_new_point))
  }
  return(list(c(cut_point0[cut_point0 != cut_new_point & cut_point0 != T + 1 & cut_point0 != 1]), X))
}

