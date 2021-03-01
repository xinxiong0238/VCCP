#' @importFrom VineCopula RVineStructureSelect
#' @importFrom VineCopula D2RVine
#' @importFrom VineCopula RVineCopSelect
#' @importFrom VineCopula RVineVuongTest
VC_OBS_Vuong <- function(X_raw, delta, CDR = "D", trunc_tree = NA, family_set = 1, pre_white = 0, ar_num = 1,
                         sig_alpha = 0.05) {
  T <- length(unique(X_raw[,1]))
  subnum <- dim(X_raw)[1]/T
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

  cut_point0 <- c(1, T+1)
  BIC_cut <- 1
  k=0
  TestMatrix <- as.data.frame(matrix(0, 1, 4))
  names(TestMatrix) <- c("t", "left Schwarz p", "right Schwarz p", "judgement")
  while (BIC_cut > 0) {
    k=k+1
    message(paste("Binary search + Vuong test for round",k,"..."))
    BIC_de <- rep(0, T)
    for (i in delta:(T-delta)) {
      if(sum(is.element((i-delta+1):(i+delta) , cut_point0)) != 0){
        next
      }else{
        t_point <- i
        t_start <- cut_point0[which.max(i <= cut_point0) - 1]
        t_end <- cut_point0[which.max(i <= cut_point0)]
        if(CDR=="R"){
          BIC_de[i] <-  RVineStructureSelect(MultiInd(X,t_start,t_end),
                                             familyset = unlist(family_set),
                                             selectioncrit = "BIC", method = "mle",
                                             indeptest = TRUE, trunclevel = trunc_tree)$BIC -
            ( RVineStructureSelect(MultiInd(X, t_start, t_point),
                                   familyset = unlist(family_set),
                                   selectioncrit = "BIC", method = "mle",
                                   indeptest = TRUE, trunclevel = trunc_tree)$BIC +
                RVineStructureSelect(MultiInd(X, t_point, t_end),
                                     familyset = unlist(family_set),
                                     selectioncrit = "BIC", method = "mle",
                                     indeptest = TRUE, trunclevel = trunc_tree)$BIC )
        }else{
          if(CDR=="C"){
            BIC_de[i] <-  RVineStructureSelect(MultiInd(X,t_start,t_end),
                                               familyset = unlist(family_set), type = 1,
                                               selectioncrit = "BIC", method = "mle",
                                               indeptest = TRUE, trunclevel = trunc_tree)$BIC -
              ( RVineStructureSelect(MultiInd(X, t_start, t_point),
                                     familyset = unlist(family_set), type = 1,
                                     selectioncrit = "BIC", method = "mle",
                                     indeptest = TRUE, trunclevel = trunc_tree)$BIC +
                  RVineStructureSelect(MultiInd(X, t_point, t_end),
                                       familyset = unlist(family_set), type = 1,
                                       selectioncrit = "BIC", method = "mle",
                                       indeptest = TRUE, trunclevel = trunc_tree)$BIC )
          }else{
            D_X = D2RVine(1:p_X, family = rep(0, p_X*(p_X-1)/2),
                          par = rep(0, p_X*(p_X-1)/2),
                          par2 = rep(0, p_X*(p_X-1)/2))
            BIC_de[i] <- RVineCopSelect(MultiInd(X, t_start, t_end), selectioncrit = "BIC",
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
          }
        }
      }
    }
    BIC_cut <- max(BIC_de)
    if(BIC_cut > 0){
      cut_new_point <- which.max(BIC_de)
      #cat(paste("Test for candidate", k, ": t =", cut_new_point), fill = TRUE)
      t_before <- cut_point0[which.max(cut_new_point <= cut_point0) - 1]
      t_after <- cut_point0[which.max(cut_new_point <= cut_point0)]
      testPoint <- Vuong_Multi_CDR_NewTestPoint(cut_new_point, t_before, t_after, X, delta, CDR, trunc_tree, family_set)
      #print(testPoint)
      if(testPoint[1] < 0 & abs(testPoint[1]) < sig_alpha &
         testPoint[2] < 0 & abs(testPoint[2]) < sig_alpha
      ){
        cut_point0 <- sort(c(cut_point0, cut_new_point))
        TestMatrix <- rbind(TestMatrix, c(cut_new_point, testPoint, "significant"))
        TestMatrix[,1] <- as.numeric(TestMatrix[,1])
        TestMatrix[,2] <- as.numeric(TestMatrix[,2])
        TestMatrix[,3] <- as.numeric(TestMatrix[,3])
      }else{
        TestMatrix <- rbind(TestMatrix, c(cut_new_point, testPoint, "not significant"))
        TestMatrix[,1] <- as.numeric(TestMatrix[,1])
        TestMatrix[,2] <- as.numeric(TestMatrix[,2])
        TestMatrix[,3] <- as.numeric(TestMatrix[,3])
        return(TestMatrix[-1,])
      }
    }else{
      return(TestMatrix[-1,])
    }
  }
}

