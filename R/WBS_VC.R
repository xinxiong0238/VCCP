MultiGenRandomInterval <- function(X_raw, delta, start, end, M, pre_white, ar_num){
  T <- length(unique(X_raw[,1]))
  subnum <- dim(X_raw)[1]/T
  X <- X_raw
  if(pre_white == 1){
    for (i in 1:subnum) {
      armodel <- stats::ar(X_raw[((i-1)*T+1):(i*T),-1], FALSE, ar_num)
      ar_resid <- armodel$resid
      ar_resid[which(is.na(ar_resid)==TRUE)] <- 0
      X[((i-1)*T+1):(i*T),-1] <- ar_resid
    }
  }
  X[,-1] <- VineCopula::pobs(X[,-1])
  t1 <- sample(start:end, M, replace = TRUE)
  t2 <- sample(start:end, M, replace = TRUE)
  t <- cbind(t1, t2)
  t_sort <- t(apply(t, 1, sort))
  t_sort <- t_sort[t_sort[,2]-t_sort[,1]>2*delta, ]
  X_subsam <- list()
  for (i in 1:dim(t_sort)[1]) {
    X_subsam1 <- c()
    for (j in 1:subnum) {
      X_subsam1 <- rbind(X_subsam1,
                         X[((j-1)*T+t_sort[i,1]):((j-1)*T+t_sort[i,2]), ])
    }
    X_subsam=c(list(X_subsam1),X_subsam)
  }
  return(X_subsam)
}

#' @importFrom VineCopula RVineStructureSelect
#' @importFrom VineCopula D2RVine
#' @importFrom VineCopula RVineCopSelect
#' @importFrom VineCopula RVineVuongTest
VC_WBS_FindPoints <- function(X_raw, delta, M=NA, CDR="D", trunc_tree=NA, family_set=1, pre_white=0, ar_num=1){
  T <- length(unique(X_raw[, 1]))
  if(is.na(M)==TRUE){
    M <- floor(9*log(T))
  }
  p_X = dim(X_raw)[2] - 1
  X_list = MultiGenRandomInterval(X_raw, delta, 1, T, M, pre_white, ar_num)
  BIC_cut_series <- cut_new_point_series <- c()
  X_list_new <- list()
  k=0
  message("Search for pseudo sample data...")
  for (j in 1:length(X_list)) {
    #message(paste("Search for pseudo sample data (",j,"/",length(X_list),") ..."),fill = TRUE)
    BIC_de <- rep(0, T)
    t_Xj <- length(unique(X_list[[j]][,1]))
    cut_point0 <- c(min(unique(X_list[[j]][,1])),max(unique(X_list[[j]][,1]))+1)
    for (i in cut_point0[1]:cut_point0[2]) {
      if(sum(is.element((i-delta):(i+delta) , cut_point0)) != 0){
        next
      }else{
        if(CDR=="R"){
          t_point <- i
          BIC_de[i] <-  RVineStructureSelect(MultiInd(X_list[[j]], 1, t_Xj+1),
                                             familyset=family_set,
                                             selectioncrit = "BIC", method = "mle",
                                             indeptest = TRUE, trunclevel = trunc_tree)$BIC -
            ( RVineStructureSelect(MultiInd(X_list[[j]], 1, i - cut_point0[1]+1),
                                   familyset=family_set,
                                   selectioncrit = "BIC", method = "mle",
                                   indeptest = TRUE, trunclevel = trunc_tree)$BIC +
                RVineStructureSelect(MultiInd(X_list[[j]],(i - cut_point0[1] + 1), t_Xj+1),
                                     familyset=family_set,
                                     selectioncrit = "BIC", method = "mle",
                                     indeptest = TRUE, trunclevel = trunc_tree)$BIC )

          #print(c(BIC_de[i],i))
        }else{
          if(CDR=="C"){
            t_point <- i
            BIC_de[i] <-  RVineStructureSelect(MultiInd(X_list[[j]], 1, t_Xj+1),
                                               familyset=family_set, type=1,
                                               selectioncrit = "BIC", method = "mle",
                                               indeptest = TRUE, trunclevel = trunc_tree)$BIC -
              ( RVineStructureSelect(MultiInd(X_list[[j]], 1, i - cut_point0[1]+1),
                                     familyset=family_set, type=1,
                                     selectioncrit = "BIC", method = "mle",
                                     indeptest = TRUE, trunclevel = trunc_tree)$BIC +
                  RVineStructureSelect(MultiInd(X_list[[j]],(i - cut_point0[1] + 1), t_Xj+1),
                                       familyset=family_set, type=1,
                                       selectioncrit = "BIC", method = "mle",
                                       indeptest = TRUE, trunclevel = trunc_tree)$BIC )

            #print(c(BIC_de[i],i))
          }else{
            t_point <- i

            D_X = D2RVine(1:p_X, family = rep(0, p_X*(p_X-1)/2),
                          par = rep(0, p_X*(p_X-1)/2),
                          par2 = rep(0, p_X*(p_X-1)/2))
            BIC_de[i] <- RVineCopSelect(MultiInd(X_list[[j]], 1, t_Xj+1), selectioncrit = "BIC",
                                        method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                                        indeptest = TRUE, trunclevel = trunc_tree)$BIC -
              (
                RVineCopSelect(MultiInd(X_list[[j]], 1, i - cut_point0[1]+1), selectioncrit = "BIC",
                               method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                               indeptest = TRUE, trunclevel = trunc_tree)$BIC +
                  RVineCopSelect(MultiInd(X_list[[j]],(i - cut_point0[1] + 1), t_Xj+1), selectioncrit = "BIC",
                                 method = "mle", Matrix = D_X$Matrix,familyset = unlist(family_set),
                                 indeptest = TRUE, trunclevel = trunc_tree)$BIC
              )
            #print(c(BIC_de[i],i))
          }
        }
      }
    }
    max_BIC_de <- max(BIC_de)
    if(max_BIC_de > 0){
      k=k+1
     # if(which.max(BIC_de) %in% cut_new_point_series){
     #    cat(paste("Find candidate", k,
     #              ": t =", which.max(BIC_de),"(duplicated) \n \n"), fill = TRUE)
     #  }else{
     #    cat(paste("Find candidate", k,
     #              ": t =", which.max(BIC_de),"\n \n"), fill = TRUE)
     #  }
      X_list_new[[k]] <- X_list[[j]]
      BIC_cut_series[k] <- max(BIC_de)
      cut_new_point_series[k] <- which.max(BIC_de)
    }else{
      # cat(paste("No candidate is found in this sample data. \n\n"), fill = TRUE)
    }
  }
  return(list(X_list_new, cut_new_point_series, BIC_cut_series))
}


Multi_Find_startend_ind <- function(rank_X_list){
  re <- matrix(0, length(rank_X_list), 2)
  for (i in 1:length(rank_X_list)) {
    re[i, 1] <- min(unique(rank_X_list[[i]][,1]))
    re[i, 2] <- max(unique(rank_X_list[[i]][,1]))
  }
  return(re)
}

TwoMulti <- function(x, y){
  z <- c()
  for (i in 1:length(x)) {
    z[i] <- x[i] * y[i]
  }
  return(z)
}


VC_WBS_Boot <- function(X_raw, delta, M=NA, CDR="D", trunc_tree=NA,
                                family_set=1, pre_white=0, ar_num=1, p=0.3, N=100, sig_alpha=0.05){
  VC_WBS_FindPoints.result <- VC_WBS_FindPoints(X_raw, delta, M, CDR, trunc_tree,
                                                family_set, pre_white, ar_num)
  X_list_new = VC_WBS_FindPoints.result[[1]]
  cut_new_point_series = VC_WBS_FindPoints.result[[2]]
  BIC_cut_series = VC_WBS_FindPoints.result[[3]]
  rank_cut_point_series <- c()
  rank_X_list <- list()
  ind <- which(!duplicated(BIC_cut_series))
  cut_new_point_series <- cut_new_point_series[!duplicated(BIC_cut_series)]
  BIC_cut_series <- BIC_cut_series[!duplicated(BIC_cut_series)]
  X_list_2 <- list()
  m=1
  tema <- as.data.frame(matrix(0, 1, 5))
  names(tema) <- c("t", "reduced BIC", "Lower_CI", "Upper_CI", "judgement")
  if(length(ind)==0){
    message("No candidate is found.")
    return(tema[-1,0])
  }else{
    for (i in 1:length(ind)) {
      X_list_2[[m]] <- X_list_new[[ind[i]]]
      m=m+1
    }
    X_list_new <- X_list_2
    ra <- floor(rank(-BIC_cut_series))
    rank_cut_point_series[ra] <- cut_new_point_series
    for (i in 1:length(cut_new_point_series)) {
      rank_X_list[[ra[i]]] <- X_list_new[[i]]
    }
    rank_ind_ma <- Multi_Find_startend_ind(rank_X_list)
    rank_BIC_cut_series <- sort(BIC_cut_series, decreasing = TRUE)
    test_I <- "significant"
    tema <- as.data.frame(matrix(0, 1, 5))
    names(tema) <- c("t", "reduced BIC", "Lower_CI", "Upper_CI", "judgement")
    k=1
    sig_cut_point <- c()
    kk <- rep(1, length(rank_ind_ma[, 1]))
    aa = 0
    message("Perform stationary bootstrap test on candidates...")
    while (test_I == "significant") {
      aa = aa + 1
      # cat(paste("Test for candidate", aa, ": t =", rank_cut_point_series[k]), fill = TRUE)
      t_start <- rank_X_list[[k]][1,1]
      t_end <- rank_X_list[[k]][dim(rank_X_list[[k]])[1],1] + 1
      test_CI <- Multi_CDR_NewTestPoint(rank_cut_point_series[k] - t_start + 1,
                                        1, t_end - t_start+1, rank_X_list[[k]][,-1],
                                        delta, CDR, trunc_tree, family_set, p, N,sig_alpha)

      if(test_CI[1] > test_CI[3]){
        test_I <- "significant"
        #cat("Significant! \n\n")
        sig_cut_point <- c(sig_cut_point, rank_cut_point_series[k])
        tema <- rbind(tema, c(rank_cut_point_series[k], test_CI, test_I))
        tema[,1] <- as.numeric(tema[,1])
        tema[,2] <- as.numeric(tema[,2])
        tema[,3] <- as.numeric(tema[,3])
        tema[,4] <- as.numeric(tema[,4])
        kk <- TwoMulti(kk ,rank_ind_ma[, 1] > rank_cut_point_series[k] |
                         rank_ind_ma[, 2] < rank_cut_point_series[k])
        if(sum(kk!=0)==0) test_I <- "not significant" else k <- which.max(kk)
      }else{
        test_I <- "not significant"
        #cat("Insignificant! \n\n")
        tema <- rbind(tema, c(rank_cut_point_series[k], test_CI, test_I))
        tema[,1] <- as.numeric(tema[,1])
        tema[,2] <- as.numeric(tema[,2])
        tema[,3] <- as.numeric(tema[,3])
        tema[,4] <- as.numeric(tema[,4])
      }
    }
    return(tema[-1,])
  }
}




VC_WBS_Vuong <- function(X_raw, delta, M=NA, CDR="D", trunc_tree=NA,
                         family_set=1, p=0.3, N=100, pre_white=0,
                         ar_num=1, sig_alpha=0.05){
  VC_WBS_FindPoints.result <- VC_WBS_FindPoints(X_raw, delta, M, CDR, trunc_tree,
                                                family_set, pre_white, ar_num)
  X_list_new = VC_WBS_FindPoints.result[[1]]
  cut_new_point_series = VC_WBS_FindPoints.result[[2]]
  BIC_cut_series = VC_WBS_FindPoints.result[[3]]
  rank_cut_point_series <- c()
  rank_X_list <- list()
  ind <- which(!duplicated(BIC_cut_series))
  cut_new_point_series <- cut_new_point_series[!duplicated(BIC_cut_series)]
  BIC_cut_series <- BIC_cut_series[!duplicated(BIC_cut_series)]
  X_list_2 <- list()
  m=1
  tema <- as.data.frame(matrix(0, 1, 4))
  names(tema) <- c("t", "left Schwarz p", "right Schwarz p", "judgement")
  if(length(ind)==0){
    message("No candidate is found.")
    return(tema[-1,0])
  }else{
    for (i in 1:length(ind)) {
      X_list_2[[m]] <- X_list_new[[ind[i]]]
      m=m+1
    }
    X_list_new <- X_list_2
    ra <- floor(rank(-BIC_cut_series))
    rank_cut_point_series[ra] <- cut_new_point_series
    for (i in 1:length(cut_new_point_series)) {
      rank_X_list[[ra[i]]] <- X_list_new[[i]]
    }
    rank_ind_ma <- Multi_Find_startend_ind(rank_X_list)
    rank_BIC_cut_series <- sort(BIC_cut_series, decreasing = TRUE)
    test_I <- "significant"
    tema <- as.data.frame(matrix(0, 1, 4))
    names(tema) <- c("t", "left Schwarz p", "right Schwarz p", "judgement")
    k=1
    sig_cut_point <- c()
    kk <- rep(1, length(rank_ind_ma[, 1]))
    aa = 0
    message("Perform Vuong test on candidates...")
    while (test_I == "significant") {
      aa = aa + 1
      # cat(paste("Test for candidate", aa, ": t =", rank_cut_point_series[k]), fill = TRUE)
      t_start <- rank_X_list[[k]][1,1]
      t_end <- rank_X_list[[k]][dim(rank_X_list[[k]])[1],1] + 1
      test_CI <- Vuong_Multi_CDR_NewTestPoint(rank_cut_point_series[k] - t_start + 1,
                                              1, t_end - t_start+1, rank_X_list[[k]][,-1],
                                              delta, CDR, trunc_tree, family_set)

      if(test_CI[1] < 0 & abs(test_CI[1]) < sig_alpha &
         test_CI[2] < 0 & abs(test_CI[2]) < sig_alpha ){
        test_I <- "significant"
        # cat("Significant! \n\n")
        sig_cut_point <- c(sig_cut_point, rank_cut_point_series[k])
        tema <- rbind(tema, c(rank_cut_point_series[k], test_CI, test_I))
        tema[,1] <- as.numeric(tema[,1])
        tema[,2] <- as.numeric(tema[,2])
        tema[,3] <- as.numeric(tema[,3])
        kk <- TwoMulti(kk ,rank_ind_ma[, 1] > rank_cut_point_series[k] |
                         rank_ind_ma[, 2] < rank_cut_point_series[k])
        if(sum(kk!=0)==0) test_I <- "not significant" else k <- which.max(kk)
      }else{
        # cat("Insignificant! \n\n")
        test_I <- "not significant"
        tema <- rbind(tema, c(rank_cut_point_series[k], test_CI, test_I))
        tema[,1] <- as.numeric(tema[,1])
        tema[,2] <- as.numeric(tema[,2])
        tema[,3] <- as.numeric(tema[,3])
      }
    }
    return(tema[-1,])
  }
}



VC_WBS <- function(X_raw, delta, M = NA, test = "V", CDR = "D", trunc_tree = NA,
                   family_set = 1, pre_white = 0, ar_num = 1,
                   p = 0.3, N = 100, sig_alpha = 0.05) {
  if (test == "V") {
    infer <- VC_WBS_Vuong(
      X_raw, delta, M, CDR, trunc_tree,
      family_set, pre_white, ar_num, sig_alpha
    )
    return(infer)
  } else {
    if (test == "B") {
      infer <- VC_WBS_Boot(
        X_raw, delta, M, CDR,
        trunc_tree, family_set,
        pre_white, ar_num, p, N, sig_alpha
      )
      return(infer)
    }
  }
}
