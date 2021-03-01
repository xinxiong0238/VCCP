VC_NBS <- function(X_raw, delta, test = "V", CDR = "D", trunc_tree = NA,
                        family_set = 1, pre_white = 0, ar_num = 1,
                         p = 0.3, N = 100, sig_alpha = 0.05) {
  result = VC.NBS.FindPoints(X_raw, delta, CDR, trunc_tree, family_set,
                             pre_white, ar_num)
  if(test=="V"){
    infer =  TestPoints.Vuong(result, X_raw, delta, CDR,
                              trunc_tree, family_set,
                              pre_white, ar_num, sig_alpha)
    return(infer)

  }else{
    if(test=="B"){
      infer = TestPoints.Boot(result, X_raw, delta, CDR,
                              trunc_tree, family_set,
                              pre_white, ar_num, p, N, sig_alpha)
      return(infer)

    }
  }
}
