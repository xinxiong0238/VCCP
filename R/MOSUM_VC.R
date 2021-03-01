VC_MOSUM <- function(X_raw, delta, G = 0.1, test = "V", CDR = "D", trunc_tree = NA,
                     family_set = 1, pre_white = 0, ar_num = 1,
                     p = 0.3, N = 100, sig_alpha = 0.05) {
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
    }
  }
}
