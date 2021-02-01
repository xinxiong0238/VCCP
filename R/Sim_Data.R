#' Simulate MVN data with 2 change points
#'
#' \code{random.mvn.simulate.2.changes} returns a multivariate normal
#'  data with 2 change points in the connection network. This function
#'  aims to display how VCCP model detects change points in examples.
#'
#' @param nobs A positive integer number equal to the length of the time series data (T).
#'  \code{nobs} must be multiple times of 3 since change points occur at T/3
#'  and 2T/3.
#'
#' @param n.ser A positive integer number indicating the dimensionality of the
#'  time series. \code{n.ser} must be larger than 8 as \code{random.mvn.simulate.2.changes}
#'  generates 3 different networks among 8 connected nodes. Other variables are generated
#'  to be independent of the 8 nodes all the time.
#'
#' @param seed A positive integer number that ensures reproducibility.
#'
#'
#' @return A \code{nobs} times \code{n.ser} matrix with 2 change points
#' at t=\code{nobs/3+1} and t=\code{nobs*2/3+1} in the network structure.
#'
#' @export
#' @examples
#' mvn <- random.mvn.simulate.2.changes(180, 8, seed = 101)
random.mvn.simulate.2.changes <- function(nobs, n.ser, seed) {
  data <- c()
  nobs <- nobs / 3

  ppp_11 <- 1
  qqq_11 <- 3
  rrr_11 <- 5
  sss_11 <- 7
  ttt_11 <- 8
  qqq_12 <- 2
  sss_12 <- 6
  ttt_12 <- 3

  qqq_21 <- 1
  rrr_21 <- 4
  sss_21 <- 5
  ttt_21 <- 7
  uuu_21 <- 3

  qqq_31 <- 2
  rrr_31 <- 3
  sss_31 <- 4
  ttt_31 <- 6
  uuu_31 <- 1
  www_31 <- 5
  sss_32 <- 1
  ttt_32 <- 7 # negatively correlated


  # =============================================================================
  # 1st Section
  set.seed(nobs * seed + 1)
  burn_in <- 100
  sigma.1 <- array(c(
    1, 0.69, 0.67, 0.76, 0.77, 0.69, 1, 0.64, 0.66, 0.69, 0.67, 0.64, 1,
    0.7, 0.72, 0.76, 0.66, 0.7, 1, 0.72, 0.77, 0.69, 0.72, 0.72, 1
  ), c(5, 5))
  data.1 <- mvtnorm::rmvnorm(nobs + burn_in, sigma = sigma.1)
  data.1 <- data.1[(burn_in + 1):(nobs + burn_in), ]


  set.seed(nobs * seed + 10001)
  burn_in <- 100
  sigma.2 <- array(c(1, 0.7, 0.65, 0.7, 1, 0.71, 0.65, 0.71, 1), c(3, 3))
  data.2 <- mvtnorm::rmvnorm(nobs + burn_in, sigma = sigma.2)
  data.2 <- data.2[(burn_in + 1):(nobs + burn_in), ]

  t <- n.ser # number of noise series
  data_1_used <- c()
  for (i in 1:t) {
    set.seed(nobs * 1 + 2 * i * 10001 + i + seed)
    data_1_used <- cbind(data_1_used, stats::rnorm(nobs))
  }

  data_1_used[, ppp_11] <- data.1[, 1]
  data_1_used[, qqq_11] <- data.1[, 2]
  data_1_used[, rrr_11] <- data.1[, 3]
  data_1_used[, sss_11] <- data.1[, 4]
  data_1_used[, ttt_11] <- data.1[, 5]
  data_1_used[, qqq_12] <- data.2[, 1]
  data_1_used[, sss_12] <- data.2[, 2]
  data_1_used[, ttt_12] <- data.2[, 3]

  # =============================================================================
  # 2nd Section
  set.seed(nobs * seed + 3 * 10001)
  burn_in <- 100
  sigma.1 <- array(c(
    1, 0.70, 0.695, 0.73, 0.67, 0.70, 1, 0.68, 0.66, 0.69, 0.695, 0.68, 1,
    0.7, 0.72, 0.73, 0.66, 0.7, 1, 0.75, 0.67, 0.69, 0.72, 0.75, 1
  ), c(5, 5))
  data.1 <- mvtnorm::rmvnorm(nobs + burn_in, sigma = sigma.1)
  data.1 <- data.1[(burn_in + 1):(nobs + burn_in), ]


  t <- n.ser # number of noise series
  data_2_used <- c()
  for (i in 1:t) {
    set.seed(nobs * seed + 5 * i * 10001 + i)
    data_2_used <- cbind(data_2_used, stats::rnorm(nobs))
  }

  data_2_used[, qqq_21] <- data.1[, 1]
  data_2_used[, rrr_21] <- data.1[, 2]
  data_2_used[, sss_21] <- data.1[, 3]
  data_2_used[, ttt_21] <- data.1[, 4]
  data_2_used[, uuu_21] <- data.1[, 5]

  # =============================================================================
  # 3rd Section
  set.seed(nobs * seed + 6 * 10001)
  burn_in <- 100
  sigma.1 <- array(c(
    1, 0.81, 0.65, 0.76, 0.68, 0.73, 0.81, 1, 0.68, 0.767, 0.69, 0.69, 0.65,
    0.68, 1, 0.7, 0.72, 0.71, 0.76, 0.767, 0.7, 1, 0.75, 0.68, 0.68, 0.69, 0.72,
    0.75, 1, 0.77, 0.73, 0.69, 0.71, 0.68, 0.77, 1
  ), c(6, 6))
  data.1 <- mvtnorm::rmvnorm(nobs + burn_in, sigma = sigma.1)
  data.1 <- data.1[(burn_in + 1):(nobs + burn_in), ]


  set.seed(nobs * seed + 8 * 10001)
  burn_in <- 100
  data.2 <- mvtnorm::rmvnorm(nobs + burn_in, sigma = array(c(1, -.8, -.8, 1), c(2, 2)))
  data.2 <- data.2[(burn_in + 1):(nobs + burn_in), ]

  t <- n.ser # number of noise series
  data_3_used <- c()
  for (i in 1:t) {
    set.seed(nobs * seed + 10 * i * 10001 + i)
    data_3_used <- cbind(data_3_used, stats::rnorm(nobs))
  }

  data_3_used[, qqq_31] <- data.1[, 1]
  data_3_used[, rrr_31] <- data.1[, 2]
  data_3_used[, sss_31] <- data.1[, 3]
  data_3_used[, ttt_31] <- data.1[, 4]
  data_3_used[, uuu_31] <- data.1[, 5]
  data_3_used[, www_31] <- data.1[, 6]
  data_3_used[, sss_32] <- data.2[, 1]
  data_3_used[, ttt_32] <- data.2[, 2]




  # Combine three data sets
  data <- rbind(data_1_used, data_2_used, data_3_used)

  return(data)
}


# =============================================================================
# =============================================================================
