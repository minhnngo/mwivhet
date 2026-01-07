#' Generate Data for Grouped Instrument Model (No Covariates)
#'
#' @description
#' Simulates data from a linear instrumental variables model with many instruments,
#' specifically for design where instruments are mutually exclusive group
#' indicators (e.g., Judges). The function allows for treatment effect heterogeneity
#' (random coefficients) and correlation between heterogeneity and selection.
#'
#'
#' @param S Numeric. Concentration parameter \eqn{\mu^2}, determining instrument strength.
#' @param Het Numeric. Heterogeneity parameter, scaling the variance of the random slope \eqn{\xi}.
#' @param sigee Numeric. Variance of the structural error \eqn{\varepsilon}.
#' @param sigvv Numeric. Variance of the first-stage error \eqn{v}.
#' @param sigexi Numeric. Covariance between structural error \eqn{\varepsilon} and heterogeneity \eqn{\xi}.
#' @param sigev Numeric. Covariance between structural error \eqn{\varepsilon} and first-stage error \eqn{v} (Endogeneity).
#' @param ConTE Logical. If \code{TRUE}, simulates Constant Treatment Effects (\eqn{\xi = 0}).
#' @param beta Numeric. The true average treatment effect (ATE).
#' @param beta0 Numeric. The null hypothesis value for \eqn{\beta} (used to compute residuals \eqn{e}).
#' @param K Integer. The number of instruments (groups) minus one. Total groups \eqn{J = K + 1}.
#' @param c Integer. The number of observations per group (balanced block size).
#'
#' @details
#' The data generating process (DGP) is:
#' \deqn{X_i = \pi_{g(i)} + v_i}
#' \deqn{Y_i = X_i (\beta + \xi_i) + \varepsilon_i}
#'
#' The instrument coefficients \eqn{\pi} are drawn from a balanced design taking values
#' \eqn{\{-s, 0, s\}}.
#'
#' If \code{ConTE = FALSE}, the error vector \eqn{(\varepsilon, \xi, v)} is drawn from a
#' trivariate normal distribution. The covariance structure varies across blocks to introduce
#' correlation between the instrument strength/direction and the correlation of \eqn{v} and \eqn{\xi}.
#'
#' @return A data frame containing:
#' \item{group}{Group identifier (1 to \eqn{J}).}
#' \item{pi}{True instrument mean for the observation.}
#' \item{eps}{Structural error \eqn{\varepsilon}.}
#' \item{xi}{Heterogeneity term \eqn{\xi}.}
#' \item{v}{First-stage error \eqn{v}.}
#' \item{X}{Endogenous regressor.}
#' \item{Y}{Outcome variable.}
#' \item{e}{Residual under the null (\eqn{Y - X\beta_0}).}
#' \item{MX}{Leverage-adjusted regressor \eqn{M X}.}
#' \item{Me}{Leverage-adjusted residual \eqn{M e}.}
#'
#' @export
GenData_nocov <- function(S = 3, Het = 3, sigee = 1, sigvv = 1, sigexi = 0, sigev = 0, ConTE = FALSE,
                          beta = 0, beta0 = 0, K = 40, c = 5) {
  # Calculations as function of parameters
  J <- K + 1
  n <- J * c
  group <- as.vector(seq(1, J) %x% rep(1, c))

  # Create Z matrix of indicators
  Z <- matrix(0, nrow = length(group), ncol = length(unique(group)))
  Z[cbind(seq_along(group), group)] <- 1

  ## Calculations related to the projection matrix
  ZZ <- t(Z) %*% Z
  ZZ_inv <- solve(ZZ)
  ZZ_inv2 <- chol(ZZ_inv)
  ZZZ_inv <- Z %*% ZZ_inv
  ZZZ_inv2 <- Z %*% ZZ_inv2 # n x k mx
  H <- rowSums(ZZZ_inv2^2) # n x 1 vector of P_ii
  P <- Z %*% ZZ_inv %*% t(Z)
  PP <- P * P

  h <- sqrt(Het / sqrt(K) / (c - 1))
  s <- sqrt(S / sqrt(K) / (c - 1))
  pi <- c(0, rep(-s, K / 2), rep(s, K / 2)) # balanced

  # Construct DGP Sigma
  sigxixi <- h / sigvv

  if (ConTE == TRUE) {
    Sigma <- matrix(c(sigee, sigev, sigev, sigvv), nrow = 2)
  } else {
    sigxiv01 <- h
    sigxiv23 <- -h
    Sigma0 <- matrix(c(sigee, sigexi, sigev, sigexi, sigxixi, 0, sigev, 0, sigvv), nrow = 3)
    # matrixcalc::is.positive.definite(Sigma0)
    Sigma01 <- matrix(c(sigee, sigexi, sigev, sigexi, sigxixi, sigxiv01, sigev, sigxiv01, sigvv), nrow = 3)
    Sigma23 <- matrix(c(sigee, sigexi, sigev, sigexi, sigxixi, sigxiv23, sigev, sigxiv23, sigvv), nrow = 3)
  }

  ## Construct dataframe
  df <- as.data.frame(group)
  df$pi <- as.vector(Z %*% pi)

  # Make random draws.
  if (ConTE == TRUE) {
    vmat <- mvtnorm::rmvnorm(n, mean = rep(0, 2), sigma = Sigma)
    df$eps <- vmat[, 1]
    df$xi <- 0
    df$v <- vmat[, 2]
  } else {
    vmat0 <- mvtnorm::rmvnorm(c, mean = rep(0, 3), sigma = Sigma0)
    vmat01 <- mvtnorm::rmvnorm(K * c / 2, mean = rep(0, 3), sigma = Sigma01)
    vmat23 <- mvtnorm::rmvnorm(K * c / 2, mean = rep(0, 3), sigma = Sigma23)
    vmat <- rbind(vmat0, vmat01[1:(K * c / 4), ], vmat23, vmat01[((K * c / 4) + 1):(K * c / 2), ])
    df$eps <- vmat[, 1]
    df$xi <- vmat[, 2]
    df$v <- vmat[, 3]
  }

  df$X <- df$pi + df$v
  df$Y <- df$X * (beta + df$xi) + df$eps
  df$e <- df$Y - df$X * beta0
  M <- diag(n) - P
  df$MX <- c(M %*% df$X)
  df$Me <- c(M %*% df$e)
  df$MY <- M %*% df$Y

  df
}

#' Generate Data with Covariates and Binary Endogenous Regressor
#'
#' @description
#' Simulates data from a grouped instrument design model that includes discrete
#' covariates and a binary endogenous regressor generated via a threshold crossing model.
#' The function incorporates essential heterogeneity, where the random treatment coefficient
#' is correlated with the first-stage selection error.
#'
#'
#' @param S Numeric. Concentration parameter \eqn{\mu^2}, scaling the instrument strength.
#' @param p1 Numeric (0 to 1). Probability parameter controlling the correlation between
#'   the first-stage error \eqn{v} and the treatment heterogeneity \eqn{\xi}.
#' @param Het Numeric. Heterogeneity parameter, scaling the magnitude of \eqn{\xi}.
#' @param sigeps Numeric. Standard deviation of the structural error \eqn{\varepsilon}.
#' @param sigev Numeric. Coefficient governing the correlation between \eqn{\varepsilon} and \eqn{v}.
#' @param beta Numeric. Average Treatment Effect (ATE).
#' @param beta0 Numeric. Null hypothesis value for \eqn{\beta}.
#' @param omega Numeric. Scaling parameter for the covariate effects \eqn{\gamma}.
#' @param K Integer. Number of instrument groups minus one.
#' @param c Integer. Number of observations per group.
#'
#' @details
#' The Data Generating Process (DGP) is structured as follows:
#'
#' \strong{Design Matrix:}
#' \itemize{
#'   \item \strong{Covariates (W):} Defined by pairs of instrument groups (e.g., strata).
#'   \item \strong{Instruments (Z):} Nested within covariates (e.g., judges within years).
#'   \item \strong{Weights:} The UJIVE weighting matrix \eqn{G} is computed as \eqn{U(P_{[Z,W]}) - U(P_W)}.
#' }
#'
#' \strong{Structural Equations:}
#' \deqn{X_i = 2 \cdot \mathbb{I}(v_i < \pi_{g(i)} + \gamma_{d(i)}) - 1}
#' \deqn{Y_i = X_i (\beta + \xi_i) + \gamma_{d(i)} + \varepsilon_i}
#'
#' \strong{Heterogeneity:}
#' The variable \eqn{X} is binary (taking values -1, 1). The random slope \eqn{\xi_i} is
#' drawn conditionally on the first-stage latent variable \eqn{v_i}, inducing correlation
#' between selection into treatment and treatment gains (Essential Heterogeneity).
#'
#' @return A data frame containing:
#' \item{group}{Overall observation index.}
#' \item{groupZ}{Instrument group identifier (the effective Z).}
#' \item{groupW}{Covariate stratum identifier.}
#' \item{pi}{Instrument effect \eqn{\pi}.}
#' \item{gammad}{Covariate effect \eqn{\gamma}.}
#' \item{v}{Latent first-stage variable (Uniform).}
#' \item{xi}{Random treatment effect heterogeneity.}
#' \item{X}{Binary endogenous regressor (-1 or 1).}
#' \item{Y}{Outcome variable.}
#' \item{e}{Residual under the null.}
#' \item{MX, Me, MY}{Variables projected onto the annihilator \eqn{M = I - P}.}
#'
#' @export
GenData_cov <- function(S = 5, p1 = 7 / 8, Het = 3, sigeps = .5, sigev = .1, beta = 0, beta0 = 0, omega = .1,
                        K = 20, c = 5) {
  # Calculations as function of parameters
  J <- 2 * (K + 1)
  n <- J * c
  group <- as.vector(seq(1, J) %x% rep(1, c))

  # Create discrete W
  groupW <- as.vector(seq(1, K + 1) %x% rep(1, 2 * c))
  groupZ <- groupW * (group %% 2)

  # Create Z matrix of indicators
  Z <- matrix(0, nrow = length(groupZ), ncol = length(unique(groupZ)) - 1)
  Z[cbind(seq_along(groupZ), groupZ)] <- 1

  ## Calculations related to Z
  ZZ <- t(Z) %*% Z
  ZZ_inv <- solve(ZZ)
  ZZ_inv2 <- chol(ZZ_inv)
  ZZZ_inv <- Z %*% ZZ_inv
  ZZZ_inv2 <- Z %*% ZZ_inv2 # n x k mx
  HZ <- Z %*% ZZ_inv %*% t(Z)

  ## Calculations related to W
  Wd <- matrix(0, nrow = length(groupW), ncol = length(unique(groupW)))
  Wd[cbind(seq_along(groupW), groupW)] <- 1
  # W <- cbind(Wd,Wcfix)
  W <- Wd
  WW <- t(W) %*% W
  WW_inv <- solve(WW)
  WW_inv2 <- chol(WW_inv)
  WWW_inv <- W %*% WW_inv
  WWW_inv2 <- W %*% WW_inv2
  HW <- W %*% WW_inv %*% t(W)
  MW <- diag(n) - HW

  ## Combine Z and W
  Q <- cbind(Z, W)
  QQ <- t(Q) %*% Q
  QQ_inv <- solve(QQ)
  QQ_inv2 <- chol(QQ_inv)
  QQQ_inv <- Q %*% QQ_inv
  QQQ_inv2 <- Q %*% QQ_inv2

  # P for full projection on Z and W
  P <- Q %*% QQ_inv %*% t(Q)
  M <- diag(n) - P

  ## G for UJIVE
  G <- solve(diag(n) - diag(diag(P))) %*% (P - diag(diag(P))) -
    solve(diag(n) - diag(diag(HW))) %*% (HW - diag(diag(HW)))

  PP_off <- G * G

  h <- sqrt(Het / sqrt(K) / (c - 1))
  s <- sqrt(S / sqrt(K) / (c - 1))
  pi <- c(0, rep(s, K / 2), rep(-s, K / 2))

  # Generate gammad
  gammad <- c(0, rep(omega, K / 2), rep(-omega, K / 2))

  ## Construct dataframe
  df <- as.data.frame(group)
  df$groupZ <- groupZ
  df$groupW <- groupW
  df$pi <- as.vector(Z %*% pi)
  df$gammad <- as.vector(Wd %*% gammad)
  dfbase <- df # This is fixed across simulations

  # Make random draws. This varies across sims
  df <- dfbase
  df$v <- runif(n, -1, 1)
  df$hetdraw <- runif(n, 0, 1)
  df$vpos <- ifelse(df$v > 0, 1, -1)
  df$ppos <- p1
  df$ppos[df$vpos < 0] <- 1 - p1

  df$eps <- sigev * df$vpos + rnorm(n, 0, sigeps)
  df$xi <- c(rep(0, 2 * c), rep(c(rep(1, c), rep(0, c), rep(0, c), rep(0, c)), K / 2)) * ifelse(df$hetdraw < df$ppos, h, -h) +
    c(rep(0, 2 * c), rep(c(rep(0, c), rep(0, c), rep(1, c), rep(0, c)), K / 2)) * ifelse(df$hetdraw < df$ppos, -h, h)
  df$X <- ifelse(df$v < df$pi + df$gammad, 1, -1)
  df$Y <- df$X * (beta + df$xi) + df$gammad + df$eps
  df$e <- df$Y - df$X * beta0
  df$MX <- c(M %*% df$X)
  df$Me <- c(M %*% df$e)
  df$MY <- c(M %*% df$Y)
  df$dM <- diag(M)

  df
}
