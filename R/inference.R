#' Title
#'
#' @param df
#' @param groupW
#' @param group
#' @param X
#' @param Y
#' @param MX
#' @param MY
#' @param q
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetCIcoef <- function(df, groupW, group, X, Y, MX, MY, q = qnorm(.975)^2, noisy = FALSE) {
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  ## Calculate Components
  C0 <- A1type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
    2 * A1type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
    A1type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = Y, lpos = MX, noisy = noisy) -
    A4type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = X, lpos = MY, noisy = noisy) -
    A4type_sum(df, group, groupW, ipos = Y, jpos = Y, kpos = X, lpos = MX, noisy = noisy)
  C1 <- -(A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
            2 * A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
            A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) -
            A4type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) -
            A4type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) +
            A1type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            2 * A1type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            A1type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy))
  C2 <- A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    2 * A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy)
  PXY <- GetLM(df, X, Y, groupW, group)
  PXX <- GetLM(df, X, X, groupW, group)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Title
#'
#' @param df
#' @param groupW
#' @param group
#' @param X
#' @param Y
#' @param MX
#' @param MY
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetSigMx <- function(df, groupW, group, X, Y, MX, MY, noisy = FALSE) {
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  ## Calculate Components
  sig11 <- 4 * A1type_sum(df, group, groupW, ipos = Y, jpos = Y, kpos = Y, lpos = MY, noisy = noisy) -
    2 * A4type_sum(df, group, groupW, ipos = Y, jpos = Y, kpos = Y, lpos = MY, noisy = noisy)
  sig22 <- A1type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
    2 * A1type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
    A1type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = Y, lpos = MX, noisy = noisy) -
    A4type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = Y, lpos = MX, noisy = noisy) -
    A4type_sum(df, group, groupW, ipos = Y, jpos = Y, kpos = X, lpos = MX, noisy = noisy)
  sig33 <- 4 * A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    2 * A4type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy)
  sig12 <- 1 * A1type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = Y, lpos = MY, noisy = noisy) +
    3 * A1type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = Y, lpos = MY, noisy = noisy) +
    1 * A4type_sum(df, group, groupW, ipos = X, jpos = Y, kpos = Y, lpos = MY, noisy = noisy) -
    3 * A4type_sum(df, group, groupW, ipos = Y, jpos = X, kpos = Y, lpos = MY, noisy = noisy)
  sig23 <- 1 * A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
    3 * A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
    1 * A4type_sum(df, group, groupW, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) -
    3 * A4type_sum(df, group, groupW, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy)
  sig13 <- 4 * A1type_sum(df, group, groupW, ipos = X, jpos = X, kpos = Y, lpos = MY, noisy = noisy) -
    2 * A4type_sum(df, group, groupW, ipos = X, jpos = X, kpos = Y, lpos = MY, noisy = noisy)

  c(sig11, sig22, sig33, sig12, sig23, sig13)
}

#' Title
#'
#' @param df
#' @param groupZ
#' @param X
#' @param Y
#' @param MX
#' @param MY
#' @param q
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetCIcoef_nocov <- function(df, groupZ, X, Y, MX, MY, q = qnorm(.975)^2, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  ## Calculate Components
  C0 <- A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
    2 * A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
    A1type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = Y, lpos = MX, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MY, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = Y, jpos = Y, kpos = X, lpos = MX, noisy = noisy)
  C1 <- -(A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
            2 * A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
            A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) +
            A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            2 * A1type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            A1type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_sum_nocov(df, groupZ, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy))
  C2 <- A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    2 * A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    A1type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_sum_nocov(df, groupZ, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy)
  PXY <- GetLM_nocov(df, X, Y, groupZ)
  PXX <- GetLM_nocov(df, X, X, groupZ)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Title
#'
#' @param df
#' @param P
#' @param G
#' @param X
#' @param Y
#' @param MX
#' @param MY
#' @param Z
#' @param W
#' @param q
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetCIcoef_iloop <- function(df, P, G, X, Y, MX, MY, Z, W, q = qnorm(.975)^2, noisy = FALSE) {
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  ## Calculate Components
  C0 <- A1type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
    2 * A2type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
    A3type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = Y, lpos = MX, noisy = noisy) -
    A4type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MY, noisy = noisy) -
    A5type_iloop_sum(df, P, G, ipos = Y, jpos = Y, kpos = X, lpos = MX, noisy = noisy)
  C1 <- -(A1type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) +
            2 * A2type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) +
            A3type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = Y, lpos = MX, noisy = noisy) -
            A4type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MY, noisy = noisy) -
            A5type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) +
            A1type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            2 * A2type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
            A3type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A4type_iloop_sum(df, P, G, ipos = X, jpos = Y, kpos = X, lpos = MX, noisy = noisy) -
            A5type_iloop_sum(df, P, G, ipos = Y, jpos = X, kpos = X, lpos = MX, noisy = noisy))
  C2 <- A1type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    2 * A2type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) +
    A3type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A4type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy) -
    A5type_iloop_sum(df, P, G, ipos = X, jpos = X, kpos = X, lpos = MX, noisy = noisy)

  ## Generate dPQ
  n <- nrow(df)
  Q <- cbind(Z, W)
  dPQ <- dPW <- rep(0, n)
  WW <- t(W) %*% W
  WWinv <- solve(WW)
  QQ <- t(Q) %*% Q
  QQinv <- solve(QQ)
  for (i in 1:n) {
    dPW[i] <- W[i, ] %*% WWinv %*% W[i, ]
    dPQ[i] <- Q[i, ] %*% QQinv %*% Q[i, ]
  }
  IdPW <- pmax(1 - dPW, .01)
  IdPQ <- pmax(1 - dPQ, .01)

  PXY <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, Y)
  PXX <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, X)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Title
#'
#' @param X
#' @param Y
#' @param P
#' @param q
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetCIcoef_iloop_nocov <- function(X, Y, P, q = qnorm(0.975)^2, noisy = FALSE) {
  n <- length(X)
  M <- diag(n) - P
  MM <- M * M

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - MM
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  MY <- M %*% Y
  MX <- M %*% X
  Poff <- P - diag(diag(P))
  XdM <- dM * X
  YdM <- dM * Y
  XX <- X * X
  XY <- X * Y
  YY <- Y * Y
  onesNMX <- tcrossprod(onesN, MX)
  onesNMY <- tcrossprod(onesN, MY)
  onesNX <- tcrossprod(onesN, X)
  onesNdM <- tcrossprod(onesN, dM)
  MXYY <- MX * YY
  MXXY <- MX * XY
  MYXX <- MY * XX
  MXX <- MX * X
  MYX <- MY * X
  MYY <- MY * Y
  MXY <- MX * Y
  XYdM <- dM * XY
  YYdM <- dM * YY
  XXdM <- dM * XX
  MXXX <- MX * XX
  MonesNX <- M * onesNX

  C0A1vec <- C0A2vec <- C0A3vec <- C0A4vec <- C0A5vec <- rep(0, n)
  C11A1vec <- C11A2vec <- C11A3vec <- C11A4vec <- C11A5vec <- rep(0, n)
  C12A1vec <- C12A2vec <- C12A3vec <- C12A4vec <- C12A5vec <- rep(0, n)

  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Pi2 <- Pi^2
    Mi2 <- Mi^2
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
    Di <- matrix(D2[, i], ncol = 1)
    D2D3i <- D2 / D3i
    D2D3i[i, ] <- 0
    D2D3i[, i] <- 0
    diag(D2D3i) <- 0
    recD3i <- 1 / D3i
    recD3i[i, ] <- 0
    recD3i[, i] <- 0
    diag(recD3i) <- 0
    recD2i <- matrix(1 / D2[, i], ncol = 1)
    recD2i[i] <- 0

    # Calculate things once
    MXi <- MX[i]
    MYi <- MY[i]
    Xi <- X[i]
    XiMii <- Xi * M[i, i]
    XPi <- X * Pi
    XMi <- X * Mi
    XMi2 <- XMi * Mi
    XPiMi <- XPi * Mi
    YPi <- Y * Pi
    YMi <- Y * Mi
    YMi2 <- YMi * Mi
    YPiMi <- YPi * Mi
    Pi2recD2i <- Pi2 * recD2i
    MrecD3i <- M * recD3i
    MMrecD3i <- recD3i * MM
    tXPiD2D3i <- crossprod(XPi, D2D3i)
    D2D3iXMi <- (D2D3i) %*% (XMi)
    Pi2XMi2 <- Pi2 * XMi2
    Pi2XMi <- Pi2 * XMi
    XMiMX <- XMi * MX
    dMXPi <- dM * XPi
    XMi2MY <- XMi2 * MY
    XMi2MX <- XMi2 * MX
    XMiMY <- XMi * MY
    XPiMX <- MX * XPi
    YPiMX <- YPi * MX
    YPidM <- YPi * dM
    XPiMY <- MY * XPi
    D2D3iMonesNX <- D2D3i * MonesNX
    Pi2XdMrecD2i <- Pi2recD2i * XdM
    Pi2YMirecD2i <- Pi2recD2i * YMi
    MiPi2recD2i <- Mi * Pi2recD2i
    Pi2YMi2recD2i <- Pi2recD2i * YMi2
    Pi2YdMrecD2i <- Pi2recD2i * YdM
    onesNdMrecD3i <- onesNdM * recD3i
    onesNMXrecD3i <- onesNMX * recD3i
    Pi2XMirecD2i <- Pi2XMi * recD2i
    Pi2YrecD2i <- Pi2recD2i * Y
    Pi2YdMMirecD2i <- Pi2YdMrecD2i * Mi
    tPi2XdMrecD2iMrecD3i <- crossprod(Pi2XdMrecD2i, MrecD3i)
    tPi2YdMrecD2iMrecD3i <- crossprod(Pi2YdMrecD2i, MrecD3i)
    tXPiMiMrecD3i <- crossprod(XPiMi, MrecD3i)
    tPi2XMi2recD2iMrecD3i <- crossprod(Pi2XMi2 * recD2i, MrecD3i)
    tMXXPiMirecD3i <- crossprod(MX * XPiMi, recD3i)
    tPi2YMirecD2iMMrecD3i <- crossprod(Pi2YMirecD2i, MMrecD3i)
    tPi2YMi2recD2iMrecD3i <- crossprod(Pi2YMi2recD2i, MrecD3i)
    tPi2YdMMirecD2ionesNdMrecD3i <- crossprod(Pi2YdMMirecD2i, onesNdMrecD3i)
    tPi2XMirecD2iMMrecD3i <- crossprod(Pi2XMirecD2i, MMrecD3i)
    tPi2XdMMirecD2irecD3i <- crossprod(Pi2XdMrecD2i * Mi, recD3i)
    tPi2recD2iXMi <- sum(Pi2recD2i * XMi)
    MrecD3iYPiMi <- MrecD3i %*% (YPiMi)
    tXPiD2D3iXPi <- tXPiD2D3i %*% (XPi)
    tPi2recD2iYMi <- sum(Pi2YMirecD2i)
    tPi2recD2iXXdM <- crossprod(Pi2recD2i, XXdM)
    tPi2YdMMirecD2irecD3i <- crossprod(Pi2YdMMirecD2i, recD3i)
    onesNMXrecD3iYPiMi <- (onesNMXrecD3i) %*% (YPiMi)
    D2D3iMonesNXMi <- (D2D3iMonesNX) %*% (Mi)
    MrecD3iXPiMi <- t(tXPiMiMrecD3i)
    Pi2recD2iMYX <- Pi2recD2i * MYX

    C0A11i <- tXPiD2D3iXPi * MYi
    C0A12i <- crossprod(MY * XPiMi, recD3i) %*% (dMXPi) -
      tXPiMiMrecD3i %*% (XPiMY)
    C0A13i <- crossprod(dMXPi, onesNMY * recD3i) %*% (XPiMi)
    C0A14i <- sum(XPiMY * MrecD3iXPiMi)
    C0A15i <- MYi * (tPi2recD2iXXdM) -
      sum(MiPi2recD2i * MYXX)

    C0A21i <- tXPiD2D3i %*% (YPi) * (MXi)
    C0A22i <- tMXXPiMirecD3i %*% (YPidM) -
      tXPiMiMrecD3i %*% (YPiMX)
    C0A23i <- crossprod(dMXPi, onesNMXrecD3iYPiMi)
    C0A24i <- crossprod(MX * XPi, MrecD3iYPiMi)
    C0A25i <- (MXi) * sum(Pi2recD2i * XYdM) -
      sum(MiPi2recD2i * MXXY)

    C0A31i <- crossprod(YPi, D2D3i) %*% (YPi) * (MXi)
    C0A32i <- crossprod(MX * YPiMi, recD3i) %*% (YPidM) -
      crossprod(YPiMi, MrecD3i) %*% (YPiMX)
    C0A33i <- crossprod(YPidM, onesNMXrecD3iYPiMi)
    C0A34i <- crossprod(YPiMX, MrecD3iYPiMi)
    C0A35i <- (MXi) * sum(Pi2recD2i * YYdM) -
      sum(MiPi2recD2i * MXYY)

    C0A41i <- crossprod(Pi2YdMrecD2i * MY, D2D3iXMi) -
      crossprod(Pi2YrecD2i * MY, D2D3iMonesNXMi)
    C0A42i <- tPi2YdMrecD2iMrecD3i %*% (XMi2) -
      tPi2YMirecD2iMMrecD3i %*% (XMi) -
      tPi2YdMMirecD2ionesNdMrecD3i %*% (XMi) +
      tPi2YMi2recD2iMrecD3i %*% (XdM)
    C0A43i <- tPi2YdMMirecD2irecD3i %*% (XMi2MY) -
      tPi2YMi2recD2iMrecD3i %*% (XMiMY) -
      M[i, i] * tPi2YdMrecD2iMrecD3i %*% (XMiMY) +
      M[i, i] * tPi2YMirecD2iMMrecD3i %*% (MYX)
    C0A44i <- XiMii * sum(Pi2recD2i * MYY) -
      Xi * (MYi) * (tPi2recD2iYMi)

    C0A51i <- crossprod(Pi2YdMrecD2i * MX, D2D3iXMi) -
      crossprod(Pi2YrecD2i * MX, D2D3iMonesNXMi)
    C0A52i <- C0A42i
    C0A53i <- tPi2YdMMirecD2irecD3i %*% (XMi2MX) -
      tPi2YMi2recD2iMrecD3i %*% (XMiMX) -
      M[i, i] * tPi2YdMrecD2iMrecD3i %*% (XMiMX) +
      M[i, i] * tPi2YMirecD2iMMrecD3i %*% (MXX)
    C0A54i <- XiMii * sum(Pi2recD2i * MXY) -
      Xi * (MXi) * (tPi2recD2iYMi)

    C11A41i <- crossprod(Pi2XdMrecD2i * MY, D2D3iXMi) -
      crossprod(Pi2recD2iMYX, D2D3iMonesNXMi)
    C11A42i <- tPi2XdMrecD2iMrecD3i %*% (XMi2) -
      tPi2XMirecD2iMMrecD3i %*% (XMi) -
      crossprod(Pi2 * XdM * Mi * recD2i, onesNdMrecD3i) %*% (XMi) +
      tPi2XMi2recD2iMrecD3i %*% (XdM)
    C11A43i <- tPi2XdMMirecD2irecD3i %*% (XMi2MY) -
      tPi2XMi2recD2iMrecD3i %*% (XMiMY) -
      M[i, i] * tPi2XdMrecD2iMrecD3i %*% (XMiMY) +
      M[i, i] * tPi2XMirecD2iMMrecD3i %*% (MYX)
    C11A44i <- XiMii * sum(Pi2recD2iMYX) -
      Xi * (MYi) * (tPi2recD2iXMi)

    C12A11i <- tXPiD2D3iXPi * (MXi)
    C12A12i <- tMXXPiMirecD3i %*% (dMXPi) -
      tXPiMiMrecD3i %*% (XPiMX)
    C12A13i <- crossprod(dMXPi, onesNMX * recD3i) %*% (XPiMi)
    C12A14i <- crossprod(XPiMX, MrecD3iXPiMi)
    C12A15i <- (MXi) * (tPi2recD2iXXdM) -
      sum(MiPi2recD2i * MXXX)

    C12A51i <- crossprod(Pi2XdMrecD2i * MX, D2D3iXMi) -
      crossprod(Pi2 * X * recD2i * MX, D2D3iMonesNXMi)
    C12A52i <- C11A42i
    C12A53i <- tPi2XdMMirecD2irecD3i %*% (XMi2MX) -
      tPi2XMi2recD2iMrecD3i %*% (XMiMX) -
      M[i, i] * tPi2XdMrecD2iMrecD3i %*% (XMiMX) +
      M[i, i] * tPi2XMirecD2iMMrecD3i %*% (MXX)
    C12A54i <- XiMii * sum(Pi2recD2i * MXX) -
      Xi * (MXi) * (tPi2recD2iXMi)

    C0A1vec[i] <- C0A11i - C0A12i - C0A13i + C0A14i + C0A15i
    C0A2vec[i] <- C0A21i - C0A22i - C0A23i + C0A24i + C0A25i
    C0A3vec[i] <- C0A31i - C0A32i - C0A33i + C0A34i + C0A35i
    C0A4vec[i] <- C0A41i + C0A42i * MYi + C0A43i + C0A44i
    C0A5vec[i] <- C0A51i + C0A52i * MXi + C0A53i + C0A54i

    C11A4vec[i] <- C11A41i + C11A42i * MYi + C11A43i + C11A44i
    C12A1vec[i] <- C12A11i - C12A12i - C12A13i + C12A14i + C12A15i
    C12A5vec[i] <- C12A51i + C12A52i * MXi + C12A53i + C12A54i

    if (noisy) {
      cat(i, "of", n, "done. ")
    }
  }

  C11A1vec <- C0A1vec
  C11A3vec <- C11A2vec <- C0A2vec
  C11A5vec <- C0A5vec
  C12A2vec <- C12A1vec
  C12A3vec <- C0A2vec
  C12A4vec <- C0A5vec
  C2A1vec <- C2A2vec <- C2A3vec <- C12A1vec
  C2A4vec <- C2A5vec <- C12A5vec

  C0 <- sum(C0A1vec * Y + 2 * C0A2vec * Y + C0A3vec * X - C0A4vec * X - C0A5vec * Y)
  C11 <- sum(C11A1vec * X + 2 * C11A2vec * X + C11A3vec * X - C11A4vec * X - C11A5vec * X)
  C12 <- sum(C12A1vec * Y + 2 * C12A2vec * Y + C12A3vec * X - C12A4vec * X - C12A5vec * Y)
  C1 <- -C11 - C12
  C2 <- sum(C2A1vec * X + 2 * C2A2vec * X + C2A3vec * X - C2A4vec * X - C12A5vec * X)

  PXY <- t(Y) %*% Poff %*% X
  PXX <- t(X) %*% Poff %*% X


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}


#' Title
#'
#' @param df
#' @param IdPW
#' @param IdPQ
#' @param dPW
#' @param dPQ
#' @param X
#' @param Y
#' @param MX
#' @param MY
#' @param q
#' @param noisyi
#' @param noisyj
#'
#' @returns
#' @export
#'
#' @examples
GetCIcoef_ijloop <- function(df, IdPW, IdPQ, dPW, dPQ, X, Y, MX, MY, q = qnorm(.975)^2, noisyi = FALSE, noisyj = FALSE) {
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  ## Calculate Components
  C0 <- A1type_ijloop_sum(df, ipos = Y, jpos = X, kpos = X, lpos = MY, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
    2 * A2type_ijloop_sum(df, ipos = Y, jpos = X, kpos = Y, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
    A3type_ijloop_sum(df, ipos = X, jpos = Y, kpos = Y, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
    A4type_ijloop_sum(df, ipos = X, jpos = Y, kpos = X, lpos = MY, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
    A5type_ijloop_sum(df, ipos = Y, jpos = Y, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj)
  C1 <- -(A1type_ijloop_sum(df, ipos = X, jpos = X, kpos = X, lpos = MY, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
            2 * A2type_ijloop_sum(df, ipos = X, jpos = X, kpos = Y, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
            A3type_ijloop_sum(df, ipos = X, jpos = X, kpos = Y, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
            A4type_ijloop_sum(df, ipos = X, jpos = X, kpos = X, lpos = MY, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
            A5type_ijloop_sum(df, ipos = X, jpos = Y, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
            A1type_ijloop_sum(df, ipos = Y, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
            2 * A2type_ijloop_sum(df, ipos = Y, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
            A3type_ijloop_sum(df, ipos = X, jpos = Y, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
            A4type_ijloop_sum(df, ipos = X, jpos = Y, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
            A5type_ijloop_sum(df, ipos = Y, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj))
  C2 <- A1type_ijloop_sum(df, ipos = X, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
    2 * A2type_ijloop_sum(df, ipos = X, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) +
    A3type_ijloop_sum(df, ipos = X, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
    A4type_ijloop_sum(df, ipos = X, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj) -
    A5type_ijloop_sum(df, ipos = X, jpos = X, kpos = X, lpos = MX, IdPQ, IdPW, noisyi = noisyi, noisyj = noisyj)
  PXY <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, Y)
  PXX <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, X)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Title
#'
#' @param X
#' @param Y
#' @param MX
#' @param MY
#' @param W
#' @param Q
#' @param WWWinv
#' @param QQQinv
#' @param IdPQ
#' @param IdPW
#' @param q
#' @param noisyi
#' @param noisyj
#'
#' @returns
#' @export
#'
#' @examples
GetCIcoef_ijloop_cov <- function(X, Y, MX, MY, W, Q, WWWinv, QQQinv, IdPQ, IdPW, q = qnorm(0.975)^2, noisyi = TRUE, noisyj = FALSE) {
  n <- length(X)
  onesN <- matrix(rep(1, n), ncol = 1)
  dPQ <- 1 - IdPQ
  dPW <- 1 - IdPW
  dM <- IdPQ
  dMW <- IdPW

  C0A1vec <- C0A2vec <- C0A3vec <- C0A4vec <- C0A5vec <- rep(0, n)
  C11A3vec <- C11A4vec <- rep(0, n)
  C12A1vec <- C12A2vec <- C12A3vec <- C12A4vec <- C12A5vec <- rep(0, n)
  C2A3vec <- C2A4vec <- rep(0, n)

  # Things that only needs to be calculated once
  MXY <- MX * Y
  MXX <- MX * X
  MYX <- MY * X
  dMX <- dM * X
  dMY <- dM * Y
  dMXX <- dMX * X
  dMXY <- dMX * Y
  MYY <- MY * Y
  dMYY <- dMY * Y
  dMYX <- dMY * X

  for (i in 1:n) {
    dPQi <- rep(0, n)
    dPQi[i] <- dPQ[i]
    dPWi <- rep(0, n)
    dPWi[i] <- dPW[i]

    Pi <- QQQinv %*% Q[i, ]
    PWi <- WWWinv %*% W[i, ]
    Girow <- (Pi - dPQi) / dM[i] - (PWi - dPWi) / dMW[i]
    Gicol <- (Pi - dPQi) / dM - (PWi - dPWi) / dMW
    Girow <- c(Girow)
    Gicol <- c(Gicol) # force vectors
    Mi <- -Pi
    Mi[i] <- 1 - Pi[i]
    D2i <- dM * dM[i] - Mi^2
    recD2i <- 1 / D2i
    recD2i[i] <- 0

    C0A11i <- C0A12i <- C0A13i <- C0A14i <- rep(0, n)
    C0A21i <- C0A22i <- C0A23i <- C0A24i <- rep(0, n)
    C0A31i <- C0A32i <- C0A33i <- C0A34i <- rep(0, n)
    C0A41i <- C0A42i <- C0A43i <- rep(0, n)
    C0A51i <- C0A52i <- C0A53i <- rep(0, n)

    C11A31i <- C11A32i <- C11A33i <- C11A34i <- rep(0, n)
    C11A41i <- C11A42i <- C11A43i <- rep(0, n)

    C12A11i <- C12A12i <- C12A13i <- C12A14i <- rep(0, n)
    C12A21i <- C12A22i <- C12A23i <- C12A24i <- rep(0, n)
    C12A31i <- C12A32i <- C12A33i <- C12A34i <- rep(0, n)
    C12A41i <- C12A42i <- C12A43i <- rep(0, n)
    C12A51i <- C12A52i <- C12A53i <- rep(0, n)

    C2A31i <- C2A32i <- C2A33i <- C2A34i <- rep(0, n)
    C2A41i <- C2A42i <- C2A43i <- rep(0, n)

    ## Save objects that must be calculated more than once
    XGirow <- X * Girow
    YGicol <- Y * Gicol
    XGicol <- X * Gicol
    MirecD2i <- Mi * recD2i
    dMXGirow <- dM * XGirow
    MYXGirow <- MY * XGirow
    MXXGirow <- MX * XGirow
    MXYGicol <- MX * YGicol
    MXXGicol <- MXX * Gicol
    XGirowMi <- XGirow * Mi
    dMYGicol <- dM * YGicol
    YGicolMi <- YGicol * Mi
    Gicol2 <- Gicol^2
    Girow2 <- Girow^2
    Gicol2Y <- Gicol2 * Y
    Gicol2X <- Gicol2 * X
    GirowGicol <- Girow * Gicol
    GirowGicolY <- GirowGicol * Y
    GirowGicolX <- GirowGicol * X
    Gicol2YMi2recD2i <- Gicol2Y * Mi^2 * recD2i
    Gicol2XMi2recD2i <- Gicol2X * Mi^2 * recD2i
    GirowGicolYMi2recD2i <- GirowGicolY * Mi^2 * recD2i
    GirowGicolXMi2recD2i <- GirowGicolX * Mi^2 * recD2i
    Gicol2XdMrecD2i <- Gicol2X * dM * recD2i
    Gicol2YMirecD2i <- Gicol2Y * MirecD2i
    Gicol2YdMrecD2i <- Gicol2Y * dM * recD2i
    MYi <- MY[i]
    MXi <- MX[i]
    Mii <- Mi[i]
    Xi <- X[i]
    MiX <- Mi * X
    MiY <- Mi * Y
    XiMii <- Xi * Mii
    XiMXi <- Xi * MXi
    XiMYi <- Xi * MYi
    Gicol2recD2i <- Gicol2 * recD2i
    Girow2MirecD2i <- Girow2 * MirecD2i
    Gicol2MirecD2i <- Gicol2 * MirecD2i
    GirowGicolrecD2i <- GirowGicol * recD2i
    Gicol2YdMMi <- Gicol2Y * dM * Mi
    GirowGicolYdMrecD2i <- GirowGicolY * dM * recD2i
    GirowGicolYdMMirecD2i <- GirowGicolYdMrecD2i * Mi
    Gicol2XMirecD2i <- Gicol2X * MirecD2i
    Gicol2XdMMirecD2i <- Gicol2XMirecD2i * dM
    Gicol2YrecD2i <- Gicol2Y * recD2i
    Gicol2XrecD2i <- Gicol2X * recD2i
    Girow2recD2i <- Girow2 * recD2i
    GirowGicolMirecD2i <- GirowGicol * MirecD2i
    tGirow2recD2idMXX <- sum(Girow2recD2i * dMXX)
    GirowGicolXdMMirecD2i <- GirowGicolX * dM * MirecD2i
    GirowGicolXrecD2iMX <- GirowGicolX * recD2i * MX
    GirowGicolMirecD2iMXX <- GirowGicolMirecD2i * MXX
    Gicol2MirecD2iMXX <- Gicol2MirecD2i * MXX
    Girow2MirecD2iX <- Girow2MirecD2i * X
    Gicol2MirecD2iMXY <- Gicol2MirecD2i * MXY
    tGicol2recD2iMiX <- sum(Gicol2recD2i * (MiX))
    tGicol2recD2iMiY <- sum(Gicol2recD2i * (MiY))



    for (j in (1:n)[-i]) {
      Pj <- QQQinv %*% Q[j, ]
      Mj <- -Pj
      Mj[j] <- 1 - Pj[j]
      D2j <- dM * dM[j] - Mj^2

      D3ij <- Mii * D2j - (dM * Mi[j]^2 + Mi^2 * dM[j] - 2 * Mj * Mi * Mi[j])
      recD3ij <- 1 / D3ij
      recD3ij[i] <- 0
      recD3ij[j] <- 0
      D2D3ij <- D2j / D3ij
      D2D3ij[i] <- 0
      D2D3ij[j] <- 0

      ## Save objects that must be calculated more than once
      MXj <- MX[j]
      MYj <- MY[j]
      Xj <- X[j]
      Yj <- Y[j]
      Mij <- Mi[j]
      dMj <- dM[j]
      MXjXj <- MXj * Xj
      dMjXj <- dMj * Xj
      Girowj <- Girow[j]
      Gicolj <- Gicol[j]
      GirowjXj <- Girowj * Xj
      GicoljYj <- Gicolj * Yj
      tXGirowD2D3ij <- sum(XGirow * D2D3ij)
      tYGicolD2D3ij <- sum(YGicol * D2D3ij)
      MjrecD3ij <- Mj * recD3ij
      recD3ijMj2 <- MjrecD3ij * Mj
      tXGirowMiMjrecD3ij <- sum(XGirowMi * MjrecD3ij)
      MijXj <- Mij * Xj
      GicoljMijXj <- Gicolj * MijXj
      GicoljMijYj <- Gicolj * Mij * Yj
      Mij2Xj <- Mij^2 * Xj
      Mij2XjMXj <- Mij2Xj * MXj
      MijXjMYj <- MijXj * MYj
      MijXjMXj <- MijXj * MXj
      onesNrecD3ij <- onesN * recD3ij
      onesNMXjrecD3ij <- onesNrecD3ij * MXj
      MirecD3ij <- Mi * recD3ij
      GirowjMijXj <- Girowj * MijXj
      D2D3ijMjonesNXj <- D2D3ij * Mj * onesN * Xj
      onesNdMjrecD3ij <- onesNrecD3ij * dMj
      MYD2D3ij <- MY * D2D3ij
      GicoljXj <- Gicolj * Xj
      GicoljXjdMj <- GicoljXj * dMj
      GicoljXjMXj <- GicoljXj * MXj
      GicoljXjMXi <- GicoljXj * MXi
      MXD2D3ijMjonesNXj <- MX * D2D3ijMjonesNXj
      MYD2D3ijMjonesNXj <- MY * D2D3ijMjonesNXj
      tGicol2YdMrecD2iMjrecD3ij <- sum(Gicol2YdMrecD2i * MjrecD3ij)
      tGicol2YMi2recD2iMjrecD3ij <- sum(Gicol2YMi2recD2i * MjrecD3ij)
      tMXXGirowMjrecD3ij <- sum(MXXGirow * MjrecD3ij)
      tMXYGicolMjrecD3ij <- sum(MXYGicol * MjrecD3ij)
      tGicol2YMirecD2irecD3ijMj2 <- sum(Gicol2YMirecD2i * recD3ijMj2)
      tdMXGirowonesNMXjrecD3ij <- sum(dMXGirow * onesNMXjrecD3ij)
      tGicol2YdMMirecD2irecD3ij <- sum(Gicol2YdMMi * recD2i * recD3ij)
      tGirowGicolXMi2recD2iMjrecD3ij <- sum(GirowGicolXMi2recD2i * MjrecD3ij)
      tGicol2YdMMirecD2ionesNdMjrecD3ij <- sum(Gicol2YdMMi * recD2i * onesNdMjrecD3ij)
      tYGicolMiMjrecD3ij <- sum(YGicolMi * MjrecD3ij)
      tGirowGicolXdMrecD2iMjrecD3ij <- sum(GirowGicolX * dM * recD2i * MjrecD3ij)
      tMXXGirowMirecD3ij <- sum(MXXGirow * MirecD3ij)
      tGirowGicolYMirecD2irecD3ijMj2 <- sum(GirowGicolY * Mi * recD2i * recD3ijMj2)
      tGirowGicolYMi2recD2iMjrecD3ij <- sum(GirowGicolYMi2recD2i * MjrecD3ij)
      tGirowGicolYdMrecD2iMjrecD3ij <- sum(GirowGicolYdMrecD2i * MjrecD3ij)
      tGicol2XdMrecD2iMjrecD3ij <- sum(Gicol2XdMrecD2i * MjrecD3ij)
      tGicol2XMirecD2irecD3ijMj2 <- sum(Gicol2XMirecD2i * recD3ijMj2)
      tGicol2XMi2recD2iMjrecD3ij <- sum(Gicol2XMi2recD2i * MjrecD3ij)
      tGirowGicolXMirecD2irecD3ijMj2 <- sum(GirowGicolX * Mi * recD2i * recD3ijMj2)
      tGicol2XdMMirecD2ionesNdMjrecD3ij <- sum(Gicol2XdMMirecD2i * onesNdMjrecD3ij)
      tdMXGicolonesNMXjrecD3ij <- sum(dM * XGicol * onesNMXjrecD3ij)
      tXGicolMiMjrecD3ij <- sum(XGicol * Mi * MjrecD3ij)
      tMXYGicolMirecD3ij <- sum(MXYGicol * MirecD3ij)
      tGicol2XdMMirecD2irecD3ij <- sum(Gicol2XdMMirecD2i * recD3ij)
      tMXXGicolMjrecD3ij <- sum(MXXGicol * MjrecD3ij)
      dMYGicolonesNMXjrecD3ij <- sum(dMYGicol * onesNMXjrecD3ij)
      tXGicolD2D3ij <- sum(XGicol * D2D3ij)
      tMXXGicolMirecD3ij <- sum(MXXGicol * MirecD3ij)

      ## C0 components
      C0A11i[j] <- (tXGirowD2D3ij) * (GirowjXj) * MYi
      C0A12i[j] <- sum(MYXGirow * MirecD3ij) * (GirowjXj * dMj) -
        tXGirowMiMjrecD3ij * (GirowjXj * MYj)
      C0A13i[j] <- sum(dMXGirow * onesN * MYj * recD3ij) * (GirowjMijXj)
      C0A14i[j] <- sum(MYXGirow * MjrecD3ij) * (GirowjMijXj)

      C0A21i[j] <- (tXGirowD2D3ij) * (GicoljYj) * MXi
      C0A22i[j] <- tMXXGirowMirecD3ij * (GicoljYj * dMj) -
        tXGirowMiMjrecD3ij * (GicoljYj * MXj)
      C0A23i[j] <- tdMXGirowonesNMXjrecD3ij * (GicoljMijYj)
      C0A24i[j] <- tMXXGirowMjrecD3ij * (GicoljMijYj)

      C0A31i[j] <- (tYGicolD2D3ij) * (GicoljYj) * MXi
      C0A32i[j] <- tMXYGicolMirecD3ij * (GicoljYj * dMj) -
        tYGicolMiMjrecD3ij * (GicoljYj * MXj)
      C0A33i[j] <- dMYGicolonesNMXjrecD3ij * (GicoljMijYj)
      C0A34i[j] <- tMXYGicolMjrecD3ij * (GicoljMijYj)

      C0A41i[j] <- sum(Gicol2YdMrecD2i * MYD2D3ij) * (MijXj) -
        sum(Gicol2YrecD2i * MYD2D3ijMjonesNXj) * Mij
      C0A42i[j] <- tGicol2YdMrecD2iMjrecD3ij * (Mij2Xj) -
        tGicol2YMirecD2irecD3ijMj2 * (MijXj) -
        tGicol2YdMMirecD2ionesNdMjrecD3ij * (MijXj) +
        tGicol2YMi2recD2iMjrecD3ij * (dMjXj)
      C0A43i[j] <- tGicol2YdMMirecD2irecD3ij * (Mij2Xj * MYj) -
        tGicol2YMi2recD2iMjrecD3ij * (MijXjMYj) -
        Mii * tGicol2YdMrecD2iMjrecD3ij * (MijXjMYj) +
        Mii * tGicol2YMirecD2irecD3ijMj2 * (MYj * Xj)

      C0A51i[j] <- sum(GirowGicolYdMrecD2i * MX * D2D3ij) * (MijXj) -
        sum(GirowGicolY * recD2i * MXD2D3ijMjonesNXj) * Mij
      C0A52i[j] <- tGirowGicolYdMrecD2iMjrecD3ij * (Mij2Xj) -
        tGirowGicolYMirecD2irecD3ijMj2 * (MijXj) -
        sum(GirowGicolYdMMirecD2i * onesNdMjrecD3ij) * (MijXj) +
        tGirowGicolYMi2recD2iMjrecD3ij * (dMjXj)
      C0A53i[j] <- sum(GirowGicolYdMMirecD2i * recD3ij) * (Mij2XjMXj) -
        tGirowGicolYMi2recD2iMjrecD3ij * (MijXjMXj) -
        Mii * tGirowGicolYdMrecD2iMjrecD3ij * (MijXjMXj) +
        Mii * tGirowGicolYMirecD2irecD3ijMj2 * (MXjXj)


      ## C11 components (note that some are already calculated and hence omitted)
      C11A31i[j] <- tXGicolD2D3ij * (GicoljYj) * MXi
      C11A32i[j] <- tMXXGicolMirecD3ij * (GicoljYj * dMj) -
        tXGicolMiMjrecD3ij * (GicoljYj * MXj)
      C11A33i[j] <- tdMXGicolonesNMXjrecD3ij * (GicoljMijYj)
      C11A34i[j] <- tMXXGicolMjrecD3ij * (GicoljMijYj)

      C11A41i[j] <- sum(Gicol2XdMrecD2i * MYD2D3ij) * (MijXj) -
        sum(Gicol2XrecD2i * MYD2D3ijMjonesNXj) * Mij
      C11A42i[j] <- tGicol2XdMrecD2iMjrecD3ij * (Mij2Xj) -
        tGicol2XMirecD2irecD3ijMj2 * (MijXj) -
        tGicol2XdMMirecD2ionesNdMjrecD3ij * (MijXj) +
        tGicol2XMi2recD2iMjrecD3ij * (dMjXj)
      C11A43i[j] <- tGicol2XdMMirecD2irecD3ij * (Mij2Xj * MYj) -
        tGicol2XMi2recD2iMjrecD3ij * (MijXjMYj) -
        Mii * tGicol2XdMrecD2iMjrecD3ij * (MijXjMYj) +
        Mii * tGicol2XMirecD2irecD3ijMj2 * (MYj * Xj)

      ## C12 components
      C12A11i[j] <- (tXGirowD2D3ij) * (GirowjXj) * MXi
      C12A12i[j] <- tMXXGirowMirecD3ij * (GirowjXj * dMj) -
        tXGirowMiMjrecD3ij * (GirowjXj * MXj)
      C12A13i[j] <- tdMXGirowonesNMXjrecD3ij * (GirowjMijXj)
      C12A14i[j] <- tMXXGirowMjrecD3ij * (GirowjMijXj)

      C12A21i[j] <- (tXGirowD2D3ij) * GicoljXjMXi
      C12A22i[j] <- tMXXGirowMirecD3ij * (GicoljXjdMj) -
        tXGirowMiMjrecD3ij * (GicoljXjMXj)
      C12A23i[j] <- tdMXGirowonesNMXjrecD3ij * (GicoljMijXj)
      C12A24i[j] <- tMXXGirowMjrecD3ij * (GicoljMijXj)

      C12A31i[j] <- (tYGicolD2D3ij) * GicoljXjMXi
      C12A32i[j] <- tMXYGicolMirecD3ij * (GicoljXjdMj) -
        tYGicolMiMjrecD3ij * (GicoljXjMXj)
      C12A33i[j] <- dMYGicolonesNMXjrecD3ij * (GicoljMijXj)
      C12A34i[j] <- tMXYGicolMjrecD3ij * (GicoljMijXj)

      C12A41i[j] <- sum(Gicol2YdMrecD2i * MX * D2D3ij) * (MijXj) -
        sum(Gicol2YrecD2i * MXD2D3ijMjonesNXj) * Mij
      C12A42i[j] <- tGicol2YdMrecD2iMjrecD3ij * (Mij2Xj) -
        tGicol2YMirecD2irecD3ijMj2 * (MijXj) -
        tGicol2YdMMirecD2ionesNdMjrecD3ij * (MijXj) +
        tGicol2YMi2recD2iMjrecD3ij * (dMjXj)
      C12A43i[j] <- tGicol2YdMMirecD2irecD3ij * (Mij2XjMXj) -
        tGicol2YMi2recD2iMjrecD3ij * (MijXjMXj) -
        Mii * tGicol2YdMrecD2iMjrecD3ij * (MijXjMXj) +
        Mii * tGicol2YMirecD2irecD3ijMj2 * (MXjXj)

      C12A51i[j] <- sum(GirowGicolXrecD2iMX * dM * D2D3ij) * (MijXj) -
        sum(GirowGicolXrecD2iMX * D2D3ijMjonesNXj) * Mij
      C12A52i[j] <- tGirowGicolXdMrecD2iMjrecD3ij * (Mij2Xj) -
        tGirowGicolXMirecD2irecD3ijMj2 * (MijXj) -
        sum(GirowGicolXdMMirecD2i * onesNdMjrecD3ij) * (MijXj) +
        tGirowGicolXMi2recD2iMjrecD3ij * (dMjXj)
      C12A53i[j] <- sum(GirowGicolXdMMirecD2i * recD3ij) * (Mij2XjMXj) -
        tGirowGicolXMi2recD2iMjrecD3ij * (MijXjMXj) -
        Mii * tGirowGicolXdMrecD2iMjrecD3ij * (MijXjMXj) +
        Mii * tGirowGicolXMirecD2irecD3ijMj2 * (MXjXj)

      ## C2 components
      C2A31i[j] <- tXGicolD2D3ij * GicoljXjMXi
      C2A32i[j] <- tMXXGicolMirecD3ij * (GicoljXjdMj) -
        tXGicolMiMjrecD3ij * (GicoljXjMXj)
      C2A33i[j] <- tdMXGicolonesNMXjrecD3ij * (GicoljMijXj)
      C2A34i[j] <- tMXXGicolMjrecD3ij * (GicoljMijXj)

      C2A41i[j] <- sum(Gicol2XdMrecD2i * MX * D2D3ij) * (MijXj) -
        sum(Gicol2XrecD2i * MXD2D3ijMjonesNXj) * Mij
      C2A42i[j] <- tGicol2XdMrecD2iMjrecD3ij * (Mij2Xj) -
        tGicol2XMirecD2irecD3ijMj2 * (MijXj) -
        tGicol2XdMMirecD2ionesNdMjrecD3ij * (MijXj) +
        tGicol2XMi2recD2iMjrecD3ij * (dMjXj)
      C2A43i[j] <- tGicol2XdMMirecD2irecD3ij * (Mij2XjMXj) -
        tGicol2XMi2recD2iMjrecD3ij * (MijXjMXj) -
        Mii * tGicol2XdMrecD2iMjrecD3ij * (MijXjMXj) +
        Mii * tGicol2XMirecD2irecD3ijMj2 * (MXjXj)

      if (noisyj == TRUE & j %% 100 == 0) cat("j:", j / n, " ")
    }

    C0A15i <- MYi * tGirow2recD2idMXX -
      sum(Girow2MirecD2iX * MYX)
    C0A25i <- MXi * sum(GirowGicolrecD2i * dMXY) -
      sum(GirowGicolMirecD2iMXX * Y)
    C0A35i <- MXi * sum(Gicol2recD2i * dMYY) -
      sum(Gicol2MirecD2iMXY * Y)
    C0A44i <- XiMii * sum(Gicol2recD2i * MYY) -
      XiMYi * tGicol2recD2iMiY
    C0A54i <- XiMii * sum(GirowGicolrecD2i * MXY) -
      XiMXi * sum(GirowGicolrecD2i * MiY)

    C11A35i <- MXi * sum(Gicol2recD2i * dMXY) -
      sum(Gicol2MirecD2iMXX * Y)
    C11A44i <- XiMii * sum(Gicol2recD2i * MYX) -
      XiMYi * tGicol2recD2iMiX

    C12A15i <- MXi * tGirow2recD2idMXX -
      sum(Girow2MirecD2iX * MXX)
    C12A25i <- MXi * sum(GirowGicolrecD2i * dMXX) -
      sum(GirowGicolMirecD2iMXX * X)
    C12A35i <- MXi * sum(Gicol2recD2i * dMYX) -
      sum(Gicol2MirecD2iMXY * X)
    C12A44i <- XiMii * sum(Gicol2recD2i * MXY) -
      XiMXi * tGicol2recD2iMiY
    C12A54i <- XiMii * sum(GirowGicolrecD2i * MXX) -
      XiMXi * sum(GirowGicolrecD2i * MiX)

    C2A35i <- MXi * sum(Gicol2recD2i * dMXX) -
      sum(Gicol2MirecD2iMXX * X)
    C2A44i <- XiMii * sum(Gicol2recD2i * MXX) -
      XiMXi * tGicol2recD2iMiX


    C0A1vec[i] <- sum(C0A11i - C0A12i - C0A13i + C0A14i) + C0A15i
    C0A2vec[i] <- sum(C0A21i - C0A22i - C0A23i + C0A24i) + C0A25i
    C0A3vec[i] <- sum(C0A31i - C0A32i - C0A33i + C0A34i) + C0A35i
    C0A4vec[i] <- sum(C0A41i + C0A43i) + sum(C0A42i) * MYi + C0A44i
    C0A5vec[i] <- sum(C0A51i + C0A53i) + sum(C0A52i) * MXi + C0A54i

    C11A3vec[i] <- sum(C11A31i - C11A32i - C11A33i + C11A34i) + C11A35i
    C11A4vec[i] <- sum(C11A41i + C11A43i) + sum(C11A42i) * MYi + C11A44i

    C12A1vec[i] <- sum(C12A11i - C12A12i - C12A13i + C12A14i) + C12A15i
    C12A2vec[i] <- sum(C12A21i - C12A22i - C12A23i + C12A24i) + C12A25i
    C12A3vec[i] <- sum(C12A31i - C12A32i - C12A33i + C12A34i) + C12A35i
    C12A4vec[i] <- sum(C12A41i + C12A43i) + sum(C12A42i) * MXi + C12A44i
    C12A5vec[i] <- sum(C12A51i + C12A53i) + sum(C12A52i) * MXi + C12A54i

    C2A3vec[i] <- sum(C2A31i - C2A32i - C2A33i + C2A34i) + C2A35i
    C2A4vec[i] <- sum(C2A41i + C2A43i) + sum(C2A42i) * MXi + C2A44i

    if (noisyi == TRUE) cat("i:", i / n, " ")
  }

  C0 <- sum(C0A1vec * Y + 2 * C0A2vec * Y + C0A3vec * X - C0A4vec * X - C0A5vec * Y)
  C11 <- sum(C0A1vec * X + 2 * C0A2vec * X + C11A3vec * X - C11A4vec * X - C0A5vec * X)
  C12 <- sum(C12A1vec * Y + 2 * C12A2vec * Y + C12A3vec * X - C12A4vec * X - C12A5vec * Y)
  C1 <- -C11 - C12
  C2 <- sum(C12A1vec * X + 2 * C12A2vec * X + C2A3vec * X - C2A4vec * X - C12A5vec * X)

  PXY <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, Y)
  PXX <- GetLM_WQ(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, X)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}


#' Title
#'
#' @param CIcoef
#'
#' @returns
#' @export
#'
#' @examples
GetCIvals <- function(CIcoef) {
  acon <- CIcoef[1]
  bcon <- CIcoef[2]
  ccon <- CIcoef[3]
  det <- CIcoef[2]^2 - 4 * CIcoef[1] * CIcoef[3]
  CILB <- (-bcon - sqrt(det)) / (2 * acon)
  CIUB <- (-bcon + sqrt(det)) / (2 * acon)
  c(CILB, CIUB)
}

#' Title
#'
#' @param CIcoef
#'
#' @returns
#' @export
#'
#' @examples
GetCItypebd <- function(CIcoef) {
  det <- CIcoef[2]^2 - 4 * CIcoef[1] * CIcoef[3]
  if (CIcoef[1] >= 0 & det >= 0) {
    CItype <- 1 # convex interval
    CIbounds <- GetCIvals(CIcoef)
  } else if (CIcoef[1] < 0 & det >= 0) {
    CItype <- 2 # donut
    CIbounds <- GetCIvals(CIcoef)
  } else if (CIcoef[1] < 0 & det < 0) {
    CItype <- 3 # accept everything
    CIbounds <- c(-100, 100)
  } else {
    CItype <- 4 # reject everything
    CIbounds <- c(NA, NA)
  }
  c(CItype, CIbounds)
}

#' Title
#'
#' @param df
#' @param groupW
#' @param groupQ
#' @param X
#' @param Y
#' @param MX
#' @param MY
#' @param q
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetL3OCIcoef_fast <- function(df, groupW, groupQ, X, Y, MX, MY, q = qnorm(.975)^2, noisy = FALSE) {
  df$X <- eval(substitute(X), df)
  df$Y <- eval(substitute(Y), df)
  df$MX <- eval(substitute(MX), df)
  df$MY <- eval(substitute(MY), df)
  df$groupW <- eval(substitute(groupW), df)
  df$groupQ <- eval(substitute(groupQ), df)

  C0A11vecs <- C0A12vecs <- C0A13vecs <- C0A14vecs <- C0A15vecs <- rep(0, max(df$groupQ))
  C0A21vecs <- C0A22vecs <- C0A23vecs <- C0A24vecs <- C0A25vecs <- rep(0, max(df$groupQ))
  C0A31vecs <- C0A32vecs <- C0A33vecs <- C0A34vecs <- C0A35vecs <- rep(0, max(df$groupQ))
  C0A41vecs <- C0A42vecs <- C0A43vecs <- C0A44vecs <- rep(0, max(df$groupQ))
  C0A51vecs <- C0A52vecs <- C0A53vecs <- C0A54vecs <- rep(0, max(df$groupQ))
  C11A11vecs <- C11A12vecs <- C11A13vecs <- C11A14vecs <- C11A15vecs <- rep(0, max(df$groupQ))
  C11A21vecs <- C11A22vecs <- C11A23vecs <- C11A24vecs <- C11A25vecs <- rep(0, max(df$groupQ))
  C11A41vecs <- C11A42vecs <- C11A43vecs <- C11A44vecs <- rep(0, max(df$groupQ))
  C11A51vecs <- C11A52vecs <- C11A53vecs <- C11A54vecs <- rep(0, max(df$groupQ))
  C12A11vecs <- C12A12vecs <- C12A13vecs <- C12A14vecs <- C12A15vecs <- rep(0, max(df$groupQ))
  C12A31vecs <- C12A32vecs <- C12A33vecs <- C12A34vecs <- C12A35vecs <- rep(0, max(df$groupQ))
  C12A51vecs <- C12A52vecs <- C12A53vecs <- C12A54vecs <- rep(0, max(df$groupQ))
  C2A11vecs <- C2A12vecs <- C2A13vecs <- C2A14vecs <- C2A15vecs <- rep(0, max(df$groupQ))
  C2A41vecs <- C2A42vecs <- C2A43vecs <- C2A44vecs <- rep(0, max(df$groupQ))

  iteration <- 1

  # outer loop
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]

    # ds <- ds %>%  group_by(groupQ) %>% mutate(numingrp = length(ds$groupQ))
    # ds <- ds[ds$numingrp>=3,]

    ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
    PW <- ZWmat / (length(ds$groupW))
    if (nrow(ds) == 0) {
      for (g in unique(ds$groupQ)) {
        C0A11vecs[g] <- C0A12vecs[g] <- C0A13vecs[g] <- C0A14vecs[g] <- C0A15vecs[g] <- 0
        C0A21vecs[g] <- C0A22vecs[g] <- C0A23vecs[g] <- C0A24vecs[g] <- C0A25vecs[g] <- 0
        C0A31vecs[g] <- C0A32vecs[g] <- C0A33vecs[g] <- C0A34vecs[g] <- C0A35vecs[g] <- 0
        C0A41vecs[g] <- C0A42vecs[g] <- C0A43vecs[g] <- C0A44vecs[g] <- 0
        C0A51vecs[g] <- C0A52vecs[g] <- C0A53vecs[g] <- C0A54vecs[g] <- 0
        C11A11vecs[g] <- C11A12vecs[g] <- C11A13vecs[g] <- C11A14vecs[g] <- C11A15vecs[g] <- 0
        C11A21vecs[g] <- C11A22vecs[g] <- C11A23vecs[g] <- C11A24vecs[g] <- C11A25vecs[g] <- 0
        C11A41vecs[g] <- C11A42vecs[g] <- C11A43vecs[g] <- C11A44vecs[g] <- 0
        C11A51vecs[g] <- C11A52vecs[g] <- C11A53vecs[g] <- C11A54vecs[g] <- 0
        C12A11vecs[g] <- C12A12vecs[g] <- C12A13vecs[g] <- C12A14vecs[g] <- C12A15vecs[g] <- 0
        C12A31vecs[g] <- C12A32vecs[g] <- C12A33vecs[g] <- C12A34vecs[g] <- C12A35vecs[g] <- 0
        C12A51vecs[g] <- C12A52vecs[g] <- C12A53vecs[g] <- C12A54vecs[g] <- 0
        C2A11vecs[g] <- C2A12vecs[g] <- C2A13vecs[g] <- C2A14vecs[g] <- C2A15vecs[g] <- 0
        C2A41vecs[g] <- C2A42vecs[g] <- C2A43vecs[g] <- C2A44vecs[g] <- 0
      }
    } else {
      ZQ <- matrix(0, nrow = length(ds$groupQ), ncol = length(unique(ds$groupQ)))
      ds$groupidx <- Getgroupindex(ds, groupQ)
      ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

      # calculate values specific to this subset
      Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))
      Ps <- PQ
      Ms <- diag(nrow(ds)) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      # Stuff to calculate once
      XMY <- ds$X * ds$MY
      XMX <- ds$X * ds$MX
      YMX <- ds$Y * ds$MX
      YMY <- ds$Y * ds$MY
      XdMs <- ds$X * dMs
      YdMs <- ds$Y * dMs
      XMYX <- XMY * ds$X
      ones <- matrix(rep(1, nrow(ds)), ncol = 1)
      onestdMs <- ones %*% t(dMs)
      onestMX <- (ones %*% t(ds$MX))
      onestMY <- (ones %*% t(ds$MY))
      onesX <- ones %x% t(ds$X)

      # inner loop
      for (g in unique(ds$groupQ)) {
        repidx <- min(which(ds$groupQ == g)) # representative index
        Pis <- matrix(Ps[, repidx], ncol = 1)
        Pgs <- ifelse(Pis == 0, 0, 1) %*% matrix(1, ncol = length(Pis), nrow = 1)
        Gis <- matrix(Gs[, repidx], ncol = 1)
        Gis[repidx, 1] <- Gis[repidx + 1, 1] # Put the G back
        Mis <- matrix(-Ps[, repidx], ncol = 1)
        recD2is <- matrix(recD2s[, repidx], ncol = 1)
        recD2is[repidx, 1] <- recD2is[repidx + 1, 1] # Put the G back
        D3is <- Ms[repidx, repidx] * D2s - (dMs %*% t(Mis)^2 + Mis^2 %*% t(dMs) - 2 * Ms * (Mis %*% t(Mis)))
        D2D3is <- D2s / D3is
        diag(D2D3is) <- 0
        recD3is <- 1 / D3is
        diag(recD3is) <- 0
        # Mes <- matrix(ds$lpos,ncol=1)

        ## Things to calculate only once
        XGis <- ds$X * Gis
        XMis <- ds$X * Mis
        YGis <- ds$Y * Gis
        YMis <- ds$Y * Mis
        YMis2 <- YMis * Mis
        XMis2 <- XMis * Mis
        Gs2 <- Gs^2
        Gis2 <- Gis^2
        Mis2 <- Mis^2
        Ms2 <- Ms^2
        GisMis <- Gis * Mis
        YGisMis <- YGis * Mis
        XGisMis <- XGis * Mis
        Mis2X <- Mis2 * ds$X
        XGisdMs <- XGis * dMs
        YGisdMs <- YGis * dMs
        XMYGis <- XMY * Gis
        XMYGisMis <- XMYGis * Mis
        Gis2YdMs <- Gis2 * YdMs
        Gis2XdMs <- Gis2 * XdMs
        recD2isMX <- ds$MX * recD2is
        YGisMX <- YGis * ds$MX
        Mis2XMX <- Mis2X * ds$MX
        XGisMX <- XGis * ds$MX
        MisXMX <- Mis * XMX
        Gis2XdMsMisrecD2is <- Gis2XdMs * Mis * recD2is
        Gis2YdMsMisrecD2is <- Gis2YdMs * Mis * recD2is
        ivectomatsY <- ivectomats(ds, ds$Y, g)
        ivectomatsX <- ivectomats(ds, ds$X, g)
        ivectomatsXdMs <- ivectomats(ds, XdMs, g)
        ivectomatsYMX <- ivectomats(ds, YMX, g)
        ivectomatsXMY <- ivectomats(ds, XMY, g)
        ivectomatsXMX <- ivectomats(ds, XMX, g)
        D2D3isivectomatsX <- D2D3is * ivectomatsX
        D2D3isivectomatsXMY <- D2D3is * ivectomatsXMY
        D2D3isivectomatsXMX <- D2D3is * ivectomatsXMX
        recD3isivectomatsY <- recD3is * ivectomatsY
        recD3isMsivectomatsY <- recD3isivectomatsY * Ms
        recD3isivectomatsX <- recD3is * ivectomatsX
        recD3isMsivectomatsX <- recD3isivectomatsX * Ms
        recD3isonestMXivectomatsX <- recD3is * onestMX * ivectomatsX
        recD2sGs2 <- recD2s * Gs2
        recD3isMs <- recD3is * Ms
        recD3isMs2 <- recD3isMs * Ms
        recD3isMsivectomatsXdMs <- recD3isMs * ivectomatsXdMs
        recD3isMsivectomatsXMY <- recD3isMs * ivectomatsXMY
        recD3isMs2ivectomatsXMY <- recD3isMs2 * ivectomatsXMY
        recD3isMsivectomatsYMX <- recD3isMs * ivectomatsYMX
        Gs2MsrecD2sPgs <- Gs2 * Ms * recD2s * Pgs
        recD2sGs2Pgs <- recD2sGs2 * Pgs
        D2D3isMsonesX <- D2D3is * Ms * (onesX)
        D2D3isMsonesXivectomatsX <- D2D3isMsonesX * ivectomatsX
        recD3isonestdMs <- recD3is * (onestdMs)
        recD3isonestdMsivectomatsXMY <- recD3isonestdMs * ivectomatsXMY
        recD3isonestdMsivectomatsXMX <- recD3isonestdMs * ivectomatsXMX

        # Calculations for every groupQ
        C0A11vecs[g] <- t(XGis) %*% (D2D3is * ivectomats(ds, YMY, g)) %*% (XGis)
        C0A12vecs[g] <- t(XMYGisMis) %*% (recD3isivectomatsY) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsY) %*% (Gis * XMY)
        C0A13vecs[g] <- t(XdMs * Gis) %*% (recD3isivectomatsY * onestMY) %*% (XGisMis)
        C0A14vecs[g] <- t(XMYGis) %*% (recD3isMsivectomatsY) %*% (XGisMis)
        C0A15vecs[g] <- t(YMY) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$Y) %*% (Gs2MsrecD2sPgs) %*% (XMYX)

        C0A21vecs[g] <- t(XGis) %*% (D2D3is * ivectomatsYMX) %*% (YGis)
        C0A22vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsY) %*% (YGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsY) %*% (YGisMX)
        C0A23vecs[g] <- t(XdMs * Gis) %*% (recD3isivectomatsY * onestMX) %*% (YGisMis)
        C0A24vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsY) %*% (YGisMis)
        C0A25vecs[g] <- t(YMX) %*% (recD2sGs2Pgs) %*% (ds$X * YdMs) -
          t(ds$Y) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$Y)

        C0A31vecs[g] <- t(YGis) %*% (D2D3isivectomatsXMX) %*% (YGis)
        C0A32vecs[g] <- t(YMX * GisMis) %*% (recD3isivectomatsX) %*% (YGisdMs) -
          t(YGisMis) %*% (recD3isMsivectomatsX) %*% (YGisMX)
        C0A33vecs[g] <- t(YdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (YGisMis)
        C0A34vecs[g] <- t(YMX * Gis) %*% (recD3isMsivectomatsX) %*% (YGisMis)
        C0A35vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$Y * YdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (YMX * ds$Y)

        C0A41vecs[g] <- t(Gis2YdMs * ds$MY * recD2is) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * YMY * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C0A42vecs[g] <- t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsXMY) %*% (Mis2X) -
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2ivectomatsXMY) %*% (XMis) -
          t(Gis2YdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMY) %*% (XMis) +
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsXMY) %*% (dMs * ds$X)
        C0A43vecs[g] <- t(Gis2YdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2 * XMY) -
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (Mis * XMY) -
          t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (Mis * XMY) +
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMY)
        C0A44vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (YMY) -
          t(ds$X * XMY) %*% (recD2sGs2Pgs) %*% (YMis)

        C0A51vecs[g] <- t(Gis2YdMs * recD2isMX) %*% (D2D3is * ivectomatsY) %*% (XMis) -
          t(Gis2 * YMX * recD2is) %*% (D2D3isMsonesX * ivectomatsY) %*% (Mis)
        C0A52vecs[g] <- t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsYMX) %*% (Mis2X) -
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsYMX) %*% (XMis) -
          t(Gis2YdMsMisrecD2is) %*% (recD3isonestdMs * ivectomatsYMX) %*% (XMis) +
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsYMX) %*% (dMs * ds$X)
        C0A53vecs[g] <- t(Gis2YdMsMisrecD2is) %*% (recD3isivectomatsY) %*% (Mis2XMX) -
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsY) %*% (MisXMX) -
          t(Gis2YdMs * recD2is) %*% (recD3isMs * ivectomats(ds, YdMs, g)) %*% (MisXMX) +
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomats(ds, YdMs, g)) %*% (XMX)
        C0A54vecs[g] <- t(ds$Y * XdMs) %*% (recD2sGs2Pgs) %*% (YMX) -
          t(ds$Y * XMX) %*% (recD2sGs2Pgs) %*% (YMis)

        C11A11vecs[g] <- t(XGis) %*% (D2D3isivectomatsXMY) %*% (XGis)
        C11A12vecs[g] <- t(XMYGisMis) %*% (recD3isivectomatsX) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsX) %*% (Gis * XMY)
        C11A13vecs[g] <- t(XdMs * Gis) %*% (recD3is * onestMY * ivectomatsX) %*% (XGisMis)
        C11A14vecs[g] <- t(XMYGis) %*% (recD3isMsivectomatsX) %*% (XGisMis)
        C11A15vecs[g] <- t(XMY) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (XMYX)

        C11A21vecs[g] <- t(XGis) %*% (D2D3isivectomatsXMX) %*% (YGis)
        C11A22vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsX) %*% (YGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsX) %*% (YGisMX)
        C11A23vecs[g] <- t(XdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (YGisMis)
        C11A24vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsX) %*% (YGisMis)
        C11A25vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$X * YdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$Y)

        C11A41vecs[g] <- t(Gis2XdMs * ds$MY * recD2is) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * XMY * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C11A42vecs[g] <- t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsXMY) %*% (Mis2X) -
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2ivectomatsXMY) %*% (XMis) -
          t(Gis2XdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMY) %*% (XMis) +
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsXMY) %*% (dMs * ds$X)
        C11A43vecs[g] <- t(Gis2XdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2 * XMY) -
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (Mis * XMY) -
          t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (Mis * XMY) +
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMY)
        C11A44vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (XMY) -
          t(ds$X * XMY) %*% (recD2sGs2Pgs) %*% (XMis)

        C11A51vecs[g] <- t(Gis2YdMs * recD2isMX) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * YMX * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C11A52vecs[g] <- t(Gis2YdMs * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (Mis2X) -
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsXMX) %*% (XMis) -
          t(Gis2YdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMX) %*% (XMis) +
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (dMs * ds$X)
        C11A53vecs[g] <- t(Gis2YdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2XMX) -
          t(Gis2 * YMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (MisXMX) -
          t(Gis2YdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (MisXMX) +
          t(Gis2 * YMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMX)
        C11A54vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (YMX) -
          t(ds$X * XMX) %*% (recD2sGs2Pgs) %*% (YMis)

        C12A11vecs[g] <- t(XGis) %*% (D2D3is * ivectomatsYMX) %*% (XGis)
        C12A12vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsY) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsY) %*% (XGisMX)
        C12A13vecs[g] <- t(XdMs * Gis) %*% (recD3is * onestMX * ivectomatsY) %*% (XGisMis)
        C12A14vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsY) %*% (XGisMis)
        C12A15vecs[g] <- t(YMX) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$Y) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$X)

        C12A31vecs[g] <- t(YGis) %*% (D2D3isivectomatsXMX) %*% (XGis)
        C12A32vecs[g] <- t(YMX * GisMis) %*% (recD3isivectomatsX) %*% (XGisdMs) -
          t(YGisMis) %*% (recD3isMsivectomatsX) %*% (XGisMX)
        C12A33vecs[g] <- t(YdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (XGisMis)
        C12A34vecs[g] <- t(YMX * Gis) %*% (recD3isMsivectomatsX) %*% (XGisMis)
        C12A35vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$Y * XdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (YMX * ds$X)

        C12A51vecs[g] <- t(Gis2XdMs * recD2isMX) %*% (D2D3is * ivectomatsY) %*% (XMis) -
          t(Gis2 * XMX * recD2is) %*% (D2D3isMsonesX * ivectomatsY) %*% (Mis)
        C12A52vecs[g] <- t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsYMX) %*% (Mis2X) -
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsYMX) %*% (XMis) -
          t(Gis2XdMsMisrecD2is) %*% (recD3isonestdMs * ivectomatsYMX) %*% (XMis) +
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsYMX) %*% (dMs * ds$X)
        C12A53vecs[g] <- t(Gis2XdMsMisrecD2is) %*% (recD3isivectomatsY) %*% (Mis2XMX) -
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsY) %*% (MisXMX) -
          t(Gis2XdMs * recD2is) %*% (recD3isMs * ivectomats(ds, YdMs, g)) %*% (MisXMX) +
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomats(ds, YdMs, g)) %*% (XMX)
        C12A54vecs[g] <- t(ds$Y * XdMs) %*% (recD2sGs2Pgs) %*% (XMX) -
          t(ds$Y * XMX) %*% (recD2sGs2Pgs) %*% (XMis)

        C2A11vecs[g] <- t(XGis) %*% (D2D3isivectomatsXMX) %*% (XGis)
        C2A12vecs[g] <- t(XMX * GisMis) %*% (recD3isivectomatsX) %*% (XGisdMs) -
          t(XGisMis) %*% (recD3isMsivectomatsX) %*% (XGisMX)
        C2A13vecs[g] <- t(XdMs * Gis) %*% (recD3isonestMXivectomatsX) %*% (XGisMis)
        C2A14vecs[g] <- t(XMX * Gis) %*% (recD3isMsivectomatsX) %*% (XGisMis)
        C2A15vecs[g] <- t(XMX) %*% (recD2sGs2Pgs) %*% (ds$X * XdMs) -
          t(ds$X) %*% (Gs2MsrecD2sPgs) %*% (XMX * ds$X)

        C2A41vecs[g] <- t(Gis2XdMs * recD2isMX) %*% (D2D3isivectomatsX) %*% (XMis) -
          t(Gis2 * XMX * recD2is) %*% (D2D3isMsonesXivectomatsX) %*% (Mis)
        C2A42vecs[g] <- t(Gis2XdMs * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (Mis2X) -
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsXMX) %*% (XMis) -
          t(Gis2XdMsMisrecD2is) %*% (recD3isonestdMsivectomatsXMX) %*% (XMis) +
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMs * ivectomatsXMX) %*% (dMs * ds$X)
        C2A43vecs[g] <- t(Gis2XdMsMisrecD2is) %*% (recD3isivectomatsX) %*% (Mis2XMX) -
          t(Gis2 * XMis2 * recD2is) %*% (recD3isMsivectomatsX) %*% (MisXMX) -
          t(Gis2XdMs * recD2is) %*% (recD3isMsivectomatsXdMs) %*% (MisXMX) +
          t(Gis2 * XMis * recD2is) %*% (recD3isMs2 * ivectomatsXdMs) %*% (XMX)
        C2A44vecs[g] <- t(ds$X * XdMs) %*% (recD2sGs2Pgs) %*% (XMX) -
          t(ds$X * XMX) %*% (recD2sGs2Pgs) %*% (XMis)
      }
    }
    if (noisy) {
      cat(iteration, "of", max(df$groupW), "done. ")
      iteration <- iteration + 1
    }
  }

  C0A1 <- C0A11vecs - C0A12vecs - C0A13vecs + C0A14vecs + C0A15vecs
  C0A2 <- C0A21vecs - C0A22vecs - C0A23vecs + C0A24vecs + C0A25vecs
  C0A3 <- C0A31vecs - C0A32vecs - C0A33vecs + C0A34vecs + C0A35vecs
  C0A4 <- C0A41vecs + C0A42vecs + C0A43vecs + C0A44vecs
  C0A5 <- C0A51vecs + C0A52vecs + C0A53vecs + C0A54vecs


  C11A1 <- C11A11vecs - C11A12vecs - C11A13vecs + C11A14vecs + C11A15vecs
  C11A2 <- C11A21vecs - C11A22vecs - C11A23vecs + C11A24vecs + C11A25vecs
  C11A4 <- C11A41vecs + C11A42vecs + C11A43vecs + C11A44vecs
  C11A5 <- C11A51vecs + C11A52vecs + C11A53vecs + C11A54vecs

  C12A1 <- C12A11vecs - C12A12vecs - C12A13vecs + C12A14vecs + C12A15vecs
  C12A3 <- C12A31vecs - C12A32vecs - C12A33vecs + C12A34vecs + C12A35vecs
  C12A4 <- C11A5
  C12A5 <- C12A51vecs + C12A52vecs + C12A53vecs + C12A54vecs

  C2A1 <- C2A11vecs - C2A12vecs - C2A13vecs + C2A14vecs + C2A15vecs
  C2A4 <- C2A41vecs + C2A42vecs + C2A43vecs + C2A44vecs

  C0 <- sum(C0A1 + 2 * C0A2 + C0A3 - C0A4 - C0A5)
  C11 <- C11A1 + 3 * C11A2 - C11A4 - C11A5
  C12 <- 3 * C12A1 + C12A3 - C12A4 - C12A5
  C1 <- -sum(C11 + C12)
  C2 <- sum(4 * C2A1 - 2 * C2A4)

  PXY <- GetLM(df, X, Y, groupW, group)
  PXX <- GetLM(df, X, X, groupW, group)


  acon <- PXX^2 - q * C2
  bcon <- -2 * PXY * PXX - q * C1
  ccon <- PXY^2 - q * C0

  c(acon, bcon, ccon)
}

#' Title
#'
#' @param df
#' @param IdPW
#' @param IdPQ
#' @param dPW
#' @param dPQ
#' @param W
#' @param Q
#' @param X
#' @param Y
#'
#' @returns
#' @export
#'
#' @examples
GetLM_WQ <- function(df, IdPW, IdPQ, dPW, dPQ, W, Q, X, Y) {
  df$Xpos <- eval(substitute(X), df)
  df$Ypos <- eval(substitute(Y), df)

  WW <- t(W) %*% W
  WWinv <- solve(WW)
  QQ <- t(Q) %*% Q
  QQinv <- solve(QQ)
  QX <- t(Q) %*% df$Xpos
  QY <- t(Q) %*% df$Ypos
  WX <- t(W) %*% df$Xpos
  WY <- t(W) %*% df$Ypos

  QXdQ <- t(Q) %*% (df$Xpos / IdPQ)
  QYdQ <- t(Q) %*% (df$Ypos / IdPQ)
  WXdW <- t(W) %*% (df$Xpos / IdPW)
  WYdW <- t(W) %*% (df$Ypos / IdPW)


  (t(QXdQ) %*% QQinv %*% QY - sum(df$Xpos * dPQ / IdPQ * df$Ypos)) -
    (t(WXdW) %*% WWinv %*% WY - sum(df$Xpos * dPW / IdPW * df$Ypos))
}

#' Title
#'
#' @param df
#' @param X
#' @param e
#' @param groupW
#' @param group
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetLM <- function(df, X, e, groupW, group, noisy = FALSE) {
  LMvecs <- rep(0, max(df$groupW))
  df$Xpos <- eval(substitute(X), df)
  df$epos <- eval(substitute(e), df)
  df$groupW <- eval(substitute(groupW), df)
  df$group <- eval(substitute(group), df)

  iteration <- 1
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]
    ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds, group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    # ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))
    ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)
    # PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)
    PW <- ZWmat / (length(ds$groupW))

    # calculate values specific to this subset
    Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
      diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))
    LMvecs[s] <- t(ds$epos) %*% Gs %*% ds$Xpos

    if (noisy) {
      cat(iteration, "of", length(unique(df$groupW)), "done. ")
      iteration <- iteration + 1
    }
  }
  sum(LMvecs)
}


#' Title
#'
#' @param df
#' @param X
#' @param e
#' @param groupZ
#' @param noisy
#'
#' @returns
#' @export
#'
#' @examples
GetLM_nocov <- function(df, X, e, groupZ, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  LMvecs <- rep(0, max(df$groupZ))
  df$Xpos <- eval(substitute(X), df)
  df$epos <- eval(substitute(e), df)

  iteration <- 1
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ == s, ]
    ZQ <- matrix(0, nrow = length(ds$groupZ), ncol = length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds, groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ)))
    LMvecs[s] <- t(ds$epos) %*% Gs %*% ds$Xpos

    if (noisy) {
      cat(iteration, "of", max(df$groupZ), "done. ")
      iteration <- iteration + 1
    }
  }
  sum(LMvecs)
}
