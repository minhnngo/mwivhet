#' Title
#'
#' @param df
#' @param P
#' @param G
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param noisy
#'
#' @returns
#' @noRd
#'
#' @examples
A1type_iloop_sum <- function(df, P, G, ipos, jpos, kpos, lpos, noisy = FALSE) {
  # A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0,max(df$group))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))

  A1vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Gi <- matrix(G[i, ], ncol = 1) # force column vector
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi^2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
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
    Poffi <- Poff[, i]
    Goffi <- Goff[i, ]
    lposmx <- matrix(df$lpos, ncol = 1)

    A11i <- t(df$jpos * Gi) %*% D2D3i %*% (Gi * df$kpos) * (df$lpos[i])
    A12i <- t(df$lpos * df$jpos * Gi * Mi) %*% recD3i %*% (Goffi * df$kpos * dM) -
      t(df$jpos * Gi * Mi) %*% (recD3i * M) %*% (Goffi * df$kpos * df$lpos)
    A13i <- t(dM * df$jpos * Gi) %*% ((onesN %x% t(lposmx)) * recD3i) %*% (Gi * Mi * df$kpos)
    A14i <- t((df$lpos) * df$jpos * Gi) %*% (M * recD3i) %*% (Gi * Mi * df$kpos)
    A15i <- (df$lpos[i]) * (t(Goffi^2 * recD2i) %*% (dM * df$jpos * df$kpos)) -
      (t(Goffi^2 * Mi * recD2i) %*% (df$lpos * df$jpos * df$kpos))


    A1vec[i] <- A11i - A12i - A13i + A14i + A15i

    if (noisy) {
      if (i %% 10 == 0) cat(i / n, " ")
    }
  }

  sum(A1vec * df$ipos)
}

#' Title
#'
#' @param df
#' @param P
#' @param G
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param noisy
#'
#' @returns
#' @noRd
#'
#' @examples
A2type_iloop_sum <- function(df, P, G, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))

  A2vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Girow <- matrix(G[i, ], ncol = 1) # force column vector
    Gicol <- matrix(G[, i], ncol = 1) # force column vector
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi^2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
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
    Poffi <- Poff[, i]
    Gofficol <- Goff[, i]
    Goffirow <- matrix(Goff[i, ], ncol = 1)
    lposmx <- matrix(df$lpos, ncol = 1)

    A21i <- t(df$jpos * Girow) %*% D2D3i %*% (Gicol * df$kpos) * (df$lpos[i])
    A22i <- t(df$lpos * df$jpos * Girow * Mi) %*% recD3i %*% (Gicol * df$kpos * dM) -
      t(df$jpos * Girow * Mi) %*% (recD3i * M) %*% (Gicol * df$kpos * df$lpos)
    A23i <- t(dM * df$jpos * Girow) %*% ((onesN %x% t(lposmx)) * recD3i) %*% (Gicol * Mi * df$kpos)
    A24i <- t((df$lpos) * df$jpos * Girow) %*% (M * recD3i) %*% (Gicol * Mi * df$kpos)
    A25i <- (df$lpos[i]) * (t(Gofficol * Goffirow * recD2i) %*% (dM * df$jpos * df$kpos)) -
      (t(Gofficol * Goffirow * Mi * recD2i) %*% (df$lpos * df$jpos * df$kpos))


    A2vec[i] <- A21i - A22i - A23i + A24i + A25i

    if (noisy) {
      if (i %% 10 == 0) cat(i / n, " ")
    }
  }

  sum(A2vec * df$ipos)
}

#' Title
#'
#' @param df
#' @param P
#' @param G
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param noisy
#'
#' @returns
#' @noRd
#'
#' @examples
A3type_iloop_sum <- function(df, P, G, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))

  A3vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Gi <- matrix(G[, i], ncol = 1) # force column vector
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi^2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
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
    Poffi <- Poff[, i]
    Goffi <- Goff[, i]
    lposmx <- matrix(df$lpos, ncol = 1)

    A31i <- t(df$jpos * Gi) %*% D2D3i %*% (Gi * df$kpos) * (df$lpos[i])
    A32i <- t(df$lpos * df$jpos * Gi * Mi) %*% recD3i %*% (Goffi * df$kpos * dM) -
      t(df$jpos * Gi * Mi) %*% (recD3i * M) %*% (Goffi * df$kpos * df$lpos)
    A33i <- t(dM * df$jpos * Gi) %*% ((onesN %x% t(lposmx)) * recD3i) %*% (Gi * Mi * df$kpos)
    A34i <- t((df$lpos) * df$jpos * Gi) %*% (M * recD3i) %*% (Gi * Mi * df$kpos)
    A35i <- (df$lpos[i]) * (t(Goffi^2 * recD2i) %*% (dM * df$jpos * df$kpos)) -
      (t(Goffi^2 * Mi * recD2i) %*% (df$lpos * df$jpos * df$kpos))


    A3vec[i] <- A31i - A32i - A33i + A34i + A35i

    if (noisy) {
      if (i %% 10 == 0) cat(i / n, " ")
    }
  }

  sum(A3vec * df$ipos)
}

#' Title
#'
#' @param df
#' @param P
#' @param G
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param noisy
#'
#' @returns
#' @noRd
#'
#' @examples
A4type_iloop_sum <- function(df, P, G, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))

  A4vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Gi <- matrix(G[, i], ncol = 1) # force column vector
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi^2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
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
    Poffi <- Poff[, i]
    Goffi <- Goff[, i]
    lposmx <- matrix(df$lpos, ncol = 1)
    kposmx <- matrix(df$kpos, ncol = 1)

    A41i <- t(Gi^2 * df$jpos * dM * recD2i * df$lpos) %*% (D2D3i) %*% (Mi * df$kpos) -
      t(Gi^2 * df$jpos * recD2i * df$lpos) %*% (D2D3i * M * (onesN %x% t(kposmx))) %*% (Mi)
    A42i <- t(Gi^2 * df$jpos * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * df$kpos) -
      t(Gi^2 * df$jpos * Mi * recD2i) %*% (recD3i * M * M) %*% (Mi * df$kpos) -
      t(Gi^2 * df$jpos * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (Mi * df$kpos) +
      t(Gi^2 * df$jpos * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * df$kpos)
    A43i <- t(Gi^2 * df$jpos * dM * Mi * recD2i) %*% (recD3i) %*% (Mi^2 * df$kpos * df$lpos) -
      t(Gi^2 * df$jpos * Mi^2 * recD2i) %*% (recD3i * M) %*% (Mi * df$kpos * df$lpos) -
      M[i, i] * t(Gi^2 * df$jpos * dM * recD2i) %*% (recD3i * M) %*% (Mi * df$kpos * df$lpos) +
      M[i, i] * t(Gi^2 * df$jpos * Mi * recD2i) %*% (recD3i * M * M) %*% (df$lpos * df$kpos)
    A44i <- df$kpos[i] * M[i, i] * (t(Goffi^2 * recD2i) %*% (df$lpos * df$jpos)) -
      df$kpos[i] * (df$lpos[i]) * (t(Goffi^2 * recD2i) %*% (Mi * df$jpos))


    A4vec[i] <- A41i + A42i * df$lpos[i] + A43i + A44i

    if (noisy) {
      if (i %% 10 == 0) cat(i / n, " ")
    }
  }
  ret <- (A4vec * df$ipos)

  sum(ret)
}

#' Title
#'
#' @param df
#' @param P
#' @param G
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param noisy
#'
#' @returns
#' @noRd
#'
#' @examples
A5type_iloop_sum <- function(df, P, G, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))

  A5vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Girow <- matrix(G[i, ], ncol = 1) # force column vector
    Gicol <- matrix(G[, i], ncol = 1) # force column vector
    D3i <- M[i, i] * D2 - (dM %*% t(Mi)^2 + Mi^2 %*% t(dM) - 2 * M * (Mi %*% t(Mi)))
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
    Poffi <- Poff[, i]
    Gofficol <- Goff[, i]
    Goffirow <- matrix(Goff[i, ], ncol = 1)
    lposmx <- matrix(df$lpos, ncol = 1)
    kposmx <- matrix(df$kpos, ncol = 1)

    A51i <- t(Girow * Gicol * df$jpos * dM * recD2i * df$lpos) %*% (D2D3i) %*% (Mi * df$kpos) -
      t(Girow * Gicol * df$jpos * recD2i * df$lpos) %*% (D2D3i * M * (onesN %x% t(kposmx))) %*% (Mi)
    A52i <- t(Girow * Gicol * df$jpos * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * df$kpos) -
      t(Girow * Gicol * df$jpos * Mi * recD2i) %*% (recD3i * M * M) %*% (Mi * df$kpos) -
      t(Girow * Gicol * df$jpos * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (Mi * df$kpos) +
      t(Girow * Gicol * df$jpos * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * df$kpos)
    A53i <- t(Girow * Gicol * df$jpos * dM * Mi * recD2i) %*% (recD3i) %*% (Mi^2 * df$kpos * df$lpos) -
      t(Girow * Gicol * df$jpos * Mi^2 * recD2i) %*% (recD3i * M) %*% (Mi * df$kpos * df$lpos) -
      M[i, i] * t(Girow * Gicol * df$jpos * dM * recD2i) %*% (recD3i * M) %*% (Mi * df$kpos * df$lpos) +
      M[i, i] * t(Girow * Gicol * df$jpos * Mi * recD2i) %*% (recD3i * M * M) %*% (df$lpos * df$kpos)
    A54i <- df$kpos[i] * M[i, i] * (t(Gofficol * Goffirow * recD2i) %*% (df$lpos * df$jpos)) -
      df$kpos[i] * (df$lpos[i]) * (t(Gofficol * Goffirow * recD2i) %*% (Mi * df$jpos))


    A5vec[i] <- A51i + A52i * df$lpos[i] + A53i + A54i

    if (noisy) {
      if (i %% 10 == 0) cat(i / n, " ")
    }
  }
  ret <- (A5vec * df$ipos)

  sum(ret)
}

#' Title
#'
#' @param df
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param IdPQ
#' @param IdPW
#' @param noisyi
#' @param noisyj
#'
#' @returns
#' @noRd
#'
#' @examples
A1type_ijloop_sum <- function(df, ipos, jpos, kpos, lpos, IdPQ, IdPW, noisyi = FALSE, noisyj = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  # M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1, n), ncol = 1)
  QQQinv <- Q %*% QQinv
  WWWinv <- W %*% WWinv


  A1vec <- rep(0, n)
  for (i in 1:n) {
    dPQi <- rep(0, n)
    dPQi[i] <- dPQ[i]
    dPWi <- rep(0, n)
    dPWi[i] <- dPW[i]

    Pi <- QQQinv %*% Q[i, ]
    Gi <- (QQQinv %*% Q[i, ] - dPQi) / dM[i] -
      (WWWinv %*% W[i, ] - dPWi) / dMW[i]
    Mi <- -Pi
    Mi[i] <- 1 - Pi[i]
    D2i <- dM * dM[i] - Mi^2
    recD2i <- 1 / D2i
    recD2i[i] <- 0

    A11i <- A12i <- A13i <- A14i <- rep(0, n)

    for (j in (1:n)[-i]) {
      dPQj <- rep(0, n)
      dPQj[j] <- dPQ[j]
      dPWj <- rep(0, n)
      dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j, ]
      Mj <- -Pj
      Mj[j] <- 1 - Pj[j]
      D2j <- dM * dM[j] - Mj^2

      D3ij <- Mi[i] * D2j - (dM * Mi[j]^2 + Mi^2 * dM[j] - 2 * Mj * Mi * Mi[j])
      recD3ij <- 1 / D3ij
      recD3ij[i] <- 0
      recD3ij[j] <- 0
      D2D3ij <- D2j / D3ij
      D2D3ij[i] <- 0
      D2D3ij[j] <- 0

      A11ij <- (t(df$jpos * Gi) %*% D2D3ij) * (Gi[j] * df$kpos[j]) * (df$lpos[i])
      A12ij <- (t(df$lpos * df$jpos * Gi * Mi) %*% recD3ij) * (Gi[j] * df$kpos[j] * dM[j]) -
        (t(df$jpos * Gi * Mi) %*% (recD3ij * Mj)) %*% (Gi[j] * df$kpos[j] * df$lpos[j])
      A13ij <- (t(dM * df$jpos * Gi) %*% (onesN * df$lpos[j] * recD3ij)) * (Gi[j] * Mi[j] * df$kpos[j])
      A14ij <- (t((df$lpos) * df$jpos * Gi) %*% (Mj * recD3ij)) * (Gi[j] * Mi[j] * df$kpos[j])

      A11i[j] <- A11ij
      A12i[j] <- A12ij
      A13i[j] <- A13ij
      A14i[j] <- A14ij

      if (noisyj == TRUE & j %% 10 == 0) cat("j:", j / n, " ")
    }

    A1vec[i] <- sum(A11i - A12i - A13i + A14i)
    A15i <- (df$lpos[i]) * (t(Gi^2 * recD2i) %*% (dM * df$jpos * df$kpos)) -
      (t(Gi^2 * Mi * recD2i) %*% (df$lpos * df$jpos * df$kpos))

    A1vec[i] <- A1vec[i] + A15i

    if (noisyi == TRUE) cat("i:", i / n, " ")
  }

  sum(A1vec * df$ipos)
}

#' Title
#'
#' @param df
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param IdPQ
#' @param IdPW
#' @param noisyi
#' @param noisyj
#'
#' @returns
#' @noRd
#'
#' @examples
A2type_ijloop_sum <- function(df, ipos, jpos, kpos, lpos, IdPQ, IdPW, noisyi = FALSE, noisyj = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  # M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1, n), ncol = 1)
  QQQinv <- Q %*% QQinv
  WWWinv <- W %*% WWinv


  A2vec <- rep(0, n)
  for (i in 1:n) {
    dPQi <- rep(0, n)
    dPQi[i] <- dPQ[i]
    dPWi <- rep(0, n)
    dPWi[i] <- dPW[i]

    Pi <- QQQinv %*% Q[i, ]
    Girow <- (QQQinv %*% Q[i, ] - dPQi) / dM[i] -
      (WWWinv %*% W[i, ] - dPWi) / dMW[i]
    Gicol <- t(QQQinv %*% Q[i, ] - dPQi) / dM -
      t(WWWinv %*% W[i, ] - dPWi) / dMW

    Mi <- -Pi
    Mi[i] <- 1 - Pi[i]
    D2i <- dM * dM[i] - Mi^2
    recD2i <- 1 / D2i
    recD2i[i] <- 0

    A21i <- A22i <- A23i <- A24i <- rep(0, n)

    for (j in (1:n)[-i]) {
      dPQj <- rep(0, n)
      dPQj[j] <- dPQ[j]
      dPWj <- rep(0, n)
      dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j, ]
      Mj <- -Pj
      Mj[j] <- 1 - Pj[j]
      D2j <- dM * dM[j] - Mj^2

      D3ij <- Mi[i] * D2j - (dM * Mi[j]^2 + Mi^2 * dM[j] - 2 * Mj * Mi * Mi[j])
      recD3ij <- 1 / D3ij
      recD3ij[i] <- 0
      recD3ij[j] <- 0
      D2D3ij <- D2j / D3ij
      D2D3ij[i] <- 0
      D2D3ij[j] <- 0

      A21ij <- (t(df$jpos * Girow) %*% D2D3ij) * (Gicol[j] * df$kpos[j]) * (df$lpos[i])
      A22ij <- (t(df$lpos * df$jpos * Girow * Mi) %*% recD3ij) * (Gicol[j] * df$kpos[j] * dM[j]) -
        (t(df$jpos * Girow * Mi) %*% (recD3ij * Mj)) %*% (Gicol[j] * df$kpos[j] * df$lpos[j])
      A23ij <- (t(dM * df$jpos * Girow) %*% (onesN * df$lpos[j] * recD3ij)) * (Gicol[j] * Mi[j] * df$kpos[j])
      A24ij <- (t((df$lpos) * df$jpos * Girow) %*% (Mj * recD3ij)) * (Gicol[j] * Mi[j] * df$kpos[j])

      A21i[j] <- A21ij
      A22i[j] <- A22ij
      A23i[j] <- A23ij
      A24i[j] <- A24ij

      if (noisyj == TRUE & j %% 10 == 0) cat("j:", j / n, " ")
    }

    A2vec[i] <- sum(A21i - A22i - A23i + A24i)
    A25i <- (df$lpos[i]) * (t(Girow * Gicol * recD2i) %*% (dM * df$jpos * df$kpos)) -
      (t(Girow * Gicol * Mi * recD2i) %*% (df$lpos * df$jpos * df$kpos))

    A2vec[i] <- A2vec[i] + A25i

    if (noisyi == TRUE) cat("i:", i / n, " ")
  }

  sum(A2vec * df$ipos)
}

#' Title
#'
#' @param df
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param IdPQ
#' @param IdPW
#' @param noisyi
#' @param noisyj
#'
#' @returns
#' @noRd
#'
#' @examples
A3type_ijloop_sum <- function(df, ipos, jpos, kpos, lpos, IdPQ, IdPW, noisyi = FALSE, noisyj = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)
  # M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1, n), ncol = 1)
  QQQinv <- Q %*% QQinv
  WWWinv <- W %*% WWinv


  A3vec <- rep(0, n)
  for (i in 1:n) {
    dPQi <- rep(0, n)
    dPQi[i] <- dPQ[i]
    dPWi <- rep(0, n)
    dPWi[i] <- dPW[i]

    Pi <- QQQinv %*% Q[i, ]
    Gi <- (QQQinv %*% Q[i, ] - dPQi) / dM -
      (WWWinv %*% W[i, ] - dPWi) / dMW
    Mi <- -Pi
    Mi[i] <- 1 - Pi[i]
    D2i <- dM * dM[i] - Mi^2
    recD2i <- 1 / D2i
    recD2i[i] <- 0

    A31i <- A32i <- A33i <- A34i <- rep(0, n)

    for (j in (1:n)[-i]) {
      dPQj <- rep(0, n)
      dPQj[j] <- dPQ[j]
      dPWj <- rep(0, n)
      dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j, ]
      Mj <- -Pj
      Mj[j] <- 1 - Pj[j]
      D2j <- dM * dM[j] - Mj^2

      D3ij <- Mi[i] * D2j - (dM * Mi[j]^2 + Mi^2 * dM[j] - 2 * Mj * Mi * Mi[j])
      recD3ij <- 1 / D3ij
      recD3ij[i] <- 0
      recD3ij[j] <- 0
      D2D3ij <- D2j / D3ij
      D2D3ij[i] <- 0
      D2D3ij[j] <- 0

      A31ij <- (t(df$jpos * Gi) %*% D2D3ij) * (Gi[j] * df$kpos[j]) * (df$lpos[i])
      A32ij <- (t(df$lpos * df$jpos * Gi * Mi) %*% recD3ij) * (Gi[j] * df$kpos[j] * dM[j]) -
        (t(df$jpos * Gi * Mi) %*% (recD3ij * Mj)) %*% (Gi[j] * df$kpos[j] * df$lpos[j])
      A33ij <- (t(dM * df$jpos * Gi) %*% (onesN * df$lpos[j] * recD3ij)) * (Gi[j] * Mi[j] * df$kpos[j])
      A34ij <- (t((df$lpos) * df$jpos * Gi) %*% (Mj * recD3ij)) * (Gi[j] * Mi[j] * df$kpos[j])

      A31i[j] <- A31ij
      A32i[j] <- A32ij
      A33i[j] <- A33ij
      A34i[j] <- A34ij

      if (noisyj == TRUE & j %% 10 == 0) cat("j:", j / n, " ")
    }

    A3vec[i] <- sum(A31i - A32i - A33i + A34i)
    A35i <- (df$lpos[i]) * (t(Gi^2 * recD2i) %*% (dM * df$jpos * df$kpos)) -
      (t(Gi^2 * Mi * recD2i) %*% (df$lpos * df$jpos * df$kpos))

    A3vec[i] <- A3vec[i] + A35i

    if (noisyi == TRUE) cat("i:", i / n, " ")
  }

  sum(A3vec * df$ipos)
}

#' Title
#'
#' @param df
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param IdPQ
#' @param IdPW
#' @param noisyi
#' @param noisyj
#'
#' @returns
#' @noRd
#'
#' @examples
A4type_ijloop_sum <- function(df, ipos, jpos, kpos, lpos, IdPQ, IdPW, noisyi = FALSE, noisyj = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)

  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1, n), ncol = 1)
  QQQinv <- Q %*% QQinv
  WWWinv <- W %*% WWinv

  A4vec <- rep(0, n)
  for (i in 1:n) {
    dPQi <- rep(0, n)
    dPQi[i] <- dPQ[i]
    dPWi <- rep(0, n)
    dPWi[i] <- dPW[i]
    Pi <- QQQinv %*% Q[i, ]
    Gi <- (QQQinv %*% Q[i, ] - dPQi) / dM -
      (WWWinv %*% W[i, ] - dPWi) / dMW
    Mi <- -Pi
    Mi[i] <- 1 - Pi[i]
    D2i <- dM * dM[i] - Mi^2
    recD2i <- 1 / D2i
    recD2i[i] <- 0

    A41i <- A42i <- A43i <- rep(0, n)

    for (j in (1:n)[-i]) {
      dPQj <- rep(0, n)
      dPQj[j] <- dPQ[j]
      dPWj <- rep(0, n)
      dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j, ]
      Mj <- -Pj
      Mj[j] <- 1 - Pj[j]
      D2j <- dM * dM[j] - Mj^2

      D3ij <- Mi[i] * D2j - (dM * Mi[j]^2 + Mi^2 * dM[j] - 2 * Mj * Mi * Mi[j])
      recD3ij <- 1 / D3ij
      recD3ij[i] <- 0
      recD3ij[j] <- 0
      D2D3ij <- D2j / D3ij
      D2D3ij[i] <- 0
      D2D3ij[j] <- 0

      A41ij <- (t(Gi^2 * df$jpos * dM * recD2i * df$lpos) %*% (D2D3ij)) * (Mi[j] * df$kpos[j]) -
        (t(Gi^2 * df$jpos * recD2i * df$lpos) %*% (D2D3ij * Mj * (onesN * df$kpos[j]))) * Mi[j]
      A42ij <- (t(Gi^2 * df$jpos * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j]^2 * df$kpos[j]) -
        (t(Gi^2 * df$jpos * Mi * recD2i) %*% (recD3ij * Mj^2)) * (Mi[j] * df$kpos[j]) -
        (t(Gi^2 * df$jpos * dM * Mi * recD2i) %*% ((onesN * dM[j]) * recD3ij)) * (Mi[j] * df$kpos[j]) +
        (t(Gi^2 * df$jpos * Mi^2 * recD2i) %*% (recD3ij * Mj)) %*% (dM[j] * df$kpos[j])
      A43ij <- (t(Gi^2 * df$jpos * dM * Mi * recD2i) %*% (recD3ij)) * (Mi[j]^2 * df$kpos[j] * df$lpos[j]) -
        (t(Gi^2 * df$jpos * Mi^2 * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * df$kpos[j] * df$lpos[j]) -
        Mi[i] * (t(Gi^2 * df$jpos * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * df$kpos[j] * df$lpos[j]) +
        Mi[i] * (t(Gi^2 * df$jpos * Mi * recD2i) %*% (recD3ij * Mj^2)) * (df$lpos[j] * df$kpos[j])

      A41i[j] <- A41ij
      A42i[j] <- A42ij
      A43i[j] <- A43ij

      if (noisyj == TRUE & j %% 10 == 0) cat("j:", j / n, " ")
    }

    A44i <- df$kpos[i] * Mi[i] * (t(Gi^2 * recD2i) %*% (df$lpos * df$jpos)) -
      df$kpos[i] * (df$lpos[i]) * (t(Gi^2 * recD2i) %*% (Mi * df$jpos))
    A4vec[i] <- sum(A41i + A42i * df$lpos + A43i) + A44i

    if (noisyi == TRUE) cat("i:", i / n, " ")
  }

  sum(A4vec * df$ipos)
}

#' Title
#'
#' @param df
#' @param ipos
#' @param jpos
#' @param kpos
#' @param lpos
#' @param IdPQ
#' @param IdPW
#' @param noisyi
#' @param noisyj
#'
#' @returns
#' @noRd
#'
#' @examples
A5type_ijloop_sum <- function(df, ipos, jpos, kpos, lpos, IdPQ, IdPW, noisyi = FALSE, noisyj = FALSE) {
  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)
  n <- nrow(df)

  # Stuff that can be calculated once
  dM <- IdPQ
  dMW <- IdPW
  onesN <- matrix(rep(1, n), ncol = 1)
  QQQinv <- Q %*% QQinv
  WWWinv <- W %*% WWinv

  A5vec <- rep(0, n)
  for (i in 1:n) {
    dPQi <- rep(0, n)
    dPQi[i] <- dPQ[i]
    dPWi <- rep(0, n)
    dPWi[i] <- dPW[i]
    Pi <- QQQinv %*% Q[i, ]
    Girow <- (QQQinv %*% Q[i, ] - dPQi) / dM[i] -
      (WWWinv %*% W[i, ] - dPWi) / dMW[i]
    Gicol <- t(QQQinv %*% Q[i, ] - dPQi) / dM -
      t(WWWinv %*% W[i, ] - dPWi) / dMW
    Mi <- -Pi
    Mi[i] <- 1 - Pi[i]
    D2i <- dM * dM[i] - Mi^2
    recD2i <- 1 / D2i
    recD2i[i] <- 0

    A51i <- A52i <- A53i <- rep(0, n)

    for (j in (1:n)[-i]) {
      dPQj <- rep(0, n)
      dPQj[j] <- dPQ[j]
      dPWj <- rep(0, n)
      dPWj[j] <- dPW[j]
      Pj <- QQQinv %*% Q[j, ]
      Mj <- -Pj
      Mj[j] <- 1 - Pj[j]
      D2j <- dM * dM[j] - Mj^2

      D3ij <- Mi[i] * D2j - (dM * Mi[j]^2 + Mi^2 * dM[j] - 2 * Mj * Mi * Mi[j])
      recD3ij <- 1 / D3ij
      recD3ij[i] <- 0
      recD3ij[j] <- 0
      D2D3ij <- D2j / D3ij
      D2D3ij[i] <- 0
      D2D3ij[j] <- 0

      A51ij <- (t(Girow * Gicol * df$jpos * dM * recD2i * df$lpos) %*% (D2D3ij)) * (Mi[j] * df$kpos[j]) -
        (t(Girow * Gicol * df$jpos * recD2i * df$lpos) %*% (D2D3ij * Mj * (onesN * df$kpos[j]))) * Mi[j]
      A52ij <- (t(Girow * Gicol * df$jpos * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j]^2 * df$kpos[j]) -
        (t(Girow * Gicol * df$jpos * Mi * recD2i) %*% (recD3ij * Mj^2)) * (Mi[j] * df$kpos[j]) -
        (t(Girow * Gicol * df$jpos * dM * Mi * recD2i) %*% ((onesN * dM[j]) * recD3ij)) * (Mi[j] * df$kpos[j]) +
        (t(Girow * Gicol * df$jpos * Mi^2 * recD2i) %*% (recD3ij * Mj)) %*% (dM[j] * df$kpos[j])
      A53ij <- (t(Girow * Gicol * df$jpos * dM * Mi * recD2i) %*% (recD3ij)) * (Mi[j]^2 * df$kpos[j] * df$lpos[j]) -
        (t(Girow * Gicol * df$jpos * Mi^2 * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * df$kpos[j] * df$lpos[j]) -
        Mi[i] * (t(Girow * Gicol * df$jpos * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * df$kpos[j] * df$lpos[j]) +
        Mi[i] * (t(Girow * Gicol * df$jpos * Mi * recD2i) %*% (recD3ij * Mj^2)) * (df$lpos[j] * df$kpos[j])

      A51i[j] <- A51ij
      A52i[j] <- A52ij
      A53i[j] <- A53ij

      if (noisyj == TRUE & j %% 10 == 0) cat("j:", j / n, " ")
    }

    A54i <- df$kpos[i] * Mi[i] * (t(Girow * Gicol * recD2i) %*% (df$lpos * df$jpos)) -
      df$kpos[i] * (df$lpos[i]) * (t(Girow * Gicol * recD2i) %*% (Mi * df$jpos))
    A5vec[i] <- sum(A51i + A52i * df$lpos + A53i) + A54i

    if (noisyi == TRUE) cat("i:", i / n, " ")
  }

  sum(A5vec * df$ipos)
}
