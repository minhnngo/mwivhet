#' Compute Generalized A1-Type Variance Component
#'
#' @description
#' Calculates a variance component structured like the \eqn{A_1} term in the L3O variance
#' estimator. This function computes the "outward" interaction \eqn{G_{ij} G_{ik}}, representing
#' the variance contribution from observation \eqn{i} acting on pairs \eqn{j, k}.
#'
#'
#' @param df Data frame. Contains the variables specified by the position arguments.
#' @param P Matrix of dimension n x n. The projection matrix of instruments.
#' @param G Matrix of dimension n x n. The UJIVE weighting matrix.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}}.
#' @param kpos Column name (unquoted). Inner summation term \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Bias correction/residual term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = \sum_{i} v_i^{(I)} \left[ v_i^{(L)} \sum_{j \neq i} \sum_{k \neq i} G_{ij} v_j^{(J)} G_{ik} v_k^{(K)} W_{ijk} - \text{BiasCorrect} \right]}
#' where \eqn{W_{ijk}} represents the Leave-Three-Out weighting derived from the annihilator matrix \eqn{M = I-P}.
#'
#' This structure is symmetric with respect to indices \eqn{j} and \eqn{k} relative to \eqn{i}.
#' It corresponds to the variance of the fitted values (signal) in the score statistic.
#'
#' @return Numeric scalar.
#'
#' @noRd
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

#' Compute Generalized A2-Type Variance Component
#'
#' @description
#' Calculates a variance component structurally equivalent to the \eqn{A_2} term.
#' This function computes the "chain" interaction \eqn{G_{ij} G_{ki}}, required for
#' covariance terms when the weighting matrix \eqn{G} is asymmetric.
#'
#' @param df Data frame. Contains the variables specified by the position arguments.
#' @param P Matrix of dimension n x n. The projection matrix of instruments.
#' @param G Matrix of dimension n x n. The UJIVE weighting matrix.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}}.
#' @param kpos Column name (unquoted). Inner summation term \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Bias correction/residual term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = \sum_{i} v_i^{(I)} \left[ v_i^{(L)} \sum_{j \neq i} \sum_{k \neq i} G_{ij} v_j^{(J)} G_{ki} v_k^{(K)} W_{ijk} - \text{BiasCorrect} \right]}
#'
#' The primary distinction from \code{A1type_iloop_sum} is the use of \eqn{G_{ki}} (column \eqn{k}, row \eqn{i})
#' in the second position. This asymmetric structure captures the interaction between
#' fitted values and residuals.
#'
#' @return Numeric scalar.
#'
#' @noRd
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

#' Compute Generalized A3-Type Variance Component
#'
#' @description
#' Calculates a variance component structurally equivalent to the \eqn{A_3} term.
#' This function computes the "inward" interaction \eqn{G_{ji} G_{ki}}, representing
#' the variance contribution from others acting on observation \eqn{i}.
#'
#' @param df Data frame. Contains the variables specified by the position arguments.
#' @param P Matrix of dimension n x n. The projection matrix of instruments.
#' @param G Matrix of dimension n x n. The UJIVE weighting matrix.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}}.
#' @param kpos Column name (unquoted). Inner summation term \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Bias correction/residual term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = \sum_{i} v_i^{(I)} \left[ v_i^{(L)} \sum_{j \neq i} \sum_{k \neq i} G_{ji} v_j^{(J)} G_{ki} v_k^{(K)} W_{ijk} - \text{BiasCorrect} \right]}
#'
#' This corresponds to the variance of the residuals projected onto the instruments,
#' typically denoted as \eqn{e' G' G e} in matrix notation (specifically the diagonal-removed part).
#'
#' @return Numeric scalar.
#'
#' @noRd
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

#' Compute Generalized A4-Type Variance Component
#'
#' @description
#' Calculates the "own-variance" bias correction component, structurally equivalent to the
#' \eqn{A_4} term. This term captures the bias arising from the squared diagonal weights
#' \eqn{G_{ji}^2}.
#'
#' @param df Data frame. Contains the variables specified by the position arguments.
#' @param P Matrix of dimension n x n. The projection matrix of instruments.
#' @param G Matrix of dimension n x n. The UJIVE weighting matrix.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}}.
#' @param kpos Column name (unquoted). Term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Residual/bias interaction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function computes:
#' \deqn{S = \sum_{i} v_i^{(I)} \sum_{j \neq i} G_{ji}^2 \widehat{Var}(v_j) W_{ij}}
#'
#' It removes the positive bias introduced by the variance of the instrument projection errors.
#'
#' @return Numeric scalar.
#'
#' @noRd
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

#' Compute Generalized A5-Type Variance Component
#'
#' @description
#' Calculates the "asymmetry" bias correction component, structurally equivalent to the
#' \eqn{A_5} term. This term captures the bias arising from the interaction of row and column
#' weights \eqn{G_{ij} G_{ji}}.
#'
#' @param df Data frame. Contains the variables specified by the position arguments.
#' @param P Matrix of dimension n x n. The projection matrix of instruments.
#' @param G Matrix of dimension n x n. The UJIVE weighting matrix.
#' @param ipos Column name (unquoted). Outer summation weight \eqn{v^{(I)}}.
#' @param jpos Column name (unquoted). Inner summation term \eqn{v^{(J)}}.
#' @param kpos Column name (unquoted). Term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Column name (unquoted). Residual/bias interaction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'   Defaults to \code{FALSE}.
#'
#' @details
#' This function computes:
#' \deqn{S = \sum_{i} v_i^{(I)} \sum_{j \neq i} G_{ij} G_{ji} \widehat{Var}(v_j) W_{ij}}
#'
#' If \eqn{G} is symmetric (e.g., \eqn{G=P}), \eqn{A_5} is identical to \eqn{A_4}.
#' For asymmetric \eqn{G} (e.g., UJIVE), this term accounts for the distinct bias channel.
#'
#' @return Numeric scalar.
#'
#' @noRd
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
