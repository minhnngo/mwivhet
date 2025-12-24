#' Leave-Three-Out Variance Estimator for LM Statistic
#'
#' @description
#' Calculates the consistent variance estimator ($\hat{V}_{LM}$) for the Lagrange Multiplier (LM)
#' test statistic. This estimator is robust to both many weak instruments and heterogeneous treatment effects.
#'
#' The function implements the "Leave-Three-Out" (L3O) approach proposed in Yap (2025),
#' which corrects for biases in variance estimation that arise when reduced-form coefficients
#' are not consistently estimable (e.g., due to many instruments).
#'
#' @param X A numeric vector of length n containing the endogenous variable.
#' @param e A numeric vector of length n containing the residuals under the null hypothesis
#'   ($e = Y - X\beta_0$).
#' @param P A numeric n x n projection matrix of the instruments (and potentially covariates).
#'   Corresponds to matrix $P$ or $H_Q$ in the paper.
#' @param G A numeric n x n weighting matrix used in the JIVE/UJIVE estimator.
#'   For standard JIVE, G is equal to P. For UJIVE with covariates, G is the adjusted
#'   matrix defined in Section 3.1.
#' @param noisy A logical indicating whether to print progress dots during the loop.
#'   Defaults to FALSE.
#'
#' @details
#' The function computes the variance estimator defined in Equation (9) of the paper:
#' \deqn{\hat{V}_{LM} = A_1 + A_2 + A_3 + A_4 + A_5}
#'
#' It iterates through each observation $i$ to compute the necessary adjustments
#' (Leave-Three-Out determinants $D_{ijk}$) and aggregates the components using
#' optimized matrix operations to handle the double sums over $j$ and $k$.
#'
#' Specifically:
#' \itemize{
#'   \item \strong{A1-A3} capture the core variance components involving interactions between
#'   the instruments, endogenous variable, and residuals.
#'   \item \strong{A4-A5} are correction terms that account for the variability from estimating
#'   the reduced-form coefficients (which cannot be treated as fixed in the many-instrument setting).
#' }
#'
#' The calculation relies on determinants $D_{ij}$ and $D_{ijk}$ derived from the annihilator
#' matrix $M = I - P$ to ensure the estimator is unbiased.
#'
#' @returns A scalar numeric value representing the estimated variance $\hat{V}_{LM}$.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity".
#'
#' @export
#'
#' @examples
L3Ovar_iloop_cov <- function(X, e, P, G, noisy = FALSE) {
  n <- length(X)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Me <- M %*% e
  MX <- M %*% X
  Poff <- P - diag(diag(P))
  Goff <- G - diag(diag(G))

  A1vec <- A2vec <- A3vec <- A4vec <- A5vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
    Gicol <- matrix(G[, i], ncol = 1) # force column vector
    Girow <- matrix(G[i, ], ncol = 1) # force column vector
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

    A11i <- t(X * Girow) %*% D2D3i %*% (Girow * X) * (t(Mi) %*% e)
    A12i <- t((Me) * X * Girow * Mi) %*% recD3i %*% (Girow * X * dM) -
      t(X * Girow * Mi) %*% (recD3i * M) %*% (Girow * X * (Me))
    A13i <- t(dM * X * Girow) %*% ((onesN %x% t(Me)) * recD3i) %*% (Girow * Mi * X)
    A14i <- t((Me) * X * Girow) %*% (M * recD3i) %*% (Gicol * Mi * X)
    A15i <- (t(Mi) %*% e) * (t(Goffirow^2 * recD2i) %*% (dM * X^2)) -
      (t(Goffirow^2 * Mi * recD2i) %*% (Me * X^2))

    A21i <- t(X * Girow) %*% D2D3i %*% (Gicol * e) * (t(Mi) %*% X)
    A22i <- t((MX) * X * Girow * Mi) %*% recD3i %*% (Gicol * e * dM) -
      t(X * Girow * Mi) %*% (recD3i * M) %*% (Gicol * e * (MX))
    A23i <- t(dM * X * Girow) %*% ((onesN %x% t(MX)) * recD3i) %*% (Gicol * Mi * e)
    A24i <- t((MX) * X * Girow) %*% (M * recD3i) %*% (Gicol * Mi * e)
    A25i <- (t(Mi) %*% X) * (t(Goffirow * Gofficol * recD2i) %*% (dM * X * e)) -
      (t(Goffirow * Gofficol * Mi * recD2i) %*% (MX * X * e))

    A31i <- t(e * Gicol) %*% D2D3i %*% (Gicol * e) * (t(Mi) %*% X)
    A32i <- t((MX) * e * Gicol * Mi) %*% recD3i %*% (Gicol * e * dM) -
      t(e * Gicol * Mi) %*% (recD3i * M) %*% (Gicol * e * (MX))
    A33i <- t(dM * e * Gicol) %*% ((onesN %x% t(MX)) * recD3i) %*% (Gicol * Mi * e)
    A34i <- t((MX) * e * Gicol) %*% (M * recD3i) %*% (Gicol * Mi * e)
    A35i <- (t(Mi) %*% X) * (t(Gofficol^2 * recD2i) %*% (dM * e * e)) -
      (t(Gofficol^2 * Mi * recD2i) %*% (MX * e * e))

    A41i <- t(Gicol^2 * e * dM * recD2i * Me) %*% (D2D3i) %*% (Mi * X) -
      t(Gicol^2 * e * recD2i * Me) %*% (D2D3i * M * (onesN %x% t(X))) %*% (Mi)
    A42i <- t(Gicol^2 * e * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * X) -
      t(Gicol^2 * e * Mi * recD2i) %*% (recD3i * M * M) %*% (Mi * X) -
      t(Gicol^2 * e * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (Mi * X) +
      t(Gicol^2 * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * X)
    A43i <- t(Gicol^2 * e * dM * Mi * recD2i) %*% (recD3i) %*% (Mi^2 * X * Me) -
      t(Gicol^2 * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (Mi * X * Me) -
      M[i, i] * t(Gicol^2 * e * dM * recD2i) %*% (recD3i * M) %*% (Mi * X * Me) +
      M[i, i] * t(Gicol^2 * e * Mi * recD2i) %*% (recD3i * M * M) %*% (Me * X)
    A44i <- X[i] * M[i, i] * (t(Gofficol^2 * recD2i) %*% (Me * e)) -
      X[i] * (t(Mi) %*% e) * (t(Gofficol^2 * recD2i) %*% (Mi * e))

    A51i <- t(Girow * Gicol * e * dM * recD2i * MX) %*% (D2D3i) %*% (Mi * X) -
      t(Girow * Gicol * e * recD2i * MX) %*% (D2D3i * M * (onesN %x% t(X))) %*% (Mi)
    A52i <- t(Girow * Gicol * e * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * X) -
      t(Girow * Gicol * e * Mi * recD2i) %*% (recD3i * M * M) %*% (Mi * X) -
      t(Girow * Gicol * e * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (Mi * X) +
      t(Girow * Gicol * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * X)
    A53i <- t(Girow * Gicol * e * dM * Mi * recD2i) %*% (recD3i) %*% (Mi^2 * X * MX) -
      t(Girow * Gicol * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (Mi * X * MX) -
      M[i, i] * t(Girow * Gicol * e * dM * recD2i) %*% (recD3i * M) %*% (Mi * X * MX) +
      M[i, i] * t(Girow * Gicol * e * Mi * recD2i) %*% (recD3i * M * M) %*% (MX * X)
    A54i <- X[i] * M[i, i] * (t(Goffirow * Gofficol * recD2i) %*% (MX * e)) -
      X[i] * (t(Mi) %*% X) * (t(Goffirow * Gofficol * recD2i) %*% (Mi * e))

    A1vec[i] <- A11i - A12i - A13i + A14i + A15i
    A2vec[i] <- A21i - A22i - A23i + A24i + A25i
    A3vec[i] <- A31i - A32i - A33i + A34i + A35i
    A4vec[i] <- A41i + A42i * (M[i, ] %*% e) + A43i + A44i
    A5vec[i] <- A51i + A52i * (M[i, ] %*% X) + A53i + A54i

    if (noisy) {
      if (i %% 10 == 0) cat(i / n, " ")
    }
  }

  sum(A1vec * e + 2 * A2vec * e + A3vec * X - A4vec * X - A5vec * e)
}



#' Leave-Three-Out Variance Estimator (Nested Loop Implementation)
#'
#' @description
#' Calculates the consistent variance estimator ($\hat{V}_{LM}$) for the Lagrange Multiplier (LM)
#' test statistic using an explicit nested loop structure.
#'
#' This function performs the same statistical estimation as \code{L3Ovar_iloop_cov} but constructs
#' the weighting matrix G dynamically from instruments (Q) and covariates (W). It iterates explicitly
#' over indices i and j, vectorizing the third summation index k. This implementation is generally
#' slower but useful for verifying results or when the full G matrix is too large to store in memory.
#'
#' @param X A numeric vector of length n containing the endogenous variable.
#' @param e A numeric vector of length n containing the residuals ($e = Y - X\beta_0$).
#' @param MX A numeric vector representing the product of the annihilator matrix M and X
#'   ($M \%*\% X$). Passed as an argument to avoid re-computation.
#' @param Me A numeric vector representing the product of the annihilator matrix M and e
#'   ($M \%*\% e$). Passed as an argument to avoid re-computation.
#' @param W A numeric matrix of covariates (controls).
#' @param Q A numeric matrix of instruments (including covariates).
#' @param WWWinv A numeric matrix representing the inverse of the gram matrix of W ($(W'W)^{-1}$).
#' @param QQQinv A numeric matrix representing the inverse of the gram matrix of Q ($(Q'Q)^{-1}$).
#' @param IdPQ A numeric vector containing the diagonal elements of the projection matrix $P_Q$.
#' @param IdPW A numeric vector containing the diagonal elements of the projection matrix $P_W$.
#' @param noisyi A logical indicating whether to print progress for the outer loop (i). Defaults to TRUE.
#' @param noisyj A logical indicating whether to print progress for the inner loop (j). Defaults to FALSE.
#'
#' @details
#' The function calculates the variance estimator defined in Equation (9) of Yap (2025):
#' \deqn{\hat{V}_{LM} = A_1 + A_2 + A_3 + A_4 + A_5}
#'
#' Unlike the vectorized version, this function computes the rows and columns of the weighting
#' matrix G on-the-fly inside the loop using the UJIVE formula:
#' \deqn{G_{ij} = \frac{P_{Q,ij}}{M_{Q,ii}} - \frac{P_{W,ij}}{M_{W,ii}}}
#'
#' The triple summation required for the variance components is handled by:
#' \enumerate{
#'   \item Outer loop over $i$.
#'   \item Inner loop over $j$.
#'   \item Vectorized operations over the vectors to handle index $k$.
#' }
#'
#' @returns A scalar numeric value representing the estimated variance $\hat{V}_{LM}$.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity".
#'
#' @export
#'
#' @examples
L3Ovar_ijloop_cov <- function(X, e, MX, Me, W, Q, WWWinv, QQQinv, IdPQ, IdPW, noisyi = TRUE, noisyj = FALSE) {
  n <- length(X)
  onesN <- matrix(rep(1, n), ncol = 1)
  dPQ <- 1 - IdPQ
  dPW <- 1 - IdPW
  dM <- IdPQ
  dMW <- IdPW

  A1vec <- rep(0, n)
  A2vec <- rep(0, n)
  A3vec <- rep(0, n)
  A4vec <- rep(0, n)
  A5vec <- rep(0, n)
  for (i in 1:n) {
    dPQi <- rep(0, n)
    dPQi[i] <- dPQ[i]
    dPWi <- rep(0, n)
    dPWi[i] <- dPW[i]

    Pi <- QQQinv %*% Q[i, ] # Should this be Pi <- Q %*% (QQQinv %*% Q[i, ])
    Girow <- (QQQinv %*% Q[i, ] - dPQi) / dM[i] -
      (WWWinv %*% W[i, ] - dPWi) / dMW[i]
    Gicol <- t(QQQinv %*% Q[i, ] - dPQi) / dM -
      t(WWWinv %*% W[i, ] - dPWi) / dMW
    Girow <- c(Girow) # force vectors
    Gicol <- c(Gicol)
    # Gicol <- t(Gicol) # Make it row vector, conformable with Gi
    Mi <- -Pi
    Mi[i] <- 1 - Pi[i]
    D2i <- dM * dM[i] - Mi^2
    recD2i <- 1 / D2i
    recD2i[i] <- 0

    A11i <- A12i <- A13i <- A14i <- rep(0, n)
    A21i <- A22i <- A23i <- A24i <- rep(0, n)
    A31i <- A32i <- A33i <- A34i <- rep(0, n)
    A41i <- A42i <- A43i <- rep(0, n)
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

      A11ij <- (t(X * Girow) %*% D2D3ij) * (Girow[j] * X[j]) * (Me[i])
      A12ij <- (t(Me * X * Girow * Mi) %*% recD3ij) * (Girow[j] * X[j] * dM[j]) -
        (t(X * Girow * Mi) %*% (recD3ij * Mj)) %*% (Girow[j] * X[j] * Me[j])
      A13ij <- (t(dM * X * Girow) %*% (onesN * Me[j] * recD3ij)) * (Girow[j] * Mi[j] * X[j])
      A14ij <- (t(Me * X * Girow) %*% (Mj * recD3ij)) * (Girow[j] * Mi[j] * X[j])

      A21ij <- (t(X * Girow) %*% D2D3ij) * (Gicol[j] * e[j]) * (MX[i])
      A22ij <- (t(MX * X * Girow * Mi) %*% recD3ij) * (Gicol[j] * e[j] * dM[j]) -
        (t(X * Girow * Mi) %*% (recD3ij * Mj)) %*% (Gicol[j] * e[j] * MX[j])
      A23ij <- (t(dM * X * Girow) %*% (onesN * MX[j] * recD3ij)) * (Gicol[j] * Mi[j] * e[j])
      A24ij <- (t(MX * X * Girow) %*% (Mj * recD3ij)) * (Gicol[j] * Mi[j] * e[j])

      A31ij <- (t(e * Gicol) %*% D2D3ij) * (Gicol[j] * e[j]) * (MX[i])
      A32ij <- (t(MX * e * Gicol * Mi) %*% recD3ij) * (Gicol[j] * e[j] * dM[j]) -
        (t(e * Gicol * Mi) %*% (recD3ij * Mj)) %*% (Gicol[j] * e[j] * MX[j])
      A33ij <- (t(dM * e * Gicol) %*% (onesN * MX[j] * recD3ij)) * (Gicol[j] * Mi[j] * e[j])
      A34ij <- (t(MX * e * Gicol) %*% (Mj * recD3ij)) * (Gicol[j] * Mi[j] * e[j])

      A41ij <- (t(Gicol^2 * e * dM * recD2i * Me) %*% (D2D3ij)) * (Mi[j] * X[j]) -
        (t(Gicol^2 * e * recD2i * Me) %*% (D2D3ij * Mj * (onesN * X[j]))) * Mi[j]
      A42ij <- (t(Gicol^2 * e * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j]^2 * X[j]) -
        (t(Gicol^2 * e * Mi * recD2i) %*% (recD3ij * Mj^2)) * (Mi[j] * X[j]) -
        (t(Gicol^2 * e * dM * Mi * recD2i) %*% ((onesN * dM[j]) * recD3ij)) * (Mi[j] * X[j]) +
        (t(Gicol^2 * e * Mi^2 * recD2i) %*% (recD3ij * Mj)) %*% (dM[j] * X[j])
      A43ij <- (t(Gicol^2 * e * dM * Mi * recD2i) %*% (recD3ij)) * (Mi[j]^2 * X[j] * Me[j]) -
        (t(Gicol^2 * e * Mi^2 * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * X[j] * Me[j]) -
        Mi[i] * (t(Gicol^2 * e * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * X[j] * Me[j]) +
        Mi[i] * (t(Gicol^2 * e * Mi * recD2i) %*% (recD3ij * Mj^2)) * (Me[j] * X[j])


      A51ij <- (t(Girow * Gicol * e * dM * recD2i * MX) %*% (D2D3ij)) * (Mi[j] * X[j]) -
        (t(Girow * Gicol * e * recD2i * MX) %*% (D2D3ij * Mj * (onesN * X[j]))) * Mi[j]
      A52ij <- (t(Girow * Gicol * e * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j]^2 * X[j]) -
        (t(Girow * Gicol * e * Mi * recD2i) %*% (recD3ij * Mj^2)) * (Mi[j] * X[j]) -
        (t(Girow * Gicol * e * dM * Mi * recD2i) %*% ((onesN * dM[j]) * recD3ij)) * (Mi[j] * X[j]) +
        (t(Girow * Gicol * e * Mi^2 * recD2i) %*% (recD3ij * Mj)) %*% (dM[j] * X[j])
      A53ij <- (t(Girow * Gicol * e * dM * Mi * recD2i) %*% (recD3ij)) * (Mi[j]^2 * X[j] * MX[j]) -
        (t(Girow * Gicol * e * Mi^2 * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * X[j] * MX[j]) -
        Mi[i] * (t(Girow * Gicol * e * dM * recD2i) %*% (recD3ij * Mj)) * (Mi[j] * X[j] * MX[j]) +
        Mi[i] * (t(Girow * Gicol * e * Mi * recD2i) %*% (recD3ij * Mj^2)) * (MX[j] * X[j])


      A11i[j] <- A11ij
      A12i[j] <- A12ij
      A13i[j] <- A13ij
      A14i[j] <- A14ij
      A21i[j] <- A21ij
      A22i[j] <- A22ij
      A23i[j] <- A23ij
      A24i[j] <- A24ij
      A31i[j] <- A31ij
      A32i[j] <- A32ij
      A33i[j] <- A33ij
      A34i[j] <- A34ij
      A41i[j] <- A41ij
      A42i[j] <- A42ij
      A43i[j] <- A43ij
      A51i[j] <- A51ij
      A52i[j] <- A52ij
      A53i[j] <- A53ij

      if (noisyj == TRUE & j %% 10 == 0) cat("j:", j / n, " ")
    }

    A1vec[i] <- sum(A11i - A12i - A13i + A14i)
    A2vec[i] <- sum(A21i - A22i - A23i + A24i)
    A3vec[i] <- sum(A31i - A32i - A33i + A34i)

    A15i <- (Me[i]) * (t(Girow^2 * recD2i) %*% (dM * X * X)) -
      (t(Girow^2 * Mi * recD2i) %*% (Me * X * X))
    A25i <- (MX[i]) * (t(Girow * Gicol * recD2i) %*% (dM * X * e)) -
      (t(Girow * Gicol * Mi * recD2i) %*% (MX * X * e))
    A35i <- (MX[i]) * (t(Gicol^2 * recD2i) %*% (dM * e * e)) -
      (t(Gicol^2 * Mi * recD2i) %*% (MX * e * e))

    A44i <- X[i] * Mi[i] * (t(Gicol^2 * recD2i) %*% (Me * e)) -
      X[i] * (Me[i]) * (t(Gicol^2 * recD2i) %*% (Mi * e))
    A54i <- X[i] * Mi[i] * (t(Girow * Gicol * recD2i) %*% (MX * e)) -
      X[i] * (MX[i]) * (t(Girow * Gicol * recD2i) %*% (Mi * e))


    A1vec[i] <- A1vec[i] + A15i
    A2vec[i] <- A2vec[i] + A25i
    A3vec[i] <- A3vec[i] + A35i
    A4vec[i] <- sum(A41i + A43i) + sum(A42i) * Me[i] + A44i
    A5vec[i] <- sum(A51i + A53i) + sum(A52i) * MX[i] + A54i


    if (noisyi == TRUE) cat("i:", i / n, " ")
  }

  sum(A1vec * e + 2 * A2vec * e + A3vec * X - A4vec * X - A5vec * e)
}


#' Title
#'
#' @param df
#' @param group
#' @param groupW
#' @param X
#' @param e
#' @param MX
#' @param Me
#'
#' @returns
#' @export
#'
#' @examples
L3Ovar_gloop_cov <- function(df, group, groupW, X, e, MX, Me) {
  df$group <- eval(substitute(group), df)
  df$groupW <- eval(substitute(groupW), df)
  df$X <- eval(substitute(X), df)
  df$e <- eval(substitute(e), df)
  df$MX <- eval(substitute(MX), df)
  df$Me <- eval(substitute(Me), df)

  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, max(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0, max(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0, max(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, max(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0, max(unique(df$group)))
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]
    ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds, group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
    ZW <- matrix(1, nrow = length(ds$groupW), ncol = length(unique(ds$groupW)))

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)
    PW <- ZW %*% solve(t(ZW) %*% ZW) %*% t(ZW)

    # calculate values specific to this subset
    Gs <- solve(diag(nrow(ds)) - diag(diag(PQ))) %*% (PQ - diag(diag(PQ))) -
      solve(diag(nrow(ds)) - diag(diag(PW))) %*% (PW - diag(diag(PW)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms), ncol = 1)
    D2s <- dMs %*% t(dMs) - Ms * Ms
    recD2s <- 1 / D2s
    diag(recD2s) <- 0

    for (g in unique(ds$group)) {
      repidx <- min(which(ds$group == g)) # representative index
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
      ones <- matrix(rep(1, nrow(ds)), ncol = 1)
      Mes <- matrix(ds$Me, ncol = 1)
      MXs <- matrix(ds$MX, ncol = 1)
      es <- matrix(ds$e, ncol = 1)

      A11vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$Me, g)) %*% (ds$X * Gis)
      A12vecs[g] <- t(ds$X * ds$Me * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Gis * ds$X * dMs) -
        t(ds$X * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Gis * ds$X * ds$Me)
      A13vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis * Mis)
      A14vecs[g] <- t(ds$X * ds$Me * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis * Mis)
      A15vecs[g] <- t(ds$Me * ds$e) %*% (recD2s * Gs^2 * Pgs) %*% (ds$X^2 * dMs) -
        t(ds$e) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$Me * ds$X^2)

      A21vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$MX, g)) %*% (ds$e * Gis)
      A22vecs[g] <- t(ds$X * ds$MX * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Gis * ds$e * dMs) -
        t(ds$X * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Gis * ds$e * ds$MX)
      A23vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis * Mis)
      A24vecs[g] <- t(ds$X * ds$MX * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis * Mis)
      A25vecs[g] <- t(ds$MX * ds$e) %*% (recD2s * Gs^2 * Pgs) %*% (ds$X * ds$e * dMs) -
        t(ds$e) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$MX * ds$X * ds$e)

      A31vecs[g] <- t(ds$e * Gis) %*% (D2D3is * ivectomats(ds, ds$X * ds$MX, g)) %*% (ds$e * Gis)
      A32vecs[g] <- t(ds$e * ds$MX * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Gis * ds$e * dMs) -
        t(ds$e * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (Gis * ds$e * ds$MX)
      A33vecs[g] <- t(ds$e * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis * Mis)
      A34vecs[g] <- t(ds$e * ds$MX * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis * Mis)
      A35vecs[g] <- t(ds$MX * ds$X) %*% (recD2s * Gs^2 * Pgs) %*% (ds$e^2 * dMs) -
        t(ds$X) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$MX * ds$e^2)

      A41vecs[g] <- t(Gis^2 * ds$e * dMs * ds$Me * recD2is) %*% (D2D3is * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * ds$Me * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$X, g)) %*% (Mis)
      A42vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis^2 * ds$X) -
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) +
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (dMs * ds$X)
      A43vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Mis^2 * ds$X * ds$Me) -
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X * ds$Me) -
        t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * dMs, g)) %*% (Mis * ds$X * ds$Me) +
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * dMs, g)) %*% (ds$X * ds$Me)
      A44vecs[g] <- t(ds$X^2 * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$Me * ds$e) -
        t(ds$X^2 * ds$Me) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$e)

      A51vecs[g] <- t(Gis^2 * ds$e * dMs * ds$MX * recD2is) %*% (D2D3is * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * ds$MX * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$e, g)) %*% (Mis)
      A52vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis^2 * ds$X) -
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) -
        t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) +
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (dMs * ds$X)
      A53vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Mis^2 * ds$X * ds$MX) -
        t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X * ds$MX) -
        t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * dMs, g)) %*% (Mis * ds$X * ds$MX) +
        t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * dMs, g)) %*% (ds$X * ds$MX)
      A54vecs[g] <- t(ds$X * ds$e * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$MX * ds$e) -
        t(ds$X * ds$e * ds$MX) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$e)

    }
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) +
    2 * (A21vecs - A22vecs - A23vecs + A24vecs + A25vecs) +
    (A31vecs - A32vecs - A33vecs + A34vecs + A35vecs) -
    (A41vecs + A42vecs + A43vecs + A44vecs) -
    (A51vecs + A52vecs + A53vecs + A54vecs)

  sum(ret)
}



#' Title
#'
#' @param df
#' @param group
#' @param X
#' @param e
#' @param MX
#' @param Me
#'
#' @returns
#' @export
#'
#' @examples
L3Ovar_gloop_nocov <- function(df, group, X, e, MX, Me) {
  df$group <- eval(substitute(group), df)
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, length(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0, length(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0, length(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, length(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0, length(unique(df$group)))
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, length(unique(df$group)))
  A21vecs <- A22vecs <- A23vecs <- A24vecs <- A25vecs <- rep(0, length(unique(df$group)))
  A31vecs <- A32vecs <- A33vecs <- A34vecs <- A35vecs <- rep(0, length(unique(df$group)))
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, length(unique(df$group)))
  A51vecs <- A52vecs <- A53vecs <- A54vecs <- rep(0, length(unique(df$group)))
  for (g in 1:length(unique(df$group))) {
    ds <- df[df$group == g, ]
    ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
    ds$groupidx <- Getgroupindex(ds, group)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- PQ
    Ps <- PQ
    Psoff <- Ps - diag(nrow(ds)) * diag(Ps)
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms), ncol = 1)
    D2s <- dMs %*% t(dMs) - Ms * Ms
    recD2s <- 1 / D2s
    diag(recD2s) <- 0
    # diag0 <- matrix(1,nrow=nrow(ds),ncol=nrow(ds)); diag(diag0) <- 0

    repidx <- min(which(ds$group == g)) # representative index
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
    ones <- matrix(rep(1, nrow(ds)), ncol = 1)
    Mes <- matrix(ds$Me, ncol = 1)
    MXs <- matrix(ds$MX, ncol = 1)
    es <- matrix(ds$e, ncol = 1)

    A11vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$Me, g)) %*% (ds$X * Gis)
    A12vecs[g] <- t(ds$X * ds$Me * Gis) %*% (Ms * Psoff * recD3is * ivectomats(ds, ds$e, g)) %*% (ds$X * dMs) -
      t(ds$X * Gis) %*% (recD3is * Ms * Ms * Psoff * ivectomats(ds, ds$e, g)) %*% (ds$X * ds$Me)
    A13vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis * Mis)
    A14vecs[g] <- t(ds$X * ds$Me * Gis) %*% (recD3is * Ms * Ms * ivectomats(ds, ds$e, g)) %*% (ds$X * Gis)
    A15vecs[g] <- t(ds$Me * ds$e) %*% (recD2s * Psoff^2) %*% (ds$X^2 * dMs) -
      t(ds$e) %*% (Ms * recD2s * Psoff^2) %*% (ds$Me * ds$X^2)

    A21vecs[g] <- t(ds$X * Gis) %*% (D2D3is * ivectomats(ds, ds$e * ds$MX, g)) %*% (ds$e * Gis)
    A22vecs[g] <- t(ds$X * ds$MX * Gis) %*% (Ms * Psoff * recD3is * ivectomats(ds, ds$e, g)) %*% (ds$e * dMs) -
      t(ds$X * Gis) %*% (recD3is * Ms * Ms * Psoff * ivectomats(ds, ds$e, g)) %*% (ds$e * ds$MX)
    A23vecs[g] <- t(ds$X * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis * Mis)
    A24vecs[g] <- t(ds$X * ds$MX * Gis) %*% (recD3is * Ms * Ms * ivectomats(ds, ds$e, g)) %*% (ds$e * Gis)
    A25vecs[g] <- t(ds$MX * ds$e) %*% (recD2s * Psoff^2) %*% (ds$X * ds$e * dMs) -
      t(ds$e) %*% (Ms * recD2s * Psoff^2) %*% (ds$MX * ds$X * ds$e)

    A31vecs[g] <- t(ds$e * Gis) %*% (D2D3is * ivectomats(ds, ds$X * ds$MX, g)) %*% (ds$e * Gis)
    A32vecs[g] <- t(ds$e * ds$MX * Gis) %*% (Ms * Psoff * recD3is * ivectomats(ds, ds$X, g)) %*% (ds$e * dMs) -
      t(ds$e * Gis) %*% (recD3is * Ms * Ms * Psoff * ivectomats(ds, ds$X, g)) %*% (ds$e * ds$MX)
    A33vecs[g] <- t(ds$e * dMs * Gis) %*% (recD3is * (ones %*% t(MXs)) * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis * Mis)
    A34vecs[g] <- t(ds$e * ds$MX * Gis) %*% (recD3is * Ms * Ms * ivectomats(ds, ds$X, g)) %*% (ds$e * Gis)
    A35vecs[g] <- t(ds$MX * ds$X) %*% (recD2s * Psoff^2) %*% (ds$e^2 * dMs) -
      t(ds$X) %*% (Ms * recD2s * Psoff^2) %*% (ds$MX * ds$e^2)

    A41vecs[g] <- t(Gis^2 * ds$e * dMs * ds$Me) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * ds$Me) %*% (recD3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$X, g)) %*% (Mis)
    A42vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis^2 * ds$X) -
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$X * ds$Me, g)) %*% (Mis * ds$X) +
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * ds$Me, g)) %*% (dMs * ds$X)
    A43vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$X, g)) %*% (Mis^2 * ds$X * ds$Me) -
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X, g)) %*% (Mis * ds$X * ds$Me) -
      t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$X * dMs, g)) %*% (Mis * ds$X * ds$Me) +
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$X * dMs, g)) %*% (ds$X * ds$Me)
    A44vecs[g] <- t(ds$X^2 * dMs) %*% (recD2s * Psoff^2) %*% (ds$Me * ds$e) -
      t(ds$X^2 * ds$Me) %*% (recD2s * Psoff^2) %*% (Mis * ds$e)

    A51vecs[g] <- t(Gis^2 * ds$e * dMs * ds$MX) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * ds$MX) %*% (recD3is * Ms * (ones %x% t(ds$X)) * ivectomats(ds, ds$e, g)) %*% (Mis)
    A52vecs[g] <- t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis^2 * ds$X) -
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) -
      t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$e * ds$MX, g)) %*% (Mis * ds$X) +
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * ds$MX, g)) %*% (dMs * ds$X)
    A53vecs[g] <- t(Gis^2 * ds$e * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$e, g)) %*% (Mis^2 * ds$X * ds$MX) -
      t(Gis^2 * ds$e * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e, g)) %*% (Mis * ds$X * ds$MX) -
      t(Gis^2 * ds$e * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$e * dMs, g)) %*% (Mis * ds$X * ds$MX) +
      t(Gis^2 * ds$e * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$e * dMs, g)) %*% (ds$X * ds$MX)
    A54vecs[g] <- t(ds$X * ds$e * dMs) %*% (recD2s * Psoff^2) %*% (ds$MX * ds$e) -
      t(ds$X * ds$e * ds$MX) %*% (recD2s * Psoff^2) %*% (Mis * ds$e)
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs) +
    2 * (A21vecs - A22vecs - A23vecs + A24vecs + A25vecs) +
    (A31vecs - A32vecs - A33vecs + A34vecs + A35vecs) -
    (A41vecs + A42vecs + A43vecs + A44vecs) -
    (A51vecs + A52vecs + A53vecs + A54vecs)

  sum(ret)
}


#' Title
#'
#' @param X
#' @param e
#' @param P
#' @param c
#'
#' @returns
#' @export
#'
#' @examples
L3Ovar_block <- function(X, e, P, c) {
  n <- length(X)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Me <- M %*% e
  MX <- M %*% X
  Poff <- P - diag(diag(P))

  D3block <- ((c - 1)^3 - (c - 1) * 3 - 2) / c^3 * (matrix(rep(1, c^2), nrow = c) - diag(c))
  recD3block <- c^3 / ((c - 1)^3 - (c - 1) * 3 - 2) * (matrix(rep(1, c^2), nrow = c) - diag(c))
  D2D3block <- c * ((c - 1)^2 - 1) / ((c - 1)^3 - (c - 1) * 3 - 2) * (matrix(rep(1, c^2), nrow = c) - diag(c))
  recD2block <- c^2 / ((c - 1)^2 - 1) * (matrix(rep(1, c^2), nrow = c) - diag(c))

  Pvec <- rep(1 / c, n)
  Mvec <- rep(-1 / c, n)

  D2D3blockrep <- as.matrix(Matrix::bdiag(replicate(J, D2D3block, simplify = FALSE)))
  recD3blockrep <- as.matrix(Matrix::bdiag(replicate(J, recD3block, simplify = FALSE)))
  recD2blockrep <- as.matrix(Matrix::bdiag(replicate(J, recD2block, simplify = FALSE)))
  blockrep <- Z %*% t(Z) - diag(n)
  ZtZ <- Z %*% t(Z)

  bindrowvec <- function(X) matrix(rep(X, n), nrow = n, byrow = TRUE)
  bindcolvec <- function(X) matrix(rep(X, n), nrow = n, byrow = FALSE)
  ivectomat <- function(X) (bindcolvec(ZtZ %*% X) * blockrep) - ((bindrowvec(X) + bindcolvec(X)) * blockrep)

  ## A1
  A11 <- t(X * Pvec) %*% (D2D3blockrep * ivectomat(Me * e)) %*% (Pvec * X)
  A12 <- t(Me * X * Pvec) %*% (recD3blockrep * Poff * M * ivectomat(e)) %*% (X * dM) -
    t(X * Pvec) %*% (recD3blockrep * M * M * Poff * ivectomat(e)) %*% (X * (Me))
  A13 <- t(dM * X * Pvec) %*% (recD3blockrep * (onesN %x% t(Me)) * M * ivectomat(e)) %*% (X * Pvec)
  A14 <- t(Me * X * Pvec) %*% (recD3blockrep * M * M * ivectomat(e)) %*% (X * Pvec)
  A15 <- t(Me * e) %*% ((Poff^2 * recD2) %*% (dM * X^2)) - t((Poff^2 * M * recD2) %*% (Me * X^2)) %*% e
  ## A2
  A21 <- t(X * Pvec) %*% (D2D3blockrep * ivectomat(MX * e)) %*% (Pvec * e)
  A22 <- t(MX * X * Pvec) %*% (M * recD3blockrep * ivectomat(e)) %*% (Pvec * e * dM) -
    t(X * Pvec) %*% (M * M * recD3blockrep * ivectomat(e)) %*% (Pvec * e * MX)
  A23 <- t(dM * X * Pvec) %*% (recD3blockrep * (onesN %x% t(MX)) * M * ivectomat(e)) %*% (e * Pvec)
  A24 <- t(MX * X * Pvec) %*% (recD3blockrep * M * M * ivectomat(e)) %*% (e * Pvec)
  A25 <- t(MX * e) %*% ((Poff^2 * recD2) %*% (dM * X * e)) - t((Poff^2 * M * recD2) %*% (MX * X * e)) %*% e
  ## A3
  A31 <- t(e * Pvec) %*% (D2D3blockrep * ivectomat(MX * X)) %*% (Pvec * e)
  A32 <- t(MX * e * Pvec) %*% (M * recD3blockrep * ivectomat(X)) %*% (Pvec * e * dM) -
    t(e * Pvec) %*% (M * M * recD3blockrep * ivectomat(X)) %*% (Pvec * e * MX)
  A33 <- t(dM * e * Pvec) %*% (recD3blockrep * (onesN %x% t(MX)) * M * ivectomat(X)) %*% (e * Pvec)
  A34 <- t(MX * e * Pvec) %*% (recD3blockrep * M * M * ivectomat(X)) %*% (e * Pvec)
  A35 <- t(MX * X) %*% ((Poff^2 * recD2) %*% (dM * e * e)) - t((Poff^2 * M * recD2) %*% (MX * e * e)) %*% X

  ## A4
  A41 <- t(Pvec^2 * e * dM * Me) %*% (recD3blockrep * ivectomat(X) * M) %*% (X) -
    t(Pvec^2 * e * Me) %*% (recD3blockrep * M^2 * (onesN %x% t(X)) * ivectomat(X)) %*% (onesN)
  A42 <- t(Pvec^2 * e * dM) %*% (recD2 * recD3blockrep * M^3 * ivectomat(X * Me + e * MX)) %*% (X) -
    t(Pvec^2 * e) %*% (recD2 * recD3blockrep * M^4 * ivectomat(X * Me + e * MX)) %*% (X) -
    t(Pvec^2 * e * dM) %*% (M * recD2 * recD3blockrep * (onesN %x% t(dM)) * M * ivectomat(X * Me + e * MX)) %*% (X) +
    t(Pvec^2 * e) %*% (recD2 * recD3blockrep * M^3 * ivectomat(X * Me + e * MX)) %*% (dM * X)
  A43 <- t(Pvec^2 * e * dM) %*% (recD2 * recD3blockrep * M^3 * ivectomat(X)) %*% (X * Me) -
    t(Pvec^2 * e) %*% (recD2 * recD3blockrep * M^4 * ivectomat(X)) %*% (X * Me) -
    t(Pvec^2 * e * dM) %*% (recD2 * recD3blockrep * M^2 * ivectomat(X * dM)) %*% (X * Me) +
    t(Pvec^2 * e) %*% (recD2 * recD3blockrep * M^3 * ivectomat(X * dM)) %*% (Me * X)
  A44 <- t(dM * X^2) %*% ((Poff^2 * recD2) %*% (Me * e)) - t((Poff^2 * recD2 * M) %*% (e)) %*% (X^2 * Me)

  ## A5
  A51 <- t(Pvec^2 * e * dM * MX) %*% (recD3blockrep * ivectomat(e) * M) %*% (X) -
    t(Pvec^2 * e * MX) %*% (recD3blockrep * M^2 * (onesN %x% t(X)) * ivectomat(e)) %*% (onesN)
  A53 <- t(Pvec^2 * e * dM) %*% (recD2 * recD3blockrep * M^3 * ivectomat(e)) %*% (X * MX) -
    t(Pvec^2 * e) %*% (recD2 * recD3blockrep * M^4 * ivectomat(e)) %*% (X * MX) -
    t(Pvec^2 * e * dM) %*% (recD2 * recD3blockrep * M^2 * ivectomat(e * dM)) %*% (X * MX) +
    t(Pvec^2 * e) %*% (recD2 * recD3blockrep * M^3 * ivectomat(e * dM)) %*% (MX * X)
  A54 <- t(dM * X * e) %*% ((Poff^2 * recD2) %*% (MX * e)) - t((Poff^2 * recD2 * M) %*% (e)) %*% (X * e * MX)


  ret <- (A11 - A12 - A13 + A14 + A15) + 2 * (A21 - A22 - A23 + A24 + A25) +
    (A31 - A32 - A33 + A34 + A35) - (A41 + A42 + A43 + A44) - (A51 + A53 + A54)

  ret
}


#' Title
#'
#' @param X
#' @param e
#' @param P
#'
#' @returns
#' @export
#'
#' @examples
L3Ovar_iloop_nocov <- function(X, e, P) {
  n <- length(X)
  M <- diag(n) - P

  # Stuff that can be calculated once
  dM <- matrix(diag(M), ncol = 1) # force column vector
  D2 <- dM %*% t(dM) - M * M
  onesN <- matrix(rep(1, n), ncol = 1)
  recD2 <- 1 / D2
  diag(recD2) <- 0
  Me <- M %*% e
  MX <- M %*% X
  Poff <- P - diag(diag(P))

  A1vec <- A2vec <- A3vec <- A4vec <- A5vec <- rep(0, n)
  for (i in 1:n) {
    # Calculation conditioned on i
    Mi <- matrix(M[, i], ncol = 1) # force column vector
    Pi <- matrix(P[, i], ncol = 1) # force column vector
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
    tMie <- t(Mi) %*% e
    tMiX <- t(Mi) %*% X
    Pie <- Pi * e
    PiX <- Pi * X
    MXe <- MX * e
    Mie <- Mi * e
    MiX <- Mi * X
    MrecD3iPiMie <- (M * recD3i) %*% (Pie * Mi)
    recD3iMMiXMX <- (recD3i * M) %*% (MiX * MX)
    recD3iMMiXMe <- (recD3i * M) %*% (MiX * Me)
    recD3iPiedM <- recD3i %*% (Pie * dM)
    Pi2eMirecD2irecD3iMM <- t(Pi^2 * e * Mi * recD2i) %*% (recD3i * M * M)
    Pi2edMMirecD2irecD3i <- t(Pi^2 * e * dM * Mi * recD2i) %*% (recD3i)
    PiXD2D3i <- t(PiX) %*% D2D3i


    A11i <- PiXD2D3i %*% (PiX) * (tMie)
    A12i <- t((Me) * PiX * Mi) %*% recD3i %*% (Poffi * X * dM) -
      t(PiX * Mi) %*% (recD3i * M) %*% (Poffi * X * (Me))
    A13i <- t(dM * PiX) %*% ((onesN %x% t(Me)) * recD3i) %*% (Pi * MiX)
    A14i <- t((Me) * PiX) %*% (M * recD3i) %*% (Pi * MiX)
    A15i <- (tMie) * (t(Poffi^2 * recD2i) %*% (dM * X^2)) -
      (t(Poffi^2 * Mi * recD2i) %*% (Me * X^2))

    A21i <- PiXD2D3i %*% (Pie) * (tMiX)
    A22i <- t((MX) * PiX * Mi) %*% recD3iPiedM -
      t(PiX * Mi) %*% (recD3i * M) %*% (Pie * (MX))
    A23i <- t(dM * PiX) %*% ((onesN %x% t(MX)) * recD3i) %*% (Pie * Mi)
    A24i <- t((MX) * PiX) %*% MrecD3iPiMie
    A25i <- (tMiX) * (t(Poffi^2 * recD2i) %*% (dM * X * e)) -
      (t(Poffi^2 * Mi * recD2i) %*% (MXe * X))

    A31i <- t(Pie) %*% D2D3i %*% (Pie) * (tMiX)
    A32i <- t((MX) * Pie * Mi) %*% recD3iPiedM -
      t(Pie * Mi) %*% (recD3i * M) %*% (Pie * (MX))
    A33i <- t(dM * Pie) %*% ((onesN %x% t(MX)) * recD3i) %*% (Pie * Mi)
    A34i <- t((MX) * Pie) %*% MrecD3iPiMie
    A35i <- (tMiX) * (t(Poffi^2 * recD2i) %*% (dM * e * e)) -
      (t(Poffi^2 * Mi * recD2i) %*% (MXe * e))

    A41i <- t(Pi^2 * e * dM * recD2i * Me) %*% ((onesN %x% t(Di)) * recD3i) %*% (MiX) -
      t(Pi^2 * e * recD2i * Me) %*% ((onesN %x% t(Di)) * recD3i * M * (onesN %x% t(X))) %*% (Mi)
    A42i <- t(Pi^2 * e * dM * recD2i) %*% (recD3i * M) %*% (Mi^2 * X) -
      Pi2eMirecD2irecD3iMM %*% (MiX) -
      t(Pi^2 * e * dM * Mi * recD2i) %*% ((onesN %x% t(dM)) * recD3i) %*% (MiX) +
      t(Pi^2 * e * Mi^2 * recD2i) %*% (recD3i * M) %*% (dM * X)
    A43i <- Pi2edMMirecD2irecD3i %*% (Mi^2 * X * Me) -
      t(Pi^2 * e * Mi^2 * recD2i) %*% recD3iMMiXMe -
      M[i, i] * t(Pi^2 * e * dM * recD2i) %*% recD3iMMiXMe +
      M[i, i] * Pi2eMirecD2irecD3iMM %*% (Me * X)
    A44i <- X[i] * M[i, i] * (t(Poffi^2 * recD2i) %*% (Me * e)) -
      X[i] * (tMie) * (t(Poffi^2 * recD2i) %*% (Mie))

    A51i <- t(Pi^2 * e * dM * recD2i * MX) %*% ((onesN %x% t(Di)) * recD3i) %*% (Mi * X) -
      t(Pi^2 * e * recD2i * MX) %*% ((onesN %x% t(Di)) * recD3i * M * (onesN %x% t(X))) %*% (Mi)
    A53i <- Pi2edMMirecD2irecD3i %*% (Mi^2 * X * MX) -
      t(Pi^2 * e * Mi^2 * recD2i) %*% recD3iMMiXMX -
      M[i, i] * t(Pi^2 * e * dM * recD2i) %*% recD3iMMiXMX +
      M[i, i] * Pi2eMirecD2irecD3iMM %*% (MX * X)
    A54i <- X[i] * M[i, i] * (t(Poffi^2 * recD2i) %*% (MXe)) -
      X[i] * (tMiX) * (t(Poffi^2 * recD2i) %*% (Mie))

    A1vec[i] <- A11i - A12i - A13i + A14i + A15i
    A2vec[i] <- A21i - A22i - A23i + A24i + A25i
    A3vec[i] <- A31i - A32i - A33i + A34i + A35i
    A4vec[i] <- A41i + A42i * (M[i, ] %*% e) + A43i + A44i
    A5vec[i] <- A51i + A42i * (M[i, ] %*% X) + A53i + A54i

    if (i %% 10 == 0) cat(i / n, " ")
  }

  sum(A1vec * e + 2 * A2vec * e + A3vec * X - A4vec * X - A5vec * e)
}
