#' Compute Generalized A1-Type Variance Component
#'
#' @description
#' Calculates a variance component structured like the \eqn{A_1} term in the L3O variance
#' estimator, allowing for arbitrary vectors in the summation positions. This generalized
#' function is used to compute \eqn{A_1}, \eqn{A_2}, and \eqn{A_3} by permuting the
#' input vectors (e.g., swapping regressors \eqn{X} and residuals \eqn{e}).
#'
#' @param df Data frame. Contains the data vectors specified by \code{ipos}, \code{jpos},
#'   \code{kpos}, and \code{lpos}.
#' @param P Matrix of dimension \eqn{n \times n}. The projection matrix of instruments.
#' @param G Matrix of dimension \eqn{n \times n}. The UJIVE weighting matrix.
#' @param ipos Name of the column in \code{df} (unquoted) corresponding to the outer summation weight \eqn{i}.
#' @param jpos Name of the column in \code{df} (unquoted) corresponding to the inner summation term \eqn{j}.
#' @param kpos Name of the column in \code{df} (unquoted) corresponding to the inner summation term \eqn{k}.
#' @param lpos Name of the column in \code{df} (unquoted) corresponding to the bias correction/residual term.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = \sum_{i} v_i^{(I)} \left[ v_i^{(L)} \sum_{j \neq i} \sum_{k \neq i} G_{ij} v_j^{(J)} G_{ik} v_k^{(K)} W_{ijk} - \text{BiasCorrect}(v^{(L)}, v^{(J)}, v^{(K)}) \right]}
#' where \eqn{v^{(I)}, v^{(J)}, v^{(K)}, v^{(L)}} correspond to the vectors specified by
#' \code{ipos}, \code{jpos}, \code{kpos}, and \code{lpos} respectively. \eqn{W_{ijk}} represents
#' the Leave-Three-Out weighting derived from the annihilator matrix \eqn{M = I-P}.
#'
#' The function implements the bias correction expansion (terms \eqn{A_{12}} through \eqn{A_{15}})
#' required for consistency under many weak instruments.
#'
#' @return A numeric scalar representing the computed sum.
#'
#' @export
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
#' Calculates a variance component structurally equivalent to the \eqn{A_2} term in the
#' L3O variance estimator. This function handles asymmetric weighting structures where
#' the summation involves terms of the form \eqn{G_{ij} G_{ki}} (chaining indices),
#' differing from the symmetric \eqn{G_{ij} G_{ik}} form calculated by \code{A1type_iloop_sum}.
#'
#' @param df Data frame. Contains the data vectors specified by \code{ipos}, \code{jpos},
#'   \code{kpos}, and \code{lpos}.
#' @param P Matrix of dimension \eqn{n \times n}. The projection matrix of instruments.
#' @param G Matrix of dimension \eqn{n \times n}. The UJIVE weighting matrix.
#' @param ipos Name of the column in \code{df} (unquoted) corresponding to the outer summation weight \eqn{i}.
#' @param jpos Name of the column in \code{df} (unquoted) corresponding to the inner summation term \eqn{j}.
#' @param kpos Name of the column in \code{df} (unquoted) corresponding to the inner summation term \eqn{k}.
#' @param lpos Name of the column in \code{df} (unquoted) corresponding to the bias correction/residual term.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = \sum_{i} v_i^{(I)} \left[ v_i^{(L)} \sum_{j \neq i} \sum_{k \neq i} G_{ij} v_j^{(J)} G_{ki} v_k^{(K)} W_{ijk} - \text{BiasCorrect} \right]}
#' where \eqn{v^{(I)}, v^{(J)}, v^{(K)}, v^{(L)}} correspond to the input vectors.
#'
#' The primary distinction from \code{A1type_iloop_sum} is the use of \eqn{G_{ki}} (the \eqn{i}-th element of column \eqn{k}, or row \eqn{k} column \eqn{i})
#' in the second position, rather than \eqn{G_{ik}}. This asymmetric structure is required
#' for interaction terms in the variance of the score statistic (e.g., \eqn{A_2} and \eqn{A_5} in Yap 2025).
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @export
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
#' Calculates a variance component structurally equivalent to the \eqn{A_3} term in the
#' L3O variance estimator. This function handles "incoming" weighting structures where
#' the summation involves weights of the form \eqn{G_{ji} G_{ki}} (where indices \eqn{j}
#' and \eqn{k} both target \eqn{i}).
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param P Matrix of dimension \eqn{n \times n}. The orthogonal projection matrix of instruments.
#' @param G Matrix of dimension \eqn{n \times n}. The UJIVE weighting matrix.
#' @param ipos Name of the column in \code{df} (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column in \code{df} (unquoted) for the inner term \eqn{v^{(J)}} weighted by \eqn{G_{ji}}.
#' @param kpos Name of the column in \code{df} (unquoted) for the inner term \eqn{v^{(K)}} weighted by \eqn{G_{ki}}.
#' @param lpos Name of the column in \code{df} (unquoted) for the scalar bias correction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'
#' @details
#' This function computes the scalar:
#' \deqn{S = \sum_{i} v_i^{(I)} \left[ v_i^{(L)} \sum_{j \neq i} \sum_{k \neq i} G_{ji} v_j^{(J)} G_{ki} v_k^{(K)} W_{ijk} - \text{BiasCorrect} \right]}
#' where \eqn{W_{ijk}} is the Leave-Three-Out weighting derived from the annihilator matrix.
#'
#' This term corresponds to \eqn{A_3} in Yap (2025). It captures the variance arising from
#' the "reverse" influence of observations \eqn{j} and \eqn{k} on observation \eqn{i}
#' through the weighting matrix. This is typically used to estimate the variance of the
#' endogenous variable \eqn{X} attributed to the residuals \eqn{e}.
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
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
#' \eqn{A_4} term in the L3O variance estimator. This term captures the bias arising from
#' the squared diagonal weights \eqn{G_{ji}^2} (or \eqn{G_{ji}G_{ij}} in symmetric cases)
#' and is used to remove the positive bias introduced by the variance of the instrument
#' projection errors.
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param P Matrix of dimension \eqn{n \times n}. The orthogonal projection matrix of instruments.
#' @param G Matrix of dimension \eqn{n \times n}. The UJIVE weighting matrix.
#' @param ipos Name of the column in \code{df} (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column in \code{df} (unquoted) for the inner term \eqn{v^{(J)}} weighted by \eqn{G_{ji}^2}.
#' @param kpos Name of the column in \code{df} (unquoted) for the term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Name of the column in \code{df} (unquoted) for the residual/bias interaction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'
#' @details
#' This function computes a scalar representing the bias correction for the "own-observation"
#' variance contributions. In the notation of Yap (2025), this corresponds to \eqn{A_4}.
#'
#' Unlike \eqn{A_1} through \eqn{A_3}, which estimate signal variances and cross-covariances,
#' \eqn{A_4} specifically estimates the quantity:
#' \deqn{S = \sum_{i} v_i^{(I)} \sum_{j \neq i} G_{ji}^2 \widehat{Var}(v_j) W_{ij}}
#' adjusted for the leverage exerted by the annihilator matrix \eqn{M} to ensure unbiasedness.
#' The term involves element-wise squaring of the weighting matrix column (\code{Gi^2}).
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
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
#' \eqn{A_5} term in the L3O variance estimator. This term captures the bias arising from
#' the interaction of row and column weights \eqn{G_{ij} G_{ji}} and is required when
#' the weighting matrix \eqn{G} is asymmetric (e.g., UJIVE).
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param P Matrix of dimension \eqn{n \times n}. The orthogonal projection matrix of instruments.
#' @param G Matrix of dimension \eqn{n \times n}. The UJIVE weighting matrix.
#' @param ipos Name of the column in \code{df} (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column in \code{df} (unquoted) for the inner term \eqn{v^{(J)}} weighted by \eqn{G_{ij}G_{ji}}.
#' @param kpos Name of the column in \code{df} (unquoted) for the term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Name of the column in \code{df} (unquoted) for the residual/bias interaction term \eqn{v^{(L)}}.
#' @param noisy Logical. If \code{TRUE}, prints progress to the console.
#'
#' @details
#' This function computes a scalar representing the bias correction for variance contributions
#' arising from the asymmetry of the weighting matrix. In the notation of Yap (2025), this
#' corresponds to \eqn{A_5}.
#'
#' The term estimates:
#' \deqn{S = \sum_{i} v_i^{(I)} \sum_{j \neq i} G_{ij} G_{ji} \widehat{Var}(v_j) W_{ij}}
#' adjusted for leverage to ensure unbiasedness.
#'
#' If \eqn{G} is symmetric (e.g., \eqn{G=P}), then \eqn{G_{ij} = G_{ji}}, making this term
#' identical to \code{A4type_iloop_sum}. For asymmetric \eqn{G}, both \eqn{A_4} and \eqn{A_5}
#' must be calculated and subtracted from the total variance.
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
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
