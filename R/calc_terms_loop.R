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

#' Compute A1-Type Variance Component (Explicit IJ-Loop)
#'
#' @description
#' Calculates the \eqn{A_1} variance component using an explicit double loop over indices
#' \eqn{i} and \eqn{j}, with vectorization over \eqn{k}. This implementation constructs
#' instrument weights row-by-row, avoiding the storage of full \eqn{N \times N} projection
#' matrices. It is primarily used for verification or memory-constrained environments.
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param ipos Name of the column (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column (unquoted) for the inner term \eqn{v^{(J)}}.
#' @param kpos Name of the column (unquoted) for the inner term \eqn{v^{(K)}}.
#' @param lpos Name of the column (unquoted) for the bias correction term \eqn{v^{(L)}}.
#' @param IdPQ Numeric vector. The diagonal elements of the full projection matrix \eqn{P} (instruments + covariates).
#' @param IdPW Numeric vector. The diagonal elements of the covariate projection matrix \eqn{P_W}.
#' @param noisyi Logical. If \code{TRUE}, prints progress for the outer loop \eqn{i}.
#' @param noisyj Logical. If \code{TRUE}, prints progress for the inner loop \eqn{j}.
#'
#' @details
#' This function computes the same scalar quantity as \code{A1type_iloop_sum} but uses a
#' different computational strategy.
#'
#' \strong{Key Differences:}
#' \itemize{
#'   \item \strong{Memory:} Does not require passing or storing the full matrices \eqn{P} or \eqn{G}.
#'   It reconstructs the rows \eqn{P_{i\cdot}} and \eqn{G_{i\cdot}} on the fly using the global
#'   objects \code{Q}, \code{QQinv}, \code{W}, and \code{WWinv}.
#'   \item \strong{Complexity:} \eqn{O(N^3)} explicitly. While mathematically equivalent,
#'   this function is generally slower in R than the fully vectorized \code{iloop} version
#'   due to the interpreted \eqn{j}-loop.
#' }
#'
#' The function relies on the existence of \code{Q}, \code{QQinv}, \code{W}, and \code{WWinv}
#' in the parent environment to construct the projection weights.
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
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

#' Compute A2-Type Variance Component (Explicit IJ-Loop)
#'
#' @description
#' Calculates the \eqn{A_2} variance component using an explicit double loop over indices
#' \eqn{i} and \eqn{j}, with vectorization over \eqn{k}. This implementation handles the
#' asymmetric weighting structure (chaining "incoming" and "outgoing" weights) by constructing
#' the necessary row and column vectors of the weighting matrix \eqn{G} on the fly.
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param ipos Name of the column (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column (unquoted) for the inner term \eqn{v^{(J)}}.
#' @param kpos Name of the column (unquoted) for the inner term \eqn{v^{(K)}}.
#' @param lpos Name of the column (unquoted) for the bias correction term \eqn{v^{(L)}}.
#' @param IdPQ Numeric vector. The diagonal elements of the full projection matrix \eqn{P} (instruments + covariates).
#' @param IdPW Numeric vector. The diagonal elements of the covariate projection matrix \eqn{P_W}.
#' @param noisyi Logical. If \code{TRUE}, prints progress for the outer loop \eqn{i}.
#' @param noisyj Logical. If \code{TRUE}, prints progress for the inner loop \eqn{j}.
#'
#' @details
#' This function calculates the same quantity as \code{A2type_iloop_sum} but is optimized
#' for memory-constrained environments where the full \eqn{N \times N} matrix \eqn{G} cannot
#' be stored.
#'
#' \strong{Asymmetry Handling:}
#' Unlike the \eqn{A_1} computation, which uses symmetric weights \eqn{G_{ik}}, this function
#' explicitly constructs:
#' \itemize{
#'   \item \code{Girow}: The \eqn{i}-th row of \eqn{G} (outgoing weights \eqn{G_{ik}}).
#'   \item \code{Gicol}: The \eqn{i}-th column of \eqn{G} (incoming weights \eqn{G_{ki}}).
#' }
#' It computes the interaction terms involving the chain \eqn{j \to i \to k} (or \eqn{G_{ji} G_{ik}})
#' and the specific bias correction \eqn{A_{25}} involving \eqn{G_{ij} G_{ji}}.
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
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

#' Compute A3-Type Variance Component (Explicit IJ-Loop)
#'
#' @description
#' Calculates the \eqn{A_3} variance component using an explicit double loop over indices
#' \eqn{i} and \eqn{j}, with vectorization over \eqn{k}. This implementation handles the
#' "incoming" weighting structure (weights \eqn{G_{ji} G_{ki}} targeting \eqn{i}) by
#' constructing the necessary column vectors of the weighting matrix \eqn{G} on the fly.
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param ipos Name of the column (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column (unquoted) for the inner term \eqn{v^{(J)}}.
#' @param kpos Name of the column (unquoted) for the inner term \eqn{v^{(K)}}.
#' @param lpos Name of the column (unquoted) for the bias correction term \eqn{v^{(L)}}.
#' @param IdPQ Numeric vector. The diagonal elements of the full projection matrix \eqn{P}.
#' @param IdPW Numeric vector. The diagonal elements of the covariate projection matrix \eqn{P_W}.
#' @param noisyi Logical. If \code{TRUE}, prints progress for the outer loop \eqn{i}.
#' @param noisyj Logical. If \code{TRUE}, prints progress for the inner loop \eqn{j}.
#'
#' @details
#' This function calculates the same quantity as \code{A3type_iloop_sum} but uses the
#' memory-efficient strategy of reconstructing weights row-by-row (or column-by-column).
#'
#' \strong{Weight Construction:}
#' The function constructs the \eqn{i}-th column of the weighting matrix \eqn{G} as:
#' \deqn{G_{\cdot i} = \frac{P_{\cdot i}}{\text{diag}(M)} - \frac{P_{W \cdot i}}{\text{diag}(M_W)}}
#' This vector \code{Gi} corresponds to the weights \eqn{G_{ki}} for all \eqn{k}, which represents
#' the influence of observation \eqn{k} on observation \eqn{i}.
#'
#' The term computed is:
#' \deqn{S = \sum_{i} v_i^{(I)} \sum_{j \neq i} \sum_{k \neq i} G_{ji} v_j^{(J)} G_{ki} v_k^{(K)} W_{ijk}}
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
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

#' Compute A4-Type Variance Component (Explicit IJ-Loop)
#'
#' @description
#' Calculates the \eqn{A_4} variance component ("own-variance" bias correction) using an
#' explicit double loop over indices \eqn{i} and \eqn{j}. This implementation handles the
#' squared incoming weights (\eqn{G_{ji}^2}) by constructing the necessary column vectors
#' of the weighting matrix \eqn{G} on the fly.
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param ipos Name of the column (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column (unquoted) for the inner term \eqn{v^{(J)}} weighted by \eqn{G_{ji}^2}.
#' @param kpos Name of the column (unquoted) for the term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Name of the column (unquoted) for the residual/bias interaction term \eqn{v^{(L)}}.
#' @param IdPQ Numeric vector. The diagonal elements of the full projection matrix \eqn{P}.
#' @param IdPW Numeric vector. The diagonal elements of the covariate projection matrix \eqn{P_W}.
#' @param noisyi Logical. If \code{TRUE}, prints progress for the outer loop \eqn{i}.
#' @param noisyj Logical. If \code{TRUE}, prints progress for the inner loop \eqn{j}.
#'
#' @details
#' This function calculates the same quantity as \code{A4type_iloop_sum} but uses a
#' memory-efficient strategy.
#'
#' \strong{Term Definition:}
#' The function estimates the bias component:
#' \deqn{S = \sum_{i} v_i^{(I)} \sum_{j \neq i} G_{ji}^2 \widehat{Var}(v_j) W_{ij}}
#' adjusted for leverage via the Sherman-Morrison expansion.
#'
#' \strong{Computational Strategy:}
#' \itemize{
#'   \item \strong{Outer Loop (i):} Constructs the column vector \eqn{G_{\cdot i}} (representing weights \eqn{G_{ki}}).
#'   Squaring this vector gives the "own-variance" weights \eqn{G_{ki}^2}.
#'   \item \strong{Inner Loop (j):} Iterates through specific observations to apply the
#'   Leave-Three-Out correction logic for the variance estimator.
#' }
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
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

#' Compute A5-Type Variance Component (Explicit IJ-Loop)
#'
#' @description
#' Calculates the \eqn{A_5} variance component ("asymmetry" bias correction) using an
#' explicit double loop over indices \eqn{i} and \eqn{j}. This implementation handles the
#' cross-product of row and column weights (\eqn{G_{ik} G_{ki}}) by constructing the
#' necessary vectors on the fly.
#'
#' @param df Data frame. Contains the data vectors specified by the position arguments.
#' @param ipos Name of the column (unquoted) for the outer summation weight \eqn{v^{(I)}}.
#' @param jpos Name of the column (unquoted) for the inner term \eqn{v^{(J)}} weighted by \eqn{G_{ik}G_{ki}}.
#' @param kpos Name of the column (unquoted) for the term interacting with leverage \eqn{v^{(K)}}.
#' @param lpos Name of the column (unquoted) for the residual/bias interaction term \eqn{v^{(L)}}.
#' @param IdPQ Numeric vector. The diagonal elements of the full projection matrix \eqn{P}.
#' @param IdPW Numeric vector. The diagonal elements of the covariate projection matrix \eqn{P_W}.
#' @param noisyi Logical. If \code{TRUE}, prints progress for the outer loop \eqn{i}.
#' @param noisyj Logical. If \code{TRUE}, prints progress for the inner loop \eqn{j}.
#'
#' @details
#' This function calculates the same quantity as \code{A5type_iloop_sum} using a
#' memory-efficient strategy.
#'
#' \strong{Asymmetry Bias:}
#' This term captures the bias arising from the interaction of "outgoing" and "incoming"
#' weights. It estimates:
#' \deqn{S = \sum_{i} v_i^{(I)} \sum_{j \neq i} G_{ij} G_{ji} \widehat{Var}(v_j) W_{ij}}
#'
#' \strong{Computational Strategy:}
#' \itemize{
#'   \item \strong{Outer Loop (i):} Constructs both \code{Girow} (\eqn{G_{i\cdot}}) and
#'   \code{Gicol} (\eqn{G_{\cdot i}}). The element-wise product of these vectors gives the
#'   weight interaction \eqn{G_{ik} G_{ki}} required for this variance component.
#'   \item \strong{Inner Loop (j):} Iterates through observations to apply the Leave-Three-Out
#'   variance correction logic.
#' }
#'
#' @return A numeric scalar.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity". Working Paper.
#'
#' @noRd
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
