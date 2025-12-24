#' Create a Matrix by Repeating a Vector Row-wise
#'
#' @description
#' A helper function that replicates a vector \code{X} into a matrix with \code{n} identical rows.
#' In linear algebra terms, this creates the rank-1 matrix \eqn{\mathbf{1}_n \otimes X'} (Kronecker product).
#'
#' @details
#' This is primarily used for "broadcasting" a vector across a matrix operation without using
#' slower `apply` loops. For example, subtracting a mean vector from every row of a data matrix.
#'
#' @param X Numeric vector. The vector to be repeated.
#' @param n Integer. The number of rows in the resulting matrix.
#'
#' @returns A numeric matrix of dimension \code{c(n, length(X))}, where every row is identical to \code{X}.
#' @noRd
bindrowvecs <- function(X, n) {
  matrix(rep(X, n), nrow = n, byrow = TRUE)
}

#' Create a Matrix by Repeating a Vector Column-wise
#'
#' @description
#' A helper function that replicates a vector \code{X} into a matrix of \code{n} rows.
#' Unlike \code{bindrowvecs}, this fills the matrix column-wise, ensuring that every column
#' is a recycled copy of \code{X}.
#'
#' @details
#' In the specific context of the MWIV package (where \code{length(X) == n}), this creates
#' a square \eqn{n \times n} matrix where every column is identical to the vector \code{X}.
#' Algebraically, this corresponds to the outer product \eqn{X \mathbf{1}_n'}.
#'
#' If \code{length(X) != n}, the function creates a matrix of dimension \eqn{n \times \text{length}(X)},
#' where each column consists of \code{X} repeated/recycled to fill \code{n} rows.
#'
#' @param X Numeric vector. The vector to be repeated (typically of length \eqn{n}).
#' @param n Integer. The number of rows in the resulting matrix.
#'
#' @returns A numeric matrix where every column contains the elements of \code{X}.
#' @noRd
bindcolvecs <- function(X, n) {
  matrix(rep(X, n), nrow = n, byrow = FALSE)
}

#' Compute Grouped Leave-Two-Out Sums
#'
#' @description
#' Calculates a matrix of leave-two-out sums for a specific group. This is a helper
#' function used to construct variance estimators that require removing the contributions
#' of observation pairs \eqn{(i,j)} from group-level totals.
#'
#' @param ds Data frame containing the grouping variable. The function assumes
#'   non-standard evaluation of a column named `group`.
#' @param X Numeric vector of length \eqn{n}. The values to be summed (e.g., instrument values).
#' @param g Scalar. The specific group identifier (level) to compute the sums for.
#'
#' @details
#' For a given group \eqn{g}, this function computes an \eqn{n \times n} matrix where the
#' \eqn{(i,j)}-th element is:
#' \deqn{S_{ij} = \sum_{k \in g} X_k - X_i \mathbb{I}(i \in g) - X_j \mathbb{I}(j \in g)}
#' This efficiently calculates the sum of \code{X} for group \code{g} while excluding
#' observations \eqn{i} and \eqn{j}.
#'
#' @return A matrix of dimension \eqn{n \times n}.
#'
#' @references
#' Yap, L. (2025). "Inference with Many Weak Instruments and Heterogeneity".
#' Working Paper.
#'
#' @noRd
ivectomats <- function(ds, X, g) {
  ds$group <- eval(substitute(group), ds)
  gidx <- which(ds$group == g)
  sumq <- sum(X[gidx])
  subtractq <- matrix(rep(0, nrow(ds)^2), nrow = nrow(ds))
  subtractq[gidx, ] <- X[gidx]
  subtractq[, gidx] <- subtractq[, gidx] + t(subtractq[gidx, ])
  out <- (matrix(rep(1, nrow(ds)^2), nrow = nrow(ds)) - diag(nrow(ds))) * sumq - subtractq
  out
}

#' Generate Integer Group Indices
#'
#' @description
#' Converts a grouping variable into a vector of contiguous integer indices ranging
#' from 1 to \eqn{G}, where \eqn{G} is the number of unique groups. This is a helper
#' function used to facilitate matrix indexing for group-level variance calculations.
#'
#' @param ds Data frame. The dataset containing the grouping variable.
#' @param group Name of the grouping variable (passed as an unquoted symbol) within \code{ds}.
#'   This variable defines the clusters (e.g., judges, examiners, or time periods).
#'
#' @details
#' The function extracts the specified \code{group} column from the data frame using
#' non-standard evaluation. It then maps each unique value of the group to an integer
#' index corresponding to its order of appearance in \code{unique(group)}.
#'
#' This transformation is functionally equivalent to \code{as.numeric(as.factor(group))},
#' but ensures explicit handling within the provided data frame context.
#'
#' @return An integer vector of length \eqn{n}, containing values in \eqn{\{1, \dots, G\}}.
#'
#' @export
#'
#' @examples
#' data <- data.frame(judge_id = c("JudgeA", "JudgeB", "JudgeA", "JudgeC"))
#' Getgroupindex(data, judge_id)
Getgroupindex <- function(ds, group) {
  ds$group <- eval(substitute(group), ds)
  ds$groupidx <- 0
  gvals <- unique(ds$group)
  for (g in unique(ds$group)) {
    ds$groupidx[ds$group == g] <- which(unique(ds$group) == g)
  }
  ds$groupidx
}

#' Title
#'
#' @param group
#' @param groupW
#' @param n
#'
#' @returns
#' @export
#'
#' @examples
GetGP <- function(group, groupW, n) {
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

  list(G = G, P = P, Z = Z, W = W)
}
