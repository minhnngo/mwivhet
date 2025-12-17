#' Title
#'
#' @param X
#' @param n
#'
#' @returns
#' @noRd
#'
#' @examples
bindrowvecs <- function(X, n) {
  matrix(rep(X, n), nrow = n, byrow = TRUE)
}

#' Title
#'
#' @param X
#' @param n
#'
#' @returns
#' @noRd
#'
#' @examples
bindcolvecs <- function(X, n) {
  matrix(rep(X, n), nrow = n, byrow = FALSE)
}

#' Title
#'
#' @param ds
#' @param X
#' @param g
#'
#' @returns
#' @noRd
#'
#' @examples
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

#' Title
#'
#' @param ds
#' @param group
#'
#' @returns
#' @noRd
#'
#' @examples
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
