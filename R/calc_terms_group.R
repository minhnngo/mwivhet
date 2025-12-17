#' Title
#'
#' @param df
#' @param group
#' @param groupW
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
A1type_sum <- function(df, group, groupW, ipos, jpos, kpos, lpos, noisy = FALSE) {
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, max(df$group))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  iteration <- 1
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]

    ds$numingrp <- 0
    for (j in unique(ds$group)) {
      ds$numingrp[ds$group == j] <- sum(ds$group == j)
    }
    ds <- ds[ds$numingrp >= 3, ]
    if (nrow(ds) == 0) {
      for (g in unique(ds$group)) {
        A11vecs[g] <- A12vecs[g] <- A13vecs[g] <- A14vecs[g] <- A15vecs[g] <- 0
      }
    } else {
      ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
      ds$groupidx <- Getgroupindex(ds, group)
      ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
      # ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))

      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)
      # PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)

      ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
      PW <- ZWmat / (length(ds$groupW))

      # calculate values specific to this subset
      Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))
      Ps <- PQ
      Ms <- diag(nrow(ds)) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      for (g in unique(ds$group)) {
        if (nrow(ds[ds$group == g, ]) <= 3) {
          A11vecs[g] <- A12vecs[g] <- A13vecs[g] <- A14vecs[g] <- A15vecs[g] <- 0
        } else {
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
          Mes <- matrix(ds$lpos, ncol = 1)

          A11vecs[g] <- t(ds$jpos * Gis) %*% (D2D3is * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (ds$kpos * Gis)
          A12vecs[g] <- t(ds$jpos * ds$lpos * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * dMs) -
            t(ds$jpos * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * ds$lpos)
          A13vecs[g] <- t(ds$jpos * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
          A14vecs[g] <- t(ds$jpos * ds$lpos * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
          A15vecs[g] <- t(ds$lpos * ds$ipos) %*% (recD2s * Gs^2 * Pgs) %*% (ds$jpos * ds$kpos * dMs) -
            t(ds$ipos) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$lpos * ds$jpos * ds$kpos)
        }
      }
    }
    if (noisy) {
      cat(iteration, "of", max(df$groupW), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs)

  sum(ret)
}

#' Title
#'
#' @param df
#' @param group
#' @param groupW
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
A4type_sum <- function(df, group, groupW, ipos, jpos, kpos, lpos, noisy = FALSE) {
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, max(df$group))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  iteration <- 1
  for (s in unique(df$groupW)) {
    ds <- df[df$groupW == s, ]
    ds$numingrp <- 0
    for (j in unique(ds$group)) {
      ds$numingrp[ds$group == j] <- sum(ds$group == j)
    }
    ds <- ds[ds$numingrp >= 3, ]
    if (nrow(ds) == 0) {
      for (g in unique(ds$group)) {
        A41vecs[g] <- A42vecs[g] <- A43vecs[g] <- A44vecs[g] <- 0
      }
    } else {
      ZQ <- matrix(0, nrow = length(ds$group), ncol = length(unique(ds$group)))
      ds$groupidx <- Getgroupindex(ds, group)
      ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1
      # ZW <- matrix(1, nrow=length(ds$groupW), ncol=length(unique(ds$groupW)))

      PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)
      # PW <- ZW %*% solve(t(ZW)%*% ZW) %*% t(ZW)

      ZWmat <- matrix(1, nrow = length(ds$groupW), ncol = length(ds$groupW))
      PW <- ZWmat / (length(ds$groupW))

      # calculate values specific to this subset
      Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ))) -
        diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PW), .99))) %*% (PW - diag(diag(PW)))
      Ps <- PQ
      Ms <- diag(nrow(ds)) - Ps
      dMs <- matrix(diag(Ms), ncol = 1)
      D2s <- dMs %*% t(dMs) - Ms * Ms
      recD2s <- 1 / D2s
      diag(recD2s) <- 0

      for (g in unique(ds$group)) {
        if (nrow(ds[ds$group == g, ]) <= 3) {
          A41vecs[g] <- A42vecs[g] <- A43vecs[g] <- A44vecs[g] <- 0
        } else {
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
          Mes <- matrix(ds$lpos, ncol = 1)

          A41vecs[g] <- t(Gis^2 * ds$jpos * dMs * ds$lpos * recD2is) %*% (D2D3is * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos) -
            t(Gis^2 * ds$jpos * ds$lpos * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$kpos)) * ivectomats(ds, ds$ipos, g)) %*% (Mis)
          A42vecs[g] <- t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis^2 * ds$kpos) -
            t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) -
            t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) +
            t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (dMs * ds$kpos)
          A43vecs[g] <- t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Mis^2 * ds$kpos * ds$lpos) -
            t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos * ds$lpos) -
            t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * dMs, g)) %*% (Mis * ds$kpos * ds$lpos) +
            t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * dMs, g)) %*% (ds$kpos * ds$lpos)
          A44vecs[g] <- t(ds$ipos * ds$kpos * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$lpos * ds$jpos) -
            t(ds$ipos * ds$kpos * ds$lpos) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$jpos)
        }
      }
    }
    if (noisy) {
      cat(iteration, "of", max(df$groupW), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A41vecs + A42vecs + A43vecs + A44vecs)

  sum(ret)
}


#' Title
#'
#' @param df
#' @param groupZ
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
A1type_sum_nocov <- function(df, groupZ, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  A11vecs <- A12vecs <- A13vecs <- A14vecs <- A15vecs <- rep(0, max(df$groupZ))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  iteration <- 1
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ == s, ]
    ZQ <- matrix(0, nrow = length(ds$groupZ), ncol = length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds, groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms), ncol = 1)
    D2s <- dMs %*% t(dMs) - Ms * Ms
    recD2s <- 1 / D2s
    diag(recD2s) <- 0

    g <- s

    repidx <- min(which(ds$groupZ == g)) # representative index
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
    Mes <- matrix(ds$lpos, ncol = 1)

    A11vecs[g] <- t(ds$jpos * Gis) %*% (D2D3is * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (ds$kpos * Gis)
    A12vecs[g] <- t(ds$jpos * ds$lpos * Gis * Mis) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * dMs) -
      t(ds$jpos * Gis * Mis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Gis * ds$kpos * ds$lpos)
    A13vecs[g] <- t(ds$jpos * dMs * Gis) %*% (recD3is * (ones %*% t(Mes)) * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
    A14vecs[g] <- t(ds$jpos * ds$lpos * Gis) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (ds$kpos * Gis * Mis)
    A15vecs[g] <- t(ds$lpos * ds$ipos) %*% (recD2s * Gs^2 * Pgs) %*% (ds$jpos * ds$kpos * dMs) -
      t(ds$ipos) %*% (Gs^2 * Ms * recD2s * Pgs) %*% (ds$lpos * ds$jpos * ds$kpos)
    if (noisy) {
      cat(iteration, "of", max(df$groupZ), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A11vecs - A12vecs - A13vecs + A14vecs + A15vecs)

  sum(ret)
}

#' Title
#'
#' @param df
#' @param groupZ
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
A4type_sum_nocov <- function(df, groupZ, ipos, jpos, kpos, lpos, noisy = FALSE) {
  df$groupZ <- eval(substitute(groupZ), df)
  A41vecs <- A42vecs <- A43vecs <- A44vecs <- rep(0, max(df$groupZ))

  df$ipos <- eval(substitute(ipos), df)
  df$jpos <- eval(substitute(jpos), df)
  df$kpos <- eval(substitute(kpos), df)
  df$lpos <- eval(substitute(lpos), df)

  iteration <- 1
  for (s in unique(df$groupZ)) {
    ds <- df[df$groupZ == s, ]
    ZQ <- matrix(0, nrow = length(ds$groupZ), ncol = length(unique(ds$groupZ)))
    ds$groupidx <- Getgroupindex(ds, groupZ)
    ZQ[cbind(seq_along(ds$groupidx), ds$groupidx)] <- 1

    PQ <- ZQ %*% solve(t(ZQ) %*% ZQ) %*% t(ZQ)

    # calculate values specific to this subset
    Gs <- diag(1 / (diag(diag(nrow(ds))) - pmin(diag(PQ), .99))) %*% (PQ - diag(diag(PQ)))
    Ps <- PQ
    Ms <- diag(nrow(ds)) - Ps
    dMs <- matrix(diag(Ms), ncol = 1)
    D2s <- dMs %*% t(dMs) - Ms * Ms
    recD2s <- 1 / D2s
    diag(recD2s) <- 0

    g <- s

    repidx <- min(which(ds$groupZ == g)) # representative index
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
    Mes <- matrix(ds$lpos, ncol = 1)

    A41vecs[g] <- t(Gis^2 * ds$jpos * dMs * ds$lpos * recD2is) %*% (D2D3is * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos) -
      t(Gis^2 * ds$jpos * ds$lpos * recD2is) %*% (D2D3is * Ms * (ones %x% t(ds$kpos)) * ivectomats(ds, ds$ipos, g)) %*% (Mis)
    A42vecs[g] <- t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis^2 * ds$kpos) -
      t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) -
      t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * (ones %*% t(dMs)) * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (Mis * ds$kpos) +
      t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * ds$lpos, g)) %*% (dMs * ds$kpos)
    A43vecs[g] <- t(Gis^2 * ds$jpos * dMs * Mis * recD2is) %*% (recD3is * ivectomats(ds, ds$ipos, g)) %*% (Mis^2 * ds$kpos * ds$lpos) -
      t(Gis^2 * ds$jpos * Mis^2 * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos, g)) %*% (Mis * ds$kpos * ds$lpos) -
      t(Gis^2 * ds$jpos * dMs * recD2is) %*% (recD3is * Ms * ivectomats(ds, ds$ipos * dMs, g)) %*% (Mis * ds$kpos * ds$lpos) +
      t(Gis^2 * ds$jpos * Mis * recD2is) %*% (recD3is * Ms^2 * ivectomats(ds, ds$ipos * dMs, g)) %*% (ds$kpos * ds$lpos)
    A44vecs[g] <- t(ds$ipos * ds$kpos * dMs) %*% (recD2s * Gs^2 * Pgs) %*% (ds$lpos * ds$jpos) -
      t(ds$ipos * ds$kpos * ds$lpos) %*% (recD2s * Gs^2 * Pgs) %*% (Mis * ds$jpos)
    if (noisy) {
      cat(iteration, "of", max(df$groupZ), "done. ")
      iteration <- iteration + 1
    }
  }

  ret <- (A41vecs + A42vecs + A43vecs + A44vecs)

  sum(ret)
}
