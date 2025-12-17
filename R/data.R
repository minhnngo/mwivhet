#' Judge data without covariates
#'
#' A simulated dataset for the judge example without covariates.
#' Generated using \code{GenData_nocov()}.
#'
#' @format A data frame with observable variables:
#' \describe{
#'   \item{group}{Group identifier}
#'   \item{X}{Treatment variable}
#'   \item{Y}{Outcome variable}
#'   \item{e}{Error term}
#'   \item{MX}{Residual of X on Z}
#'   \item{Me}{Residual of e on Z}
#'   \item{MY}{Residual of Y on Z}
#' }
"dnc"

#' Interacted data with covariates
#'
#' A simulated dataset for the QOB example with covariates.
#' Generated using \code{GenData_cov()}.
#'
#' @format A data frame with observable variables:
#' \describe{
#'   \item{group}{Group identifier}
#'   \item{groupW}{Covariate group identifier (e.g., State)}
#'   \item{X}{Treatment variable}
#'   \item{Y}{Outcome variable}
#'   \item{e}{Error term}
#'   \item{MX}{Residuals}
#'   \item{Me}{Residuals}
#'   \item{MY}{Residuals}
#' }
"dc"
