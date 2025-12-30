#' @keywords internal
"_PACKAGE"

#' @importFrom stats rnorm runif qnorm df
NULL

## Fix for Note 2: Tell R these variables are safe to use
utils::globalVariables(c(
  "Q", "QQinv", "W", "WWinv", "dPQ", "dPW",
  "J", "Z", "df", "group", "X", "Y", "e",
  "MX", "MY", "Me", "n"
))
