#' Birthday Function
#'
#' @param x A numeric value representing the number of people in the group
#'
#' @returns The probability (a numeric value between 0 and 1) that at least two people in the group of size x share a birthday
#' @export
#'
#' @examples
#' birthday(10)
birthday <- function(x) {
  1-exp(lchoose(365,x) + lfactorial(x) - x*log(365))
}
