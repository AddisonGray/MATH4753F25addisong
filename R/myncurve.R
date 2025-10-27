#' myncurve
#'
#' @param mu mean
#' @param sigma standard deviation
#' @param a upper bound for shaded area
#'
#' @returns a plot with normal curve and shaded area that goes from negative infinity to x=a
#' @export
#'
#' @examples
#' myncurve(mu=0, sigma=1, a=2.5)
myncurve = function(mu, sigma,a){
  curve(dnorm(x,mean=mu,sd=sigma),
        xlim = c(mu-3*sigma, mu+3*sigma))
  list(mu = mu, sigma = sigma)
  xcurve = seq(mu-3*sigma, a, length=1000)
  ycurve = dnorm(xcurve, mean=mu, sd=sigma)
  polygon(c(xcurve, a, mu-3*sigma), c(ycurve, 0, 0))
  area <- pnorm(a, mean=mu, sd=sigma)
  list(mu = mu, sigma = sigma, area = area)
}
