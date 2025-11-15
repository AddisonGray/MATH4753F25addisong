#' ntickets for Overbooking Problem
#'
#' @param N the number of seats on the flight
#' @param gamma the probability that overbooking will occur
#' @param p the probability that a passenger will be a "show" for the flight
#' @param n_extra how many extra tickets to check beyond N
#'
#' @returns Named list with N, p, gamma, nd, and nc as well as two plots that show how many tickets should be sold (discrete and continuous)
#' @export
#'
#' @examples
#' ntickets(N=400, gamma=0.02, p=0.95)
ntickets <- function(N, gamma, p, n_extra = 25) {
  #N = number of seats
  #gamma = probability of overbooking
  #p = probability of a "show" for the flight
  #n_extra = how many extra tickets to check beyond N

  # candidate ticket numbers
  n_vals <- seq(N, N + n_extra, 1)

  #Binomial (Discrete)
  f_disc <- 1 - gamma - pbinom(N, n_vals, p)
  nd <- n_vals[which.min(abs(f_disc))] #find the first n where f_disc >= 0 (like we did in class)

  #Normal Approximation (Continuous)
  f_cont <- function(n) {
    1 - gamma - pnorm(N+0.5, mean=n*p, sd=sqrt(n*p*(1-p))) #1-gamma-pnorm()=0
  }
  nc <- uniroot(f_cont, lower=N, upper=N+n_extra)$root #find the root numerically where f_cont(n)=0

  #Plotting
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  par(mfrow=c(2, 1), mar=c(4, 4, 3, 2))

  #Discrete Plot
  plot(n_vals, f_disc, type="b", pch=19, cex=0.7,
       xlab="n (tickets sold)", ylab="Objective",
       main=paste("Objective function vs n (Discrete), N=400, n=", nd))
  abline(h=0, col="Red", lwd=2)
  abline(v=nd, col="Red", lwd=2)

  #Continuous Plot
  f_vals_cont <- sapply(n_vals, f_cont)
  plot(n_vals, f_vals_cont, type="l", lwd=2,
       xlab="n (tickets sold)", ylab="Objective",
       main=paste("Objective function vs n (Continuous), N=400, n=", round(nc, 3)))
  abline(h=0, col="Blue", lwd=2)
  abline(v=nc, col="Blue", lwd=2)

  #Return a named list
  result <- list(N = N, p = p, gamma = gamma, nd = nd, nc = nc)
  print(result)
}
