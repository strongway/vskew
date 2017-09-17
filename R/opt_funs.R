# function related to option calculation
# Strongway - 27 Sept. 2016
# revision 1st Feb. 2017

#' interpolate U and R
#'
#' @import dplyr
#' @import tidyr
#' @import cowplot
#'
#' Interpolate Underlying U and interest rate R at given tau and price
#' @param param SVI parameter output table, including U, R etc.
#' @param tau_interp tau of interpolation
#' @param price Spot price
#'
#' @return data.frame of interpolated U, R, mean Forward Price mF, tau, price
#' @export
#' Note: tau must be scalar, price could be an array
UR.interp <- function(param, tau_interp, price) {
  # first get difference between the price and spot
  dp <- price - param$spot[1]
  # interpolate R
  if (nrow(param) == 1){
    R = param$R
    mF1 = param$mF
  } else {
    R <-  approxExtrap(param$maturity, param$R, xout = tau_interp)$y
    # approximate the forward price
    mF1 <-  approxExtrap(param$maturity, param$mF, xout = tau_interp)$y
  }
  mF = mF1 + dp
  U <-  mF*exp(-R*tau_interp)
  UR <-  data.frame(U=U, R=R, mF = mF , tau = tau_interp, price = price)
  UR
}

#' plot skews from a surface
#' @param surf surface of a given chain
#' @param dtes DTEs for skews
#' @param strikes range of strikes (or logmoneyness)
#' @export
#' @return estimated ivs (including skews, pdf etc. )
surface.skew <- function(surf, dtes, strikes = seq(-0.2,0.1, length.out = 300)){
  ivs <- do.call(rbind, lapply(dtes, vs.options, surface = surf, x = strikes))
  ivs$dte <- factor(ivs$dte)
  print(ggplot(ivs, aes(strike, skew, color = dte)) + geom_line() )
  return(ivs)
}

#' estimate combo prices
#'
#' Estimate combo (e.g., butterfly) price, delta, theta based on a given
#' surface.
#' @param surf a given surface
#' @param dte Days to expiration
#' @param strikes strikes of the fly
#' @param positions positions of the strikes
fly.price <- function(surf, dte, strikes, positions = c(1,-2,1), spot = NULL, type = 'put',trueDelta = TRUE){
  ivs <- vs.options(surf, dte, x = strikes , spot = spot, type = type)
  price <- sum(ivs$value*positions)*100
  delta <- sum(ivs$delta*positions*100)
  theta <- sum(ivs$theta*positions*100)
  if (trueDelta){
    # calculate true delta move
    if (is.null(spot))
      spot0 = surf$spot[1]
    else
      spot0 = spot
    # move 1 dollar
    ivs1 <- vs.options(surf, dte, x = strikes , spot = spot0+1, type = type)
    price1 <- sum(ivs1$value*positions)*100
    delta_true = price1 - price
  } else {
    delta_true = NA
  }
  return(data.frame(dte,price, delta, delta_true,theta ))
}

#' estimate calendar prices
#'
#' Estimate calendar combo price, delta, theta based on a given
#' surface.
#' @param surf a given surface
#' @param dte Days to expiration
#' @param strikes strikes of the calendar
#' @param positions positions of the strikes
cal.price <- function(surf, dte, strike, positions = c(-1,1), spot = NULL, type = 'put',trueDelta = TRUE){
  ivs = data.frame()
  for (id in dte){
    ivs0 <- vs.options(surf, id, x = strike , spot = spot, type = type)
    ivs <- rbind(ivs,ivs0)
  }
  price <- sum(ivs$value*positions)*100
  delta <- sum(ivs$delta*positions*100)
  theta <- sum(ivs$theta*positions*100)
  if (trueDelta){
    # calculate true delta move
    if (is.null(spot))
      spot0 = surf$spot[1]
    else
      spot0 = spot
    # move 1 dollar
    ivs1 = data.frame()
    for (id in dte){
      ivs0 <- vs.options(surf, id, x = strike , spot = spot0+1, type = type)
      ivs1 <- rbind(ivs1,ivs0)
    }
    price1 <- sum(ivs1$value*positions)*100
    delta_true = price1 - price
  } else {
    delta_true = NA
  }
  return(data.frame(price, delta, delta_true,theta ))
}

#' plot raw chain
#'
#' Plot the raw chain
plot.rawchain <- function(chainvol, type = 'strike', drange = c(14,100)){
  pchain <- chainvol %>% filter(dte < drange[2] & dte > drange[1]) %>% mutate(dte = round(dte))
  if (type == 'strike') {
    fig <-  ggplot( pchain, aes(strike, iv*sqrt(maturity), color = as.factor(dte)))
  } else { #logmoneyness
    fig <-  ggplot( pchain, aes(logmoneyness, iv*sqrt(maturity), color = as.factor(dte)))
  }
  fig + geom_line() + geom_point(size = 1) + facet_wrap(~as.factor(date))
}

#' Obtain greeks from vertical skew
#'
#' Obtain greeks from a surface with a given tau
#' Rework on vs.options with array input
vs.greeks <- function(surf, tau, callput = 'p', k = seq(-0.5, 0.3, length.out = 1200)){
  # change to own algo, so it is faster than fOption
  para = vs.interParameter(surf, tau)
  mF = para$mF # underlying
  R = para$R # interest rate
  if (max(k) < 1) {
    X = exp(k) * mF # strike
  }
  else {
    X = k # strike
    k = log(X/mF)
  }
  ivs = vs.ivs(k, para)
  # calculate greeks
  d1 = (-k + ( ivs * ivs/2) * tau)/(ivs* sqrt(tau))
  d2 = d1 - ivs*sqrt(tau)
  Theta1 = -(mF * dnorm(d1) * ivs)/(2 * sqrt(tau))
  if (callput == "c") {
    deltas = pnorm(d1)
    thetas = (Theta1  - R * X * pnorm(+d2))* exp(-R * tau)
    prices = (mF * pnorm(d1) - X  * pnorm(d2))* exp(-R*tau)
  } else {# put
    deltas = pnorm(d1) - 1
    thetas = (Theta1  + R * X  * pnorm(-d2)) * exp(-R * tau)
    prices = (X * pnorm(-d2) - mF * pnorm(-d1))*exp(-R*tau)
  }
  gammas = dnorm(d1)/(mF * ivs * sqrt(tau))
  vegas = mF * dnorm(d1) * sqrt(tau) * exp(-R*tau)

  xt = k/sqrt(para$maturity)
  param = c(para$a, para$b, para$m, para$rho, para$sigma)
  svi <- trans.svi_fun(xt, param)
  dsvi = para$b * (para$rho + (xt - para$m)/sqrt((xt - para$m)^2 +
                                                   para$sigma^2))
  slope = para$H * svi^(para$H - 1) * dsvi/sqrt(para$maturity)
  strikeDeltas = deltas - vegas * slope/para$U[1]
  data.table(iv = ivs, strike = X, tau = tau, delta = deltas,
             gamma = gammas, theta = thetas, vega = vegas, price = prices,
             date = surf$date[1], mF =mF, spot = para$spot[1],
             slope = slope, strikeDelta = strikeDeltas)

}

#' Delta Strikes
#'
#' From Deltas to calculate the strikes from a surface and a given tau.
#' Note, by default the input Delta is a normal BS delta. It is also possible to
#' calculate strike-delta (delta - vega * skew slope) by set strikeDelta = TRUE.
#' by default, input delta are positive (even for Put side)
#' @param deltas input deltas (either BS delta or true delta)
#' @param surf surf parameters
#' @param tau maturity
#' @param callput call/put, default 'p'
#' @param minGap minimum interval in the option chain, by default is 5
#' @param strikeDelta logical option for use strikeDelta
vs.deltaStrikes <- function(deltas, surf, tau, callput = 'p', minGap = 5, strikeDelta = FALSE){
  # make sure get greeks that cover the deltas
  gk <- as.data.table(vs.greeks(surf, tau, callput))
  # add minus if it is put (default is positive)
  if (tolower(callput) == 'p' & deltas[1] > 0)
    deltas = - deltas

  output = gk[0,] # initialization
  for (idelta in deltas){
    # find closest delta
    if (strikeDelta == FALSE)
      idx = which.min(abs(idelta-gk$delta))
    else
      idx = which.min(abs(idelta-gk$strikeDelta))

    strike = gk$strike[idx]
    # find closest strike with given minGap
    strike = round(strike/5)*5
    # then interpolate
    idx = which.min(abs(strike-gk$strike))
    # find two closes taus
    if (gk$strike[idx] <= strike)
      idx = idx + 1
    if (idx>nrow(gk))
      idx = nrow(gk)
    # interpolate alpha
    alpha = (gk$strike[idx] - strike)/(gk$strike[idx]-gk$strike[idx-1])
    if (alpha>1)      alpha = 1
    if (alpha<0)      alpha = 0
    # interpolate the rest
#    curDate = gk[idx,'date']
#    gk[,c('date'):=NULL] # remove date
    inter_para = gk[idx,] - alpha *(gk[idx,]-gk[idx-1,])
#    inter_para$date = curDate
    output = rbind(output, inter_para)
  }
  output
}

#' Obtain strikes price
#'
#' Calculate strike prices and greeks from surface
#'
vs.strikes <- function(strikes, surf, tau, callput ='p'){
  para = vs.interParameter(surf, tau)
  k = log(strikes/para$U)
  gk <- vs.greeks(surf, tau, callput, k)
  gk
}

