#' Log-normal vol by SABR
#'
#' Calculate log-normal vol based on sabr model parameters
#' @param f Forward
#' @param k strike
#' @param alpha alpha parameter
#' @param beta beta parameter
#' @param rho rho parameter
#' @param volvol volvol parameter
#' @export
sabrVol <- function(f,k,t,alpha,beta,rho,volvol){
  eps = 1e-05
  logfk = log(f / k)
  fkbeta = (f*k)**(1 - beta)
  a = (1 - beta)^2 * alpha^2 / (24 * fkbeta)
  b = 0.25 * rho * beta * volvol * alpha / fkbeta^0.5
  c = (2 - 3*rho^2) * volvol**2 / 24
  d = fkbeta^0.5
  v = (1 - beta)^2 * logfk^2 / 24
  w = (1 - beta)^4 * logfk^4 / 1920
  z = volvol * fkbeta^0.5 * logfk / alpha
  # if |z| > eps
  q =  (1 - 2*rho*z + z^2)^.5 + z - rho
  p = suppressWarnings(log(q/(1-rho)))
  vz =  alpha * z * (1 + (a + b + c) * t)/(d * (1 + v + w) * p)
  # if |z| <= eps
  v0 = alpha * (1 + (a + b + c) * t) / (d * (1 + v + w))
  ifelse(abs(z)>eps, vz, v0)
}

# find Alpha parameter
findAlpha <- function(f,t,atmvol,b,r,v){
  # complete solution based on West (2005)
  c0 <- -atmvol*(f^(1-b))
  c1 <- (1+((2-3*(r^2))*(v^2)*t/24))
  c2 <- (r*b*v*t)/(4*f^(1-b))
  c3 <- (((1-b)^2)*t)/(24*f^(2-2*b))
  roots <- polyroot(c(c0,c1,c2,c3))
  # remove any complext or negative roots
  roots <- Re(roots[which(Im(roots)==0 & Re(roots)>0)])
  if(length(roots)==0) roots <- 0
  # pick the minimum value
  return(min(roots))
  # apporximation from Hagan (2003) 3.1a, but very fast and stable
  #a <- atmvol*(f^(1-b))
  # return(a)
}

# optimization function
squared_error <- function(params,strikes,IV,fwd,t,a, b){
  N <- NROW(strikes)
  err <- rep(0,N)
  f <- fwd
  r <- params[1]
  v <- params[2]
  err = (sabrVol(f,strikes,t,a,b,r,v) - IV)^2
  y <- sum(err)

  # penalities for non-admissible values of r and v
  if(abs(r)>1) y <- Inf
  if(v<0) y <- Inf

  return(y)
}

#' Calibrate SABR model
#'
#' Calibrate SABR model. Note the parameters t and b must be scalar.
#' @param fwd forward
#' @param strikes strikes
#' @param ivs implied volatilities
#' @param t tau
#' @param b beta parameter from SABR
#' @export
calibrate.sabrvol <- function(fwd, strikes, ivs, t, b, v0){

  params <- c(-0.5,0.5)
  # apporximation from Hagan (2003) 3.1a, but very fast and stable
  #a <- atmvol*(f^(1-b))
  a <- v0*(fwd^(1-b))
  fit <- stats::optim(params,fn=squared_error,
                      strikes=strikes,IV=ivs,f=fwd,t=t,b=b,a = a)
  err = sqrt(fit$value)/NROW(strikes)
  if(err > 1e-03){
    rms <- round(err,6)
    warning(paste("Average error:",rms,". Goodness of fit is below threshold."))
  }
  r = fit$par[1]
  v=fit$par[2]
  vol <- data.frame( a = a, b=b,r=r,v=v, t = t, err = err)
  return(vol)
}

#' Calibrate SABR
#'
#' calibrate SABR model and return the SABR parameters.
#' @param df a data.frame containing strikes, ivs, forwards, tau
#' @param b beta
#' @export
calibrate.sabr <- function(df,b=0.5){
  # using relative strikes,given that sig(X,F) = sig(S/F,1)
  # with relative strikes, parameter alpha -> vol_atm
  # estimate atm vol
  fit0 = lm(iv ~ logmoneyness + I(logmoneyness^2) + I(logmoneyness^3), data = df)
  v0 = as.numeric(fit0$coefficients[1]) # atm vol

  k = df$strike/df$mF[1]
  ivs = df$iv
  t = df$maturity[1]

  para = calibrate.sabrvol(1,k,ivs,t,b,v0)
  return(para)
}

#' volatility surface
#'
#' Volatility surface based on SABR model
#' @param sdf data.frame of sabr parameters
#' @param kr relative strike k range
#' @param tr tau range
#' @export
plot.surface <- function(sdf, kr = c(0.7, 1.3), tr = c(10,180)){
  ks = seq(kr[1],kr[2],by = 0.005)
  ts = unique(sdf$dte)
  ts = ts[ts>tr[1] & ts <=tr[2]]
  ivs = matrix(NROW(ts)*NROW(ks),nrow = NROW(ts),ncol = NROW(ks))
  for (i in 1:NROW(ts)){
    v = sdf[i,]
    ivs[i,] = sabrVol(1,ks,v$t,v$a,v$b,v$r,v$v)
  }
  # plot surface using plot_ly
  plot_ly(x = ks, y = ts, z = ivs) %>% add_surface() %>%
    plotly::layout(title = 'Volatility Surface', scene = list(xaxis = list(title ='K/F'),
           yaxis = list(title = 'Tau'), zaxis = list(title = 'Vol')))
}

#' Term structure of vol surface
#'
#' Visualize the term structure (Vol ATM) from SABR model
#' @param sdf data.frame of sabr model
#' @export
plot.termstructure <- function(sdf, tr = c(14,180)){
  sdf %>% filter(dte >=tr[1], dte <=tr[2]) -> sdf
  # fittd curve
  fita = lm(a ~ dte + I(dte^2) + I(dte^0.5), data = sdf)
  fitr = lm(log(1.1+r) ~ dte + I(dte^2) , data = sdf)
  fitv = lm(log(v) ~ dte + I(dte^2), data = sdf)
  ndf = data.frame(dte = seq(tr[1],tr[2])) # new data
  na = predict(fita,ndf)
  nr = exp(predict(fitr,ndf)) - 1.1
  nv = exp(predict(fitv,ndf))

  figa = sdf %>% plot_ly(x = ~dte, y = ~a, type = 'scatter',
                        mode = 'markers', alpha = 0.3, name = 'V0') %>%
    plotly::layout(xaxis = list(title = 'Tau'),
           yaxis = list(title = 'Volatility')) %>%
      add_trace(x = ndf$dte,  y = na,
                           mode = 'lines', name = 'V0',alpha = 1 )

  # Add SABR other parameters rho, volvol
  figr = sdf %>% plot_ly(x = ~dte, y = ~r, type = 'scatter',
                        mode = 'markers', alpha = 0.3, name = 'rho') %>%
    plotly::layout(xaxis = list(title = 'Tau'),
           yaxis = list(title = 'rho'))%>%
    add_trace(x = ndf$dte,  y = nr,
              mode = 'lines', name = 'rho',alpha = 1 )
  figv = sdf %>% plot_ly(x = ~dte, y = ~v, type = 'scatter',
                         mode = 'markers', alpha = 0.3, name = 'volvol') %>%
    plotly::layout(xaxis = list(title = 'Tau'),
           yaxis = list(title = 'Volvol'))  %>%
    add_trace(x = ndf$dte,  y = nv,
              mode = 'lines', name = 'volvol',alpha = 1 )
  subplot(figa,subplot(figr,figv,nrows = 2, titleX = TRUE,
                       titleY = TRUE))
}

#' Fit Term structure over tau
#'
#' fit non-linear relation of v0(t), rho(t), and volvol(t)
#' @param sdf data.frame of sabr model
#' @export
fit.termstructure <- function(sdf, tr = c(14,180), b = 0.5){
  sdf %>% filter( dte >=tr[1], dte <=tr[2]) -> sdf
  # fittd curve
  fita = lm(a ~ dte + I(dte^2) + I(dte^0.5), data = sdf)
  fitr = lm(log(1.1+r) ~ dte + I(dte^2) , data = sdf)
  fitv = lm(log(v) ~ dte + I(dte^2), data = sdf)
  ndf = data.frame(dte = seq(tr[1],tr[2])) # new data
  ndf$a = predict(fita,ndf)
  ndf$r = exp(predict(fitr,ndf)) - 1.1
  ndf$v = exp(predict(fitv,ndf))
  ndf$b = b
  ndf$t = ndf$dte/365
  return(ndf)
}

#' Visualize SABR fit
#'
#' Visualize SABR goodness of fit
#' @param sdf data.frame of sabr model
#' @param chain raw chain
#' @param rt range of tau
plot.sabrfit <- function(sdf, chain, tr = c(0.1, 1)){
  ks = seq(0.6,1.2,by = 0.005)
  sdf %>% filter(t>tr[1], t<=tr[2]) -> sdf
  ts = unique(sdf$t)
  for (i in 1:NROW(ts)){
    v = sdf[i,]
    fwd = v$U*exp(v$R*v$t)
    strike = ks*fwd
    ivs = sabrVol(1,ks,v$t,v$a,v$b,v$r,v$v)
    tb = tibble(k = ks, iv = ivs, t = v$t, strike = strike)
    if (i==1)
      sa = tb
    else
      sa = rbind(sa, tb)
  }
  chain %>% filter(dte/365>tr[1], dte/365<=tr[2]) %>%
    unnest() -> ch
  ggplot(sa, aes(strike, iv, color = t, group = t)) + geom_line() +
    geom_point(data = ch, aes(strike, iv, color = maturity, group = maturity))
}

