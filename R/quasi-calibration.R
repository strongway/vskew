# Using quai-explicit calibration model for SVI model
# strongway 23.02.2017

#' svi curve function
#'
#' Modified SVI function with fraction exponent H
#' modified version: svi^H, additional parameter
#' @param k log moneyness
#' @param param 6 parameters: a, b,m, rho, sigma, H
#' @return estimated values according to the curve
#' @export
svi_fun <- function(k,param){
  a = param[1];
  b = param[2];
  m = param[3];
  rho = param[4];
  sigma = param[5];
  H = param[6];
  result = (a + b*(rho*(k-m) + sqrt((k-m)^2 + sigma^2)))^H
}

#' Inner optimization for parameter: c,d,a
#'
#' optimization of the inner parameters
#' let y(k) = (k - m)/sig
#' svi w = a + b*sig*(rho*y(k) + sqrt(y(x)^2+1))
#'        = a + dy(x) + cz(x)
#' d = roh*b*sig, c = b*sig, z = sqrt(y(x)^2+1)
#' Constrains:
#'     0 <= c <= 4sig
#'     |d| <= c
#'     |d| <= 4sig -c
#'     0 <= a <= max_i(w_i)
#' @param y (k-m)/sigma
#' @param w  real time-weighted iv data, usually iv*sqrt(tau)
#' @param sig sigma
#' @param H fractional exponent
#' @return a,b,rho parameters
#' @export
fit_inner <- function(y, w, sig, H){

  # constrain optimization for (a, d, c, H)
  fitfun <- function(x) {
    yhat = (x[1] + x[2]*y + x[3]*sqrt(y^2+1))^H
    sum(( yhat - w )^2)
  }

  # gradient of fitfun
#  grr <- function(x){ # gradient of fitfun
#    y1 = (x[1]+x[2]*y + x[3]*sqrt(y^2+1))
#    y2 = y1^(H-1)
#    y3 = (y1^(H) - w)/4
#    c(mean(H*y2*y3), mean(H*y2*y*y3),
#      mean(H*y2*sqrt(y^2+1)*y3))
#  }

  # constrains (see above)
  ui = rbind(c(0,0,1),c(0,0,-1),
             c(0,-1,1), c(0, 1, 1),
             c(1, 0, 0), c(-1, 0, 0))
  ci = c(0, -4*sig, 0, 0, 0, -max(w))

  # initial guess
  x0 = c(max(w)/2, -sig, 2*sig)

  # constrained optimization
#  op <- constrOptim(x0,fitfun, grr, ui, ci, method = "BFGS")
  op <- constrOptim(x0,fitfun, NULL, ui, ci)
  parameters <- op$par

  a = parameters[1]
  b = parameters[3]/sig
  rho = parameters[2]/parameters[3]
#  H = parameters[4]
  # return raw parameters
  return(c(a,b,rho))

}


#' Estimated svi parameters with a given H
#'
#' Estimate modified svi parameters with a given H.
#' @param k log moneyness
#' @param iv iv
#' @param tau tau
#' @param H exponent H for the modified SVI function
#' @export
fit_quasi_svih <- function(k, iv, tau, H) {
  # two-step optimization:
  # 1. obtain inner parametplers: a,b, rho
  # 2. estimate m and sigma
  w = iv*sqrt(tau[1])
  ofitfun <- function(x) { # m, sigma, exponent
    y = (k-x[1])/x[2] # m, sigma
    par <- fit_inner(y, w, x[2],H) # a, b, rho, H
    param = c(par[1],par[2], x[1], par[3], x[2], H) #a, b, m, rho, sigma, H
    raw <- svi_fun(k, param)
    sum((raw - w)^2)
  }

  # intial guess
  x0 = c( 0.0, 0.02)

  ui = rbind(c(1,0),c(-1,0), c(0,1), c(0,-1))
  ci = c( min(k), -max(k), 0.0001, -5)

  op <- constrOptim(x0, ofitfun, NULL, ui, ci)
  opar = op$par
  y = (k-opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2],H)
  # parameters: a,b, m, rho, sigma, H
  parameters = data.frame(maturity =tau[1], a=par[1],b=par[2],
         m =opar[1], rho =par[3], sigma=opar[2],H= H)
}

#' Estimate svi parameers with GA algorithm
#'
#' This require ga package, GA may give a robust result
#' @export
#'
#' @param k log moneyness
#' @param iv iv
#' @param tau tau
fit_ga_svi <- function(k, iv, tau){
  # two-step optimization:
  # 1. obtain inner parameters: a,b, rho
  # 2. estimate m and sigma, H using GA
  w = iv*sqrt(tau[1])

  ofitfun <- function(x) { # m, sigma, exponent
    y = (k-x[1])/x[2] # m, sigma
    par <- fit_inner(y, w, x[2],x[3]) # a, b, rho, H
    param = c(par[1],par[2], x[1], par[3], x[2], x[3]) #a, b, m, rho, sigma, H
    raw <- svi_fun(k, param)
    value = 1/mean((raw - w)^2)
    return(value)
  }

  op <- GA::ga("real-valued", fitness = ofitfun,
           min = c(min(k),0.00001, 0.4),
           max = c(max(k),10, 1.1),
           maxiter = 100, run =100, optim = TRUE,
           monitor = TRUE)
  opar = op@solution[1,]

  y = (k-opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2], opar[3])
  # parameters: a,b, m, rho, sigma, H
  parameters = data.frame(maturity =tau[1], a=par[1],b=par[2],
                          m =opar[1], rho =par[3], sigma=opar[2],H= opar[3])

}

#' estimate svi parameter with know exponent H
#'
#' Using GA method to estimate parameters
#' @param k log moneyness
#' @param iv iv
#' @param tau tau
#' @export
fit_ga_svih <- function(k, iv, tau, H) {
  # two-step optimization:
  # 1. obtain inner parameters: a,b, rho
  # 2. estimate m and sigma
  w = iv*sqrt(tau[1])
  ofitfun <- function(x) { # m, sigma, exponent
    y = (k-x[1])/x[2] # m, sigma
    par <- fit_inner(y, w, x[2],H) # a, b, rho, H
    param = c(par[1],par[2], x[1], par[3], x[2], H) #a, b, m, rho, sigma, H
    raw <- svi_fun(k, param)
    1/mean((raw - w)^2)
  }

  op <- ga("real-valued", fitness = ofitfun,
           min = c(min(k),0.00001),
           max = c(max(k),10),
           maxiter = 100, run =100, optim = TRUE,
           monitor = TRUE)
  opar = op@solution[1,]

  y = (k-opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2],H)
  # parameters: a,b, m, rho, sigma, H
  parameters = data.frame(maturity =tau[1], a=par[1],b=par[2],
                          m =opar[1], rho =par[3], sigma=opar[2],H= H)
}


#' fit modifed SVI function
#'
#' The function will fit the 6 parameters, including exponent H, for the SVI function
#' Due to the exponent H, the fitting time is slow. It need to be improved.
#'
#' @param k log moneyness
#' @param iv iv
#' @param tau tau
#' @return svi 6 parameters
#' @export
fit_quasi_svi <- function(k, iv, tau) {
  # two-step optimization:
  # 1. obtain inner parameters: a,b, rho
  # 2. estimate m and sigma, H
  w = iv*sqrt(tau[1])
  ofitfun <- function(x) { # m, sigma, exponent
    y = (k-x[1])/x[2] # m, sigma
    par <- fit_inner(y, w, x[2],x[3]) # a, b, rho, H
    param = c(par[1],par[2], x[1], par[3], x[2], x[3]) #a, b, m, rho, sigma, H
    raw <- svi_fun(k, param)
    sum((raw - w)^2)
  }

  # intial guess
  x0 = c( 0.001, 0.02,0.7)

  ui = rbind(c(1,0,0),c(-1,0,0), c(0,1,0), c(0,-1,0),
             c(0,0,1),c(0,0,-1))
  ci = c( min(k), -max(k), 0.0001, -10, 0.4, -1.0)

  op <- constrOptim(x0, ofitfun, NULL, ui, ci)
  opar = op$par
  y = (k-opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2], opar[3])
  # parameters: a,b, m, rho, sigma, H
  parameters = data.frame(maturity =tau[1], a=par[1],b=par[2],
                          m =opar[1], rho =par[3], sigma=opar[2],H= opar[3])
}

#' Fit Chain SVI parameters with quasi-calibration approach
#'
#' Fit iv surface using svi approach with a given day chain. Due to slowness, the function
#' first calculate H for one-month chain, and use this H for the whole surface.
#'
#' @param chain an iv option chain containing iv, maturity, logmoneyness
#' Note: chain should be from one single day
#' @return surface parameter (a,b,m, rho, sigma, H)
#'
#' @export
svi.quasiparam <- function(chain){
  # fit surface, and return data.frame parameters
  # to reduce complexity and slowness of parameter calibration
  # first to fit one dte curve to find out the H parameter
  # then use the same H parameter for the rest

  taus = unique(chain$maturity)
  # find 1 month
  idx = which.min(abs(taus-0.1))
  slice = chain %>% filter(maturity == taus[idx])
  param = fit_quasi_svi(slice$logmoneyness, slice$iv, taus[idx])
  H = param[,'H']

  surface = data.frame()
  for (tau in taus){
    c = chain %>% filter(maturity == tau)
    param = fit_quasi_svih(c$logmoneyness,c$iv, tau,H)
    surface = rbind(surface, param)
  }
  ch <- chain %>% distinct(dte, maturity, mF, R, U, spot)

  left_join(surface, ch, by = c('maturity'))
}

#' interpolation of svi parameters
#'
#' @param para surface parameters
#' @param tau_interp to-be-interpolated tau
#' @return
#' interpoplated parameter at tau_interp
#' Note: Interpolation works only within the maturity range.
#' In this version, interpolation also for the linear part
#' @return values contains 5 of svi parameters, 3 of linear model
#' @export
quasi.interparam <- function(para, tau_interp) {

  #ensure scalar input
  if (!is.null(dim(tau_interp))) stop('tau_interp has to be scalar')
  if (nrow(para) == 1){
    inter_para = as.numeric(para[1,c('a','b','m','rho','sigma','H')])
  } else {
    # interpolation
    a = approxExtrap(para$maturity,para$a, xout = tau_interp)$y
    b = approxExtrap(para$maturity,para$b, xout = tau_interp)$y
    m = approxExtrap(para$maturity,para$m, xout = tau_interp)$y
    rho = approxExtrap(para$maturity,para$rho, xout = tau_interp)$y
    sigma = approxExtrap(para$maturity,para$sigma, xout = tau_interp)$y
    H = approxExtrap(para$maturity,para$H, xout = tau_interp)$y
    inter_para = c(a,b,m,rho,sigma,H)
  }
  return(inter_para)
}

#' interpolation of IVs with surface
#'
#' @param surf surface parameters
#' @param k logmoneyness
#' @param tau_interp to-be-interpolated tau
#' @return
#' interpoplated IVs at tau_interp
#' Note: Interpolation works only within the maturity range.
#' @export
quasi.interIV <- function(surf, k, tau_interp) {

  #ensure scalar input
  if (!is.null(dim(tau_interp))) stop('tau_interp has to be scalar')
  if (nrow(surf) == 1){
    inter_iv = svi_fun(k,surf)/sqrt(tau_interp)
  } else { # multiple DTEs
    # interpolation
    idx = which.min(abs(tau_interp-surf$maturity))
    # find two closes taus
    if (surf$maturity[idx] <= tau_interp)
      idx = idx+1
    if (idx == 1) # tau < min_tau
      idx = idx + 1
    if (idx > nrow(surf)) # tau > max_tau
      idx = nrow(surf)

    alpha = (surf$maturity[idx] - tau_interp)/(surf$maturity[idx]-surf$maturity[idx-1])

    # estimate iv curves with svi parameterization
    ivs1 = svi_fun(k, as.numeric(surf[idx-1,2:7]))/sqrt(surf[idx-1,1])
    ivs2 = svi_fun(k, as.numeric(surf[idx,2:7]))/sqrt(surf[idx,1])

    # interpolation
    inter_iv = ivs2 - alpha*(ivs2-ivs1)
  }
  return(inter_iv)
}


#' estimate iv curve and skew at a given DTE using modified SVI function
#'
#' This function will estimate a given iv curve with a specified dte
#' @param surface a surface param, estimated from SVI
#' @param dte Days to expiration
#' @param x a sequence of logmoneyness or strike. If max(x)<1, it will be regarded as log moneyness.
#' It return a table with iv, skew
#' @param n a smoothing factor for skew curve
#' @export
quasi.ivs <- function(surface, dte, spot = NULL, type = 'put', x = seq(-0.2,0.1,length.out = 600)){
  tau = dte/365
  if (is.null(spot))  spot = surface$spot[1]
  UR <- UR.interp(surface, tau, spot)
  if (max(x) <1 ) {
    # log moneyness
    lm = x
    strike = round(exp(lm)*UR$mF*10)/10
  } else {
    # strike input
    strike = x
    lm = log(strike/UR$mF)
  }
  # interpolate parameter - no longer linear (be cautious)
#  raw_para <- quasi.interparam(surface, tau)
  # raw: a,b,m, rho, sigma, H
  # estimate iv curve with svi parameterization
#  ivs = svi_fun(lm, raw_para)/sqrt(tau)
  ivs = quasi.interIV(surface, lm, tau)
  # skew of svi:
#  v1 = svi_fun(lm, c(raw_para[1:5],raw_para[6]-1))
#  b = raw_para[2]
#  m = raw_para[3]
#  rho = raw_para[4]
#  sigma = raw_para[5]
  skew <- c(NA, diff(ivs)/diff(lm))
#  skew <- raw_para[6] *v1 * b * (rho + (lm-m)/sqrt((lm-m)^2+sigma^2))
  # atm iv
  #  iv0 = sqrt(param_interp[9]/tau)
#  iv0 =svi_fun(0,raw_para)/sqrt(tau)

  # calculate Put price across K space
  opts = NULL
  for (idx in c(1:length(lm))){
    opt <- EuropeanOption(type, UR$U, strike[idx], 0, UR$R, tau, ivs[idx])
    opts = rbind(opts, do.call(cbind, opt))
  }
  iv_chain <- data.frame(lm = lm, strike = strike, iv = ivs, skew = skew,
                         spot = spot)

  iv_chain <- cbind(iv_chain,opts)
  # calculate delta true, or risk neutral probability (sign diff)
  if (type == 'put')
    iv_chain$delta_true <- -(iv_chain$delta - iv_chain$vega*iv_chain$skew/spot)
  else
    iv_chain$delta_true <- iv_chain$delta - iv_chain$vega*iv_chain$skew/spot

  iv_chain$delta_true <-  exp(UR$tau*UR$R) * iv_chain$delta_true
  iv_chain$theta <- iv_chain$theta/365 #convert to 1-day theta
  if (length(lm)>100){ # smooth the pdf
    iv_chain$pdf <-  c(NA, diff(iv_chain$delta_true)) # no meaning if x is not continuous
    n = round(length(lm)/30)
    iv_chain$pdf <- as.numeric(stats::filter(iv_chain$pdf, rep(1/n,n), sides = 2))
  }
  iv_chain$dte <- dte

  return(iv_chain)
}

#' plot iv surface with fitted model
#'
#' @param ivchain iv chain
#' @param para svi parameters
#' @param range dte range (default c(7,90))
#' @param type default ('iv') plot normal IV figure, or 'ivt', plot iv/sqrt(t), or 'var', plot total
#' @param x default ('strike'), alternative 'logmoneyness'
#' implied variance (i.e., iv^2*tau)
#' @export
quasi.plotsurface <- function(ivchain, para, range =c(7,90), type = 'iv', x = 'strike'){
  iv <- ivchain %>% dplyr::filter(dte >= range[1] & dte <=range[2])
  iv$svi.iv = 0
  for (id in unique(iv$dte)){
    svi.iv = quasi.ivs(para, id, x = iv$logmoneyness[iv$dte==id])
    iv$svi.iv[iv$dte==id] = svi.iv$iv
  }
  iv$err = iv$iv - iv$svi.iv

  # plot fitted data
  if (x == 'strike')
    spot = iv$spot[1]
  else
    spot = 0

  iv$expiration = factor(iv$expiration)
  iv$dte  = factor(round(iv$dte))
  switch(type,
         'iv' = {
           figs <- ggplot(iv, aes_string(x, "iv", color = "dte")) +
             geom_line(aes( y = svi.iv)) + geom_point() +
             geom_vline(xintercept = spot, color = 'gray') +
             xlab(x) + ylab('IV')
         },
         'ivt' = { # plot variance
           figs <- ggplot(iv, aes_string(x, y = "iv*sqrt(maturity)", color = "dte")) +
             geom_line(aes( y = svi.iv*sqrt(maturity))) + geom_point() +
             geom_vline(xintercept = spot, color = 'gray') +
             xlab(x) + ylab('IV * sqrt(T)')
         },
         'var' = { # plot variance
           figs <- ggplot(iv, aes_string(x, "iv*iv*maturity", color = "dte")) +
             geom_line(aes( y = svi.iv*svi.iv*maturity)) + geom_point() +
             geom_vline(xintercept = spot, color = 'gray') +
             xlab(x) + ylab('Implied Variance * T')
         })
  print(figs)
  print(c(sd(iv$err),max(iv$err)))
}

#' estimate skews across DTES
#'
#' estimate skews from the modified SVI surface
#' @export
#' @param surface  svi surface parameters
#' @param dtes dtes
#' @param x log moneyness or strikes
#'
quasi.skew <- function(surface, dtes = c(10,20,30,40,50,60,70), x = seq(-0.1,0.1, length.out = 300), plot = TRUE){
  ivs = data.frame()
  for (dte in dtes){
    iv <- quasi.ivs(surface, dte, x = x)
    ivs = rbind(ivs, iv)
  }
  ivs$dte <- factor(ivs$dte)
  if (plot)
    print(ggplot(ivs, aes(strike, skew, color = dte)) + geom_line())
  return(ivs)
}

#' Plot term structure
#'
#' Based on modified svi surface, estimate ATM IVs across DTEs
#' @export
#' @param surface surface parameters
#' @param plot plot the graph or not
quasi.termstructure <- function(surface, plot = TRUE){
  dtes = seq(round(min(surface$dte)), round(max(surface$dte)))
  ivs = data.frame()
  for (dte in dtes){
    iv <- quasi.ivs(surface, dte, x = 0)
    ivs = rbind(ivs, iv)
  }
  if (plot)
    print(ggplot(ivs, aes(dte, iv0)) + geom_line() + ylab('ATM IV'))
  return(ivs[,c('dte','iv0')])
}

#' estimate surface from a chain file
#'
#' Estimate surface with modified SVI model using a chain file
#' The chain file name should be 'rds'. This function uses parallel
#' computing.
#'
#' Return the whole chain.
quasi.estimateSurface <- function(filename, range = c(7,90)){
  r_date = gsub('.*_(.*).rds', '\\1', filename)
  chain = readRDS(filename)
  chain <- chain %>% dplyr::filter(dte >= range[1] & dte <= range[2])
  lchain <- split(chain, chain$date)
  # use parallel computing
  no_cores <- detectCores()-1
  # Initiate cluster
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {library(dplyr)})
  surfs <- parLapply(cl, lchain, svi.quasiparam)
  stopCluster(cl)
  surface = do.call(rbind, surfs);
  return(surface)
}

#' estimate surface error
#'
#' Provide a chain and surface parameters, show the error
#' By default, the graph output is off.
#' @export
#'
quasi.err <- function (ivchain, para, range = c(7, 90), show = FALSE)
{
  iv <- ivchain %>% dplyr::filter(dte >= range[1] & dte <=
                                    range[2])
  iv$svi.iv = 0
  for (id in unique(iv$dte)) {
    svi.iv = quasi.ivs(para, id, x = iv$logmoneyness[iv$dte ==
                                                       id])
    iv$svi.iv[iv$dte == id] = svi.iv$iv
  }
  iv$err = iv$iv - iv$svi.iv
  iv$dte = factor(round(iv$dte))
  if (show)
    print(ggplot(iv, aes(strike, err/iv * 100, color = dte))
          + geom_line() + ylab('Error (%)'))
  err = data.frame(date = iv$date[1], sd = sd(iv$err), max = max(iv$err))
  return(err)
}

#' Fit Chain SVI parameters with quasi-calibration approach
#'
#' Fit iv surface using svi approach with a given day chain. Due to slowness, the function
#' first calculate H for one-month chain, and use this H for the whole surface.
#' This version uses parallel computing.
#'
#' @param chain an iv option chain containing iv, maturity, logmoneyness
#' Note: chain should be from one single day
#' @return surface parameter (a,b,m, rho, sigma, H)
#'
#' @export
svi.quasiparam.p <- function(chain, H = NULL){
  # fit surface, and return data.frame parameters
  # to reduce complexity and slowness of parameter calibration
  # first to fit one dte curve to find out the H parameter
  # then use the same H parameter for the rest

  # use parallel computing
  no_cores <- detectCores() -1
  # Initiate cluster
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {library(dplyr)})

  # first estimate H
  if (is.null(H)){
    taus = unique(chain$maturity)
    # find 1 month, and get H parameter for the whole surface
    idx = which.min(abs(taus-0.1))
    slice = chain %>% filter(maturity == taus[idx])
    param = fit_quasi_svi(slice$logmoneyness, slice$iv, taus[idx])
    H = param[,'H']
  }
  fit_curve <- function(c, h){
    param = fit_quasi_svih(c$logmoneyness,c$iv, c$maturity[1],h)
    return(param)
  }

  lchain <- split(chain, chain$maturity)
  paras <- parLapply(cl, lchain,fit_curve, h = H)
  surface = do.call(rbind, paras)
  ch <- chain %>% distinct(dte, maturity, mF, R, U, spot)
  surf <- left_join(surface, ch, by = c('maturity'))

  stopCluster(cl)
  return(surf)
}

#' plot surface with a given hour and minutes in a day
#'
#' This function to plot furface with a given hour and minutes, such as "10:30".
#' It will then extract the chain and surface and plot it.
#' @export
#' @param chain chains
#' @param hourtime time, e.g., "10:30"
#'
quasi.surface.at <- function(chain,surf, hourtime, type = 'ivt', range = c(7,100)){
    c <- chain %>% filter(substr(date,12,16) == hourtime)
    s <- surf %>% filter(substr(date,12,16) == hourtime)
    print(quasi.plotsurface(c,s, type = type, range = range))
    return(list(chain = c, surface = s))
}
