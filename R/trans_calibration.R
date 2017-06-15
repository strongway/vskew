# using translation method to parameterize the surface
# idea: implied vol curve change over time in a regular way.
# This method is to find their transform and translation across DTEs, and this would give insights of the surface


#' svi original function
#'
#' @param k log moneyness
#' @param param 5 parameters: a, b,m, rho, sigma
#' @return estimated values according to the curve
#' @export
trans.svi_fun <- function(k,param){
  a = param[1];
  b = param[2];
  m = param[3];
  rho = param[4];
  sigma = param[5];
  result = a + b*(rho*(k-m) + sqrt((k-m)^2 + sigma^2))
}

#' Inner optimzation of the SVI function
#'
#' Here given m and sigma (in x), k, and total iv w, estimate a, b, rho
#' @param x two parameters (m, sigma)
#' @param k log moneyness or alike
#' @param w total implied vol, in the format of iv^(1/H). When H = 1/2, it is normal total variance
#' @export
#' @return estimated 5 SVI parameters (a, b, m, rho, sigma)
trans.fit_inner <- function(x, k,w){
  y = (k-x[1])/x[2] # m, sigma
  z = sqrt(y^2+1)
  lr <- lm(w ~ y + z)
  par <- lr$coefficient
  rho = par[2]/par[3]
  # constrains on rho
  if (abs(rho) > 1){ # -1< rho <1
    rho <- -1 # in most case, left skew
    y2 <- z - y # rho = -1
    lr <- lm(w ~ y2)
    par <- lr$coefficient
    a <- par[1]
    b <- par[2]/x[2] # b*sigma
  } else {
    a <-  par[1]
    b <-  par[3]/x[2] # c' = b*sigma (x[2])
  }
  # return svi param
  param = c(a,b, x[1], rho, x[2]) #a, b, m, rho, sigma
  return(param)
}

#' Inner optimzation of the SVI function with known rho
#'
#' Here given m and sigma and rho (in x), k, and total iv w, estimate a, b
#' @param x two parameters (m, rho,sigma)
#' @param k log moneyness or alike
#' @param w total implied vol, in the format of iv^(1/H). When H = 1/2, it is normal total variance
#' @export
#' @return estimated 5 SVI parameters (a, b, m, rho, sigma)
trans.fit_ab <- function(x, k,w){
  y = (k-x[1])/x[3] # m, sigma
  z = x[3]*(y*x[2] + sqrt(y^2+1)) # sigma*[rho*y + sqrt(y^2+1)]
  lr <- lm(w ~ z)
  par <- lr$coefficient
  a <- par[1]
  b <- par[2]
  # return svi param
  param = c(par[1],par[2], x[1], x[2], x[3]) #a, b, m, rho, sigma
  return(param)
}

#' fit SVI function
#' Here we assume iv^H = SVI_Fun(Param)
#' The function will fit the 6 parameters, including exponent H, for the SVI function
#'
#' @param k log moneyness
#' @param iv iv
#' @return svi 6 parameters
#' @export
trans.fit_svi <- function(k, iv) {
  # two-step optimization:
  # 1. obtain inner parameters: a,b, rho
  # 2. estimate m and sigma, H
  ofitfun <- function(x, k, iv) { # m, sigma, exponent
    w = iv^(1/x[3]) # total iv: iv^(1/H)
    param = trans.fit_inner(x, k, w) # estimate inner
    raw <- trans.svi_fun(k, param) # estimate raw
    sum((raw^x[3] - iv)^2)
  }

  # intial guess
  x0 = c( mean(k), 0.02,0.7)

  # constrained boxes
  ui = rbind(c(1,0,0),c(-1,0,0), c(0,1,0), c(0,-1,0),
             c(0,0,1),c(0,0,-1))
  ci = c( min(k), -max(k), 0.0001, -10, 0.4, -1.0)

  op <- constrOptim(x0, ofitfun, NULL, ui, ci, k=k, iv= iv)
  opar = op$par # outer optimization parameters: m, sigma, H
  # get parameters of SVI function
  p = trans.fit_inner(opar,k, iv^(1/opar[3]))
  # parameters: a,b, m, rho, sigma, H
  return(data.frame(a = p[1],b = p[2],m = p[3], rho = p[4], sigma = p[5], H = opar[3])) #a, b, m, rho, sigma, H
}

#' fit SVI function with a known H and a known rho
#' Here we assume iv^H = SVI_Fun(Param)
#' The function will fit the 5 parameters, including exponent H, for the SVI function
#'
#' @param k log moneyness
#' @param iv iv
#' @param H H exponent
#' @param rho rho parameter in svi function
#' @param weights weights for individual rows (weighted optimization)
#' @return svi 6 parameters
#' @export
trans.fit_svih <- function(k, iv, H, rho, weights = 1) {
  # known H and rho
  w = iv^(1/H) # iv^(1/H_t)
  # two-step optimization:
  # 1. obtain inner parameters: a,b
  # 2. estimate m and sigma
  ofitfunh <- function(x, k, w, H, rho) { # x -> m, sigma
    xi = c(x[1],rho,x[2])
    param = trans.fit_ab(xi, k, w)
    raw <- trans.svi_fun(k, param)
    # avoid negative iv
    raw[raw<=0] =0
    sum((raw^H - iv)^2*weights)
  }

  # intial guess of m, sigma
  x0 = c( max(0,max(k)/2), 0.02)

  ui = rbind(c(1,0),c(-1,0), c(0,1), c(0,-1))
  ci = c( -0.1, -max(0.1,max(k)), 0.0001, -10)

  op <- constrOptim(x0, ofitfunh, NULL, ui, ci, k=k, w= w, H = H, rho = rho)
  opar = op$par
  p = trans.fit_ab(c(opar[1],rho,opar[2]),k, w)
  # parameters: a,b, m, rho, sigma, H
  return(data.frame(a = p[1],b = p[2],m = p[3], rho = p[4], sigma = p[5], H = H)) #a, b, m, rho, sigma, H
}

#' Fit Vol Surface with the transform method
#'
#' Estimate surface by transforming and translating the surface to a uniform surface
#' Using this to find rho and H. The whole surface uses the same rho and H.
#' Here we assume the curvation of the iv against k/sqrt(tau) is fixed (H), and
#' the smile rotation is fixed as well
#' @param ochain input chain
#' @return estimated parameters
#' @export
#'
vs.fitSurface <- function(ochain, parallel = TRUE){
  # step 1: roughly estimate weights based on log moneyness
#  ochain$weight = dnorm(ochain$logmoneyness, mean = mean(ochain$logmoneyness), sd = sd(ochain$logmoneyness)*1.5)
  # centered on the ATM
  ochain$weight = dnorm(ochain$logmoneyness, mean = 0, sd = sd(ochain$logmoneyness)*1.5)

  # step 1.1: Estimate H_t Hurst component
  # it is a scaling factor for b, m, and sigma, so it can be skipped here
  #op <- optim(0.5,trans.fit_H, NULL, method = 'Brent',
  #            lower = 0.25, upper = 1, chain = ochain)
  #H <- op$par

  # x = K/T^H
  #ochain$kth = ochain$logmoneyness/ochain$maturity^0.5
  #ochain$kt2 = ochain$logmoneyness/sqrt(ochain$maturity)
  # use standard one
  ochain$kth = ochain$logmoneyness/sqrt(ochain$maturity)

  # Step 2: Estimate IV0 using linear regression, and
  #    translation scalex, and scaley
  ochain %>% filter(abs(kth) < 0.1) %>%
    group_by(date, maturity) %>%
    do(tidy(lm(iv ~  poly(logmoneyness,2, raw = TRUE), data = .))) %>%
    select(date, maturity, term, estimate) %>%
    spread(term, estimate) %>% rename(iv0 = `(Intercept)` ) %>%
    select(date,maturity, iv0) -> lr

  ochain <- left_join(ochain, lr, by = c('date','maturity'))
  # find best translation for chains within dte (30,70)
  # this is the best chains to estimate rho
  midchain <- ochain %>% filter(kth > -0.4 & dte > 30 & dte < 70)
  op <- optim(c(1.2, 0.7),trans.fit_curve, NULL,  method = 'L-BFGS-B',
              lower = c(0.01,0.01), upper = c(10,10),
              chain = midchain)
  scalex <- op$par[1]
  scaley <- op$par[2]
  # ggplot(ochain, aes(kth - iv0*scalex, iv - iv0*scaley, color = dte)) + geom_point()

  # step 4: estimate a general curve, then estimate H and rho

  midchain <- midchain %>% mutate(iv2 = iv - iv0*scaley,
                              kth2 = kth - iv0*scalex)
  # avoid negative iv
  midchain$iv2 = midchain$iv2 - min(midchain$iv2) + 0.1
  #  ggplot(ochain, aes(kth2, iv2, color = as.factor(dte))) + geom_point() + geom_line()
  #  lr2 <- lm(iv2 ~ poly(kth2,3, raw = TRUE), data = ungroup(ochain), weights = weight)
  #  p <- lr2$coefficients #a,b,c,d
  p <- trans.fit_svi(midchain$kth2, midchain$iv2)
  # --- p$H and p$rho ----

  # step 5: calibrate each DTE
  # estimate each with kth and ivs
  if (parallel) {
    # use parallel computing
    no_cores <- detectCores() -1
    # Initiate cluster
    cl <- makeCluster(no_cores)
    clusterEvalQ(cl, {library(dplyr)})
    clusterEvalQ(cl, {library(vskew)})

    fit_curve <- function(c, par){
      param = trans.fit_svih(c$kth,c$iv, par$H[1], par$rho[1], c$weight )
      param$date = c$date[1]
      param$maturity = c$maturity[1]
      return(param)
    }

    lchain <- split(ochain, paste(ochain$date,ochain$maturity))
    paras <- parLapply(cl, lchain, fit_curve, par = p)
    lr3 <-  do.call(rbind, paras)

  #  stopCluster(cl)
  }
  else {

  # now we use kt2 (k/sqrt(tau)) to estimate parameters
  lr3 <- ochain %>% group_by(date, maturity) %>%
      do(trans.fit_svih(.$kth,.$iv, p$H[1], p$rho[1],.$weight))
  }
  spara <- left_join(lr,lr3, by=c('date','maturity'))

  ru <- ochain %>% ungroup(.) %>% distinct(date, maturity, spot, R, U, mF)
  spara <- left_join(spara, ru, by = c('date','maturity'))
  return(spara)
}


# Optimization function for fitting Hust component H_t
#
# @param H Hurst component
# @param chain chain
# @param range range of k for linear regression
# @export
#trans.fit_H <- function(H, chain, range = c(-0.2,0)){
  # estimate slopes for each DTEs,
  # see if the slopes are similar
#  lr <-   chain %>% mutate(y = iv*maturity^H) %>%
#    filter(logmoneyness > range[1] & logmoneyness < range[2]) %>%
#    group_by(date, dte) %>% do(tidy(lm(y ~ logmoneyness, data = .))) %>%
#    select(date, dte, term, estimate) %>% spread(term, estimate) %>%
#    lm(logmoneyness ~ 1, data = .) # constant?
#  return(sum(lr$residuals^2))
#}

#' find optimal translation x, y
#'
trans.fit_curve <- function(x, chain){
  scalex = x[1]
  scaley = x[2]
  ungroup(chain) %>% mutate(iv2 = iv - iv0*scaley,
                            kth2 = kth - iv0*scalex) %>%
    lm(iv2 ~ poly(kth2,3, raw = TRUE), data = ., weights = weight) -> lr2
  return(sum(lr2$residuals^2*lr2$weight))
}

#' interpolation of parameters
#'
#' @param surf surface parameters
#' @param tau_interp to-be-interpolated tau
#' @return interpolated parameters
#' Note: Interpolation works only within the maturity range.
#' To interpolate the parameters, one has to convert raw parameters to
#' jump wing parameters, given that jump wing parameters have clear meanings.
#' That is, iv0, skew, left/right slopes, minimum. Those parameters can be
#' easily interpolated.
#' @export
vs.interParameter <- function(surf, tau_interp) {
  #ensure scalar input
  if (!is.null(dim(tau_interp))) stop('tau_interp has to be scalar')
  if (nrow(surf) == 1){
    surf$maturity = tau_interp
    inter_para <- surf[,-1]
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

    # convert parameters to jump wing
    para1 = surf[idx-1, c('a','b','m','rho','sigma')]
    para2 = surf[idx, c('a','b','m','rho','sigma')]
    svi_j1 = svi_convertparameters(as.numeric(para1),'raw','jumpwing')
    svi_j2 = svi_convertparameters(as.numeric(para2),'raw','jumpwing')
    inter_jump = svi_j2 - alpha*(svi_j2-svi_j1)
    inter_raw = svi_convertparameters(inter_jump, 'jumpwing', 'raw')
    # interpolate the rest parameters
    inter_para <- surf[idx,]
    inter_para$maturity <- tau_interp
    inter_para$iv0 <- surf[idx,'iv0'] -alpha * (surf[idx,'iv0']-surf[idx-1,'iv0'])
    # replace the 5 svi parameters.
    inter_para$a = inter_raw[1]
    inter_para$b = inter_raw[2]
    inter_para$m = inter_raw[3]
    inter_para$rho = inter_raw[4]
    inter_para$sigma = inter_raw[5]
  }
  return(inter_para)
}

#' Estimate iv from trans.method parameters
#' @param k log moneyness
#' @param para parameters include maturity, iv0, a,b,m, rho, sigma H, Ht, Sx,Sy
vs.ivs <- function(k,para){
  if (!(nrow(para) == 1)) stop('Please provide one single row parameters')
#  para <- para %>% mutate(d_iv = iv0*Sy, d_k = iv0*Sx, tau_h = maturity^Ht)
  kth = k/sqrt(para$maturity) #- para$d_k
  param = c(para$a,para$b,para$m, para$rho,para$sigma)
  w = trans.svi_fun(kth,param)
  iv =  w^para$H
}

#' plot iv surface with fitted model
#'
#' @param ivchain iv chain
#' @param para svi parameters
#' @param range dte range (default c(7,90))
#' @param type default ('iv') plot normal IV figure, or 'ivt', plot iv/sqrt(t), or 'var', plot total
#' @param x default ('strike'), alternative 'logmoneyness'
#' @export
vs.plotSurface <- function(ivchain, para, range =c(10,90), type = 'iv', x = 'strike'){
  iv <- ivchain %>% dplyr::filter(dte > range[1] & dte <=range[2])
  iv$svi.iv = 0
  for (id in unique(iv$maturity)){
    par <- para %>% filter(maturity == id)
    svi.iv = vs.ivs(iv$logmoneyness[iv$maturity==id],par)
    iv$svi.iv[iv$maturity==id] = svi.iv
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
         })
  print(figs)
  print(c(sd(iv$err),max(iv$err)))
}

#' calculate skew slope curve based on numeric approach
#'
#'
vs.slop0 <- function(k, para){
  skews = data.frame()
  for (id in unique(para$maturity)){
    par <- para %>% filter(maturity == id)
    iv = vs.ivs(k,par)
    skew = diff(iv)/diff(k)
    df.skew = data.frame(dte = round(id*365.25), k = k,
                         iv = iv, skew = c(NA, skew))
    skews = rbind(skews,df.skew)
  }
  skews$dte = factor(skews$dte)
  return(skews)
}

#' Skew Slope
#'
#' Esimate Skew slope based derative on equation
#'
vs.slope <- function(k, s){
  skews = data.frame()
  for (id in unique(s$maturity)){
    par <- s %>% filter(maturity == id)
    x = k/sqrt(par$maturity)
    param = c(par$a, par$b, par$m, par$rho, par$sigma)
    svi <- trans.svi_fun(x, param)
    dsvi = par$b * (par$rho + (x-par$m)/sqrt((x-par$m)^2 + par$sigma^2))
    # using chain rule
    skew = par$H*svi^(par$H - 1)* dsvi /sqrt(id)
    iv = vs.ivs(k, par)
    df.skew = data.frame(dte = round(id*365.25), k = k,
                         iv = iv, skew = skew)
    skews = rbind(skews,df.skew)
  }
  skews$dte = factor(skews$dte)
  return(skews)
}


#' estimate price, skew and nrd at a given DTE
#'
#' This function will estimate options prices etc at a given iv
#' curve with a specified dte
#' @param surface a surface param, estimated from SVI
#' @param dte Days to expiration
#' @param x a sequence of logmoneyness or strike. If max(x)<1, it will be regarded as log moneyness.
#' @export
vs.options <- function(surface, dte, spot = NULL, type = 'put', x = seq(-0.2,0.1,length.out = 600)){
  tau = dte/365.25
  para = vs.interParameter(surface, tau)
  if (is.null(spot))  spot = para$spot[1]
  dp = spot - para$spot[1] # price change
  if (max(x) <1 ) {
    # log moneyness
    k = x
    strike = round(exp(k)*para$mF*10)/10
  } else {
    # strike input
    strike = x
    mF = para$mF + dp # adjust forward price
    k = log(strike/mF)
  }
  df <- vs.slope(k,para)

  # calculate Put price across K space
  opts = NULL
  for (idx in c(1:length(k))){
    U = para$U + dp
    opt <- EuropeanOption(type, U, strike[idx], 0, para$R, tau, df$iv[idx])
    opts = rbind(opts, do.call(cbind, opt))
  }
  iv_chain <- data.frame(k = k, strike = strike, iv = df$iv, skew = df$skew,
                         spot = spot)

  iv_chain <- cbind(iv_chain,opts)
  # calculate strike-delta and risk neutral probability
  # strike-delta = delta + vega * skew
  iv_chain$strike_delta <- iv_chain$delta + iv_chain$vega*iv_chain$skew/iv_chain$strike
  iv_chain$strike_delta <-  exp(para$maturity*para$R) * iv_chain$strike_delta
  iv_chain$theta <- iv_chain$theta/365 #convert to 1-day theta
  if (length(k)>100){ # smooth the pdf
    iv_chain$rnd <-  c(NA, abs(diff(iv_chain$strike_delta))) # no meaning if x is not continuous
    n = round(length(k)/30)
    iv_chain$rnd <- as.numeric(stats::filter(iv_chain$rnd, rep(1/n,n), sides = 2))
  }
  iv_chain$dte <- dte
  return(iv_chain)
}

#' Plot term structure
#'
#' Based on modified svi surface, estimate ATM IVs across DTEs
#' @export
#' @param surface surface parameters
vs.termstructure <- function(surface){
  print(ggplot(surface, aes(maturity, iv0)) +
          geom_line() + ylab('IV0') + geom_point() + xlab('tau') +
    geom_smooth(method = 'lm', formula = y ~ poly(x,2)) )
}

#' Esimate Vertical Skew
#'
#' Esimate vertical skew curve with given x, surface, and dte
#' @param x Strikes or log(K/F)
#' @param surf surface parameters
#' @param dte Day-to-expiration
#' @param spot spot value for sticky-delta calculation
#'
#' When the spot is not specify, vertical estimation is soly dependent
#' on the surf, so it is sticky-strike. When the spot is specify,
#' the surface will move and allocate to the spot (i.e., sticky-delta)
#' @export
vs.vskew <- function(x, surf, dte, spot=NULL){
  if (is.null(spot)) spot = surf$spot[1]
  # estimate vertical skew parameter at DTE
  para0 <- vs.interParameter(surf, dte/365.25)
  # if spot is different
  mF = para0$mF + (spot - para0$spot[1]) # rough approx.

  if (max(x)<1){ # log strike
    k = x
    strikes = mF * exp(k)
  } else {
    strikes = x
    k = log(strikes/mF)
  }
  # calculate ivs
  ivs <-   vs.ivs(k, para0)
  data.frame(strike = strikes, k = k, iv = ivs)
}
