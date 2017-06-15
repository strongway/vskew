#' estimate iv curve and skew at a given DTE
#'
#' This function will estimate a given iv curve with a specified dte
#' @param param a surface param, estimated from SVI
#' @param dte Days to expiration
#' @param x a sequence of logmoneyness or strike. If max(x)<1, it will be regarded as log moneyness.
#' It return a table with iv, skew
#' @param n a smoothing factor for skew curve
#' @export
stiching.curve <- function(param, dte, x = seq(-0.2,0.1,length.out = 600)){
  tau = dte/365.25
  UR <- UR.interp(param, tau, param$spot[1])
  if (max(x) <1 ) {
    # log moneyness
    lm = x
    strike = round(exp(lm)*UR$mF*10)/10
  } else {
    # strike input
    strike = x
    lm = log(strike/UR$mF)
  }
  # interpolate parameter
  param_interp <- svi.interparam(param, tau)
  # raw: a,b,m, rho, sigma
  raw_para <- svi_convertparameters(param_interp[1:5], 'jumpwing','raw', tau)
  # estimate iv curve with svi parameterization
  iv_interp <- svi_jumpwing(lm, param_interp[1:5], tau)
  iv <- iv_interp$impliedvolatility
  # skew of svi: 0.5*b*(rho+(k-m)/sqrt((k-m)^2)+sigma^2)/sqrt(v*tau)
  skew <- 0.5*raw_para[2] *(raw_para[4]+
                              (lm-raw_para[3])/(sqrt((lm-raw_para[3])^2+raw_para[5]^2)))/
    sqrt(param_interp[1]*tau)
  # adjust left side with hybrid model
  if (length(param_interp) > 5) { # left side adjustment
    # adjust the left linear wing using linear model
    rep_idx <-  lm <= param$boundary[1]
    iv[rep_idx] <- (param_interp[6] + param_interp[7]*lm[rep_idx]+
                      param_interp[8]*lm[rep_idx]^2 )/sqrt(tau)
    skew[rep_idx] <- param_interp[7] + 2*param_interp[8]*lm[rep_idx]
    # transition between two region
    tran_idx <-  lm > param$boundary[1] & lm <= 0
    iv_l <- (param_interp[6] + param_interp[7]*lm[tran_idx]+
               param_interp[8]*lm[tran_idx]^2)/sqrt(tau)
    iv_r <- iv[tran_idx]
    r_weight = seq(1,length(iv_r))/length(iv_r)
    #weighted transition of iv
    iv[tran_idx] <- iv_r * r_weight + iv_l*(1-r_weight)

    skew_l <- param_interp[7] + 2*param_interp[8]*lm[tran_idx]
    skew_r <- skew[tran_idx]
    # weighted transition of skew
    skew[tran_idx] <- skew_l * (1-r_weight) + skew_r*r_weight
  }
  # atm iv
  iv0 = sqrt(param_interp[9]/tau)
  #  skew = as.numeric(stats::filter(skew, rep(1/n,n), sides = 2))
  iv.df <- data.frame(lm = lm, strike = strike, iv = iv, skew = skew, iv0 = iv0)
  return(iv.df)
}

#' option price based on svi-estimated iv surface
#'
#' Inputs: type: 'call' or 'put'
#'         strike, s(pot), tau (maturity),
#'         param: iv surface parameter data.frame
#' Return: EuropeanOption (value, delta, gamma, vega, theta, rho, divRho)
#'§
#'
#'  optionPrice is now accept array of strike, spot, and type. Three must be the same length
#' @export
stiching.optionPrice <- function(strike, param, tau, s = NULL, type = "put") {
  if (is.null(s)) {# use the spot from the surface
    s = rep(param$spot[1], length(strike))
  }
  if (length(strike) > length(type)){
    type = rep(type, length.out = length(strike))
  }
  UR <- UR.interp(param, tau, s )
  # jump-wing parameter at given tau
  param_interp <- svi.interparam(param, tau) # svi parameters at tau
  # we assume the volatility surface remains the same, and estimate the iv from the surface
  # in reality iv changes. It must consider the price movement alter the whole surface in a certain way.
  # this step requires modeling.
  lm = log(strike/UR$mF)
  iv_interp <- stiching.curve(param, tau*365.25, lm)
  iv <- iv_interp$iv
  # calculate EuropeanOption parameters
  opts = NULL
  for (idx in c(1:nrow(UR))){
    opt <- EuropeanOption(as.character(type[idx]), UR$U[idx], strike[idx], 0, UR$R[idx], tau, iv[idx])
    opts = rbind(opts, do.call(cbind, opt))
  }
  opts
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
plot.stiching_surface <- function(ivchain, para, range =c(7,90), type = 'iv', x = 'strike'){
  iv <- ivchain %>% dplyr::filter(dte >= range[1] & dte <=range[2])
  iv$svi.iv = 0
  for (id in unique(iv$dte)){
    svi.iv = stiching.curve(para, id, iv$logmoneyness[iv$dte==id])
    iv$svi.iv[iv$dte==id] = svi.iv$iv
  }

  # plot fitted data
  if (x == 'strike')
    spot = iv$spot[1]
  else
    spot = 0

  iv$expiration = factor(iv$expiration)
  iv$dte  = factor(round(iv$dte))
  switch(type,
         'iv' = {
           ggplot(iv, aes_string(x, "iv", color = "dte")) +
             geom_line(aes( y = svi.iv)) + geom_point() +
             geom_vline(xintercept = spot, color = 'gray') +
             xlab(x) + ylab('IV')
         },
         'ivt' = { # plot variance
           ggplot(iv, aes_string(x, y = "iv*sqrt(maturity)", color = "dte")) +
             geom_line(aes( y = svi.iv*sqrt(maturity))) + geom_point() +
             geom_vline(xintercept = spot, color = 'gray') +
             xlab(x) + ylab('IV * sqrt(T)')
         },
         'var' = { # plot variance
           ggplot(iv, aes_string(x, "iv*iv*maturity", color = "dte")) +
             geom_line(aes( y = svi.iv*svi.iv*maturity)) + geom_point() +
             geom_vline(xintercept = spot, color = 'gray') +
             xlab(x) + ylab('Implied Variance * T')
         })
}

#' Risk neutral probability based on the iv
#'
#' This function estimated a risk neutral probability /density function at a given DTE
#' @param param svi surfance paramters
#' @param dte Days to expiration
#' @param sopt default is the spot price of the surface that obtained. It can be different.
#' @param type put or call
#' @param x x-axis (log moneyness or strike)
#' @param n smoothing factor for skew and pdf
#'
#' @export
stiching.rnd <- function(param, dte, spot = NULL, type = 'put',x = seq(-0.2,0.1,length.out = 600)){
  #estimated iv curve
  tau = dte/365.25
  if (is.null(spot))  spot = param$spot[1]
  iv1 <- stiching.curve(param, dte, x)
  UR <- UR.interp(param, tau, spot)
  # calculate Put price across K space
  opts = NULL
  for (idx in c(1:nrow(iv1))){
    opt <- EuropeanOption(type, UR$U, iv1$strike[idx], 0, UR$R, tau, iv1$iv[idx])
    opts = rbind(opts, do.call(cbind, opt))
  }
  iv_chain <- cbind(iv1,opts)
  # calculate delta true, or risk neutral probability (sign diff)
  if (type == 'put')
    iv_chain$delta_true <- -(iv_chain$delta - iv_chain$vega*iv_chain$skew/spot)
  else
    iv_chain$delta_true <- iv_chain$delta - iv_chain$vega*iv_chain$skew/spot

  pdf <- c(0, diff(iv_chain$delta_true))

  iv_chain$pdf <- as.numeric(stats::filter(pdf, rep(1/n,n), sides = 2))

  iv_chain$delta_true <-  exp(UR$tau*UR$R) * iv_chain$delta_true
  iv_chain$pdf <- exp(UR$tau*UR$R) * iv_chain$pdf
  iv_chain$dte <- dte
  return(iv_chain)
}


#' obtain optio price history from the chain
#'
#' If the chain or strike is not available, based on svi interpolation
#' Input: start date, expiration date, strike combo, whole option chain
#' @export
stiching.option.History <- function(start_date, exp_date, strike_combo, wholechain) {
  # chain contains options that we are interested in.
  max_dte = as.numeric(as.POSIXct(exp_date) - as.POSIXct(start_date) , units = 'days')
  chain <- wholechain %>% filter( date >= as.POSIXct(start_date) &  date <= as.POSIXct(expiration))
  # dates we have in the chain
  dates <-  unique(chain$date)
  opt <- data.frame(date = rep(dates, each = length(strike_combo)), strike = strike_combo, expiration = as.POSIXct(exp_date))
  opt <- left_join(opt, chain, by = c('date','strike', 'expiration')) %>% select( date, expiration, strike, C, P, dte, iv)

  # check if there are missing values
  na.indx <- which(is.na(opt$iv))
  for (indx in na.indx){
    # using vol surface to estimate
    cur_dte = as.numeric(as.POSIXct(exp_date) - opt$date[indx], units = 'days')
    tau = cur_dte / 365.25
    # refine the dte to 5 to max_dte that containing cur_dte
    cur_chain <- chain %>% filter(date == opt$date[indx] & dte <= max_dte & dte >= min(cur_dte, 5))
    param = svi.hybridparam(cur_chain)
    P = stiching.optionPrice(opt$strike[indx], param, tau, 'put')
    C = stiching.optionPrice(opt$strike[indx], param, tau, 'call')

    UR <- UR.interp(param, tau, opt$strike[indx])
    param_interp <- svi.interparam(param, tau) # svi parameters at tau
    iv_interp <- svi_jumpwing(log(opt$strike[indx]/UR$mF), param_interp, tau)

    # now fill the NA
    opt$C[indx] = C[1,"value"]
    opt$P[indx] = P[1,"value"]
    opt$dte[indx] = cur_dte
    opt$iv[indx] = iv_interp$impliedvolatility
  }
  opt
}


#' Estimate P/L curves
#'
#' Using obtained svi surface to estimate pl curves
#' @param start The intial date
#' @param exp_date expiration date
#' @param combo a combo option strategy, including three components:
#'     strike, position size, and type (call/put)
#' @param prange price range that to be evaluated
#' @param chain option chain of the intial date
#' @param t days from the intial date >=0.
#' @export
stiching.plcurves <- function(start, exp_date, combo, prange, chain, t = 0){
  start = as.POSIXct(start)
  exp_date=as.POSIXct(exp_date)
  # estimate svi surface parameter
  param = svi.hybridparam(chain)
  #svi.plotskew(chain, param)
  # expand combo with the price range (prange)
  pcombo = combo[rep(seq_len(nrow(combo)),length(prange) ),]
  pcombo$price = as.numeric(rep(prange, each = nrow(combo)))

  max_dte = as.numeric(exp_date - start, units = 'days')
  opts = NULL
  combos = NULL
  # walk through t curves
  for (day in t){
    dte = max_dte - day
    opt <- stiching.optionPrice(pcombo$strike, param,  dte/365.25, pcombo$price, as.character(pcombo$type) )
    cur_pcombo = cbind(pcombo, opt)
    cur_pcombo$dte = dte
    cur_pcombo$t = day
    combos = rbind(combos, cur_pcombo)
  }
  pl.curves <- combos %>% group_by(price, t) %>%
    summarise(pl = 100*sum(value*position), theta = sum(theta*position),
              delta = 100*sum(delta*position), gamma = 100*sum(gamma*position))
  pl.curves
}
