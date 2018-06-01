# Black-Scholes model and implied volatility

#' BS total IV approximation (tight lower bound)
#'
#' Black-Scholes implied volatility approximation, adopted from Stefanica and Radoicic (2017)
#' An explicit implied volatility formula, the approximation error < 2%, maximum error < 10%
#' @param c normalized option price C/D/F (call), or P/D/F (put). D = e^(-rt) discount, F forward price
#' @param k log(K/F)
#' @param cp call/put indicator, call = 1, put = -1
#' @export
l_bs <- function(c,k, cp = 1){
  # change k to k>0, put to call
  para = parity_change(c,k, cp)
  c = para$c
  k = para$k
  cp = para$cp

  R = 2*c* exp(-k) - exp(-k) + 1
  A = (exp((1-2/pi)*k) - exp(-(1-2/pi)*k))^2
  B = 4*(exp(2/pi*k)+ exp(-2/pi*k)) - 2*exp(k)*
    (exp((1-2/pi)*k) + exp(-(1-2/pi)*k)) * (exp(-2*k) + 1 -R^2)
  C = exp(2*k)*(R^2 - (exp(-k)-1)^2)*((exp(-k)+1)^2-R^2)
  gamma = -pi/2 * log(2*C/(B + sqrt(B^2 + 4*A*C)))

  C0 = 0.5 - exp(k)*apx(-sqrt(2*k))
  sigs = sqrt(gamma + k) + (-1)^(c<=C0) * sqrt(gamma -k)
  sigs[is.na(sigs)] = 0.05 # replace NA
  sigs
}

#' BS total IV approximation (upper bound)
#'
#' The upper bound of total IV, given by Gatheral et al., 2017
#' @param c normalized option price C/D/F (call), or P/D/F (put). D = e^(-rt) discount, F forward price
#' @param k log(K/F)
#' @param cp call/put indicator, call = 1, put = -1
#' @export
u_bs <- function(c,k, cp = 1){
  # change k to k>0, put to call
  para = parity_change(c,k, cp)
  c = para$c
  k = para$k
  cp = para$cp

  a = (1-c)/(1+exp(k)) # eq. 3.4
  phi = qnorm(a, 0,1) # eq. 3.6
  theta = - dnorm(phi,0,1)/phi * (exp(k)-1)/(exp(2*k)+1)

  sigs = - qnorm(a + theta * exp(k)) - qnorm(a - theta)  # eq. 3.3
  sigs[is.na(sigs)] = 0.99 # replace NA
  sigs
}

# Polya approximation function A
apx <- function(x){
  0.5 + 0.5*(-1)^(x<0)*sqrt(1-exp(-2*x^2/pi))
}

# put-call parity: c = p + 1 - exp(k)
# change put to call, k<0 to k>0
parity_change <- function(c,k, cp){
  # make all parameters same length
  n = max(NROW(c),NROW(k),NROW(cp))
  c = rep(c, length.out = n)
  k = rep(k, length.out = n)
  cp = rep(cp, length.out = n)

  idx = cp != 1 # select put price
  c[idx] = c[idx] + 1 - exp(k[idx]) # change put to call price

  # if k < 0, use call-put parity, change to k > 0
  idx = k < 0
  c[idx] = exp(-k[idx])*c[idx] + 1 - exp(-k[idx])
  k[idx] = -k[idx]
  list(c = c, k = k, cp = cp)
}

# normalized call price (i.e., C/D/F)
cbs <- function(k, sig){
  dp = - k / sig + sig/2
  dn = - k / sig - sig/2
  pnorm(dp) - exp(k)* pnorm(dn)
}

#' Using bisection method to find the iv
#'
#' Vectoring IV estimation with the bisection method
#' @param cm  option market mid price
#' @param k log(K/F)
#' @param cp call/put flag (call 1, put -1)
#' @param vl lower bound of total IV
#' @param vh upper bound of total IV
#' @export
bisection <- function(cm, k,cp = 1, vl, vh){
  # change k to k>0, put to call
  para = parity_change(cm,k, cp)
  cm = para$c
  k = para$k
  cp = para$cp

  eps = 1e-6
  flags = rep(FALSE, NROW(cm)) # flag for the exit
  cl = rep(0, NROW(cm))
  ch = rep(0, NROW(cm))
  vi = rep(0, NROW(cm))
  ct = rep(0, NROW(cm))
  # small number of iteration
  for (i in 1:20) { # maximum 20 times
    idx = flags == FALSE # not yet converged
    cl[idx] = cbs(k[idx], vl[idx])
    ch[idx] = cbs(k[idx], vh[idx])
    # avoid NA
    dc = ch[idx] - cl[idx]
    dc = ifelse(dc<1e-5, 1e-5, dc)
    vi[idx] = vl[idx] + (cm[idx] - cl[idx])/dc * (vh[idx]-vl[idx])
    ct[idx] = cbs(k[idx], vi[idx])
    flags[idx] = abs(cm[idx]-ct[idx]) < eps
    if (sum(flags) == NROW(flags)) break
    vl[idx] = ifelse(ct[idx]<cm[idx], vi[idx],vl[idx])
    vh[idx] = ifelse(ct[idx]>cm[idx], vi[idx],vh[idx])
  }
  #  print(i)
  vi
}

#' Calculate BS option price
#' @param forward Forward price
#' @param strike strike
#' @param tau maturity
#' @param iv  implied volatility
#' @param cp call/put flag (call 1, put -1)
#' @param R interest rate
#' @export
bsPrice <- function(forward, strike, tau, iv, cp = 1, R =0.02){
  if (NROW(strike) > NROW(cp))  cp = rep(cp, length.out = NROW(strike))
  k = log(strike/forward) # log moneyness
  sig = iv*sqrt(tau) # total variance

  call = exp(-tau*R)*forward*cbs(k,sig)
  put = call - exp(-tau*R)*(forward - strike)
  ifelse(cp==1, call, put)
}

