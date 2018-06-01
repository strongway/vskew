
usePackage <- function(pk){
  for (p in pk) {
    if (!is.element(p, installed.packages()[,1])) install.packages(p, dependencies =TRUE)
    library(p, character.only = TRUE)
  }
}

#' load necessary packages for vskew
#'
#' @export
loadSkewPackages <- function(){
  usePackage(c('data.table', 'tidyverse', 'plotly',
               'quantmod', 'broom','parallel',
               'doParallel','foreach',
               'anytime','lubridate'))
}

#' Read ivolatility option chain
#' read option chain raw data downloaded from ivolatility
#' read ivolatility cvs file.
#' @export
read.ivolchain <- function(filename){
  idat <- fread(filename)
  idat$date <-  as.POSIXct(paste(idat$date, '16:30:00'),
                           format = "%m/%d/%Y %H:%M:%S")
  # split option symbol
  idat <- idat %>% separate(., col = `option symbol`,into = c('root','opsymbol'))
  # based on 'root', we know the settle time
  # expiration monthly and weekly are different
  idat$exp_time <- ifelse(idat$root == "SPX","09:30:00", "16:30:00")
  # concatenate settle time to the expiration date
  idat$expiration <-  as.POSIXct(paste(idat$expiration,idat$exp_time),
                                 format = "%m/%d/%Y %H:%M:%S")

  if ("stock price for iv" %in% colnames(idat)) {#option with iv data
    ndat <- idat %>%
      rename(callput = `call/put`, spot = `stock price for iv`, oi = `open interest`,
             Bid = bid, Asked = ask) %>%
      dplyr::filter(Asked > spot *0.0002) %>% # remove tiny
      mutate(mPrice = (Bid+Asked)/2, dte = as.numeric(expiration - date, units = 'days'),
             maturity = dte/365)
  }  else {# option data
    ndat <- idat %>%
      rename(callput = `call/put`, spot = `unadjusted`, oi = `open interest`,
             Bid = bid, Asked = ask) %>%
      dplyr::filter(Asked > spot *0.0002) %>% # remove tiny
      mutate(mPrice = (Bid+Asked)/2, dte = as.numeric(expiration - date, units = 'days'),
             maturity = dte/365)
  }
  # convert to CP paired chain
  chainP <- toChainPairs(ndat) %>% dplyr::filter(dte > 3)
  dates <- as.Date(unique(chainP$date))
  # get Libor interest rate
  libor <- getSymbols('USD12MD156N',src = 'FRED', auto.assign = F)
  lb <- libor[dates] #
  if (length(lb) < 1 )
    print('Error in retrieving Libor rate')
  # fill with those NA data
  lb <- na.locf(lb)
  lb[1,] <- lb[2,]
  rates <- lb[,'USD12MD156N']/100
  # calculate chain iv

  sc <- split(chainP, chainP$date)
  lchains <- do.call(rbind, Map(estimate.iv, sc, rates))

  ofile <- gsub('.csv','_chain.rds',filename)
  saveRDS(lchains, file = ofile)

  lchains

}

#' Read csv data from zip file
#'
#' Read a large csv file from a zip file
#' If the zip file contains multiple csv files, all files will be read and concatenated.
#' @export
#'
#' @param zipfile filename of the zip file
#' @param pattern the pattern that match the file name, such as '.csv'
#'
read.csv.zip <- function(zipfile, pattern="\\.csv$", ...){

  # Create a name for the dir where we'll unzip
  zipdir <- tempfile()

  # Create the dir using that name
  dir.create(zipdir)

  # Unzip the file into the dir
  unzip(zipfile, exdir=zipdir)

  # Get a list of csv files in the dir
  files <- list.files(zipdir, rec=TRUE, pattern=pattern)

  # Create a data.table of the imported csv files
  csv.data <- do.call(rbind, lapply(files,
                                    function(f){
                                      fp <- file.path(zipdir, f)
#                                      dat <- fread(fp,...) # problem with 'Could not find first good line start '
                                      dat <- read.csv(fp, stringsAsFactors = FALSE)
                                      return(dat)
                                    }) )


  # Use csv names to name list elements
  #names(csv.data) <- basename(files)

  # Return data
  return(csv.data)
}

#' Read livevol zip csv chain
#'
#' Import csv raw chain, and calculate iv based on bid/ask
#' return a refined chain
#' @export
read.livevol <- function(filename, LR = 0.017){
  raw <- read.csv.zip(filename)
  raw <- raw[,c(1:6,12:17)] # discard the last few duplicated columns
  setnames(raw, c("option_type","bid","ask","underlying_symbol","quote_datetime"),
           c("callput","Bid","Asked","symbol","date"))
  # change time format
  raw <- raw %>% mutate(date = anytime(date),
                 expiration = anytime(expiration) + ifelse(root == 'SPX', as.ITime("9:30"), as.ITime("16:30")),
                 spot = (underlying_bid + underlying_ask)/2,
                 dte = as.numeric(expiration - date, units = 'days'),
                 maturity = dte/365)
  # remove shortest chain
  chainP <- toChainPairs(raw) %>% filter(dte > 3)
  # calculate chain iv
  chain <- estimate.iv(chainP, LR)
  chain
}

# convert chain to chain call/put pairs
# remove unpaired chains, average Bid/Ask as mid price
toChainPairs <- function(chain){
  chain %>% dplyr::filter(callput == 'C') %>% mutate(C = (Bid+Asked)/2) %>%
    select(symbol, date, expiration, strike, C) %>%
    left_join(.,
        chain %>% dplyr::filter(callput == 'P') %>% mutate(P = (Bid+Asked)/2) %>%
          select(symbol, date, expiration, strike, P, spot, dte, maturity),
        by = c('symbol','date','expiration','strike')) %>%
    arrange(symbol, date, expiration, strike) %>%
    dplyr::filter(!is.na(P) & C > spot*0.0002 & P > spot*0.0002 ) -> chainpair
  as.data.table(chainpair)
}

# estimate forward price with libor rate
# calculate forward price for the CP-pair chain with known interest rates
# Input: paired chain, interest rate
chainForwardPrice <- function(pairChain, LR=0.0169){
 # convert libor rate LR to continuously compounded (CC) rate R
  # estimate Forward price F based on Call/put parity
  fChain <- pairChain %>% mutate(R = log(1+LR*maturity)/maturity) %>%
    mutate( F = strike + exp(R*maturity)*(C-P), cp = C-P) %>%
    arrange(symbol, date, expiration, cp)

  # clean chain outliers (+/- 1.5 sd)
  fChain %>% group_by(symbol, date, expiration) %>%
    mutate(mF = mean(F), sdF = sd(F)) %>%
    dplyr::filter((F > mF-1.5*sdF) & (F<mF+1.5*sdF)  ) %>%
    select(-mF, -sdF)-> fChain
  # balance the negative and positive parities
  fChain %>% group_by(symbol, date, expiration) %>%
    mutate(netParity = cumsum(cp)) %>%
    dplyr::filter(netParity < 0.01*spot) %>% # remove those two far put side
    group_by(symbol, date, expiration) %>%
    summarise(mF = mean(F))-> FPrice
  # join in backforward price, calculate underlying price
  left_join(fChain, FPrice, by = c('symbol','date','expiration')) %>%
    mutate(U = exp(-R*maturity)*mF) %>% dplyr::filter(!is.na(U)) %>% # combined mean F, and U
    mutate(logmoneyness = log(strike/mF))     -> chain2
  chain2
}

#' Estimate chain ivs, F, U
#'
#' Using bisection method and tighter bounds proposed by Gatheral et al. 2017
#' @param chain call/put paired option chain
#' @param LR liborrate
#'
#' @export
estimate.iv <- function(chain, LR=0.01){

  chain <- chainForwardPrice(chain, LR)
  # normalized call and put prices
  rc = chain$C/chain$U
  rp = chain$P/chain$U
  k = chain$logmoneyness
  lb_c = l_bs(rc,k) #lower bound
  ub_c = u_bs(rc,k) # upper bound
  lb_p = l_bs(rp,k, cp = -1)
  ub_p = u_bs(rp,k, cp = -1)
  tiv_c = bisection(rc,k, 1, lb_c, ub_c)
  tiv_p = bisection(rp,k, -1, lb_p, ub_p)

  chain$ivc <- tiv_c/sqrt(chain$maturity)
  chain$ivp <- tiv_p/sqrt(chain$maturity)

  # use otm iv, and nest the chain data (2017-08-25)
  chain <- chain %>%mutate(iv = ifelse(logmoneyness > 0, ivc,ivp)) %>%
    select(-c(cp, F)) %>%
    group_by(symbol, date, expiration, dte, spot, R, U) %>%
    nest()
  chain
}

#' Read trade log text file

#' @return stock transaction,
#'    option transaction,
#'    stock position,
#'    option position
#'    to the global environment
#' @export
readTradeLog <- function(logfile){
  lines <- readLines(logfile, n =-1)
  sdat <- split(lines, cumsum(nchar(lines)==0))
  sdat[1] <- NULL # remove account info
  sdat[length(sdat)] <- NULL # remove last 'eof'
  lapply(sdat , function(lins) {
    good<- lins[nchar(lins)>0]
    if (length(good) > 0) {
      df  <-  read.table(text=good[-1], sep = '|')
      if (ncol(df)==12) { # position
        names(df) <- c('type','account','symbol','fullSymbol','currency',
                       'd','t','position','quantity','price','value','v1')
      } else if (ncol(df)==16) {# transaction
        names(df) <- c('type','tr_no','symbol','fullSymbol','location',
                       'sb','co','d','t','currency','position',
                       'quantity','price','value','fee','v1')
      }
      # convert date, time columns
      df$date  <-  as.POSIXct(paste(as.character(df$d), df$t), format = "%Y%m%d %H:%M:%S")
      assign(make.names(good[1]),   #name
             df , envir=.GlobalEnv)
      }
    }
    )
  return('loaded into global Environment')
}

#' select option transactions from trade log table
#'
#' @param op_df option table
#' @param symbol target symbol, e.g., 'SPX'
#' @return option transaction with expiration, strike, callput, ...
#' @export
optTransactions <- function(op_df,sym) {
  df <- dplyr::filter(op_df, grepl(sym, symbol)) %>%
    separate(fullSymbol,c('stk','expiration','strike','callput'),sep = ' ') %>%
    select(one_of(c('stk','date','expiration','strike','callput','sb','position',
                    'quantity','price','value','fee'))) %>%
    arrange(date)
  df
}

#' Read Libor data from FRED
#'
#' 6, and 12 months libor rates, return as a table.
getLibors <- function(){
  getSymbols('USD6MTD156N',src = 'FRED')
  getSymbols('USD12MD156N',src = 'FRED')
  merge(USD6MTD156N,USD12MD156N)
}




