#' load package
#'
#' load package with package array
#' @export
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
  usePackage(c('data.table', 'dplyr','tidyr','ggplot2','tidyr','Hmisc',
               'RQuantLib','quantmod', 'broom','parallel',
               'doParallel','foreach',
               'anytime','lubridate'))
}

#' Read ivolatility option chain
#' read option chain raw data downloaded from ivolatility
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @import quantmod
#' @import RQuantLib
#' read ivolatility cvs file.
#' @export
read.ivolchain <- function(filename){
  idat <- fread(filename)
  idat$date <-  as.POSIXct(paste(idat$date, '16:00:00'),
                           format = "%m/%d/%y %H:%M:%S")
  # split option symbol
  idat <- idat %>% separate(., col = `option symbol`,into = c('root','opsymbol'))
  # based on 'root', we know the settle time
  # expiration monthly and weekly are different
  idat$exp_time <- ifelse(idat$root == "SPX","09:30:00", "16:30:00")
  # concatenate settle time to the expiration date
  idat$expiration <-  as.POSIXct(paste(idat$expiration,idat$exp_time),
                                 format = "%m/%d/%y %H:%M:%S")

  if ("stock price for iv" %in% colnames(idat)) {#option with iv data
    ndat <- idat %>%
      rename(callput = `call/put`, spot = `stock price for iv`, oi = `open interest`,
             Bid = bid, Asked = ask) %>%
      dplyr::filter(Asked > spot *0.0002) %>% # remove tiny
      mutate(mPrice = (Bid+Asked)/2, dte = as.numeric(expiration - date, units = 'days'),
             maturity = dte/365.25)
  }  else {# option data
    ndat <- idat %>%
      rename(callput = `call/put`, spot = `unadjusted`, oi = `open interest`,
             Bid = bid, Asked = ask) %>%
      dplyr::filter(Asked > spot *0.0002) %>% # remove tiny
      mutate(mPrice = (Bid+Asked)/2, dte = as.numeric(expiration - date, units = 'days'),
             maturity = dte/365.25)
  }
  # convert to CP paired chain
  chainP <- toChainPairs(ndat) %>% dplyr::filter(dte > 3)
  dates <- as.Date(unique(chainP$date))
  # get Libor interest rate
  getLibors() -> libor
  lb <- libor[dates] #
  if (length(lb) < 1 )
    print('Error in retrieving Libor rate')
  # fill with those NA data
  lb <- na.locf(lb)
  lb[1,] <- lb[2,]
  rates <- lb[,'USD12MD156N']/100
  # calculate chain iv

  sc <- split(chainP, chainP$date)
  lchains <- do.call(rbind, Map(chain.iv, sc, rates))

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
                                      dat <- fread(fp,...)
#                                      dat <- read.csv(fp, ...)
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
                 maturity = dte/365.25)
  #  raw$date <- as.POSIXct(raw$date, format = "%Y-%m-%d %H:%M:%S")
  # expiration monthly and weekly are different
  #raw$exp_time <- ifelse(raw$root == "SPX","09:30:00", "16:30:00")
  #raw$expiration <- as.POSIXct(paste(raw$expiration, raw$exp_time), format = "%Y-%m-%d %H:%M:%S")
  #  raw <- raw %>% mutate(spot = (underlying_bid + underlying_ask)/2,
  #                      dte = as.numeric(expiration - date, units = 'days'),
  #                      maturity = dte/365.25)
  # --- change the above to fast performance
#  raw[,`:=`(date = anytime(date),
#            expiration = anytime(expiration) + ifelse(root == 'SPX', as.ITime("9:30"), as.ITime("16:30")),
#            spot = (underlying_bid + underlying_ask)/2)]
#  raw[, dte := as.numeric(expiration - date, units = 'days')]
#  raw[,date := anytime(date)]
#  raw[,expiration := anytime(expiration) +
#        ifelse(root == 'SPX', as.ITime("9:30"), as.ITime("16:30"))]
#  raw[,spot := (underlying_bid + underlying_ask)/2]
#  raw[,dte := as.numeric(expiration - date, units = 'days')]
#  raw[, maturity := dte/365.25]
  # remove shortest chain
  chainP <- toChainPairs(raw)
  # calculate chain iv
  chain <- chain.iv(chainP, LR)
  chain
}

#' Read optionVue option chain file
#'
#' read optionVue exported cvs file. It should contain:
#' Th.Price, Date, Exp Date, Strike Price,
#' Mid IV, Call/Put, Item Type, Market Price
#' @export
read.vuechain <- function(filename){
  vuedat <- fread(filename)
  vuedat$Date <-  as.POSIXct(paste(as.character(vuedat$Date), vuedat$Time), format = "%y%m%d %H:%M:%S")
  vuedat$Exp.Date <- as.POSIXct(paste(as.character(vuedat$Exp.Date), vuedat$Time), format = "%y%m%d %H:%M:%S")
  # separate Symbol
  ndat <- vuedat %>% separate(Symbol, c('symbol','OpName')) %>%
    rename(meanprice = Th.Price, date = Date, expiration = Exp.Date,
           strike = `Strike Price`, miv = `Mid IV`, callput = `Call/Put`) %>%
    mutate(dte = as.numeric(expiration - date, units = 'days'), maturity = dte/365.25)

  # separate spot price
  spot <- ndat %>% dplyr::filter(`Item Type` == 'I' | `Item Type` == 'S') %>%
    select(symbol, date, `Market Price`) %>% rename(spot = `Market Price`)

  # join the spotPrice
  ndat <- ndat %>% dplyr::filter(`Item Type` == 'O') %>%
    dplyr::filter(!is.na(date) & !is.na(miv)) %>% # filter out extreme values
    left_join(., spot, by = c('symbol','date')) %>%
    select(symbol, date, expiration, strike, callput, meanprice, Bid, Asked, miv, spot, dte,
           maturity)
  ndat
}

#' convert chain to chain call/put pairs
#'
#' remove unpaired chains, average Bid/Ask as mid price
toChainPairs <- function(chain){
  chain %>% dplyr::filter(callput == 'C') %>% mutate(C = (Bid+Asked)/2) %>%
    select(symbol, date, expiration, strike, C) %>%
    left_join(.,
        chain %>% dplyr::filter(callput == 'P') %>% mutate(P = (Bid+Asked)/2) %>%
          select(symbol, date, expiration, strike, P, spot, dte, maturity),
        by = c('symbol','date','expiration','strike')) %>%
    arrange(symbol, date, expiration, strike) %>%
    dplyr::filter(!is.na(P) & C > spot*0.0002 & P > spot*0.0002 & dte > 3) -> chainpair
  as.data.table(chainpair)
}

#' estimate forward price with libor rate
#'
#' calculate forward price for the CP-pair chain with known interest rates
#' Input: paired chain, interest rate
#'
chainForwardPrice <- function(pairChain, LR=0.0169){
 # convert libor rate LR to CC rate R
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
    mutate(logmoneyness = log(strike/mF)) %>% # calculate log moneyness
    dplyr::filter((strike -P <= mF ) & #refine chain
                    (strike + C >= mF) )    -> chain2
  chain2
}

#' estimate iv from optionVue chain
#'
#' @param filename optionVue chain
#' @param LR libor rate
#'
#' @return refined chain with estimated IV
#'
#' @export
vue.iv <- function(filename, LR){
  rawChain <- toChainPairs(readVueChain(filename))
  ivChain <- chain.iv(rawChain, LR)
  ivChain
}

#' Calculate chain iv, F, U
#'
#' Base on BS model
#' @param chain call/put paired option chain
#' @param LR liborrate
#'
#' @export
chain.iv <- function(chain, LR=0.01){
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
  clusterEvalQ(cl, {library(RQuantLib)})

  chain <- chainForwardPrice(chain, LR)
  ivc = rep(0,nrow(chain))
  ivp = rep(0,nrow(chain))
  foreach(idx = 1:nrow(chain)) %dopar% {
#  for (idx in 1:nrow(chain)){
    t = chain[idx,]
    ivc[idx] = EuropeanOptionImpliedVolatility('call', value = t$C, t$U, t$strike, 0,t$R,t$maturity,0.1)
    ivp[idx] = EuropeanOptionImpliedVolatility('put', value = t$P, t$U, t$strike, 0,t$R,t$maturity,0.1)
  }
  #stop cluster
  stopCluster(cl)

  chain$ivc = ivc
  chain$ivp = ivp
  chain$iv = (chain$ivc+chain$ivp)/2
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

#' Import livevol zip file and generate chain and surface
#'
#' @param filename filename of an input zip file
#' @param rate  Interest rate at that day
#' @return a list contain the chain and surface, and save two to files
#' @export
livevol.chain <- function(filename, rate = 0.017){
  # read rawdate from zip file
  r_date = gsub('.*_(.*).zip', '\\1', filename)
  print(r_date)
  # convert to chain
  chain <- read.livevol(filename, LR = rate)
  saveRDS(chain, file = paste0('spxchain_',r_date,'.rds'))

  # estimate surface
#  surface = data.table();
#  for (d in unique(chain$date)){
#    print(as.POSIXct(d, origin = '1970-01-01'))
#    c1 <- chain %>% filter(date == d & dte > range[1] & dte < range[2])
#    ptm = proc.time()
#    s1 <- svi.quasiparam.p(c1)
#    print(proc.time()-ptm)
#    if (show == TRUE)
#        print(quasi.plotsurface(c1, s1, type = 'ivt'), range = range)
#    surface = rbind(surface,s1)
#  }
#  saveRDS(surface, file = paste0('spxsurf_',r_date,'.rds'))
  return(chain)
}

#' Import optionVue csv file and generate chain and surface
#'
#' @param filename filename of an input zip file
#' @param rate  Interest rate at that day
#' @return a list contain the chain and surface, and save two to files
#' @export
vue.chain <- function(filename, rate = 0.017){
  rawchain = toChainPairs(read.vuechain(filename))
  chain <- chain.iv(rawchain, rate)
  saveRDS(chain, file = paste0('chain_',filename,'.rds'))
#  surface <- svi.quasiparam.p(chain)

#  quasi.plotsurface(chain, surface, range = range, type = 'ivt')
#  saveRDS(surface, file = paste0('surface_',filename, '.rds'))
#  return(list(chain = chain, surface = surface))
  return(chain)
}

