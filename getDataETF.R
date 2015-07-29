#####################################################################################
#Get the price data. Here I grab as many tickers as I can find using different online 
#sources.   
#I also grab the grab 3-month daily US Tbill rates for calculating the risk-free 
#rates of return.
#Calculate monthly and daily returns.
#####################################################################################

rm(list=ls())

suppressMessages(library(xts))
suppressMessages(library(quantmod))

##############################
#user inputs here
##############################

#list of etf symbols we may want to include in portfolio
usable = c("VTI", "VTV", "VOE", "VBR", "VEA", "VWO", "LQD", 
           "AGG", "TIP", "VNQ", "DBC")

# length of rolling time window for estimating rolling CAPM betas
window.t <- 60 

# source the required functions
source("allFunctions.R")

##############################
# download ETF data from Yahoo Finance
# calculate daily and monthly total returns
##############################

# first 10 are the ETFs in the Betterment portfolio
# the rest of the ETFs are included to:
# achieve similiar risk exposures over longer histories or
# to augment the overall portfolio (e.g., with DBC and VNQ)

symbols = c("VTI", "VTV", "VOE", "VBR", "VEA", "VWO", "LQD", "BND", "BNDX", "VWOB",
            "IJH", # longer history than VOE but without value tilt
            "AGG", # longer history than BND 
            "TIP", # for exposure to TIPS markets (inflation factor)
            "VNQ", # REIT for exposure to US real estate markets
            "DBC"  # broad commodities exposure just for kicks (expensive!)
            )

getSymbols(Symbols=symbols, from="2001-05-31", adjust=TRUE, auto.assign=TRUE)

adjClose <- getAdj(symbols)
# extract month-end data
adjClose.M <- toMonthly(adjClose)
  
# daily total returns (arithmetic) of ETFs
etf.Rt <- xts(apply(adjClose, 2, Delt), order.by = index(adjClose)) 

# monthly total returns (arithmetic) of ETFs
etf.Rt.M <- xts(apply(adjClose.M, 2, Delt), order.by = index(adjClose.M)) 

# daily
# if VOE return is NA, then replace with IJH return
# note that monthly returns have high correlation of 0.97
# achieves exposure to US Mid Cap stocks with VOE having a value tilt not present in 
# IJH
na.idx <- is.na(etf.Rt[ ,"VOE"])
etf.Rt[na.idx, "VOE"] <- etf.Rt[na.idx, "IJH"]
etf.Rt <- etf.Rt[, usable]

# monthly
# if VOE return is NA, then replace with IJH return
# note that monthly returns have high correlation of 0.97
# achieves exposure to US Mid Cap stocks with VOE having a value tilt not present in 
# IJH
na.idx <- is.na(etf.Rt.M[ ,"VOE"])
etf.Rt.M[na.idx, "VOE"] <- etf.Rt.M[na.idx, "IJH"]
etf.Rt.M <- etf.Rt.M[, usable]

##############################
# Download 3M Treasury Bill Rate data from FRED
# calculate daily and monthly total returns
##############################

getSymbols('DTB3', src='FRED', auto.assign=TRUE)
dtb3.xts <- window(DTB3, start = "2001-05-31")
# calculate daily returns based on actual/365 day count convention
riskfree <- lag(dtb3.xts)/100 * c(NA, diff(index(dtb3.xts)))/365
colnames(riskfree) <- "riskfree"

# extract month-end data
mtb3.xts <- toMonthly(dtb3.xts)
# calculate monthly returns based on actual/365 day count convention
riskfree.M <- lag(mtb3.xts)/100 * c(NA, diff(index(mtb3.xts)))/365
colnames(riskfree.M) <- "riskfree"


##############################
# put final data set together for daily and monthly returns
##############################

dailyReturnsData <- na.omit(merge(etf.Rt, riskfree))
monthlyReturnsData <- na.omit(merge(etf.Rt.M, riskfree))

save(file = "returnsData.RData", dailyReturnsData, monthlyReturnsData, window.t)
