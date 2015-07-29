#####################################################################################
# laundry list of all functions used
#####################################################################################

##############################
#for extracting price data from quantmod package
##############################
getAdj <- function(symbols){
  require(quantmod)
  # function for extracting adjusted close prices given data extracted from 
  # Yahoo Finance using quantmod package
  for (sym in symbols){
    
    if (sym==symbols[1]) prices = Ad(get(sym))
    else
      prices = merge(prices, Ad(get(sym)))
    
  }
  colnames(prices) <- symbols
  return(prices)
} 

##############################
#functions required for going from daily data to monthly data
##############################

lastWkDayNextMon <- function(date) {
  # function for finding the last weekday of the next month
  x <- as.POSIXlt(date)
  next.mon = as.POSIXlt(seq(x, length=2, by = "1 month")[2])
  last.day = ISOdate(next.mon$year + 1900, next.mon$mon + 1, next.mon$mday)
  if (next.mon$wday == 0) { 
    last.day = last.day + 1
  } else {
    if (next.mon$wday == 6) {
      last.day = last.day - 1  
    }    
  }
  return(as.Date(last.day))
}

toMonthly <- function(dxts) {
  # extract month-end observations from daily xts object
  # with month-end being the last weekday of each month
  require(xts)
  mxts <- dxts[endpoints(dxts)]
  # if last date in mxts is not last weekday of the month, only take the series
  # through the end of the previous month
  t.end <- length(index(mxts))
  if (lastWkDayNextMon(index(mxts)[t.end-1]) != index(mxts)[t.end]){ 
    mxts <- mxts[1:(t.end-1)]
  }
  return(mxts)
}

##############################
#for portfolio optimization
##############################

momentsReturn <- function(R, tscale=1) {
  # this function allows us to calculate mu and sigma over a given rolling window 
  out <- list()
  out$mu <- tscale * apply(R, 2, mean)
  out$sigma <- tscale * cov(R)
  return(out)
}

momentsReturnBL <- function(R, Betas, market, riskfree, tscale=1) {
  # find the expected returns based on the Black-Litterman model
  out <- list()
  rf <- mean(riskfree)
  out$mu <- tscale * (rf + (Betas * (mean(R[, market]) - rf)))
  out$sigma <- tscale*cov(R)
  return(out)
}


##############################
#Summary performance measures for portfolios
##############################

sharpeRatio <- function(r.ts, rf.ts, tscale=1){
  #r.ts and rf.ts must be univariate zoo or xts objects
  #with weekly frequency
  mg <- na.omit(merge(r.ts, rf.ts))
  ret.e <- mg[,1] - mg[,2]
  
  sr <- mean(ret.e) / sd(mg[,1])  
  sr <- sqrt(tscale)*sr
  
  return(sr)
}

SR.StdErr <- function(r.ts, rf.ts, tscale=1) {
  #r.ts and rf.ts must be univariate zoo or xts objects
  r <- coredata(r.ts)
  n <- nrow(r)
  sk <- skewness(r)
  krt <- kurtosis(r)
  sr.w <- sharpeRatio(r.ts, rf.ts, tscale)
  
  sr.var <- tscale*(1 + 0.25*(krt+2)*sr.w^2 - sr.w*sk)/n
  sqrt(sr.var)
  
}

etl <- function(r, alpha = 0.05) {
  #function for calculating non-parametric ETL as in section 4.2 of Ch 4 lecture notes
  r <- sort(r) #r must be a numeric vector
  #this generates the index of the quantile corresponding to the alpha level
  n.tail <- ifelse( alpha == 0, 1, ceiling(alpha*length(r)))
  #ETL expressed as a positive number
  -1/n.tail * sum(r[which((1:length(r)) <= n.tail)])
}

starrRatio <- function(r.ts, rf.ts, alpha = 0.05, tscale=1) {
  #r.ts and rf.ts must be univariate zoo or xts objects
  mg <- na.omit(merge(r.ts, rf.ts))
  ret.e <- mg[,1] - mg[,2]
  
  starr <- mean(ret.e) / etl(as.numeric(mg[,1], alpha=alpha)) 
  starr <- sqrt(tscale)*starr
  
  return(starr)
  
}

geomMean <- function(r, tscale=1) {
  #r must be a vector or univariate xts (or zoo) object
  g.m <- (prod(1+r))^(tscale/length(r)) - 1  
  return(g.m)
}

#compute drawdowns of portfolio
drawdowns <- function(r) {
  out <- list()
  ret.cumul <- cumprod(c(1,(1+r)))
  out$dd <- cummax(ret.cumul) - ret.cumul
  out$rd <- out$dd / cummax(ret.cumul) 
  out$max <- max(out$rd) #maximum relative drawdown
  out$avg <- mean(out$rd) # average relative drawdown
  return(out)
}

TO <- function(weights) {
  # calculate turnover of portfolio over time
  dates = index(weights)
  weights=coredata(weights)
  n.asset=ncol(weights)
  n.dates=nrow(weights)
  if(n.dates<2){
    print("Less than 2 balancing dates!")
    return()
  }
  TurnOver=rep(0,n.dates-1)
  for(i in 1:length(TurnOver)){
    TurnOver[i]=sum(abs(weights[i+1,]-weights[i,]))
  }
  dates=dates[-1]
  res=zoo(TurnOver,order.by = dates)
  res
}

DIV <- function(weights){
  # quantify degree to which portfolio is diversified
  # 1 - Herfindahl Hirschman Index (HHI) 
  if (length(weights) == 1) { 
    # this is the case for "1 point in time" weights
    return(1 - sum(weights^2))
  }
  # code below for when weights is an xts or zoo  object
  n.dates <- nrow(weights)
  if(n.dates<1){
    print("empty data set")
    return()
  }
  diversification=rep(0,n.dates)
  for(i in 1:n.dates){
    diversification[i]=1-sum(weights[i,]^2)
  }
  dates=index(weights)
  Div=zoo(diversification,dates)
  return(Div)
}


performSummary <- function(r.ts, rf.ts, alpha = 0.05, tscale=1) {
  #r.ts and rf.ts must be univariate zoo or xts objects
  mu.ann <- geomMean(r.ts, tscale)
  SR <- sharpeRatio(r.ts, rf.ts, tscale)
  SR.se <- SR.StdErr(r.ts, rf.ts, tscale)
  starr <- starrRatio(r.ts, rf.ts, alpha, tscale)
  maxdd <- drawdowns(r.ts)$max
  avgdd <- drawdowns(r.ts)$avg
  PerfTable <- data.frame(GeomMean=100*mu.ann, 
                          SR=SR, SR.se=SR.se, 
                          STARR=starr, 
                          Max.Drawdown = 100*maxdd, 
                          Avg.Drawdown = 100*avgdd
                          )
  return(t(PerfTable))

}

##############################
# functions needed for gauging robustness of portfolio optimization approach
##############################

sd.shocks <- function(Rt, n.sd=1) {
  # randomly shock 1/2 of observations in a univariate xts object (Rt) by 
  # +/- n.sd standard deviations 
  N.t <- length(Rt)
  half <- round(N.t/2)
  idx <- sample(seq(1,N.t))
  pRt <- Rt
  pRt[idx[1:half]] <- Rt[idx[1:half]] + sign(rnorm(half))*n.sd*sd(Rt)
  return(pRt)
}


