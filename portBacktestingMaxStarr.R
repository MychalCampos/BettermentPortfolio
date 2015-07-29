#####################################################################################
#dynamic portfolio analysis: Black-Litterman, Max Return per Unit ETL
#Does not perform as well as max QU
#####################################################################################

rm(list=ls())

#load necessary Libraries
suppressMessages(library(xts))
suppressMessages(library(quantmod))
suppressMessages(library(PortfolioAnalytics))
suppressMessages(library(mpo))
suppressMessages(library(foreach))
suppressMessages(library(iterators))
suppressMessages(library(ROI))
suppressMessages(library(ROI.plugin.quadprog)) 
suppressMessages(library(ROI.plugin.glpk))

load(file = "returnsData.RData")

# source the required functions
source("allFunctions.R")

##############################
#user inputs
##############################
#usable <- c("VTI", "VTV", "VOE", "VBR", "VEA", "VWO", "LQD", "AGG", "TIP") 
#let's focus on equity exposure for now just to make our main points about
#diversification
usable <- c("VTI", "VTV", "VOE", "VBR", "VEA", "VWO", "LQD", "AGG", "TIP")

##############################
#Portfolio backtesting: Black Litterman vs. Mean-variance
#only 2 constraints: 1) full investment and 2) long only
#Max return per unit ETL 
##############################

# assume that on our month-end rebalancing period we only have the daily returns
# through the previous day
daily.R <- na.omit(lag(dailyReturnsData[, usable]))
riskfree <- na.omit(lag(dailyReturnsData[, "riskfree"]))

# initialize the portfolio objects
funds <- colnames(daily.R)[colnames(daily.R) != "riskfree"]
init.portf <- portfolio.spec(assets=funds)

# # for group constraints (ignore, as doesn't lead to great performance)
# asset.class = c(rep("stocks",6), rep("bonds",3))
# init.portf <- portfolio.spec(assets=funds, category_labels=asset.class)

pspec.fi <- add.constraint(init.portf, type="full_investment")
pspec.fi.lo <- add.constraint(pspec.fi, type="long_only")
pspec.maxStarrLo = pspec.fi.lo
pspec.maxStarrLo = add.objective(portfolio=pspec.maxStarrLo, type="return", 
                                 name="mean")
pspec.maxStarrLo = add.objective(portfolio=pspec.maxStarrLo, type="risk", 
                                 name="ES",
                                 arguments=list(p=0.925))

# # for group constraints (ignore, as doesn't lead to great performance)
# #x% stocks, (1-x)% bonds 
# pspec.maxStarrLo <- add.constraint(portfolio=pspec.maxStarrLo, type="group",
#                          groups=init.portf$category_labels,
#                          group_min=c(0.7, 0.1),
#                          group_max=c(0.9, 0.3))

opt.wts.BL.maxStarr <- matrix(0, (nrow(daily.R) - window.t + 1), ncol(daily.R))
colnames(opt.wts.BL.maxStarr) <- colnames(daily.R)

for (i in 1:(nrow(daily.R) - window.t + 1)) {
  
  t.win <- i:(i+window.t-1)
  
  #estimate betas
  stocks <- daily.R[t.win,]
  rf <- riskfree[t.win]
  R.e <- stocks - matrix(rf, nrow(stocks), ncol(stocks))
  betas <- lm(R.e[ , colnames(R.e) != "VTI"] ~ R.e[, "VTI"])$coef[2,]  
  betas <- c(1, betas)
  names(betas) <- c("VTI", colnames(R.e[,colnames(R.e) != "VTI"]))
  
  port.BL.maxStarr <- optimize.portfolio(R=stocks, pspec.maxStarrLo, 
                                         optimize_method="ROI", 
                        trace=TRUE, momentFUN = "momentsReturnBL", 
                        Betas=betas, market="VTI", riskfree=rf, tscale = 20)
  opt.wts.BL.maxStarr[i, ] <- extractWeights(port.BL.maxStarr)
  
}

opt.wts.BL.maxStarr <- xts(opt.wts.BL.maxStarr, 
                           order.by = index(daily.R)[window.t:nrow(daily.R)])
opt.wts.BL.maxStarr.M <- na.approx(toMonthly(opt.wts.BL.maxStarr))

Diversity.BL.maxStarr <- DIV(opt.wts.BL.maxStarr.M)

monthly.R <- monthlyReturnsData[,usable]
ret.BL.maxStarr <- Return.rebalancing(monthly.R, opt.wts.BL.maxStarr.M) 

ret.Naive <- 0.7*monthly.R[,"VTI"] + 0.3*monthly.R[,"AGG"]
ret.EqW <- xts(apply(monthly.R, 1, mean), order.by=index(monthly.R))
colnames(ret.EqW) <- "Equal-Weighted"

ret.all <- merge(ret.BL.maxStarr, monthly.R[,"VTI"], ret.Naive, ret.EqW)
colnames(ret.all) <- c("BL MaxSTARR", "Market", "Naive", "Equal-Weighted")

png(file="PerformanceSummary.png", width=7, height=5, units="in", res=300)
charts.PerformanceSummary(ret.all, wealth.index = TRUE,
                          lty = c(1,2,3,4), colorset = c("blue","black", "green", 
                                                         "slategray"))
dev.off()


