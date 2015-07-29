#####################################################################################
#dynamic portfolio analysis: Black-Littermen, Max QU 
#####################################################################################

rm(list=ls())

#load necessary Libraries
suppressMessages(library(xts))
suppressMessages(library(quantmod))
suppressMessages(library(PortfolioAnalytics))
suppressMessages(library(mpo))
suppressMessages(library(PerformanceAnalytics))
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

# need this parameter as we are maximizing quadratic utility in our optimizations
# risk aversion of the market (i.e., typical participant in the equity market
# based on calculations in Ang(2014), section 3.3 (location 1201 in kindle version)
risk.aversion <- 3 

##############################
#Portfolio backtesting: Black Litterman vs. Mean-variance
#only 2 constraints: 1) full investment and 2) long only
#maximize quadratic utility with typical investor's risk aversion
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

pspec.fi.lo.QU <- add.objective(portfolio=pspec.fi.lo, type="quadratic_utility",
                                risk_aversion=risk.aversion)

# # for group constraints (ignore, as doesn't lead to great performance)
# #x% stocks, (1-x)% bonds 
# pspec.fi.lo.QU <- add.constraint(portfolio=pspec.fi.lo.QU, type="group",
#                          groups=init.portf$category_labels,
#                          group_min=c(0.7, 0.1),
#                          group_max=c(0.9, 0.3))

opt.wts.BL <- matrix(0, (nrow(daily.R) - window.t + 1), ncol(daily.R))
colnames(opt.wts.BL) <- colnames(daily.R)

opt.wts.MeanVar <- opt.wts.BL

for (i in 1:(nrow(daily.R) - window.t + 1)) {
  
  t.win <- i:(i+window.t-1)
  
  #estimate betas
  stocks <- daily.R[t.win,]
  rf <- riskfree[t.win]
  R.e <- stocks - matrix(rf, nrow(stocks), ncol(stocks))
  betas <- lm(R.e[ , colnames(R.e) != "VTI"] ~ R.e[, "VTI"])$coef[2,]  
  betas <- c(1, betas)
  names(betas) <- c("VTI", colnames(R.e[,colnames(R.e) != "VTI"]))
  
  port.BL <- optimize.portfolio(R=stocks, pspec.fi.lo.QU, optimize_method="ROI", 
                     trace=TRUE, momentFUN = "momentsReturnBL", 
                     Betas=betas, market="VTI", riskfree=rf, tscale = 20)
  opt.wts.BL[i, ] <- extractWeights(port.BL)
  
  port.MeanVar <- optimize.portfolio(R=stocks, pspec.fi.lo.QU, optimize_method="ROI", 
                                     trace=TRUE, momentFUN = "momentsReturn", 
                                     tscale = 20)
  opt.wts.MeanVar[i, ] <- extractWeights(port.MeanVar)
  
}

opt.wts.BL <- xts(opt.wts.BL, order.by = index(daily.R)[window.t:nrow(daily.R)])
opt.wts.BL.M <- na.approx(toMonthly(opt.wts.BL))

opt.wts.MeanVar <- xts(opt.wts.MeanVar, 
                       order.by = index(daily.R)[window.t:nrow(daily.R)])
opt.wts.MeanVar.M <- na.approx(toMonthly(opt.wts.MeanVar))


##############################
# Performance Summary
##############################

monthly.R <- monthlyReturnsData[,usable]
ret.BL <- Return.rebalancing(monthly.R, opt.wts.BL.M) 
ret.MeanVar <- Return.rebalancing(monthly.R, opt.wts.MeanVar.M) 
ret.Naive <- 0.7*monthly.R[,"VTI"] + 0.3*monthly.R[,"AGG"]
ret.EqW <- xts(apply(monthly.R, 1, mean), order.by=index(monthly.R))
colnames(ret.EqW) <- "Equal-Weighted"

ret.all <- merge(ret.BL, ret.MeanVar, monthly.R[,"VTI"], ret.Naive, ret.EqW)
colnames(ret.all) <- c("Black-Litterman", "Mean-Variance", "Market", "Naive (70/30)",
                       "Equal-Weighted")

# summary plot
png(file="PerformanceSummary.png", width=7, height=5, units="in", res=300)
charts.PerformanceSummary(ret.all, wealth.index = TRUE,
                          lty = c(1,2,3,4,5), colorset = c("blue","black", "green", 
                                                         "slategray", "skyblue"),
                          main = "Performance: Black Litterman vs. 4 Benchmarks")
dev.off()

# create summary table
for (i in 1:ncol(ret.all)) {
  
  if (i == 1) {
    perf.tbl <- performSummary(na.omit(ret.all[,i]), 
                               monthlyReturnsData[ , "riskfree"], 
                               alpha=0.05,
                               tscale=12)
  } else {
    perf.tbl <- cbind(perf.tbl,
                      performSummary(na.omit(ret.all[,i]), 
                                     monthlyReturnsData[ , "riskfree"], 
                                     alpha=0.05,
                                     tscale=12))
  }
}
perf.tbl <- round(perf.tbl, 2)
colnames(perf.tbl) <- colnames(ret.all)


##############################
# Diversification analysis
# Comparing Black-Litterman to Mean-Variance
##############################

Diversity <- merge(DIV(opt.wts.BL.M), DIV(opt.wts.MeanVar.M))
colnames(Diversity) <- c("Black-Litterman", "Mean-Variance")

#find and plot a good example of the discrepancies between the 2 methods
idx <- which(Diversity[,1] == max(Diversity[,1]))

png(file="wts.png", width=7, height=5, units="in", res=300)

par(mfrow=c(1,2))

barplot(100*opt.wts.BL.M[idx], beside=TRUE, 
        names.arg=funds, las=2, cex.names=0.95,
        ylab="optimal portfolio weights (%)",
        ylim = c(0,100), 
        main = paste("Black Litterman Portfolio",
                     paste("\n", index(opt.wts.BL.M)[idx])),
        cex.main=0.9,
        col="skyblue3")

barplot(100*opt.wts.MeanVar.M[idx], beside=TRUE, 
        names.arg=funds, las=2, cex.names=0.95,  
        ylab="optimal portfolio weights (%)",
        ylim = c(0,100), 
        main = paste("Mean Variance Portfolio",
                     paste("\n", index(opt.wts.BL.M)[idx])),
        cex.main=0.9,
        col="skyblue3")

par(mfrow = c(1,1))

dev.off()


# box plots for plotting (1 - HHI(w)) over time
png(file="Diversification.png", width=7, height=5, units="in", res=300)
chart.Boxplot(Diversity, main= paste("Diversification of Portfolios", 
                                     paste("\n(", paste(start(Diversity), "to", 
                                           end(Diversity), sep = " "), ")", sep="")),
                                     xlab = "1 - HHI(w)")
dev.off()



