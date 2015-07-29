# BettermentPortfolio
Betterment employs a Black-Litterman (BL) based approach to portfolio optimization. I wanted to focus on getting a better understanding of the benefits of applying BL portfolio optimization to the specific ETFs in the Betterment portfolio. You can find my R code and results for the analysis in this repository. Please note that this is still very much a work in progress.

getDataETF.R: Extracts the required ETF data from Yahoo Finance using the quantmod package and stores the Adjusted Close returns in returnsData.RData (returns are calculated at daily and monthly frequency)

allFunctions.R: All the required custom functions that are sourced by all other R scripts can be found here.  

portBacktestingQU.R: This is where most of the action happens. You can find the code for portfolio backtesting, diversification analysis, and performance analysis here.

portBacktestingMaxStarr.R: Please ignore. Portfolio backtesting in this R script has nothing to do with Black-Litterman

BettermentPortAnalysis.Rmd: This is the rmarkdown file where I write up the analysis. This is still very much in draft form and some sections are not yet complete. BettermentPortAnalysis.pdf is the resulting document. 