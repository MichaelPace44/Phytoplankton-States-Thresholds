# ADF (Augmented Dickey Fuller) stationarity test
# and ARIMA for Markov lag
# 

rm(list = ls())
graphics.off()

library(stats)
library(tseries)
library(forecast)

#Load data based on DLM results
#Lines provided are for Peter Lake DOsat
#load(file="DLM_DOsat_median_Peter13.Rdata")
#load(file="DLM_DOsat_median_Peter14.Rdata")
load(file="DLM_DOsat_median_Peter15.Rdata")
#load(file="DLM_DOsat_median_Peter19.Rdata")
#load(file="DLM_DOsat_median_Peter24.Rdata")

#Assign variate for analysis and set-up time
Xvar0 = X0
nx = length(Xvar0)
Tstep0 = Tstep[1:nx]

#Plot time series
windows()
plot(Tstep0,Xvar0,type='l',lwd=1,col='blue')

print('tails of Xvar0',quote=F)
print(range(Xvar0))

# find optimal AR lags using auto.arima from forecast library
arfit = auto.arima(Xvar0) # use defaults

lagopt = arimaorder(arfit)
print('optimal order using autoarima() and arimaorder()',quote=F)
print(lagopt)
aropt = unname(lagopt[1])  # save optimal AR order
# if optimal lag is 0 then data are uncorrelated, use original data
aropt = ifelse(aropt==0,1,aropt)  

# subsample Xvar0 according to lagopt
nx2 = length(Xvar0)
ikeep = seq(1,nx2,by=aropt)
Xvar0 = Xvar0[ikeep]
Tstep = Tstep0[ikeep]

# check AR order of thin data
arfit1 = auto.arima(Xvar0)
lagopt1 = arimaorder(arfit1)
print('',quote=F)
print('optimal order of thinned data using autoarima() and arimaorder()',quote=F)
print(lagopt1)

# test stationarity
ADF.result = adf.test(Xvar0,alternative='stationary')
pvalue = ADF.result$p.value

print('ADF for thinned series',quote=F)
print('p value for HO: series NOT stationary',quote=F)
print(pvalue)
