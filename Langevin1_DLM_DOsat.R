# Program runs a dynamic linear model on chosen lake, variable, year combination
# Results are saved as an input file for B-spline analysis
# When using a regular variable like DOsat for B-spline, results from this program are not used
# not further except to provide input file for ADF test and B-spline programs 
#  


rm(list = ls())
graphics.off()

#Function needed for DLM 
source('ODLMAR_NoBoot_2018-10-20.R')

# LOAD DATA
# Select data file with lake, year, variable needed
# Lines included here are for daily DOsat from each year of Peter Lake
# Data files are the output from 'DailyDataStats_DOsat' program

#load(file="Daily_R_DOsat2013.Rdata")
#load(file="Daily_R_DOsat2014.Rdata")
load(file="Daily_R_DOsat2015.Rdata")
#load(file="Daily_R_DOsat2019.Rdata")
#load(file="Daily_R_DOsat2024.Rdata")

#Select file name for output; lines included here are for Peter Lake DOsat 
#Fname = c('DLM_DOsat_median_Peter13.Rdata')
#Fname = c('DLM_DOsat_median_Peter14.Rdata')
Fname = c('DLM_DOsat_median_Peter15.Rdata')
#Fname = c('DLM_DOsat_median_Peter19.Rdata')
#Fname = c('DLM_DOsat_median_Peter24.Rdata')

#Title for Plots
title = c('Dissolved Oxygen Saturation')

# select data; 'DOx' is from the loaded data file
X0 = DOx$DOmed
nx = length(X0)
Tstep = DOx$doy

windows()
plot(Tstep,X0,type='l',lwd=1,col='seagreen',main='original data')

# transform - optional
X.dlm = (X0 - mean(X0))/sd(X0)  # Z-score
print(c('range X.dlm = ',range(X.dlm,na.rm=T)),quote=F)

# Start DLM
windows(width=12,height=6)
plot(Tstep,X.dlm,type='l',col='forestgreen',xlab='DoY index',ylab='X.dlm',
     main=title)
grid()

# Set up DLM
nobs = length(X.dlm)
nl = 1 # number of lags
print('**************',quote=F)
print(' ',quote=F)
print('**************',quote=F)
delta = 0.97 # 0<delta<1; see advice in functions

# Run DLM
ODL.out = ODLMAR(nl,delta,X.dlm,Tstep,title)

# Output matrices are stored sideways, like MARSS
Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]  
B.sd = ODL.out[[4]]
errvar = ODL.out[[5]] # updated error variance

# Post process DLM -----------------------------------------------

# Calculate moving equilibrium
X.eq = B.ests[1,]/(1 - B.ests[2,])
# Calculate its variance
deno = (1 - B.ests[2,])
Vterm1 = (B.sd[1,]*B.sd[1,] + errvar)/(deno*deno)
Vterm2 = ((B.ests[1,]*B.sd[2,])/(deno*deno))*((B.ests[1,]*B.sd[2,])/(deno*deno))
SD.eq = sqrt(Vterm1+Vterm2)
  
# Z score
Z.eq = X.eq/SD.eq

# Time steps start at 2
Nstep = length(Tstep)

# Plot components of steady-state estimate
windows(width=6,height=12)
par(mfrow=c(3,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],X.eq,type='l',col='blue',ylim=c(-10,10),
     ylab='Steady State',xlab='DoY',
     main='Local Steady-State estimate, sd, and ratio')
grid()
plot(Tstep[2:Nstep],SD.eq,type='l',col='red',ylim=c(0,10),
     ylab='S.D.',xlab='DoY')
grid()
plot(Tstep[2:Nstep],Z.eq,type='l',col='purple',
     ylab='Z score',xlab='DoY')
grid()

# Calculate level estimates
level = B.ests[1,]
levelsd = B.sd[1,]
stdlevel = B.ests[1,]/B.sd[1,]

# Plot components of level estimate
windows(width=12,height=9)
par(mfrow=c(3,1),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tstep[2:Nstep],level,type='l',col='deepskyblue',#ylim=c(-10,10),
     ylab='level',xlab='DoY',
     main='Level and Std level estimate')
grid()
plot(Tstep[2:Nstep],levelsd,type='l',col='blue',#ylim=c(-10,10),
          ylab='level s.d',xlab='DoY')
plot(Tstep[2:Nstep],stdlevel,type='l',col='red',#ylim=c(0,10),
     ylab='Std Level',xlab='DoY')
grid()

# Density plots
dens.Xdlm = density(X.dlm,bw='SJ',window="epanechnikov",n=512,na.rm='T')
dens.X0 = density(X0,bw='SJ',window="epanechnikov",n=512,na.rm='T')
dens.slev = density(stdlevel,bw='SJ',window="epanechnikov",n=512,na.rm='T')

windows(width=10,height=6)
par(mfrow=c(1,2),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(dens.Xdlm$x,dens.Xdlm$y,type='l',lwd=2,col='forestgreen',xlab='Input series',
     ylab='density')
plot(dens.slev$x,dens.slev$y,type='l',lwd=2,col='forestgreen',xlab='Standardized Level',
     ylab='density')

save(Tstep,X0,X.dlm,level,levelsd,stdlevel,Yyhat,B.ests,B.sd,errvar,nl,delta,file=Fname)
print(Fname,quote=F)


