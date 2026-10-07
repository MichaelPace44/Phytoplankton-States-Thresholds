#Bootstrap Analysis
#Part one executes the B-spline analysis 
#Be sure to use same combination of knots and trimming as original
#Part two is the bootstrap; results are saved under Fname

rm(list = ls())
graphics.off()

source('EPFunction+EQ_ssML_nKnot.R')  # use ss() from npreg


library(stats)
library(cubature)
library(numDeriv)
library(splines)
library(npreg)
library(tictoc)

#Load data
load(file="DLM_DOsat_median_Peter15.Rdata")
#Set file name for bootstrap output; number equal nboot, number of cycles
Fname = c('Boot_10_Spline_DOsat_Peter2015.Rdata')

#Part One B-spline*************************************************************
# SET OPTIONS FOR THE ANALYSIS =========================================3=============
# CHOOSE VARIATE, SET NUMBERS OF KNOTS AND SPLINE ORDER, TRIM OUTLIERS IF NEEDED

# select data 
Xvar0 = X0
nx = length(Xvar0)
Tstep0 = Tstep[1:nx]

# number of knots for mu (deterministic core) and sigma (sqrt(conditional variance))
nk.mu = 5
nk.sig = 7
# order of polynomial splin5e (usually 3 for cubic spline; 2 may dampen fluctuations)
npoly.mu = 3
npoly.sig = 3
# Should data be trimmed off to remove outliers?
Trim = 0  # 0 selects original data, 1 selects trimmed data

# Number of bootstrap cycles
nboot = 10 #Use a low number like 10 to test and 1000 for final version

# END OF OPT1ONS ===================================================================

windows()
plot(Tstep0,Xvar0,type='l',lwd=1,col='seagreen',main='original data')

# find optimal AR lags using auto.arima from forecast library 
aropt=1  # enter aropt from diagnostics script; use 1 given prior trimming

# DT
DT = aropt  # daily data

# subsample Xvar0 according to Markov lag
ikeep = seq(1,nx,by=aropt)
Xvar = Xvar0[ikeep]
Tstep = Tstep0[ikeep]
nx = length(Xvar)

windows()
plot(Tstep,Xvar,type='l',lwd=1,col='blue',main='Markov-thinned data')

tstart = Sys.time() # start clock --------------------------------------------

# ======================================================================

# Make data set with first and second moments

#Use these lines if when ARIMA autocorrelation is 0 or 1
x0 = Xvar0[1:(nx-1)]
x1 = Xvar0[2:nx]

#Use these lines when ARIMA autocorrelation is 2 
#x0 = Xvar0[1:(nx-2)]
#x1 = Xvar0[3:nx]

dx = (x1-x0)/DT   # first moment
dx2 = 0.5*(x1-x0)^2/DT  # second moment
dsd = sqrt(2*dx2)  # sigma

# Make data frame and sort by x0
DF0 = as.data.frame(cbind(x0,x1,dx,dx2,dsd))
DF0s = DF0[order(DF0$x0),]  # sort by x0

windows()
par(mfrow=c(2,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(DF0s$x0,DF0s$dx,type='b',pch=20,col='blue')
plot(DF0s$x0,DF0s$dx2,type='b',pch=20,col='red')

# check extremes and trim if needed
#
xtr = quantile(X0,probs=c(0.01,0.02,0.03,0.04,0.96,0.97,0.98,0.99))
#xtr = quantile(X0,probs=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08))
print(xtr)

print('tails of x0',quote=F)
print(range(DF0s$x0))
#

# If data is to be trimmed choose a statement for degree of trimming 
#DF1 = subset(DF0s,subset=(x0 >= xtr[1] & x0 <= xtr[8])) # remove 1st and 99th percentiles
DF1 = subset(DF0s,subset=(x0 >= xtr[4] & x0 <= xtr[5])) # remove through 4th and above 96th percentiles
#DF1 = subset(DF0s,subset=(x0 >= xtr[4])) #remove lower percentiles
#DF1 = subset(DF0s,subset=(x0 >= xtr[8])) # remove bottom 8 percentiles
#DF1 = subset(DF0s, subset=(x0 <= xtr[5]))
#DF1 = subset(DF0s, subset=(x0 <= xtr[5] & x0 > 88)) #remove top 4 percentiles and all values below )

# compare number of data points in original vs trimmed data
print(c('dimension of original data',dim(DF0s)),quote=F)
print(c('dimension of trimmed data',dim(DF1)),quote=F)

# if Trim is 0 then analyze original data; otherwise analyze trimmed data
if(Trim == 0) {DF1 = DF0s}

windows()
par(mfrow=c(2,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(DF1$x0,DF1$dx,type='b',pch=20,col='blue')
plot(DF1$x0,DF1$dx2,type='b',pch=20,col='red')

# Quantiles for knots
xqt.mu = quantile(DF1$x0,probs=seq(0,1,length.out=nk.mu),na.rm=T)
print('quantiles',quote=F)
print(xqt.mu,quote=F)

# Internal knots based on quantiles
kmid.mu = xqt.mu[2:(nk.mu-1)]

# Fit D1 
# Defaults for bs() are order = 3, boundary knots are limits of data
basisD1 = bs(DF1$x0,knots=kmid.mu,degree=npoly.mu,intercept=F,Boundary.knots = range(DF1$x0))
lmD1 = lm(dx ~ basisD1 -1,data=DF1,x=T,y=T)  # least-squares fit of the B-spline
D1mat = unname(lmD1$x) # remove row and col names
D1y = unname(lmD1$y) # dx conforming to spline design

# Solve spline for weights with normal equations
D1D1 = t(D1mat)%*%D1mat
iD1D1 = solve(D1D1)
# as.vector strips the row names from the result
wD1 = as.vector(iD1D1%*%t(D1mat)%*%D1y) # weights for D1 spline
print('Solution complete for D1 spline',quote=F)

print('summary of linear spline fit',quote=F)
print(summary(lmD1))
print('dimension of D1mat',quote=F)
print(dim(D1mat))

# Find errors of dx model
lmD1err = lmD1$residuals  # from lm()
# We want to predict x1 = x0 + D1*DT which has error x1 - x0 - D1*DT
eP = DF1$x1 - DF1$x0 - (D1mat%*%wD1)*DT  # which should match lmD1err!
#windows()  # plot shows that errors are identical by lm() & normal equations
#plot(lmD1err,eP,type='p',pch='+')

# sigma vector from eP is sqrt(2*D2) 
DF1$ep2 = 0.5*(eP^2)
DF1$sep = sqrt(2*DF1$ep2)

windows()
par(mfrow=c(2,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(DF1$x0,DF1$dx,type='b',pch=20,col='blue')
plot(DF1$x0,DF1$ep2,type='b',pch=20,col='red')

# Build spline for sqrt(error^2), i.e. sqrt(D2 corrected for mean) for effective potential
# Quantiles
xqt.sig = quantile(DF1$x0,probs=seq(0,1,length.out=nk.sig),na.rm=T)
print('quantiles',quote=F)
print(xqt.sig,quote=F)

# Internal knots based on quantiles
kmid.sig = xqt.sig[2:(nk.sig-1)]

# Fit sigma 
# Defaults for bs() are order = 3, boundary knots are limits of data
basisig = bs(DF1$x0,knots=kmid.sig,degree=npoly.sig,intercept=F,Boundary.knots = range(DF1$x0))
lmsig = lm(sep ~ basisig -1,data=DF1,x=T,y=T)
sigmat = unname(lmsig$x) # remove row and col names
print('dimensions of design matrix for log(sigma)',quote=F)
print(dim(sigmat))
sigy = unname(lmsig$y) # sigma conforming to spline design
# Solve spline of sigma for weights using normal equations
sig2 = t(sigmat)%*%sigmat
isig2 = solve(sig2)
# as.vector strips the row names from the result
wsig = as.vector(isig2%*%t(sigmat)%*%sigy) # weights for sigma spline
print('Solution complete for sigma spline',quote=F)

# plot LS estimates
LSD1 = D1mat%*%wD1
LSsigma = sigmat%*%wsig 

windows(width=10,height=10)
par(mfrow=c(2,2),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(DF1$x0,LSD1,type='l',lwd=2,col='blue',
     xlab='x0',ylab='D1 by LS')
abline(h=0,lty=3,lwd=2)
plot(DF1$x0,LSsigma,type='l',lwd=2,col='red',
     xlab='x0',ylab='sigma by LS')

# equilibria from spline fit of D1
sdrift = sign(LSD1)
dsdrift = c(0,-diff(sdrift))
xeq.D1 = DF1$x0[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('',quote=F)
print('equilibria of D1 on x0 axis',quote=F)
print(xeq.D1,quote=F)
print(ixeq,quote=F)  

# compute and plot effective potential
nk.ss = max(nk.mu,nk.sig) #sets smooth spline knots

EPout = EPFEQ(DF1$x0,LSD1,LSsigma,nk.ss) 
EPX = EPout[[1]]
EPF = EPout[[2]]
dEPdx1 = EPout[[3]]
xeq.epf = EPout[[4]]
NEPF = length(EPF)

# Print equilibria from EPF
print('',quote=F)
print('Equilibria from fitted effective potential, least squares',quote=F)
print(xeq.epf,quote=F)
print('',quote=F)

# Plot EPF
#windows(width=8,height=4)
#par(mfrow=c(1,2),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(EPX[2:100],EPF,type='l',lwd=2,col='blue',xlab='X',ylab='Effective Potential',
     main='Least Squares')
grid()
# sign of derivative was corrected in EPfunction
plot(EPX,dEPdx1,type='l',lwd=2,col='blue',xlab='X',ylab='-d(EP)/dx')
abline(h=0,lty=3,lwd=2,col='black')
grid()

#Part 2 Bootstrap***************************************************************
# Bootstrap calculation ========================================================
# Remember that D2 is calculated from the residuals of D1
print(c('Number of boot cycles ',nboot),quote=F)
rboot = length(LSD1)
D1boot = matrix(0,nr=rboot,nc=nboot)  # save bootstraps of D1
sigboot = matrix(0,nr=rboot,nc=nboot)  # save bootstraps of sigma
booteqD1 = as.vector(rep(0,nboot))  # save bootstrapped number of eq of D1
EPXboot = matrix(0,nr=(NEPF+1),nc=nboot)
EPFboot = matrix(0,nr=NEPF,nc=nboot) # save bootstrapped EPF
booteqepf = as.vector(rep(0,nboot))  # save bootstrapped number of eq of EPF

# rationale:
# (1) The nominal D1 estimate LSD1 = D1mat%*%wD1
# (2) The spline for D1 predicts x1 = x0 + D1*DT which has error x1 - x0 - D1*DT
# (3) Therefore we calculate x1 pseudodata as
#       x1.pseudo =  DF1$x0 + (LSD1 + randomized(lmD1err))*DT
#           where the lmD1err values are randomized without replacement (all are used)
#       dx.pseudo = x1.pseudo - DF1$x0 
# (4) and fit 1 column of D1boot using DF1$x0 and dx.pseudo
# Then repeat for nboot cycles
# Because x0 stays the same the spline basis basisD1 
#   and the design matrix D1mat are unchanged and can be recycled
# 
nD1err = length(lmD1err)
for(i in 1:nboot) {      # start bootstrap cycle
  ranD1err = sample(lmD1err,size=nD1err,replace=T)
  ranx1 = DF1$x0 + (LSD1 + ranD1err)*DT
  randx = ranx1 - DF1$x0
  ranD1y = randx  
  ranwD1 = as.vector(iD1D1%*%t(D1mat)%*%ranD1y) # random weights for D1 spline
  D1boot[,i] = D1mat%*%ranwD1
  # find D1 equilibria
  sdrift = sign(D1boot[,i])
  dsdrift = c(0,-diff(sdrift))
  raneq = DF1$x0[which(!dsdrift == 0)]
  raneq.4 = round(raneq,4)
  #print(raneq)
  neq = length(raneq)
  booteqD1[i] = neq
  # find sigma equilibria
  eP = ranx1 - DF1$x0 - D1boot[,i]*DT  # which should match ranD1err!
  DF1$ep2 = 0.5*(eP^2)  # D2
  DF1$sep = sqrt(2*DF1$ep2)
  sigy = DF1$sep   # dependent variate for design matrix of sigma
  wsig = as.vector(isig2%*%t(sigmat)%*%sigy) # weights for sigma spline
  sigboot[,i] = sigmat%*%wsig 
  # attempt to integrate EPF with EPFEQ
  #EPFEQ = function(X,D1,sigma,nknots) 
  #outlist = list(X.ep,EPF,dEPdx,xeq,D12)
  epf.out = EPFEQ(DF1$x0,D1boot[,i],sigboot[,i],nk.ss)
  EPXboot[,i] = epf.out[[1]]   
  EPFboot[,i] = epf.out[[2]]
  eqepf = epf.out[[4]]
  eqepf.4 = round(eqepf,4)
  neq = length(eqepf)
  booteqepf[i] = neq
  print('',quote=F)
  print(c('complete cycles: ',i,' ----------------------------------------'),quote=F)
  print(c('D1 eq ',raneq.4),quote=F)
  print(c('EPF eq ',eqepf.4),quote=F)
}

# -----------------------------------------------------------------------------
tstop = Sys.time()
print('----------------------------------------------------------',quote=F)
runtime = difftime(tstop,tstart,units='mins')
print(c('runtime for bootstrapping, minutes ',runtime),quote=F)
print('----------------------------------------------------------',quote=F)
print('',quote=F)
print('frequencies of numbers of deterministic D1 equilibria',quote=F)
freq.table = table(booteqD1)
print('frequency table for number of equilibria',quote=F)
print('first row is number of equilibria ',quote=F)
print(c('second row is number of bootstrap cycles of total ',nboot),quote=F)
print(freq.table)
print(c('median number of equilibria = ',median(booteqD1)),quote=F)
print('frequencies of numbers of stochastic EPF equilibria',quote=F)
freq.table = table(booteqepf)
print('frequency table for number of equilibria',quote=F)
print('first row is number of equilibria ',quote=F)
print(c('second row is number of bootstrap cycles of total ',nboot),quote=F)
print(freq.table)
print(c('median number of equilibria = ',median(booteqepf)),quote=F)

D1range = range(c(LSD1,D1boot),na.rm=T)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5) 
plot(DF1$x0,LSD1,type='l',ylim=D1range,lwd=3,col='darkblue',xlab='X0',ylab='D1')
for(i in 1:nboot) {
  points(DF1$x0,D1boot[,i],type='l',lwd=1,col='skyblue2')
}
abline(h=0,lty=3,lwd=3,col='darkred')
grid()

sigrange = range(LSsigma,sigboot)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5) 
plot(DF1$x0,LSsigma,type='l',ylim=sigrange,lwd=3,col='darkred',xlab='X0',ylab='sigma')
for(i in 1:nboot) {
  points(DF1$x0,sigboot[,i],type='l',lwd=1,col='hotpink')
}
grid()

EPFrange = range(EPF,EPFboot,rm=T)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(EPX[1:NEPF],EPF,ylim=EPFrange,log='y',type='l',lwd=3,col='darkgreen',xlab='X',ylab='EPF')
for(i in 1:nboot)  {
  points(EPXboot[1:NEPF,i],EPFboot[,i],type='l',lwd=1,col='limegreen')
}
grid()


# nboot is number of bootstrap cycles, DF1 is data frame of data to bootstrap,
#  DF1$x0 is the X vector for plots of Langevin bootstraps,
#  D1boot is D1 bootstraps, booteqD1 is number of equilibria of D1 bootstraps, 
#  sigboot is sigma bootstraps, EPXboot is X axes of bootstrapped EPFboot,
#  booteqepf is number of equilibria of EPF bootstraps.
#  Nominal estimates from data:  LSD1, LSsigma, EPX, EPF, xeq.D1, xeq.epf
save(nboot,DF1,D1boot,booteqD1,sigboot,EPXboot,EPFboot,booteqepf,NEPF,
     LSD1,LSsigma,xeq.D1,EPX,EPF,xeq.epf,file=Fname)
print(c('results saved in ',Fname),quote=F)

