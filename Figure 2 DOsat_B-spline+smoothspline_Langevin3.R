
rm(list = ls())
graphics.off()

source('EPFunction+EQ_ssML_nKnot.R')  # use ss() from npreg
#source('EPFunction+EQ_simple+deriv_sm_spl_nK.R')  # use smooth.spline() from stats

library(stats)
library(cubature)
library(numDeriv)
library(splines)
library(npreg)


load(file="DLM_DOsat_median_Peter15.Rdata")
Fname = c('Bspline_DOsat_median_Peter15.Rdata')


# SET OPTIONS FOR THE ANALYSIS =========================================3=============
# CHOOSE VARIATE, SET NUMBERS OF KNOTS AND SPLINE ORDER, TRIM OUTLIERS IF NEEDED

# select data 
Xvar0 = X0
nx = length(Xvar0)
Tstep0 = Tstep[1:nx]

# number of knots for mu (deterministic core) and sigma (sqrt(conditional variance))
nk.mu = 5
nk.sig = 7
# order of polynomial spline (usually 3 for cubic spline; 2 may dampen fluctuations)
npoly.mu = 3
npoly.sig = 3
# Should data be trimmed off to remove outliers?
Trim = 0  # 0 selects original data, 1 selects trimmed data

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
#xtr = quantile(X0,probs=c(0.01,0.02,0.03,0.04,0.96,0.97,0.98,0.99))
xtr = quantile(X0,probs=c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08))
print(xtr)

print('tails of x0',quote=F)
print(range(DF0s$x0))
#

# If data is to be trimmed choose a statement for degree of trimming 
#DF1 = subset(DF0s,subset=(x0 >= xtr[1] & x0 <= xtr[8])) # remove 1st and 99th percentiles
#DF1 = subset(DF0s,subset=(x0 >= xtr[4] & x0 <= xtr[5])) # remove through 4th and above 96th percentiles
DF1 = subset(DF0s,subset=(x0 >= xtr[4])) 
#DF1 = subset(DF0s,subset=(x0 >= xtr[8])) # remove bottom 8 percentiles

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
xeq = DF1$x0[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('',quote=F)
print('equilibria of D1 on x0 axis',quote=F)
print(xeq,quote=F)
print(ixeq,quote=F)  

# compute and plot effective potential
nk = max(nk.mu,nk.sig)
run1 = EPFEQ(DF1$x0,LSD1,LSsigma,nk) # last entry is nKnots for EPFEQ function
#outlist = list(xvec.ep,EPF,dEPdx,xeq)
xvec.ep1 = run1[[1]]
epf1 = run1[[2]]
dEPdx1 = run1[[3]]
xeq1 = run1[[4]]

# Print equilibria from EPF
print('',quote=F)
print('Equilibria from fitted effective potential, least squares',quote=F)
print(xeq1,quote=F)
print('',quote=F)

# Plot EPF
#windows(width=8,height=4)
#par(mfrow=c(1,2),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(xvec.ep1[2:100],epf1,type='l',lwd=2,col='blue',xlab='X',ylab='Effective Potential',
     main='Least Squares')
grid()
# sign of derivative was corrected in EPfunction
plot(xvec.ep1,dEPdx1,type='l',lwd=2,col='blue',xlab='X',ylab='-d(EP)/dx')
abline(h=0,lty=3,lwd=2,col='black')
grid()

# save result
#save(DF1,nk.mu,nk.sig,npoly.mu,npoly.sig,Trim,LSD1,LSsigma,xeq,xeq1,file=Fname)
#print(c('output file name ',Fname),quote=F)

###Plot for Paper Figure 2
windows(width=10,height=10)
par(mfrow=c(2,2),mar=c(4.5, 5, 2, 1.5) 
    + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(DF1$x0,LSD1, type='l',lwd=3, col='blue',
     yaxt = 'n',
     xlab='DO (% sat)', ylab='')
axis(2, las = 1, at = c(-1.2, -0.8, -0.4, 0.0))
mtext("D1", cex=1.5, side=2, line=3.5)
mtext("(a)", side=3, line=0, at=90, cex=1.2)
abline(h=0,lty=3,lwd=2)

plot(DF1$x0,LSsigma,type='l',lwd=3,col='orange',
     yaxt = 'n',
     xlab='DO (% sat)', ylab='')
axis(2, las = 1)
mtext("(b)", side=3, line=0, at=90, cex=1.2)
mtext("D2", cex=1.5, side=2, line=3.5)

plot(xvec.ep1[2:100],epf1,type='l',lwd=3,col='red',
     yaxt = 'n',
     xlab='DO (% sat)',ylab='')
axis(2, las = 1)
mtext("(c)", side=3, line=0, at=90, cex=1.2)
mtext("Effective Potential", cex=1.5, side=2, line=3.5)

plot(xvec.ep1,dEPdx1,type='l',lwd=3,col='darkgreen',
     yaxt = 'n',
     xlab='DO (% sat)', ylab='')
axis(2, las = 1, at = c(-0.2, 0.0, 0.2))
mtext("(d)", side=3, line=0, at=90, cex=1.2)
mtext("-d(EP)/dx", cex=1.5, side=2, line=3.5)
abline(h=0,lty=3,lwd=2,col='black')