#Program to plot and analyze results of bootstraps


rm(list = ls())
graphics.off()

library(stats)
library(cubature)
library(numDeriv)
library(splines)
library(npreg)


# load bootstrap result; file options 10 or 1000 cycles for Peter 2015 DOsat
# nboot is number of bootstrap cycles, DF1 is data frame of data to bootstrap,
#  DF1$x0 is the X vector for plots of Langevin bootstraps,
#  D1boot is D1 bootstraps, booteqD1 is number of equilibria of D1 bootstraps, 
#  sigboot is sigma bootstraps, EPXboot is X axes of bootstrapped EPFboot,
#  booteqepf is number of equilibria of EPF bootstraps.
#  Nominal estimates from data:  LSD1, LSsigma, EPX, EPF, xeq.D1, xeq.epf
#save(nboot,DF1,D1boot,booteqD1,sigboot,EPXboot,EPFboot,booteqepf,NEPF,
#     LSD1,LSsigma,xeq.D1,EPX,EPF,xeq.epf,file=Fname)

#load(file='Boot_10_Spline_DOsat_Peter2015.Rdata')
load(file='Boot_1000_Spline_DOsat_Peter2015.Rdata')

# Analysis & plots from the bootstrap program

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

windows()
hist(booteqepf, main='', col=c("blue"), xaxt='n', xlim= c(0,11), ylim=c(0,400),
     xlab="", ylab="", las=1)
abline(h=0)
mtext ("1       2      3      4      5      6      7      8     9     10      11", 
       side=1, line=0, adj=0.65, cex=1.2)
mtext("Equilibria (Number)", side=1, line=2, cex=1.4)
mtext("Frequency", side=2, line=2.7, cex=1.4)

D1range = range(c(LSD1,D1boot),na.rm=T)

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5) 
plot(DF1$x0,LSD1,type='l',ylim=D1range,lwd=3,col='darkblue',xlab='X0',ylab='D1')
for(i in 1:nboot) {
  points(DF1$x0,D1boot[,i],type='l',lwd=1,col='skyblue2')
}
points(DF1$x0,LSD1,type='l',ylim=D1range,lwd=3,col='darkblue')
abline(h=0,lty=3,lwd=3,col='darkred')
grid()

sigrange = range(LSsigma,sigboot)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5) 
plot(DF1$x0,LSsigma,type='l',ylim=sigrange,lwd=3,col='darkred',xlab='X0',ylab='sigma')
for(i in 1:nboot) {
  points(DF1$x0,sigboot[,i],type='l',lwd=1,col='hotpink')
}
points(DF1$x0,LSsigma,type='l',ylim=sigrange,lwd=3,col='darkred')
grid()

EPFrange = range(EPF,EPFboot,rm=T)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(EPX[1:NEPF],EPF,ylim=EPFrange,log='y',type='l',lwd=3,col='darkgreen',xlab='X',ylab='EPF')
for(i in 1:nboot)  {
  points(EPXboot[1:NEPF,i],EPFboot[,i],type='l',lwd=1,col='limegreen')
}
points(EPX[1:NEPF],EPF,ylim=EPFrange,log='y',type='l',lwd=3,col='darkgreen')
grid()

# Percentile plots
dimEPF = dim(EPFboot)
nx = dimEPF[1]
pvec = c(0.1,0.25,0.5,0.75,0.9)
np = length(pvec)
EPFpct = matrix(0,nr=nx,nc=np)
for(i in 1:nx)  {
  qvals = quantile(EPFboot[i,],probs=pvec,na.rm=T)
  EPFpct[i,] = qvals[1:5]
}

EPFrange = range(EPF,EPFpct,rm=T)
windows()
par(mfrow=c(1,1),mar=c(4.5, 5, 2, 2) + 0.1,cex.axis=1.3,cex.lab=1.5)
# EPF of data
plot(EPX[1:NEPF],EPF,ylim=EPFrange,log='y',type='l',lwd=3,col='darkblue',
     xlab='Dissolved Oxygen (% sat)',ylab='',las=1) # data
mtext("Effective Potential", side=2, line=3.4, cex=1.5)
points(EPX[1:NEPF],EPFpct[,1],type='l',lwd=3,col='lightseagreen') # 10%
points(EPX[1:NEPF],EPFpct[,5],type='l',lwd=3,col='lightseagreen') # 90%
points(EPX[1:NEPF],EPFpct[,2],type='l',lwd=3,col='skyblue') # 25%
points(EPX[1:NEPF],EPFpct[,4],type='l',lwd=3,col='skyblue') # 75%
points(EPX[1:NEPF],EPFpct[,3],type='l',lwd=3,col='blue') # 50%
#points(EPX[1:NEPF],EPF,ylim=EPFrange,log='y',type='l',lwd=3,col='darkblue')  # replot data if needed
#grid()
legend('topleft',legend=c('Data','50%','25% & 75%','10% & 90%'),lwd=c(3,3,3,3),
       col=c('darkblue','blue','skyblue','lightseagreen'),cex=1,bty='n')
       #title='EPF of Data and Bootstrap Percentiles')

# find equilibria of percentiles 25, 50, 75
epfeq=function(epx,epf,nknots)  {
  # convert epf to a function
  xEPqt = quantile(epx,probs=seq(0,1,length.out=nknots),na.rm=T) # quantiles of X used for knots
  EPspline = ss(x=epx,y=epf,method='ML',m=2,knots=xEPqt)
  EPfun = function(x)  {
    yhat = predict(EPspline,x)$y
    return(yhat)
  }
  # take first derivative and find the roots
  dEPdx = grad(EPfun,epx,'Richardson')
  # Potential is a negative integral; reverse sign of derivative for plots
  # find roots of first derivative of effective potential
  sdrift = sign(dEPdx)
  dsdrift = c(0,-diff(sdrift))
  xeq = epx[which(!dsdrift == 0)]
  return(xeq)
}

print('',quote=F)
print('equilibria from original data',quote=F)
print(c('deterministic, D1:  ',round(xeq.D1,4)),quote=F)
print(c('stochastic, EPF:  ',round(xeq.epf,4)),quote=F)
print('',quote=F)
print('equilibria of percentiles 25, 50, 75',quote=F)
print('percentile 25',quote=F)
xeq25 = epfeq(EPX[1:NEPF],EPFpct[,2],5)  
print(xeq25)
print('percentile 50',quote=F)
xeq50 = epfeq(EPX[1:NEPF],EPFpct[,3],9)
print(xeq50)
print('percentile 75',quote=F)
xeq75 = epfeq(EPX[1:NEPF],EPFpct[,4],9)
print(xeq75)
print('',quote=F)
print('=============================================================',quote=F)
print('Bootstrap bias estimate for equilibria at the median',quote=F)
print(round( (xeq.epf - xeq50),3))  # note E&T reverse the order as xeq50 - xeq.EPF
print('Relative bootstrap bias for median equilibria',quote=F)
RBB = (xeq.epf - xeq50)/xeq.epf
print(round(RBB,3))
print('',quote=F)
print('bias-corrected equilibria and errors',quote=F)
BB = (xeq50 - xeq.epf)  # bootstrap bias by E&T text between eqs 10.40 & 10.41 p. 138
xeq.epf.adj = 2*xeq.epf - xeq50 # identical to xeq.epf - BB eq. 10.40
xeq25.adj = xeq25 - BB
xeq75.adj = xeq75 - BB
print(c('lower quartile',round(xeq25.adj,3)),quote=F)
print(c('median adjusted',round(xeq.epf.adj,3)),quote=F)
print(c('upper quartile',round(xeq75.adj,3)),quote=F)
print('',quote=F)
print('bias-corrected equilibria and errors using RBB',quote=F)
print('caution, this fixup may not be correct',quote=F)
RBB = abs(BB/xeq.epf) 
xeq.epf.adj = xeq.epf - (RBB*xeq.epf) 
xeq25.adj = xeq25 - (RBB*xeq25)
xeq75.adj = xeq75 - (RBB*xeq75)
print(c('lower quartile',round(xeq25.adj,3)),quote=F)
print(c('median adjusted',round(xeq.epf.adj,3)),quote=F)
print(c('upper quartile',round(xeq75.adj,3)),quote=F)

# Density plots for cases of 3 equilibria
print('',quote=F)
print('Density plots for cases of 3 equilibria',quote=F)
tstart = Sys.time() # start clock --------------------------------------------
print(freq.table[4])
N3 = unname(freq.table[4])
print(c('number of bootstrap EPF with 3 eq',N3),quote=F)
# if there were 3 EPF equilibria then find them
EQ3mat = matrix(0,nr=nboot,nc=3)  # matrix to hold cases of 3 EPF equilibria
Ncase = 0  # counter for cases of 3 equilibria
i3eq = c(0)  # save row numbers with 3 equilibia
for(i in 1:nboot) {
  if(booteqepf[i] != 3) {next}  # if there are not 3 equilibria then skip
  else{
    # epfeq=function(epx,epf,nknots)
    xeqboot = epfeq(EPXboot[(1:NEPF),i],EPFboot[,i],7) 
    #EPX = EPout[[1]]
    #EPF = EPout[[2]]
    #xeq.epf = EPout[[4]]
    #NEPF = length(EPF)
  }
  if(length(xeqboot) != 3){next}
  EQ3mat[i,] = xeqboot
  Ncase = Ncase + 1
  i3eq = c(i3eq,i)
  # print(c('row ',i,' done, total ',Ncase),quote=F) 
}

# -----------------------------------------------------------------------------
tstop = Sys.time()
print('----------------------------------------------------------',quote=F)
runtime = difftime(tstop,tstart,units='mins')
print(c('runtime for EPF equilibria, minutes ',runtime),quote=F)
print(c('confirmed cases of 3 stochastic (EPF) equilibria ',Ncase),quote=F)
print('----------------------------------------------------------',quote=F)
print('',quote=F)

# Select cases for density plots, then make the plots
EQ3case = EQ3mat[i3eq,]  # save only the cases with 3 verified equilibria

eqden1 = density(EQ3case[,1],bw='SJ',window="epanechnikov",n=512,na.rm='T')
eqden2 = density(EQ3case[,2],bw='SJ',window="epanechnikov",n=512,na.rm='T')
eqden3 = density(EQ3case[,3],bw='SJ',window="epanechnikov",n=512,na.rm='T')

windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(eqden1$x,eqden1$y,type='l',lwd=2,col='blue',xlab='Low Stable Eq',ylab='Density')
abline(v=xeq.epf[1],lty=3,lwd=3,col='black')
plot(eqden2$x,eqden2$y,type='l',lwd=2,col='red',xlab='Unstable Eq',ylab='Density')
abline(v=xeq.epf[2],lty=3,lwd=3,col='black')
plot(eqden3$x,eqden3$y,type='l',lwd=2,col='blue',xlab='High Stable Eq',ylab='Density')
abline(v=xeq.epf[3],lty=3,lwd=3,col='black')

print('Direct calculation of mean & sd for cycles w 3 equilibria',quote=F)
print(colMeans(EQ3case))
print(apply(EQ3case,2,sd))
print('means do not match xeq.epf',quote=F)
print('use deviation of means from nominal eq to correct bias via densities',quote=F)

# Bias-correct densities to match nominal equilibria  
# If x = equilibrium axis and y(x) = density then
# mean(x) = sum(y(x)*x*dx)/sum(y(x)*dx) 
# Note dx cancels in the ratio.  The densities have constant dx.
# Delta = nominal - mean(x)
# x_corrected = x + delta
#
# sum(x) = sum(y(x)*x*dx) because y(x) is proportional to N at each x
# sum(x^2) = sum(y(x)*x*x*dx)
# SS(x) = sum(x^2) - sum(x)^2/N; N is sum(y(x)dx) when y is density
# SS(x) = sum(y(x)*x*x*dx) - sum(y(x)*x*dx)^2/sum((y(x)dx))
# var_population(x) = SS(x)/sum(y(x)*dx) # Note that dx cancels in this step
# var_sample(x) = SS(x)/(sum(y(x)*dx) - 1)
# std dev is sqrt(var) for population or sample
# 
# Note that mean & sd of the densities will not be identical to mean and sd 
#  of the data because the density curve is smoothed with constant dx. 
#  The data are not smoothed and dx is not constant.  Nonetheless the approximation  
#  is close enough for visualization, the point of plotting the densities.

print('',quote=F)
print('++++++++++++++++++++++++++++++++++++++++++++++++++++',quote=F)
print('RESULTS FOR BOOTSTRAPPED EQUILIBRIA',quote=F)

# test variance and s.d. calculations
# Low equilibrium
sx1 = sum(eqden1$y * eqden1$x,na.rm=T)
sy1 = sum(eqden1$y,na.rm=T)  # N for density data
s2x1 = sum(eqden1$y * eqden1$x * eqden1$x,na.rm=T)
SSx1 = s2x1 - (sx1^2)/sy1
varsam1 = SSx1/(sy1 - 1)
sdsam1 = sqrt(varsam1)
print('',quote=F)
print('Nominal equilibria from data',quote=F)
print(xeq.epf)
print('',quote=F)
print('sample variance & standard deviation for low stable equilibrium',quote=F)
print(c(varsam1,sdsam1))
# Unstable equilibrium
sx2 = sum(eqden2$y * eqden2$x,na.rm=T)
s2x2 = sum(eqden2$y * eqden2$x * eqden2$x,na.rm=T)
sy2 = sum(eqden2$y,na.rm=T)
SSx2 = s2x2 - (sx2^2)/sy2
varsam2 = SSx2/(sy2 - 1)
sdsam2 = sqrt(varsam2)
print('',quote=F)
print('sample variance & standard deviation for unstable equilibrium',quote=F)
print(c(varsam2,sdsam2))
#
sx3 = sum(eqden3$y * eqden3$x,na.rm=T)
s2x3 = sum(eqden3$y * eqden3$x * eqden3$x,na.rm=T)
sy3 = sum(eqden3$y,na.rm=T)
SSx3 = s2x3 - sx3^2/sy3
varsam3 = SSx3/(sy3 - 1)
sdsam3 = sqrt(varsam3)
print('',quote=F)
print('sample variance & standard deviation for high stable equilibrium',quote=F)
print(c(varsam3,sdsam3))

# shift distributions to match means
xbar1 = sum(eqden1$y * eqden1$x)/sum(eqden1$y)
xbar2 = sum(eqden2$y * eqden2$x)/sum(eqden2$y)
xbar3 = sum(eqden3$y * eqden3$x)/sum(eqden3$y)

del1 = xeq.epf[1] - xbar1
del2 = xeq.epf[2] - xbar2
del3 = xeq.epf[3] - xbar3

print('',quote=F)
print('nominal equilibria, means of distributions, delta',quote=F)
print(xeq.epf)
print(c(xbar1,xbar2,xbar3))
print(c(del1,del2,del3))
print('sd of distributions',quote=F)
print(c(sdsam1,sdsam2,sdsam3))
#
x1adj = eqden1$x + del1
x2adj = eqden2$x + del2
x3adj = eqden3$x + del3

windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(5, 5, 3, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(x1adj,eqden1$y,type='l',lwd=2,col='blue',xlab='Dissolved Oxygen (% sat)',ylab='Density')
abline(v=xeq.epf[1],lty=3,lwd=3,col='black')
mtext("(a) Low Stable Equilibrium", side=3, line=0.2, adj=0.5, cex=1.1)
plot(x2adj,eqden2$y,type='l',lwd=2,col='red',xlab='Dissolved Oxygen (% sat)',ylab='Density')
abline(v=xeq.epf[2],lty=3,lwd=3,col='black')
mtext("(b) Middle Unstable Equilibrium", side=3, line=0.2, adj=0.5, cex=1.1)
plot(eqden3$x,eqden3$y,type='l',lwd=2,col='blue',xlab='Dissolved Oxygen (% sat)',ylab='Density')
abline(v=xeq.epf[3],lty=3,lwd=3,col='black')
mtext("(c) High Stable Equilibrium", side=3, line=0.2, adj=0.5, cex=1.1)

#Option to save results of 3 equilibria cases
#save(EQ3mat,EQ3case,i3eq,Ncase,eqden1,eqden2,eqden3,xeq.epf,
#     file='Densities_3eq_Peter2015.Rdata')
