# Calculates Daily Median, 5%, 95% values, and 5th to 95th range
# of Dissolved Oxygen %Saturation
# Modify program to make similar calculations for pH and BGA (defined below)
# Data sets are organized by project
# Project 1 = years 2013, 2014, 2015; Hydrolab primary sonde
# Project 2 = year 2019; Hydrolab and YSI sonde used
# Project 3 = year 2024; YSI primary sonde
# Analyst must select:
#       project data file to load, 
#       the lake and sonde combination,
#       the year.
# Analyst must comment out the sonde code not used.
# Analyst must set output files to be saved with lake code, variable, year.

# Clear lists and graphics
rm(list = ls())
graphics.off()

library(dplyr)

# Function to calculate the daily statistics 
dayDOsat = function(dat1) {
  dat2 = dat1  # optionally, select part of dat1
  # option to select part of a day
  #dat1$dec = dat1$doy - floor(dat1$doy)
  #dat2 = subset(dat1,subset=(dec >= 0.375 & dec <= 0.625))  # keep 0900 - 1500
  # stats by day
  udoy = unique(trunc(dat2$DoY))
  ndoy = length(udoy)
  dat3 = matrix(0,nrow=ndoy,ncol=6) 
  for(i in 1:ndoy) {
    daykeep = subset(dat2,subset=(trunc(dat2$DoY) == udoy[i]))
    qDOsat = quantile(daykeep$DOsat_calc,probs=c(0.05,0.5,0.95))
    rDOsat = qDOsat[3] - qDOsat[1]  # daily range of DOsat
    dat3[i,]=c(daykeep$Year[1],udoy[i],qDOsat,rDOsat)
  }
  dat4 = as.data.frame(dat3)
  colnames(dat4) = c('year','doy','DO5','DOmed','DO95','DOrange')
  outlist = list(dat2,dat4)
  return(dat4) 
}

# LOAD DATA
# Data frames for Project 1 are Peter, Paul, and Tuesday lakes, years: 2013-2015
# Data frame for Project 2 is Peter and Paul lakes for 1 year (2019)
# Data frames for Project 3 are by Lake either Peter or Paul (select which)

# Columns are:  'Year', 
#               'Lake', as one letter code: L=Paul, R=Peter, T=Tuesday
#               'DoY', day of year
#               'BGA_HYLB', blue green algae, HYLB is sonde, other name is YSI
#               'BGA_logged_HYLB', log is base 10
#               'DOsat_calc_HYLB', dissolved oxygen percent saturation
#               'PH_HYLB', pH

#Data to calculate daily means from sonde data, see text for details
#on which sonde type (Hydro or YSI) was used
#Select one line to load data; comment out other lines
#For BGA data in 2024; Peter transformed to Hydrolab units using 2019 (see text)
#BGA data in 2024 for Paul not transformed; original YSI units

dat <- load(file="Project1_HYLB.Rdata") #Hydrolab all variables 2013-2015
#dat <- load(file="Project2_YSI.Rdata")
#dat <- load(file="Project2_HYLB.Rdata")
#dat <- load(file="Project3_YSI_Peter.Rdata")
#dat <- load(file="Project3_YSI_Paul.Rdata")

# Use commands below to assign uselake for Projects 1 and 2
# Select one command comment out others

uselake = Pete_HYLB # lake choice
#uselake = Pete_YSI # lake choice
#uselake - Paul_HYLB # lake choice
#uselake = Paul_YSI # lake choice

#Observations with NA are omitted
dat0 = na.omit(uselake)
print('dim before na.omit: ',quote=F)
print(dim(uselake))
print('dim after na.omit: ',quote=F)
print(dim(dat0))

# Choose Year depending on project loaded and year needed
# Project 1 has three years, 2013, 2014, 2015, Project 2 2019, Project 3 2024
dat0c = subset(dat0,subset=(Year==2015))

#Renane columns to remove "YSI" or "HYLB"
#Comment out lines for sonde not used
dat0c <- dat0c%>%
  rename(BGA = BGA_HYLB, BGA_logged = BGA_logged_HYLB, 
         DOsat_calc = DOsat_calc_HYLB, PH = PH_HYLB)
# rename(BGA = BGA_YSI, BGA_logged = BGA_logged_YSI, 
#        DOsat_calc = DOsat_calc_YSI, PH = PH_YSI)


# simple rate plot for DOsat
nDO = length(dat0c$DOsat_calc)
dDO = dat0c$DOsat_calc[2:nDO ] - dat0c$DOsat_calc[(1:(nDO-1))]
windows()
par(mfrow=c(1,1),mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(dat0c$DOsat_calc[(1:(nDO-1))],dDO,type='p',pch=20,cex=0.5,xlab='median',ylab='delta DOsat')

# Check QQ plot for DOsat
windows()
par(mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
qqnorm(dat0c$DOsat_calc,main='DOsat Normal QQ plot')
qqline(dat0c$DOsat_calc,lwd=2,col='blue')

# Get daily stats for DOsat using the function dayDOsat at start of program
DOx = dayDOsat(dat0c)   

# save daily DOsat; 
#col names are c('year','doy','DO5','DOmed','DO95','DOrange')
#adjust file name for lake, variable, and year

save(DOx,file='Daily_R_DOsat2015.Rdata')
write.csv(DOx,file='Daily_R_DOsat2015.csv')

# Plot stats for daily DOsat
windows()
par(mfrow=c(4,1),mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(DOx$doy,DOx$DO5,type='l',lwd=2,col='blue',xlab='Day of Year',
     ylab='5th pctile',main='DOsat')
plot(DOx$doy,DOx$DO95,type='l',lwd=2,col='red',xlab='Day of Year',
     ylab='95th pctile')
plot(DOx$doy,DOx$DOmed,type='l',lwd=2,col='purple',xlab='Day of Year',
     ylab='median')
plot(DOx$doy,DOx$DOrange,type='l',lwd=2,col='black',xlab='Day of Year',
     ylab='range')

# simple rate plot of daily medians for DOsat
DOmed = DOx$DOmed
nDO = length(DOmed)
dDO = DOmed[2:nDO ] - DOmed[(1:(nDO-1))]
windows()
par(mfrow=c(1,1),mar=c(2.5,4.5,1,2)+0.1,cex.axis=1.6,cex.lab=1.6)
plot(DOmed[(1:(nDO-1))],dDO,type='p',pch=20,xlab='median',ylab='delta DOsat')
abline(h=0,lty=3,lwd=2)

windows()
#plot of daily DO median versus DOY; nicer plot
plot(DOx$doy,DOx$DOmed, type="b", pch=16, col='darkorange3', 
     xlab= "Day of Year", cex.lab=1.4,
  ylab= "DO (Percent Saturation)", cex.lab=1.4)
#use abline to show fertilization period or comment out
#fertilization start end Day of Years were:
  #2013: 154 238
  #2014  153 240
  #2015  152 180
  #2019  161 237
  #2024  162 233
abline (v=152, col="black", lwd=2, lty=3)
abline (v=180, col="black", lwd=2, lty=3)


# Select data; make a series of transformations
# This examples uses median DOsat
X0 = DOx$DOmed
Tstep = DOx$doy
# try a few Box-Cox transforms
X0.recip = 1/X0
X0.recipsqrt = 1/sqrt(X0)
X0.log10 = log10(X0)
X0.sqrt = sqrt(X0)
#
nx = length(X0)

tforms = list(X0.recip,X0.recipsqrt,X0.log10,X0.sqrt,X0)
tnames = c('X0.recip','X0.recipsqrt','X0.log10','X0.sqrt','X0')
for(i in 1:5) {
  X0test = tforms[[i]]
  windows()
  par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
  qqnorm(X0test,pch=20,col='red',xlab='Quantiles, sd',ylab='DOsat',
         main=tnames[i])
  qqline(X0test,lwd=2,col='blue')
  grid()
}

