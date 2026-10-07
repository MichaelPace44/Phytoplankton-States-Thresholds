#Plots of DOsat for year 2013, 14, 19, 24
#Figure 3
#Data are daily medians

rm(list = ls())
graphics.off()

#Load data
load(file='Daily_R_DOsat2013.Rdata')
dat1 <- DOx
load(file='Daily_R_DOsat2014.Rdata')
dat2 <- DOx
load(file='Daily_R_DOsat2015.Rdata')
dat3 <- DOx
load(file='Daily_R_DOsat2019.Rdata')
dat4 <- DOx
load(file='Daily_R_DOsat2024.Rdata')
dat5 <- DOx

load(file='Daily_L_DOsat2013.Rdata')
L_dat1 <- DOx
load(file='Daily_L_DOsat2014.Rdata')
L_dat2 <- DOx
load(file='Daily_L_DOsat2015.Rdata')
L_dat3 <- DOx
load(file='Daily_L_DOsat2019.Rdata')
L_dat4 <- DOx
load(file='Daily_L_DOsat2024.Rdata')
L_dat5 <- DOx



#Plot 4 years separately on one page
#Do not include 2015 because already plotted in Figure 1
windows()
par(mfrow=c(2,2),cex.axis=1.2,cex.lab=1.4)

plot(dat1$doy,dat1$DOmed,type='l',col='blue',
     xlab= 'Day of Year',ylab='DO (%sat)',ylim=c(70,130),
     las=1)
points(L_dat1$doy,L_dat1$DOmed, type='l', col='lightslategray')
text(180,125,"2013",cex=1.2)
abline(v=c(154,238),col="black",lty="dashed")

plot(dat2$doy,dat2$DOmed,type='l',col='blue',
     xlab= 'Day of Year',ylab='DO (%sat)',
     ylim=c(70,130),
     las=1)
points(L_dat2$doy,L_dat2$DOmed,type='l',col='lightslategray')
text(180,125,"2014",cex=1.2)
abline(v= c(153,240) ,col="black",lty="dashed")

plot(dat4$doy,dat4$DOmed,type='l',col='blue',
xlab= 'Day of Year',ylab='DO (%sat)',ylim=c(70,130),
las=1)
points(L_dat4$doy,L_dat4$DOmed,type='l',col='lightslategray')
text(180,125,"2019",cex=1.2)
abline(v=c(161,237),col="black",lty="dashed")

plot(dat5$doy,dat5$DOmed,type='l',col='blue',
     xlab= 'Day of Year',ylab='DO (%sat)',ylim=c(70,130),
     las=1)
points(L_dat5$doy,L_dat5$DOmed,type='l',col='lightslategray')
text(180,125,"2024",cex=1.2)
abline(v= c(162, 233), col="black",lty="dashed")

#Plot 4 years separately on one page
#Same plot but with one x-axis and one y-axis
#Do not include 2015 because already plotted in Figure 1
windows()
par(mfrow=c(2,2),cex.axis=1.2,cex.lab=1.4, 
    mar=c(4.5, 5, 2, 1.5) + 0.1)

plot(dat1$doy,dat1$DOmed,type='l',col='blue',lwd=3,
     xlab= '',ylab='',ylim=c(70,130),
     las=1)
points(L_dat1$doy,L_dat1$DOmed, type='l', col='lightslategray',lwd=3)
#text(180,125,"2013",cex=1.2)
mtext("(a) 2013", side=3, line=0, at=140, cex=1.2)
abline(v=c(154,238),col="black",lty="dashed")
legend("top", legend= c("Peter Lake", "Paul Lake"), lwd=2,
       col=c("blue","lightslategray"))

plot(dat2$doy,dat2$DOmed,type='l',col='blue', lwd=3,
     xlab= '',ylab='',ylim=c(70,130),
     las=1)
points(L_dat2$doy,L_dat2$DOmed,type='l',col='lightslategray',lwd=3)
#text(180,125,"2014",cex=1.2)
mtext("(b) 2014", side=3, line=0, at=140, cex=1.2)
abline(v= c(153,240) ,col="black",lty="dashed")

plot(dat4$doy,dat4$DOmed,type='l',col='blue', lwd=3,
     xlab= '',ylab='',ylim=c(70,130),
     las=1)
points(L_dat4$doy,L_dat4$DOmed,type='l',col='lightslategray',lwd=3)
#text(180,125,"2019",cex=1.2)
mtext("(c) 2019", side=3, line=0, at=140, cex=1.2)
abline(v=c(161,237),col="black",lty="dashed")

plot(dat5$doy,dat5$DOmed,type='l',col='blue', lwd=3,
     xlab= '',ylab='',ylim=c(70,130),
     las=1)
points(L_dat5$doy,L_dat5$DOmed,type='l',col='lightslategray',lwd=3)
#text(180,125,"2024",cex=1.2)
mtext("(d) 2024", side=3, line=0, at=140, cex=1.2)
abline(v= c(162, 233), col="black",lty="dashed")

mtext("Day of Year", side = 1, outer = TRUE, line = -2, cex = 1.6)
mtext("DO (% saturation)", side = 2, outer = TRUE, line = -2, cex = 1.6)
     