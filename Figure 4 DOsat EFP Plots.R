#Plots of DOsat Effective Potential Curves  
#Results are from B-spline program 

rm(list = ls())
graphics.off()

#Load data
#Each Rdata file has x by lake and year code x_R13
#Each Rdata file has EFP by lake and year cod epl_R13
load(file='EP_R13.Rdata')
load(file='EP_R14.Rdata')
load(file='EP_R15.Rdata')
load(file='EP_R19.Rdata')
load(file='EP_R24.Rdata')
 

#Plot 4 years separately on one page
windows()
par(mfrow=c(2,2),mar=c(5, 5, 1.5, 1.5) 
    + 0.1, cex.axis=1.2,cex.lab=1.4)
plot(x_R13,epl_R13,type='l',lwd=3, col='blue',
     xlab= '',ylab='',
     yaxt='n',
     ylim=c(0,6))
axis (2, las=1, at=c(2, 4, 6))
mtext("(a)", side=3, line=0, at=90, cex=1.2)
text(92,5.5,"2013", cex=1.2)

plot(x_R14,epl_R14,type='l', lwd=3, col='blue',
     xlab= '',ylab='',
     yaxt='n',
     ylim=c(0,6))
axis (2, las=1, at=c(2, 4, 6))
mtext("(b)", side=3, line=0, at=90, cex=1.2)
text(95,5.5,"2014", cex=1.2)

#plot(x_R15,epl_R15,type='l',col='blue',
#     xlab= 'DO (%sat)',ylab='Effective Potential')
#text(90,2,"2015")

plot(x_R19,epl_R19,type='l', lwd=3, col='blue',
     xlab= '',ylab='',
     yaxt='n',
     ylim=c(0,6))
axis (2, las=1, at=c(2, 4, 6))
mtext("(c)", side=3, line=0, at=90, cex=1.2)
text(96,5.5,"2019", cex=1.2)

plot(x_R24,epl_R24,type='l',lwd=3, col='blue',
     xlab= '',ylab='',
     yaxt='n',
     ylim=c(0,5))
axis (2, las=1, at=c(1, 3, 5))
mtext("(d)", side=3, line=0, at=90, cex=1.2)
text(95,5.5,"2024", cex=1.2)

mtext("DO (% saturation)", side=1, outer=TRUE, line=-2, cex=1.6)
mtext("Effective Potential", side=2, outer=TRUE, line=-2, cex=1.6)
