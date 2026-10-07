#Plots of Effective Potential Curves
#Variable are phycocyanin, chlorphyll, pH
#Years are 2015, 2019, 2024
#Results are from B-spline program 

rm(list = ls())
graphics.off()
#Load data
#Each Rdata file is the effective potential curve by lake, by year, by variablehas
#Lake: R = Peter, Year = 15 or 19 or 24, Variable = log Chl, log phycocyanin (BGA), pH


load(file='EP_R15_LogChl.Rdata')
load(file='EP_R19_LogChl.Rdata')
load(file='EP_R24_LogChl.Rdata')
load(file='EP_R15_LogBGA.Rdata')
load(file='EP_R19_LogBGA.Rdata')
load(file='EP_R24_LogBGA.Rdata')
load(file='EP_R15_pH.Rdata')
load(file='EP_R19_pH.Rdata')
load(file='EP_R24_pH.Rdata')

#Plot 3 variables (columns) by 3 years (rows)
windows(width=10,height=10)
par(mfrow=c(3,3),mar=c(4.5, 5, 2, 1.5) 
    + 0.1,cex.axis=1.5)

plot(x_R15_LogChl,epl_R15_LogChl,type='l',lwd=3, col='forestgreen',
     yaxt='n', xlab='', ylab='')
axis(2, las = 1, at = c(0, 1, 2, 3))
mtext("2015", side=3, line=0.2, cex=1.6)
mtext("(a)", side=3, line=0.2, at=0.32, cex=1.2)

plot(x_R19_LogChl,epl_R19_LogChl,type='l', lwd=3, col='forestgreen',
     yaxt='n', xlab='', ylab='')
axis(2, las = 1, at = c(0, 1, 2, 3))
mtext (expression(paste("Log Chlorophyll (", mu,"g L"^"-1",")")),
        side=1, line=3.4, cex=1.5)
mtext ("2019", side=3, line=0.2, cex=1.6)
mtext("(b)", side=3, line=0.2, at=0.34, cex=1.2)

plot(x_R24_LogChl,epl_R24_LogChl,type='l', lwd=3, col='forestgreen',
     yaxt='n', xlab ='', ylab= '', ylim=c(0,16))
axis(2, las=1, at = c(5, 10, 15))
mtext("2024", side=3, line=0.2, cex=1.6)
mtext("(c)", side=3, line=0.2, at=0.42, cex=1.2)

plot(x_R15_LogBGA,epl_R15_LogBGA,type='l',lwd=3, col='turquoise3',
     yaxt='n', xlab='', ylab= '', ylim=c(0,11))
axis(2, las=1, at=c(3, 6, 9))
mtext("(d)", side=3, line=0.2, at=2.3, cex=1.2)
 
plot(x_R19_LogBGA,epl_R19_LogBGA,type='l', lwd=3, col='turquoise3',
     yaxt= 'n', xlab ='', ylab='', ylim=c(0,5))
axis(2, las=1, at=c(1, 3, 5))
mtext("(e)", side=3, line=0.2, at=2.5, cex=1.2)
mtext(expression("Log Cyanobacteria (cells ml"^{-1}*")"),
      side=1, line=3.8, cex=1.5)

plot(x_R24_LogBGA,epl_R24_LogBGA,type='l', lwd=3, col='turquoise3',
     yaxt= 'n', xlab='', ylab='', ylim = c(0,9))
axis(2, las=1, at=c(3, 6, 9))
mtext("(f)", side=3, line=0.2, at=2.3, cex=1.2)

plot(x_R15_pH,epl_R15_pH,type='l',lwd=3, col='brown',
     yaxt= 'n',
     xlab= '',ylab='')
axis (2, las=1, at=c(2, 4, 6))
mtext("(g)", side=3,line=0.2, at=6.4, cex=1.2)

plot(x_R19_pH,epl_R19_pH,type='l', lwd=3, col='brown',
     yaxt='n',
     xlab= '',ylab='')
axis (2, las=1, at=c(2, 6, 10))
mtext("(h)", side=3,line=0.2, at=6.6, cex=1.2)
mtext("pH", side=1,line=3, cex=1.5)

plot(x_R24_pH,epl_R24_pH,type='l', lwd=3, col='brown',
     yaxt= 'n',
     xlab='',ylab='', ylim=c(0,9.5))
axis (2, las=1, at=c(3, 6, 9))
mtext("(i)", side=3, line=0.2, at=6.5, cex=1.2)

mtext("Effective Potential", side=2, outer=TRUE, line=-2, cex=1.6)
