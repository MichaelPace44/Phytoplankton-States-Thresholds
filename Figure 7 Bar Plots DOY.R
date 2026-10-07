rm(list = ls())
graphics.off()

# Data compiled from individual runs of B-spline script to generate 
# estimate of the first day a threshold was crossed for a given year and varible
DOY2015 <- c(157, 171, 166, 168)
DOY2019 <- c(202, 192, 206, 203)
DOY2024 <- c(192, 197, 165, 193)
categories <- c("DO sat", "Cyano" , "Chl", "pH")
#categories <- c("","","","")

windows(width=6,height=10)
par(mfrow=c(3,1),mar=c(2.5, 5, 2.5, 1.5) 
    + 0.1,cex.axis=1.3,cex.lab=1.5)

#bar plot 2015 threshold crossing 
bp1 <- barplot(DOY2015, names.arg = categories, 
        col = c("blue","aquamarine2", "green1","brown"),
        main = "2015", 
        yaxt='n', 
        ylim = c(100,220),
#        xlab = "Variable", 
       ylab = "Day of Year",xpd=F)
axis(2, las=1, c(100, 140, 180, 220))
abline(h=100)
mtext("(a)", side=3, line=0.2, at=0.2, cex=1.1)
text(x=bp1, y=DOY2015, labels=DOY2015, pos=3, cex=1.2, col ="black")
abline(h=152, lty="dashed")

#legend("topright",legend = c("DO saturation", "Cyanobacteria",
#        "Chlorophyll", "pH"), fill = c("blue", "aquamarine2",
#        "green1","brown"))

#bar plot 2019 threshold crossing
bp2 <- barplot(DOY2019, names.arg = categories, 
        col = c("blue","aquamarine2", "green1","brown"), main = "2019",
        yaxt='n',
        ylim = c(100,220), 
#        xlab = "Variable", 
        ylab = "Day of Year",xpd=F)
axis(2, las=1, c(100, 140, 180, 220))
abline(h=100)
abline(h=161, lty="dashed")
text(x=bp2, y=DOY2019, labels=DOY2019, pos=3, cex=1.2, col ="black")
mtext("(b)", side=3, line=0.2, at=0.2, cex=1.1)
#mtext("Day of Year", side=2, line=1, cex=1.6)           
           
#bar plot 2024 threshold crossing

bp3 <-  barplot(DOY2024, names.arg = categories, 
        col = c("blue","aquamarine2", "green1","brown"), main = "2024",
        yaxt='n',
        ylim = c(100,220), 
#        xlab = "Variable",
        ylab = "Day of Year",xpd=F)
axis(2, las=1, c(100, 140, 180, 220))
abline(h=100)
abline (h=162, lty="dashed")
text(x=bp3, y=DOY2024, labels=DOY2024, pos=3, cex=1.2, col ="black")
mtext("(c)", side=3, line=0.2, at=0.2, cex=1.1)
#mtext("Day of Year", side=2, outer=TRUE, line=-2, cex=1.6)  