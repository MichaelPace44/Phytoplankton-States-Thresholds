

rm(list = ls())
graphics.off()

data1 <- read.csv(file="DOsat_Equilibria.csv")
str(data1)
data2 <- read.csv(file="DOsat_EQ_Mean+SD.csv")
str(data2)

data2$lower <- data2$Mean - data2$SD
data2$upper <- data2$Mean + data2$SD
point_colors_nominal <- c("blue", "blue", "blue", "red", "green3",
                  "blue", "red", "green3", "blue", "red", "green3")
point_colors_mean <- c("blue", "blue","blue","red", "green3",
                       "blue", "red", "green3", "blue","red",
                       "green3")
arrow_colors <- c("blue","blue","blue","red","green3","blue","red","green3","blue","red","green3")

windows()
par(mar=c(5.1,5.0,4.1,4.0))
plot(data1$x,data1$Equil,xaxt='n',ylim=c(85,130),pch=16, cex=1.4,
     col=point_colors_nominal, xlab="", ylab="", cex.lab = 1.8, las=1)
mtext("Dissolved Oxygen (% saturation)", side=2, line=3.0, cex=1.4)

mtext(("2013  2014        2015             2019                2024"),
      side=1, line=0.8, adj=0, cex=1.2)
mtext("Nutrient Manipulation Year", side=1, line=2, cex=1.4)
points(data2$x,data2$Mean,pch=17, cex=1.4, col=point_colors_mean)
arrows(data2$x, data2$lower, data2$x, data2$upper, length=0.05, angle=90, code=3,
       col=arrow_colors)
abline (v=3, lty=3, col="gray")
abline (v=6, lty=3, col="gray")
abline (v=13, lty=3, col="gray")
abline(v=20, lty=3, col="gray")
legend("topleft",
       legend = c("LS Data", "LS Mean", "MU Data", "MU Mean", "US Data", "US Mean"),
       pch = c(16,17,16,17,16,17),
       col = c("blue","blue", "red", "red","green3","green3"),
       pt.cex=1.2, bty='n')                   


#Version 2 with only blue (stable) and red (unstable) points

point_colors_nominal <- c("blue", "blue", "blue", "red", "blue",
                          "blue", "red", "blue", "blue", "red", "blue")
point_colors_mean <- c("blue", "blue","blue","red", "blue",
                       "blue", "red", "blue", "blue","red",
                       "blue")
arrow_colors <- c("blue","blue","blue","red","blue","blue","red" , "blue","blue","red","blue")

windows()
par(mar=c(5.1,5.0,4.1,4.0))
plot(data1$x,data1$Equil,xaxt='n',ylim=c(85,130),pch=16, cex=1.4,
     col=point_colors_nominal, xlab="", ylab="", cex.lab = 1.8, las=1)
mtext("Dissolved Oxygen (% saturation)", side=2, line=3.0, cex=1.4)

mtext(("2013  2014        2015             2019                2024"),
      side=1, line=0.8, adj=0, cex=1.2)
mtext("Nutrient Manipulation Year", side=1, line=2, cex=1.4)
points(data2$x,data2$Mean,pch=17, cex=1.4, col=point_colors_mean)
arrows(data2$x, data2$lower, data2$x, data2$upper, length=0.05, angle=90, code=3,
       col=arrow_colors)
abline (v=3, lty=3, col="gray")
abline (v=6, lty=3, col="gray")
abline (v=13, lty=3, col="gray")
abline(v=20, lty=3, col="gray")
legend("topleft",
       legend = c("Stable Data", "Stable Bootstrap", "Unstable Data", "Unstable Bootstrap"),
       pch = c(16,17,16,17,16,17),
       col = c("blue","blue", "red", "red"),
       pt.cex=1.2, bty='n')                   

