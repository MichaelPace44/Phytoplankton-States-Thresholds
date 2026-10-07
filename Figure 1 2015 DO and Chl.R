# Plot for manuscript Figure 1
# Chlorophyll and DO saturation in 2015 nutrient addition
# Daily time series

rm(list = ls())
graphics.off()

load(file="Daily_R_DOsat2015.Rdata")
load(file="Daily_R_Chl2015.Rdata")


#Drop columns Lake and Year from Chlorophyll data
#Add rows to chlorophyll data make same length as oxygen data
dat1 <- subset(dat0c, select = -c(Lake, Year))
new_rows <- data.frame(
  DOY = c(129:144),
  Manual_Chl = NA,
  Log_Chl = NA)
dat1 <- rbind(new_rows,dat1)
new_row <- data.frame(
  DOY = c(247), Manual_Chl=NA, Log_Chl=NA)
dat1 <- rbind(dat1,new_row)

#Change column name to DOY and merge two dataframes by DOY
colnames(DOx)[colnames(DOx) == "doy"] <- "DOY"
dat2 <- merge(DOx,dat1, by = "DOY")

#Plot oxygen and chlorophyll on same graph with different
#y-axes

dat3 <- load(file="Daily_L_DOsat2015.Rdata")

windows()
par(mar=c(5.1,5.0,4.1,4.0))
#First plot for DO data
plot(dat2$DOY,dat2$DOmed, type = 'b', pch=16, col='blue', lwd=2, cex=1.1,
     ylim=c(83,135),
     xlab = "Day of Year", ylab = '', cex.lab=1.4, las=1)
points(DOx$doy,DOx$DOmed, type ='b', pch=16, col='lightslategray', lwd=2, cex=1.2)
mtext("Dissolved Oxygen (% saturation)", side=2, line=3.0, cex=1.4)
#Indicate period of nutrient addition
abline(v=152, col='black', lty="dashed")
abline(v=170, col='black', lty="dashed")

# Allow a second plot on the same graph
par(new=TRUE)

# Plot the second plot and put axis scale on right
plot(dat2$DOY,dat2$Manual_Chl, type = 'b', pch=16, lwd=2,  
     col = 'forestgreen', cex=1.1, xlab="", ylab="",
     axes=FALSE)
#mtext("Chlorophyll (ug L-1)",side=4,line=2.5, cex=1.4) 
mtext(expression(paste("Chlorophyll ( ", mu,"g L"^"-1",")")), 
      side=4,line=2.8, cex=1.4)
axis(4, ylim=c(0,40),las=1)
legend(x = "topright",
       legend = c("Peter Lake DO", "Peter Lake Chlorophyll", "Paul Lake DO"),
       lty = c(1,1),
       col = c("blue","forestgreen","lightslategray"),
       lwd=2)