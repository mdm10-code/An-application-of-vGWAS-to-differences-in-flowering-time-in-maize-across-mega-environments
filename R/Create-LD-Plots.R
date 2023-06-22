###################################################################
################## LD visualization script #######################
##################################################################


#Remove all objects (i.e., data sets) from R
rm(list = ls())

#load in libraries
library(lattice)
library(Hmisc)

#Remove all objects (i.e., data sets) from R 
rm(list = ls())

#################################################################################################################################################################
##First, let me set the working directory to visuals
setwd(wd$visuals)

##Let me order everything by positions
pre.plot.BFT.GDDanthesis.narrowed <- pre.plot.BFT.GDDanthesis.narrowed[order(pre.plot.BFT.GDDanthesis.narrowed$Position1), ]
pre.plot.DGLM.GDDanthesis.narrowed <- pre.plot.DGLM.GDDanthesis.narrowed[order(pre.plot.DGLM.GDDanthesis.narrowed$Position1), ]
pre.plot.BFT.GDDsilk.narrowed <- pre.plot.BFT.GDDsilk.narrowed[order(pre.plot.BFT.GDDsilk.narrowed$Position1), ]
pre.plot.DGLM.GDDsilk.narrowed <- pre.plot.DGLM.GDDsilk.narrowed[order(pre.plot.DGLM.GDDsilk.narrowed$Position1), ]


#Identify what the maximum of the -log10 palue was for the graph using the results without the covariate. 
max.log.y.first.graph <- max(y)

#before I can run this script, I will need to convert position as an integer
pre.plot.BFT.GDDanthesis.narrowed$Position1 <- as.integer(pre.plot.BFT.GDDanthesis.narrowed$Position1)

colnames(pre.plot.BFT.GDDanthesis.narrowed)[2] <- "R.2"

#Create a vector of symbols to differentiate which symbol corresponds to 4K, 55K, and GBS, SNPs
#GBS SNPs
pch.vector <- rep(17, nrow(pre.plot.BFT.GDDanthesis.narrowed))
col.vector <- rep("Navyblue", nrow(pre.plot.BFT.GDDanthesis.narrowed))
cex.vector <- rep(1.0, nrow(pre.plot.BFT.GDDanthesis.narrowed))

#Identify the most statistically significant SNP as purple
col.vector[grep(154942763, pre.plot.BFT.GDDanthesis.narrowed$Position1)] = "Orange"
pch.vector[grep(154942763, pre.plot.BFT.GDDanthesis.narrowed$Position1)] = 25
cex.vector[grep(154942763, pre.plot.BFT.GDDanthesis.narrowed$Position1)] = 1.5

#Here, I will attach the pre plot info to the graph display
attach(pre.plot.BFT.GDDanthesis.narrowed)

min <- (min(pre.plot.BFT.GDDanthesis.narrowed[, 1]) - 1000) / 10^6
max <- (max(pre.plot.BFT.GDDanthesis.narrowed[, 1]) + 1000) / 10^6
x <- (pre.plot.BFT.GDDanthesis.narrowed[, 1] / (10^6))
y <- -log10(pre.plot.BFT.GDDanthesis.narrowed[, 4])

#Identify what the maximum of the -log10 palue was for the graph using the results without the covariate. 
max.log.y.first.graph <- max(y)

#Distance from significant SNP to end position of ca5p7 - CCAAT-HAP5-transcription factor 57
114344506 - 113424601
#Distance from signficant SNP to end position of UNUSUAL FLORAL ORGANS
114344506 - 114157054

#Now, I will create an increment for the window
increment <- (max - min) / 4

ticks.x.axis <- c(min, (min + increment), (min + (2 * increment)), (min + (3 * increment)), max)

setwd(wd$visuals)
pdf("LD_plot_Chr_9_BFT_difference_in_GDD_to_Anthesis.pdf", width = 10)
  par(mar = c(5, 5, 5, 5))
  plot(y ~ x, col = "darkgrey", type = "h", cex.lab = 1.5, xlim = c(min, max), ylim = c(0, max(y)), axes = FALSE, xlab = "", ylab = "", lwd = 1.5)
  axis(1, at = ticks.x.axis, labels = round(ticks.x.axis, 2), cex.axis = 1.5, tick = F)
  mtext(expression(Position ~ (RefGen_v4) ~ x ~ 10^6), side = 1, line = 3, cex = 1.8)
  axis(2, cex.axis = 1.5, tick = F)
  mtext(expression(-log[10](italic(p))), side = 2, line = 3, cex = 2.0)
  abline(h = min.most.signficant, col = "black", lty = 2)
  par(new = T)
  plot(pre.plot.BFT.GDDanthesis.narrowed$R.2 ~ x, pch = 17, col = col.vector, type = "p", axes = F, cex = cex.vector, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(4, at = c(0, 0.25, 0.50, 0.75, 1), cex.axis = 1.5, labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), tick = F)
  mtext(expression(r^2), side = 4, line = 3, cex = 2.0)
  abline(v = c(154923673 / 10^6, 154975225 / 10^6), col = "Purple", lty = 1)
  abline(h = 6.79, col = "black", lty = 2)
  box()
dev.off()

####Now, I will do this with UvBFT for difference in GDD to Anthesis
#before I can run this script, I will need to convert position as an integer
pre.plot.DGLM.GDDanthesis.narrowed$Position1 <- as.integer(pre.plot.DGLM.GDDanthesis.narrowed$Position1)

colnames(pre.plot.DGLM.GDDanthesis.narrowed)[2] <- "R.2"

#Create a vector of symbols to differentiate which symbol corresponds to 4K, 55K, and GBS, SNPs
#GBS SNPs
pch.vectorb <- rep(17, nrow(pre.plot.DGLM.GDDanthesis.narrowed))
col.vectorb <- rep("Navyblue", nrow(pre.plot.DGLM.GDDanthesis.narrowed))
cex.vectorb <- rep(1.0, nrow(pre.plot.DGLM.GDDanthesis.narrowed))

#Identify the most statistically significant SNP as purple
col.vectorb[grep(154942763, pre.plot.DGLM.GDDanthesis.narrowed$Position1)] = "Orange"
pch.vectorb[grep(154942763, pre.plot.DGLM.GDDanthesis.narrowed$Position1)] = 25
cex.vectorb[grep(154942763, pre.plot.DGLM.GDDanthesis.narrowed$Position1)] = 1.5

#Here, I will attach the pre plot info to the graph display
attach(pre.plot.DGLM.GDDanthesis.narrowed)

minb <- (min(pre.plot.DGLM.GDDanthesis.narrowed[, 1]) - 1000) / 10^6
maxb <- (max(pre.plot.DGLM.GDDanthesis.narrowed[, 1]) + 1000) / 10^6
xb <- (pre.plot.DGLM.GDDanthesis.narrowed[, 1] / (10^6))
yb <- -log10(pre.plot.DGLM.GDDanthesis.narrowed[, 4])

#Identify what the maximum of the -log10 palue was for the graph using the results without the covariate. 
max.log.y.first.graphb <- max(yb)

#Now, I will create an increment for the window
incrementb <- (maxb - minb) / 4

ticks.x.axisb <- c(minb, (minb + incrementb), (minb + (2 * incrementb)), (minb + (3 * incrementb)), maxb)

setwd(wd$visuals)
pdf("LD_plot_Chr_9_DGLM_difference_in_GDD_to_Anthesis.pdf", width = 10)
  par(mar = c(5, 5, 5, 5))
  plot(yb ~ xb, col = "darkgrey", type = "h", cex.lab = 1.5, xlim = c(minb, maxb), ylim = c(0, max(yb)), axes = FALSE, xlab = "", ylab = "", lwd = 1.5)
  axis(1, at = ticks.x.axisb, labels = round(ticks.x.axisb, 2), cex.axis = 1.5, tick = F)
  mtext(expression(Position ~ (RefGen_v4) ~ xb ~ 10^6), side = 1, line = 3, cex = 1.8)
  axis(2, cex.axis = 1.5, tick = F)
  mtext(expression(-log[10](italic(p))), side = 2, line = 3, cex = 2.0)
  abline(h = min.most.signficant, col = "black", lty = 2)
  par(new = T)
  plot(pre.plot.DGLM.GDDanthesis.narrowed$R.2 ~ xb, pch = 17, col = col.vectorb, type = "p", axes = F, cex = cex.vectorb, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(4, at = c(0, 0.25, 0.50, 0.75, 1), cex.axis = 1.5, labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), tick = F)
  mtext(expression(r^2), side = 4, line = 3, cex = 2.0)
  abline(v = c(154923673 / 10^6, 154975225 / 10^6), col = "Purple", lty = 1)
  abline(h = 6.79, col = "black", lty = 2)
  box()
dev.off()

###
pre.plot.DGLM.GDDsilk.narrowed$Position1 <- as.integer(pre.plot.DGLM.GDDsilk.narrowed$Position1)

colnames(pre.plot.DGLM.GDDsilk.narrowed)[2] <- "R.2"

#Create a vector of symbols to differentiate which symbol corresponds to 4K, 55K, and GBS, SNPs
#GBS SNPs
pch.vectorc <- rep(17, nrow(pre.plot.DGLM.GDDsilk.narrowed))
col.vectorc <- rep("Navyblue", nrow(pre.plot.DGLM.GDDsilk.narrowed))
cex.vectorc <- rep(1.0, nrow(pre.plot.DGLM.GDDsilk.narrowed))

#Identify the most statistically significant SNP as purple
col.vectorc[grep(154942763, pre.plot.DGLM.GDDsilk.narrowed$Position1)] = "Orange"
pch.vectorc[grep(154942763, pre.plot.DGLM.GDDsilk.narrowed$Position1)] = 25
cex.vectorc[grep(154942763, pre.plot.DGLM.GDDsilk.narrowed$Position1)] = 1.5

#Here, I will attach the pre plot info to the graph display
attach(pre.plot.DGLM.GDDsilk.narrowed)

minc <- (min(pre.plot.DGLM.GDDsilk.narrowed[, 1]) - 1000) / 10^6
maxc <- (max(pre.plot.DGLM.GDDsilk.narrowed[, 1]) + 1000) / 10^6
xc <- (pre.plot.DGLM.GDDsilk.narrowed[, 1] / (10^6))
yc <- -log10(pre.plot.DGLM.GDDsilk.narrowed[, 4])

#Identify what the maximum of the -log10 palue was for the graph using the results without the covariate. 
max.log.y.first.graphc <- max(yc)

#Now, I will create an increment for the window
incrementc <- (maxc - minc) / 4

ticks.x.axisc <- c(minc, (minc + incrementc), (minc + (2 * incrementc)), (minc + (3 * incrementc)), maxc)

setwd(wd$visuals)
pdf("LD_plot_Chr_9_DGLM_difference_in_GDD_to_Silk.pdf", width = 10)
  par(mar = c(5, 5, 5, 5))
  plot(yc ~ xc, col = "darkgrey", type = "h", cex.lab = 1.5, xlim = c(minc, maxc), ylim = c(0, max(yc)), axes = FALSE, xlab = "", ylab = "", lwd = 1.5)
  axis(1, at = ticks.x.axisc, labels = round(ticks.x.axisc, 2), cex.axis = 1.5, tick = F)
  mtext(expression(Position ~ (RefGen_v4) ~ xc ~ 10^6), side = 1, line = 3, cex = 1.8)
  axis(2, cex.axis = 1.5, tick = F)
  mtext(expression(-log[10](italic(p))), side = 2, line = 3, cex = 2.0)
  abline(h = min.most.signficant, col = "black", lty = 2)
  par(new = T)
  plot(pre.plot.DGLM.GDDsilk.narrowed$R.2 ~ xc, pch = 17, col = col.vectorc, type = "p", axes = F, cex = cex.vectorc, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(4, at = c(0, 0.25, 0.50, 0.75, 1), cex.axis = 1.5, labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), tick = F)
  mtext(expression(r^2), side = 4, line = 3, cex = 2.0)
  abline(v = c(154923673 / 10^6, 154975225 / 10^6), col = "Purple", lty = 1)
  abline(h = 6.79, col = "black", lty = 2)
  box()
dev.off()

####
pre.plot.BFT.GDDsilk.narrowed$Position1 <- as.integer(pre.plot.BFT.GDDsilk.narrowed$Position1)

colnames(pre.plot.BFT.GDDsilk.narrowed)[2] <- "R.2"

#Create a vector of symbols to differentiate which symbol corresponds to 4K, 55K, and GBS, SNPs
#GBS SNPs
pch.vectord <- rep(17, nrow(pre.plot.BFT.GDDsilk.narrowed))
col.vectord <- rep("Navyblue", nrow(pre.plot.BFT.GDDsilk.narrowed))
cex.vectord <- rep(1.0, nrow(pre.plot.BFT.GDDsilk.narrowed))

#Identify the most statistically significant SNP as purple
col.vectord[grep(154942763, pre.plot.BFT.GDDsilk.narrowed$Position1)] = "Orange"
pch.vectord[grep(154942763, pre.plot.BFT.GDDsilk.narrowed$Position1)] = 25
cex.vectord[grep(154942763, pre.plot.BFT.GDDsilk.narrowed$Position1)] = 1.5

#Here, I will attach the pre plot info to the graph display
attach(pre.plot.BFT.GDDsilk.narrowed)

mind <- (min(pre.plot.BFT.GDDsilk.narrowed[, 1]) - 1000) / 10^6
maxd <- (max(pre.plot.BFT.GDDsilk.narrowed[, 1]) + 1000) / 10^6
xd <- (pre.plot.BFT.GDDsilk.narrowed[, 1] / (10^6))
yd <- -log10(pre.plot.BFT.GDDsilk.narrowed[, 4])

#Identify what the maximum of the -log10 palue was for the graph using the results without the covariate. 
max.log.y.first.graphd <- max(yd)

#Now, I will create an increment for the window
incrementd <- (maxd - mind) / 4

ticks.x.axisd <- c(mind, (mind + incrementd), (mind + (2 * incrementd)), (mind + (3 * incrementd)), maxd)

setwd(wd$visuals)
pdf("LD_plot_Chr_9_BFT_difference_in_GDD_to_Silk.pdf", width = 10)
  par(mar = c(5, 5, 5, 5))
  plot(yd ~ xd, col = "darkgrey", type = "h", cex.lab = 1.5, xlim = c(mind, maxd), ylim = c(0, max(yd)), axes = FALSE, xlab = "", ylab = "", lwd = 1.5)
  axis(1, at = ticks.x.axisd, labels = round(ticks.x.axisd, 2), cex.axis = 1.5, tick = F)
  mtext(expression(Position ~ (RefGen_v4) ~ xd ~ 10^6), side = 1, line = 3, cex = 1.8)
  axis(2, cex.axis = 1.5, tick = F)
  mtext(expression(-log[10](italic(p))), side = 2, line = 3, cex = 2.0)
  abline(h = min.most.signficant, col = "black", lty = 2)
  par(new = T)
  plot(pre.plot.BFT.GDDsilk.narrowed$R.2 ~ xd, pch = 17, col = col.vectord, type = "p", axes = F, cex = cex.vectord, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  axis(4, at = c(0, 0.25, 0.50, 0.75, 1), cex.axis = 1.5, labels = c("0.00", "0.25", "0.50", "0.75", "1.00"), tick = F)
  mtext(expression(r^2), side = 4, line = 3, cex = 2.0)
  abline(v = c(154923673 / 10^6, 154975225 / 10^6), col = "Purple", lty = 1)
  abline(h = 6.79, col = "black", lty = 2)
  box()
dev.off()

########################### END ###################################