####################################################################
####################    Dot Chart Script   #########################
####################################################################

#Let's create a dotchart for the variance componenets
####Dot chart instead

silk.var.comp <- silk.var[["fullRandVC"]]
tassel.var.comp <- tassel.var[["fullRandVC"]]

#Before I begin plotting, I need to change a few things
rownames(silk.var.comp) <- c("mega-environment", "mega-environment:environment interaction", "genotype", "genotype:mega-environment interaction", "residuals")
rownames(tassel.var.comp) <- c("mega-environment", "mega-environment:environment interaction", "genotype", "genotype:mega-environment interaction", "residuals")

# Dot chart of a single numeric vector
pdf("Variance-components-for-GDDaystoSilking.pdf", width = 10)
dotchart(x = silk.var.comp$vcovPerc, labels = row.names(silk.var.comp),
         cex = 0.70, xlab = "Variance", main = "Variance Components for GDDaystoSilk")
dev.off()


#Dotchart for GDDaystoAnthesis
pdf("Variance-components-for-GDDaystoAnthesis.pdf", width = 10)
dotchart(x = tassel.var.comp$vcovPerc, labels = row.names(tassel.var.comp),
         cex = 0.70, xlab = "Variance", main = "Variance Components for GDDaystoAnthesis")
dev.off()
