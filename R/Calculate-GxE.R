###################################################
### GxE calculations ###
###################################################

#### Section #1: Load required libraries and data

#load dataset
load(file.choose())

#set working directory

#Next, we load the required libraries
pkgs.for.GxE <- c("statgenGxE", "lme4")

inst <- lapply(pkgs.for.GxE, library, character.only = T) #load them

#### Section 2: Mega-environment determination #####
mega.GDDsilk$year <- NA
mega.GDDtassel$year <- NA

sixers <- grep("*6", mega.GDDsilk$Env)

seveners <- grep("*7", mega.GDDsilk$Env)

#populating for 2006
mega.GDDsilk[c(sixers), 6] <- "2006"
mega.GDDsilk[c(seveners), 6] <- "2007"

#
sixers <- grep("*6", mega.GDDtassel$Env)

seveners <- grep("*7", mega.GDDtassel$Env)

#populating for 2006
mega.GDDtassel[c(sixers), 6] <- "2006"
mega.GDDtassel[c(seveners), 6] <- "2007"

colnames(mega.GDDsilk)[4] <- "region"
colnames(mega.GDDtassel)[4] <- "region"


#In this section, I am going to determine the number the mega-environments
TD.GDDDaysSilk.woNAM <- statgenSTA::createTD(data = mega.GDDsilk, genotype = "TAXA", trial = "Env")
TD.GDDDaysTassel.woNAM <- statgenSTA::createTD(data = mega.GDDtassel, genotype = "TAXA", trial = "Env")

#variance components
silk.var <- gxeVarComp(TD = TD.GDDDaysSilk.woNAM, trait = "GDDDaystoSilk", nestingFactor = "region")
tassel.var <- gxeVarComp(TD = TD.GDDDaysTassel.woNAM, trait = "GDDDaystoTassel", nestingFactor = "region")

write.csv(silk.var[["fullRandVC"]], "var.comp.GDDSilk.csv", sep = ",")
write.csv(tassel.var[["fullRandVC"]], "var.comp.GDDtassel.csv", sep = ",")

###############################################
save.image("chaptII.RData")
save.image("chaptII_copy.RData")

