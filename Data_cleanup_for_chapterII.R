############################################################################
########### BLUE/BLUP Script ##########
###########################################################################

##Load workspace
load("chaptII.RData")

#Here, I am going to upload my R packages as according to Efficient Programming in R

#Now we can install the readxl R package via pkgs

#Here, I am going to add an already uploaded library
pkgs <- c("data.table", "tidyverse", "psych", "dplyr", "lme4", "stats", "dplyr", "tibble")

inst <- lapply(pkgs, library, character.only = T) #load them

##Now, I will create my working directory shortcuts
wd <- list()

wd$chaptII <- "~/projects/dissertation/ChaptII-RI/"
wd$data <- paste0(wd$chaptII, "0_data/")
wd$datmanual <- paste0(wd$data, "Manual/")
wd$out <- paste0(wd$chaptII, "2_pipeline/out/")

#Let's see what files I need from the manual data folder
list.files(wd$datmanual)

############### Section 3: Load required data ###############
#For this script, I need to load two things: the genotypic and phenotypic
#datasets for the maize US Buckler Goodman Diversity populations. To do this, I 
#am going to use readr since it's much quicker

#I will read both sheets from the "2006-2009_Maize_Data-MMURPHY-7-12-22.xlsx" file
phenotypic.data <- fread(file = paste0(wd$datmanual, "2006-2009_Maize_Data-MMURPHY-7-12-22.csv"))

#I will need to convert this to a dataframe
phenotypic.data <- as.data.frame(phenotypic.data) #3-28-2023

#Now, I do the same for the genotypic dataset. This will also load the genotypic and map datasets
load("~/projects/dissertation/ChaptII/0_data/External/Genotypic_data/Maize/GBSv4282numericdata.Rdata")

#Now, I am going to do some quick quality control steps. I am going to use
#str to see if the variables are read-in correctly
str(phenotypic.data) 
str(genotypic.data)

#To see how many individuals are present for each trait,
#I am going to run the describe function from pysch
described.pheno <- describe(phenotypic.data[-c(1), -c(1)])
View(described.pheno)

################ Section 4: Data manipulation #################

#For this section, I am going to manipulate the genotypic and phenotypic
#datasets for downstream analyses. I will first work on the phenotypic dataset

#Here, I am going to extract info for only flowering time data
flowering.pheno <- phenotypic.data[, c(1, grep("^GDDDays*", colnames(phenotypic.data)))]
colnames(flowering.pheno) <- paste0(as.character(flowering.pheno[1, ]), "_", colnames(flowering.pheno))
colnames(flowering.pheno) <- gsub("\\..*", "", colnames(flowering.pheno))

#Now, let's see how many individuals does my dataset have
described.flowering <- describe(flowering.pheno[-c(1), -c(1)])

flowering.pheno2 <- flowering.pheno[-c(1), -c(1)]
rownames(flowering.pheno2) <- flowering.pheno[-c(1), c(1)]

flowering.pheno3 <- as.matrix(flowering.pheno2)
flowering.pheno3 <- matrix(as.numeric(flowering.pheno3), ncol = ncol(flowering.pheno3))

#As a quality control measure, I need to make sure that 
#everything lines up correctly. 
flowering.pheno3[1:5, 1:5]
flowering.pheno2[1:5, 1:5]
flowering.pheno[1:5, 1:5]

#Now, with confidence, I can add everything back in
rownames(flowering.pheno3) <- rownames(flowering.pheno2)
colnames(flowering.pheno3) <- colnames(flowering.pheno2)

#For this next step, I am going to see which environments have more than 50% missing data or do not have a corresponding environment
described.BG <- describe(flowering.pheno3[1:286, ])

#I am going to remove columns 10, 11, 21
flowering.pheno4 <- flowering.pheno3[, -c(10, 11, 21)]

#Now, let's create a TAXA columns
flowering.pheno4 <- rownames_to_column(as.data.frame(flowering.pheno4), var = "TAXA")

#Then, I will create a trait vector
trait.vec <- c("GDDDaystoSilk", "GDDDaystoTassel")

pheno.list <- vector(mode = "list", length = length(trait.vec))
names(pheno.list)[[1]] <- trait.vec[1]
names(pheno.list)[[2]] <- trait.vec[2]

#now, I am going to create a for loop that iterates through flowering.pheno 4
for (i in seq_along(trait.vec)) {
  print(names(pheno.list)[[i]])
  trait.cols <- grep(trait.vec[i], colnames(flowering.pheno4))
  trait.data <- flowering.pheno4[, c(1, trait.cols)] #This includes the TAXA column
  colnames(trait.data) <- gsub("\\_.*", "", colnames(trait.data)) #This is to keep the environmental codes
  trait.data.long <- pivot_longer(as.data.frame(trait.data), cols = 2:ncol(trait.data))
  colnames(trait.data.long)[2:3] <- c("Env", trait.vec[i])
  pheno.list[[i]] <- trait.data.long
}



################ Section 6: End of Session ##################
#Before we end this session, let's do some project quality control.

#First, I am going to save the workspace for this script

save.image("C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Three/Chapter3.RData")

#Before I push to my github repository, I will lint this script
lint("C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Three/Initial/Chapter-3/1_code/Perform-data-read-in-and-manipulations-Script.R")

############################ END ###########################