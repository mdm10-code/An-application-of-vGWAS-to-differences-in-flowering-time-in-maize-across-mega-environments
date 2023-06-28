#####################################################################
######################### GAPIT for Chapter 2 #######################
#####################################################################

#First, we will load our R workspace
load(file = file.choose())

#First, we need to upload GAPIT
source("http://zzlab.net/GAPIT/GAPIT.library.R") 
source("http://zzlab.net/GAPIT/gapit_functions.txt")

#To make reading in the Hapmap file much quicker, let's using the readr
#R package

# Installing
install.packages("readr")

# Loading
pkgs <- c("readr", "dplyr")

inst <- lapply(pkgs, library, character.only = T) #load them

##In this section, I am going to load the required files
myG <- read_delim(file = file.choose(), col_names = F)
myG <- as.data.frame(myG)

##Now, read in kinship
myKI <- read_delim(file = file.choose(), col_names = T)
myKI <- as.data.frame(myKI)

##Here, I need to do some slight modications with some of the inbreds
myKI[1, 2] <- "4226"
myKI[1, 3] <- "4722"

myKI[, 2] <- as.numeric(myKI[,2])
myKI[, 3] <- as.numeric(myKI[,3])

###Loading the PCAs to control for pop structure
myPC <- read_delim(file = file.choose(), col_names = T)
myPC <- as.data.frame(myPC)

###I have to modify the colnames of myG to match phenotypes
###Now, I am going to modify the colnames of ames_geno (which contains the inbred names)
inbred_vect <- myG[1, 12:292]
#
modinbred_vect <- gsub(":.*", "", inbred_vect)
#
length(grep(T, blups.data$Taxa %in% modinbred_vect))
#
#Now, let me check where the falses are coming from
#grep(F, blups.data$taxa %in% modinbred_vect)

modinbred_vect[33]
modinbred_vect[159]
modinbred_vect[278]

modinbred_vect[33] <- "B2"
modinbred_vect[159] <- blups.data[160, 1]
modinbred_vect[278] <- blups.data[279, 1]

myG[1, 12:292] <- modinbred_vect

#Now, let's save the modified 
write_delim(myG, "BG282.Hapmap.mod.txt", col_names = T, delim = "\t")

blups.data <- as.data.frame(blups.data)

#working directory for pre-vGWAS
setwd("~/projects/dissertation/ChaptII-RI/2_pipeline/out/GAPITIII")

colnames(blups.data)[1] <- "Taxa"

#Now, time to run GAPIT
GAPIT(
  Y = blups.data, 
  G = myG,
  CV = myPC,
  KI = myKI,
  SNP.effect = "Add",
  Model.selection = T,
  model = "MLM"
)

###outputting

#first, let's read-in myG with colnames set to true
myG2 <- read_delim(file = file.choose(), col_names = T)
myG2 <- as.data.frame(myG2)

#change the colnames
myG[1, 12:292] <- modinbred_vect
colnames(myG2)[12:292] <- modinbred_vect

write_delim(myG2, "BG282.Hapmap.mod.txt", col_names = T, delim = " ")

write_delim(blups.data, "BLUPs.chaptII.txt", delim = "\t")

