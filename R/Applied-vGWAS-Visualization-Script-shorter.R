###############################################################
############# Applied vGWAS Visualization Script #############
###############################################################


install.packages("CMplot")
library(CMplot)
library(dplyr)

#Next, we set the working directory to where the results files are located
setwd("~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Two/Revision-1/2_pipeline/out/vGWAS-results")

#Now, I will set a working directory to each folder
wd <- list()
wd$vGWAS.results <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Two/Revision-1/2_pipeline/out/vGWAS-results/"

list.dirs(wd$vGWAS.results)

#each mode
wd$BFT <- paste0(wd$vGWAS.results, "BFT/")
wd$DGLM <- paste0(wd$vGWAS.results, "DGLM/")
wd$DGLM.mQTNs <- paste0(wd$vGWAS.results, "DGLM/DGLM-results-w-mQTNs-as-fixed-effects/")

length(list.files(wd$BFT))
list.files(wd$DGLM)
list.files(wd$DGLM.mQTNs) #Still need south GDD.silk

#Now, for all the results, I am going to create separate lists
BFT.list <- vector(mode = "list", length = length(list.files(wd$BFT)))
DGLM.list <- vector(mode = "list", length = length(grep(T, grepl(".csv", list.files(wd$DGLM)))))
DGLM.nomQTNs.list <- vector(mode = "list", length = length(list.files(wd$DGLM.mQTNs)))


length(grep(T, grepl(".csv", list.files(wd$DGLM))))

##Read in results files
###For this next section, I am going to read-in all the results for each
#list from a for loop
BFT.vec <- list.files(wd$BFT)
DGLM.vec <- list.files(wd$DGLM)
DGLM.vec <- DGLM.vec[2:9]
DGLM.nomQTN.vec <- list.files(wd$DGLM.mQTNs)

feedme_results <- function(results.vec = NULL, wd = NULL, output.list = NULL) {
  for (file in seq_along(results.vec)) {
    print(file)
    r.file <- read_csv(file = paste0(wd, results.vec[file]), col_names = T)
    r.file <- as.data.frame(r.file)
    output.list[[file]] <- r.file
    names(output.list)[[file]] <- results.vec[file]
  }
  return(output.list)
}

##chatGPT suggestion
feedme_results <- function(results.vec = NULL, wd = NULL) {
  output.list <- lapply(results.vec, function(file) {
    r.file <- read_csv(file = paste0(wd, file), col_names = T)
    r.file <- as.data.frame(r.file)
    names(r.file) <- file
    return(list(file = r.file))
  }) 
}



###
BFT.list <- feedme_results(results.vec = BFT.vec, wd = wd$BFT, output.list = BFT.list)
DGLM.list <- feedme_results(results.vec = DGLM.vec, wd = wd$DGLM, output.list = DGLM.list)
DGLM.nomQTNs.list <- feedme_results(results.vec = DGLM.nomQTN.vec, wd = wd$DGLM.mQTNs, output.list = DGLM.nomQTNs.list)

#Now, I am going to modify the list names
names(BFT.list) <- gsub(".csv", "", names(BFT.list))
names(BFT.list) <- gsub("BFT_results", "", names(BFT.list))

names(DGLM.list) <- gsub(".csv", "", names(DGLM.list))
names(DGLM.list) <- gsub("DGLM_results", "", names(DGLM.list))

names(DGLM.nomQTNs.list) <- gsub(".csv", "", names(DGLM.nomQTNs.list))
names(DGLM.nomQTNs.list) <- gsub("DGLM_results", "", names(DGLM.nomQTNs.list))

########################### Section #: Manhattan Plots ########

#This is just for the results to present
BFT.list.sub <- vector(mode = "list", length = 2)
DGLM.list.sub <- vector(mode = "list", length = 2)

BFT.list.sub[[1]] <- BFT.list[[5]]
names(BFT.list.sub)[1] <- names(BFT.list[5])

BFT.list.sub[[2]] <- BFT.list[[6]]
names(BFT.list.sub)[2] <- names(BFT.list[6])

#Now, doing the same for DGLM
DGLM.list.sub[[1]] <- DGLM.list[[5]]
names(DGLM.list.sub)[1] <- names(DGLM.list[5])

DGLM.list.sub[[2]] <- DGLM.list[[6]]
names(DGLM.list.sub)[2] <- names(DGLM.list[6])

###Before I begin plotting, I wil lneed to create a switch for the visuals
wd$output <- "~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Two/Revision-1/3_output/"

#Now, I am going to create the appropriate directories
dir.create(paste0(wd$output, "visuals/"))

wd$visuals <- paste0(wd$output, "visuals/")

#Now create folders for manhattan plots and qqplots
dir.create(paste0(wd$visuals, "manhattan/"))
dir.create(paste0(wd$visuals, "qqplots/"))

#Now, create the switches
wd$manhattan <- paste0(wd$visuals, "manhattan/")
wd$qqplot <- paste0(wd$visuals, "qqplots/")

##Create a function for manhattan plots
setwd(wd$manhattan)

ondeck <- BFT.list.sub[[1]]
ondeck <- ondeck[, c(1:3, 8)]
ondeck <- na.omit(ondeck)
colnames(ondeck)[4] <- names(BFT.list.sub)[1]

create_manhattan <- function(wd = NULL, results.list = NULL, ymax = NULL, testname = NULL) {
  for (i in seq_along(results.list)) {
    setwd(wd)
    print(i)
    ondeck <- results.list[[i]]
    ondeck <- ondeck[, c(1:3, 8)]
    ondeck <- na.omit(ondeck)
    colnames(ondeck)[4] <- paste0(testname, "_", names(results.list)[i])
    CMplot(ondeck, ylim = c(0, ymax), col = c("darkblue", "darkorange"), plot.type = c("m"), multracks = F, threshold = c(0.05 / dim(ondeck)[1]), threshold.lty = c(1),
           threshold.lwd = c(1, 1), threshold.col = c("black", "grey"), amplify = T, bin.size = 1e6,
           chr.den.col = NULL, signal.col = c("magenta3"),
           signal.cex = 1, file = "jpg", memo = "", dpi = 300, file.output = T, verbose = T,
           highlight = NULL, highlight.text = NULL, highlight.text.cex = 1.4, LOG10 = T)
  }
}

#manhattan plots
create_manhattan(wd = wd$manhattan, results.list = BFT.list.sub, ymax = 16, testname = "BFT")
create_manhattan(wd = wd$manhattan, results.list = DGLM.list.sub, ymax = 16, testname = "DGLM")

#just mQTNs
create_manhattan(wd = wd$manhattan, results.list = DGLM.nomQTNs.list, ymax = 16, testname = "DGLM")


#For the requested analyses
create_manhattan(wd$manhattan, results.list = BFT.list, ymax = 18, testname = "BFT")
create_manhattan(wd$manhattan, results.list = DGLM.list, ymax = 18, testname = "DGLM")


###Now, I can am going to create a similar function for qqplots
create_qqplot <- function(wd = NULL, results.list = NULL, ymax = NULL, testname = NULL) {
  for (i in seq_along(results.list)) {
    setwd(wd)
    print(i)
    ondeck <- results.list[[i]]
    ondeck <- ondeck[, c(1:3, 8)]
    ondeck <- na.omit(ondeck)
    colnames(ondeck)[4] <- paste0(testname, "_", names(results.list)[i])
    CMplot(ondeck, ylim = c(0, ymax), col = c("darkblue"), plot.type = c("q"), multracks = F, threshold = c(0.05 / dim(ondeck)[1]), threshold.lty = c(1),
           threshold.lwd = c(1, 1), threshold.col = c("black", "grey"), amplify = T, bin.size = 1e6,
           chr.den.col = NULL, signal.col = c("darkorange3"),
           signal.cex = 1, file = "tiff", memo = "", dpi = 400, file.output = T, verbose = T,
           highlight = NULL, highlight.text = NULL, highlight.text.cex = 1.4, LOG10 = T)

  }
}

#
create_qqplot(wd = wd$qqplot, results.list = BFT.list.sub, ymax = 16, testname = "BFT")
create_qqplot(wd = wd$qqplot, results.list = DGLM.list.sub, ymax = 16, testname = "DGLM")

#mQTNs
create_qqplot(wd = wd$qqplot, results.list = DGLM.nomQTNs.list, ymax = 16, testname = "DGLM")

#rest of the traits
create_qqplot(wd = wd$qqplot, results.list = BFT.list, ymax = 18, testname = "BFT")
create_qqplot(wd = wd$qqplot, results.list = DGLM.list, ymax = 18, testname = "DGLM")

######################
save.image("ChaptII.visuals.RData")
