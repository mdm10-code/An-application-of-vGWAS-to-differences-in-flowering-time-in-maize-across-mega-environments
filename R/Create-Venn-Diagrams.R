####################################################################################
################### Create Venn Diagrams Script ####################################
####################################################################################

#Outline of this script
#1.) Load Workspace
#2.) Load required libraries
#3.) Read-in required results files
##a.) data.table R package
##b.) for loop into a list
#4.) Filter results files by MAF
##a.) Dplyr and/or tidyr
##b.) for loop through this
#5.) Further subset files for significant snps
##a.) Dplyr and/or tidyr
##b.) for loop through this
#6.) Get snps as characters from objects
##a.) input into a new list for VennDiagram
#7.) Feed into Venn Diagram
#8.) Get intersection from Venn-Diagram


#Section 1: Load workspace
load("chaptII.RData")

#As always, I am going to install any required R packages if not done so
install.packages("ggVennDiagram")

#Know, I will load the R package
library(data.table)
library(tidyr)
library(dplyr)
library("ggVennDiagram")
library(ggplot2)

#I am also going to set up my working directory shortcuts
wd <- list()
wd$vGWAS <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Two/Initial-Sub/2_pipeline/out/Applied-vGWAS/"

#Here, I am going to create a list
BFT.sigsnp.list <- vector(mode = "list", length = 8)
DGLM.sigsnp.list <- vector(mode = "list", length = 8)
DGLM.nomQTNS.sigsnp.list <- vector(mode = "list", length = 2)

#Let's also create another list for all the information with the significant SNPs
BFT.sigsnpallinfo.list <- vector(mode = "list", length = 8)
DGLM.sigsnpallinfo.list <- vector(mode = "list", length = 8)
DGLM.nomQTNs.sigsnpallinfo.list <- vector(mode = "list", length = 2)

Extract_sigsnps <- function(result.list = NULL, sigsnp.list = NULL, allinfo.list = NULL, test = NULL) {
  for (i in seq_along(sigsnp.list)) {
    print(i)
    names(sigsnp.list)[i] <- names(result.list)[i]
    names(allinfo.list)[i] <- names(result.list)[i]
    result_file <- result.list[[i]]
    result_file <- na.omit(result_file)
    if(test == "DGLM") {
      sigsnps <- result_file[result_file$P.disp <= (0.05 / 327056), ]
    } else {
      if(test == "BFT") {
        sigsnps <- result_file[result_file$p.value <= (0.05 / 327056), ]
      }
    }
    allinfo.list[[i]] <- sigsnps
    sigsnp.list[[i]] <- as.character(sigsnps$snp)
    print(sigsnp.list[[i]])
  }
  return(sigsnp.list)
}

DGLM.sigsnps <- Extract_sigsnps(result.list = DGLM.list, sigsnp.list = DGLM.sigsnp.list, allinfo.list = DGLM.sigsnpallinfo.list, test = "DGLM")
BFT.sigsnps <- Extract_sigsnps(result.list = BFT.list, sigsnp.list = BFT.sigsnp.list, allinfo.list = BFT.sigsnpallinfo.list, test = "BFT")

DGLM.mQTN.sigsnps <- Extract_sigsnps(result.list = DGLM.nomQTNs.list, sigsnp.list = DGLM.nomQTNS.sigsnp.list, allinfo.list = DGLM.nomQTNs.sigsnpallinfo.list, test = "DGLM")

#Now, I will Create venn diagrams specifically for each trait. For this, I will
#create new lists. Before I do this, I will rename the list names
names(sigsnp_list) <- c("UvDGLM Silk", "UvDGLM GDD Anthesis", "MvBFT", "UvBFT GDD Silk", "UvBFT GDD Anthesis")

Anthesis_list <- vector(mode = "list", length = 2)
Anthesis_list[[1]] <- BFT.sigsnps$diff.tassel
Anthesis_list[[2]] <- DGLM.sigsnps$diff.tassel
names(Anthesis_list) <- c("BFT", "DGLM")

  
Silking_list <- vector(mode = "list", length = 2)
Silking_list[[1]] <- BFT.sigsnps$diff.silk
Silking_list[[2]] <- DGLM.sigsnps$diff.silk.nomQTN
names(Silking_list) <- c("BFT", "DGLM")

#Silking list, but with mQTN
Silking_mQTN_list <- vector(mode = "list", length = 2)
Silking_mQTN_list[[1]] <- BFT.sigsnps$diff.silk
Silking_mQTN_list[[2]] <- DGLM.mQTN.sigsnps$diff.silk.wMQTN
names(Silking_mQTN_list) <- c("BFT", "DGLM")



#Now, I will create the venn diagrams for this publication
setwd(wd$visuals)

tiff(file = "Difference-in-GDD-to-Anthesis-Venn.tiff",
     width = 800, height = 400, res = 150)
ggVennDiagram(Anthesis_list, set_color = "Black", label_alpha = 0,
              show_intersect = F, label = "count", label_color = "white") +
  scale_fill_gradient(low = "blue", high = "orange")
dev.off()

tiff(file = "Difference-in-GDD-to-Silking-Venn.tiff",
     width = 800, height = 400, res = 150)
ggVennDiagram(Silking_list, set_color = "Black", label_alpha = 0,
              show_intersect = F, label = "count", label_color = "white") +
  scale_fill_gradient(low = "blue", high = "orange")
dev.off()

##DGLM mQTNs
tiff(file = "Difference-in-GDD-to-Silking-mQTN-Venn.tiff",
     width = 800, height = 400, res = 150)
ggVennDiagram(Silking_mQTN_list, set_color = "Black", label_alpha = 0,
              show_intersect = F, label = "count", label_color = "white") +
  scale_fill_gradient(low = "blue", high = "orange")
dev.off()

################################## END #################################
