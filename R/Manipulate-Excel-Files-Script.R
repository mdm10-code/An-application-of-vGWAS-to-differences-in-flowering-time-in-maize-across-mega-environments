##################################################################
################ Excel Manipulation Script #######################
##################################################################

#The purpose of this script is to manipulate the LD results
#from TASSEL in order to be plotted. 

#Load R workspace
load("chaptII.RData")

#For dealing with excel files, I am going to use the rio R package
install.packages("rio")

#Let's load our required libraries
library(dplyr)
library(tidyr)
library(rio)

#Now that I have my required R packages, I am going to now make
#a path to the required excel file. As always, I am going to create
#switches

wd <- list()
wd$pipeline <- "C:/Users/mdm10/Documents/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Two/Revision-1/2_pipeline/"

#Now that I have created my path variable, I am going to
#list the files to see which one is the excel file
list.files(wd$pipeline)

#from here, I will know import the data needed
Excel.results <- import_list(paste0(wd$pipeline, "LD.results.ALL.xlsx"))

#To see which sheets I need, I will print the names of the data objct
names(Excel.results)

#before I join the datasets, I am going to filter each vGWAS
BFT.GDDanthesis.filt <- Excel.results[[5]]
BFT.GDDanthesis.filt <- BFT.GDDanthesis.filt[BFT.GDDanthesis.filt$chr == 9, ]

#
BFT.GDDsilk.filt <- Excel.results[[6]]
BFT.GDDsilk.filt <- BFT.GDDsilk.filt[BFT.GDDsilk.filt$chr == 9, ]

#
DGLM.GDDanthesis.filt <- Excel.results[[8]]
DGLM.GDDanthesis.filt <- DGLM.GDDanthesis.filt[DGLM.GDDanthesis.filt$chr == 9, ]

#
DGLM.GDDsilk.filt <- Excel.results[[7]]
DGLM.GDDsilk.filt <- DGLM.GDDsilk.filt[DGLM.GDDsilk.filt$chr == 9, ]



#Now, I can join sheets for each vGWAS with the input file
pre.plot.BFT.GDDanthesis <- inner_join(BFT.GDDanthesis.filt, Excel.results[[9]], by = "Position2")
pre.plot.BFT.GDDsilk <- inner_join(BFT.GDDsilk.filt, Excel.results[[9]], by = "Position2")
pre.plot.DGLM.GDDanthesis <- inner_join(DGLM.GDDanthesis.filt, Excel.results[[9]], by = "Position2")
pre.plot.DGLM.GDDsilk <- inner_join(DGLM.GDDsilk.filt, Excel.results[[9]], by = "Position2")

#To subset this for the LD plot, I need to do two things.
#first, I am going to need to subset to four columns. Second, 
##I am going to further subset based on LD and basepair distance
colnames(pre.plot.BFT.GDDanthesis)
colnames(pre.plot.BFT.GDDsilk)
colnames(pre.plot.DGLM.GDDanthesis)
colnames(pre.plot.DGLM.GDDsilk)

#Now, I am going to subset everything
pre.plot.BFT.GDDanthesis <- pre.plot.BFT.GDDanthesis[, c(3, 9, 10, 8)]
pre.plot.BFT.GDDsilk <- pre.plot.BFT.GDDsilk[, c(3, 9, 10, 8)]
pre.plot.DGLM.GDDanthesis <- pre.plot.DGLM.GDDanthesis[, c(3, 9, 10, 8)]
pre.plot.DGLM.GDDsilk <- pre.plot.DGLM.GDDsilk[, c(3, 9, 10, 8)]

#let's change the colnames of these objects
colnames(pre.plot.BFT.GDDanthesis) <- c("Position1", "R^2", "SNP", "P.value")
colnames(pre.plot.BFT.GDDsilk) <- c("Position1", "R^2", "SNP", "P.value")
colnames(pre.plot.DGLM.GDDanthesis) <- c("Position1", "R^2", "SNP", "P.value")
colnames(pre.plot.DGLM.GDDsilk) <- c("Position1", "R^2", "SNP", "P.value")


#Now, I will further subset this based on 3 million bp on each side of the
##most significant SNP
pre.plot.BFT.GDDanthesis.narrowed <- pre.plot.BFT.GDDanthesis[pre.plot.BFT.GDDanthesis$Position1 >= 151942763, ]
pre.plot.BFT.GDDanthesis.narrowed <- pre.plot.BFT.GDDanthesis.narrowed[pre.plot.BFT.GDDanthesis.narrowed$Position1 <= 157942763, ]

#
pre.plot.BFT.GDDsilk.narrowed <- pre.plot.BFT.GDDsilk[pre.plot.BFT.GDDsilk$Position1 >= 151942763, ]
pre.plot.BFT.GDDsilk.narrowed <- pre.plot.BFT.GDDsilk.narrowed[pre.plot.BFT.GDDsilk.narrowed$Position1 <= 157942763, ]

#
pre.plot.DGLM.GDDanthesis.narrowed <- pre.plot.DGLM.GDDanthesis[pre.plot.DGLM.GDDanthesis$Position1 >= 151942763, ]
pre.plot.DGLM.GDDanthesis.narrowed <- pre.plot.DGLM.GDDanthesis.narrowed[pre.plot.DGLM.GDDanthesis.narrowed$Position1 <= 157942763, ]

#
pre.plot.DGLM.GDDsilk.narrowed <- pre.plot.DGLM.GDDsilk[pre.plot.DGLM.GDDsilk$Position1 >= 151942763, ]
pre.plot.DGLM.GDDsilk.narrowed <- pre.plot.DGLM.GDDsilk.narrowed[pre.plot.DGLM.GDDsilk.narrowed$Position1 <= 157942763, ]


#Finally, Let's save our progress
save.image("chaptII.RData")

###################################### END ###################################
