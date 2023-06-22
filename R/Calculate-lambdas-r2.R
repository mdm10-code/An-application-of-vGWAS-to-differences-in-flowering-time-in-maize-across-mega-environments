############################################################################
################ Genomic Inflation Factor Script ################
############################################################################

#Load workspace
load("~/Ph.D._Dissertation/Projects/Dissertation_research/Chapter_Two/Revision-1/ChaptII.visuals.RData")

###Here, I am going to create a function that processes my result lists

calculate_inflation <- function(results.list = NULL, test = NULL) {
  lambda.matrix <- matrix(data = NA, nrow = length(results.list), ncol = 2)
  colnames(lambda.matrix) <- c("test", "inflation.factor")
  for (i in seq_along(results.list)) {
    print(i)
    ondeck <- results.list[[i]]
    ondeck <- na.omit(ondeck)
    if (test == "BFT") {
      lambdaBFT <- median(ondeck$F.value / qf(0.5, ondeck$DF.genotype, ondeck$DF.individ))
      lambda.matrix[i, 1] <- names(results.list)[i]
      lambda.matrix[i, 2] <- as.numeric(lambdaBFT)
    } else {
      if(test == "DGLM") {
        teststatDGLM <- qchisq(ondeck$P.disp, df = 1, lower.tail = FALSE)
        lambdaDGLM <- (median(teststatDGLM) / qchisq(0.5, 1))
        lambda.matrix[i, 1] <- names(results.list)[i]
        lambda.matrix[i, 2] <- lambdaDGLM
      }
    }
  }
  return(lambda.matrix)
}

#calculating lambdas
BFT.lambdas <- calculate_inflation(results.list = BFT.list, test = "BFT")
DGLM.lambdas <- calculate_inflation(results.list = DGLM.list, test = "DGLM")
DGLM.nomQTN.lambdas <- calculate_inflation(results.list = DGLM.nomQTNs.list, test = "DGLM")

#To save progress, I will save the workspace
save.image("~/chaptII.RData")