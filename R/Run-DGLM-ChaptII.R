###########################################################################
######################### Run DGLM ########################################
###########################################################################

#first, load libraries 
library(dglm)
library(data.table)
library(readr)
library(tictoc)

####################### Section #2 ############################

#These next few lines are for testing purposes
Map <- jus.map.list[[5]] 
mymap <- Map[1, ]
p <- residuals[, c(1, 3)]
p <- Phenos[, c(1, cT)]
geno <- jus.geno.list[[5]]
g <- geno[, c(1, 2)]
mQTN <- geno[, c(1, 1000)]

temp_123456789 <- merge(p, g, by = "Taxa")
temp_123456789 <- merge(temp_123456789, mQTN, by = "Taxa")
colnames(temp_123456789) <- c("Taxa", "y", "snp", "mQTN")
rm(p, g)

###Defining the DGLM function
my.pdglm <- function(temp_123456789 = NULL, mQTN = NULL, Map = NULL) {
  #print(paste("--------- Fitting DGLM model for SNP out of ", dim(Geno)[2], "  ----------", sep = ""))
  if (is.null(mQTN) == T) {
    #colnames(temp_123456789) <- c("Taxa", "y", "snp")
    #temp_123456789 <- na.omit(temp_123456789)
    model <-
      dglm(
        formula = y ~ snp,
        ~ snp, data = temp_123456789,
        family = gaussian(link = "identity")
      )
  } else {
    if (mQTN == "present") {
    #colnames(temp_123456789) <- c("Taxa", "y", "snp", "mQTN")
    #temp_abcd <- na.omit(temp_abcd)
    model <-
      dglm(
        formula = as.formula(paste("y ~ snp + mQTN")),
        ~ snp, data = temp_123456789,
        family = gaussian(link = "identity")
      )
  }
  }
  Map_info <- Map[i, ]
  P.mean <- summary(model)$coef["snp", "Pr(>|t|)"]  # Extarct p values for mean part
  an <- anova(model)
  P.disp <- pchisq(q = an["Dispersion model", "Adj.Chisq"], df = an["Dispersion model", "DF"], lower.tail = FALSE)
  df <- an["Dispersion model", "DF"]
  Teststat = an["Dispersion model", "Adj.Chisq"]
  s.model <- summary(model$dispersion.fit)
  beta <- s.model$coef["snp", "Estimate"]  # Extarct cofficients
  se <- s.model$coef["snp", "Std. Error"]  # Extract standard errors
  out <- data.frame(snp = as.character(mymap$snp), chr = as.character(mymap$chr), pos = as.character(mymap$pos), maf = as.numeric(mymap$MAF), DF = df, Teststat = Teststat, Beta = beta, SE = se, P.mean = P.mean, P.disp = P.disp, stringsAsFactors = FALSE)  # Save all the extracted variables in data frame out
  rm(list = "temp_123456789", envir = .GlobalEnv)
  return(out)
}

#Test run
my.pdglm(temp_123456789 = temp_123456789, Map = mymap, mQTN = NULL)

##This is for all residuals but south silk
residuals2 <- residuals[, c(1:2, 4:9)]

for (cT in 2:ncol(residuals)) {
  p <- residuals[, c(1, cT)]
  TF.list <- vector(mode = "list", length(seq_along(jus.geno.list)))
  for (j in seq_along(jus.geno.list)) {
    geno <- jus.geno.list[[j]]
    map <- jus.map.list[[j]]
    TF <- matrix(NA, nrow = dim(map)[1], ncol = 10)
    TF <- as.data.frame(TF)
    colnames(TF) <- c("snp", "chr", "pos", "maf", "DF", "Teststat", "Beta", "SE", "P.mean", "P.disp")
    tic()
    for (k in 2:dim(geno)[2]) {
      g <- geno[, c(1, k)]
      mymap <- map[k - 1, ]
      temp_123456789 <- merge(p, g, by = "Taxa")
      colnames(temp_123456789) <- c("Taxa", "y", "snp")
      try({
        outm <- my.pdglm(temp_123456789 = temp_123456789, Map = mymap, mQTN = NULL)
        TF[k - 1, ] <- outm
        if (k %% 1000 == 0) {
          print(TF[k - 1, ])
        }
      }, silent = T)
    }
    toc()
    TF.list[[j]] <- TF
    #write_csv(TF, file = paste0("BFT_results", colnames(P[2]), j, ".csv"))
  }
  dglm_results <- rbindlist(TF.list)
  dglm_results <- dglm_results[order(P.disp), ]
  write_csv(dglm_results, file = paste0("DGLM_results", colnames(residuals[cT]), ".csv"))
}

#Now for traits with a mQTN
##First, I need to find the statistically significant mQTN
grep(T, colnames(jus.geno.list[[6]]) == "S5_9802633")
mQTN <- jus.geno.list[[6]][, c(1, 5180)]

residuals2 <- residuals[, c(1, 3)]

for (cT in 2:ncol(residuals2)) {
  p <- residuals2[, c(1, cT)]
  TF.list <- vector(mode = "list", length(seq_along(jus.geno.list)))
  for (j in seq_along(jus.geno.list)) {
    geno <- jus.geno.list[[j]]
    map <- jus.map.list[[j]]
    TF <- matrix(NA, nrow = dim(map)[1], ncol = 10)
    TF <- as.data.frame(TF)
    colnames(TF) <- c("snp", "chr", "pos", "maf", "DF", "Teststat", "Beta", "SE", "P.mean", "P.disp")
    tic()
    for (k in 2:dim(geno)[2]) {
      g <- geno[, c(1, k)]
      mymap <- map[k - 1, ]
      temp_123456789 <- merge(p, g, by = "Taxa")
      temp_123456789 <- merge(temp_123456789, mQTN, by = "Taxa")
      colnames(temp_123456789) <- c("Taxa", "y", "snp", "mQTN")
      try({
        outm <- my.pdglm(temp_123456789 = temp_123456789, Map = mymap, mQTN = "present")
        TF[k - 1, ] <- outm
        if (k %% 1000 == 0) {
          print(TF[k - 1, ])
        }
      }, silent = T)
    }
    toc()
    TF.list[[j]] <- TF
    #write_csv(TF, file = paste0("BFT_results", colnames(P[2]), j, ".csv"))
  }
  dglm_results <- rbindlist(TF.list)
  dglm_results <- dglm_results[order(P.disp), ]
  write_csv(dglm_results, file = paste0("DGLM_results", colnames(residuals2[cT]), ".csv"))
}

###Now with a differenty mQTN for a different trait
grep(T, colnames(jus.geno.list[[4]]) == "S3_159384497")
mQTN2 <- jus.geno.list[[4]][, c(1, 18583)]

residuals3 <- residuals[, c(1, 6)]

for (cT in 2:ncol(residuals3)) {
  p <- residuals3[, c(1, cT)]
  TF.list <- vector(mode = "list", length(seq_along(jus.geno.list)))
  for (j in seq_along(jus.geno.list)) {
    geno <- jus.geno.list[[j]]
    map <- jus.map.list[[j]]
    TF <- matrix(NA, nrow = dim(map)[1], ncol = 8)
    TF <- as.data.frame(TF)
    colnames(TF) <- c("snp", "chr", "pos", "maf", "Beta", "SE", "P.mean", "P.disp")
    tic()
    for (k in 2:dim(geno)[2]) {
      g <- geno[, c(1, k)]
      mymap <- map[k - 1, ]
      temp_123456789 <- merge(p, g, by = "Taxa")
      temp_123456789 <- merge(temp_123456789, mQTN2, by = "Taxa")
      colnames(temp_123456789) <- c("Taxa", "y", "snp", "mQTN")
      try({
        outm <- my.pdglm(temp_123456789 = temp_123456789, Map = mymap, mQTN = "present")
        TF[k - 1, ] <- outm
        if (k %% 1000 == 0) {
          print(TF[k - 1, ])
        }
      }, silent = T)
    }
    toc()
    TF.list[[j]] <- TF
    #write_csv(TF, file = paste0("BFT_results", colnames(P[2]), j, ".csv"))
  }
  dglm_results <- rbindlist(TF.list)
  dglm_results <- dglm_results[order(P.disp), ]
  write_csv(dglm_results, file = paste0("DGLM_results", colnames(residuals3[cT]), ".csv"))
}

########################## Section
save.image("chaptII.RData")
save.image("chaptII_copy.RData")
