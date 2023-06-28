###################################################################
######### Brown-Forsythe Test: Chapter 2 ##############
##################################################################

#First, let's the required libraries
library(car)
library(data.table)
library(readr)

#Next, let's read in the required files for these analyses
##Residuals
residuals <- read_delim(file.choose())
residuals <- as.data.frame(residuals)

##Numeric genotype file
myG.filt <- read_delim(file = file.choose(), col_names = T)
myG.filt <- as.data.frame(myG.filt)

#to get to 0,1,2 numericalization, I need to add ones
one.matrix <- matrix(1, nrow = nrow(myG.filt), ncol = length(colnames(myG.filt)[6:286]))

myG.filt[, 6:286] <- myG.filt[, 6:286] + one.matrix

###Now, I have to convert the taxa names
inbred_vect <- myGD$taxa

#
modinbred_vect <- gsub(":.*", "", inbred_vect)

#Quality control
length(grep(T, residuals$Taxa %in% modinbred_vect))

#Now, let me check where the falses are coming from
grep(F, blups.data$taxa %in% modinbred_vect)

modinbred_vect[33]
modinbred_vect[159]
modinbred_vect[278]

modinbred_vect[33] <- "B2"
modinbred_vect[159] <- blups.data[160, 1]
modinbred_vect[278] <- blups.data[279, 1]

myGD$taxa <- modinbred_vect


#Now, let me calculate MAF
#Before I run the BFT, I need to calculate MAFs
zm_MAF <- apply(myGD[, -c(1)], 2, mean)
zm_MAF <- matrix(zm_MAF, nrow = 1)
zm_MAF <- apply(zm_MAF, 2, function(x) min(1 - x / 2, x / 2))

#Now, that I have MAFs, I will modify my genotypic dataframe to reflect this change
#myG.filt <- myG.filt[, c(1, 3, 4, 6:286)]
myGM2 <- data.frame(myGM, zm_MAF)

myGD2 <- myGD[, -c(1)]
rownames(myGD2) <- myGD$taxa
myGD2 <- t(myGD2)

myGDM <- data.frame(myGM2, myGD2)

#Now, I need to rename some of the columns
colnames(myGDM)[4:8] <- c("MAF", "4226", "4722", "33-16", "38-11")
colnames(myGDM)[1:3] <- c("snp", "chr", "pos")

####################### Section #2 #####################
###create a geno list
geno.list <- split(myGDM, myGDM$chr)

#Here, I need to create lists for the genotypes and map files
jus.geno.list <- vector(mode = "list", length = 10)
jus.map.list <- vector(mode = "list", length = 10)

##Now, I need to format a few thngs slightly
tic()
for (chr in seq_along(geno.list)) {
  print(chr)
  gm <- geno.list[[chr]]
  map <- gm[, 1:4]
  geno <- gm[, -c(1:4)]
  geno <- t(geno)
  colnames(geno) <- map$snp
  geno <- rownames_to_column(as.data.frame(geno), var = "Taxa")
  geno <- as.data.frame(geno)
  jus.geno.list[[chr]] <- geno
  jus.map.list[[chr]] <- map
  rm(geno)
  rm(map)
}
toc()

### vGWAS via BFT ###
##Now, let's define the BFT function

geno <- as.data.frame(geno)
geno1 <- rownames_to_column(geno, var = "Taxa") #This is for aligning the genotype and phenotype files correctly
my.bft <- function(cT = NULL, phenos = NULL, i = NULL, geno = NULL, map = NULL) {
  mymap <- map[i - 1, ] 
  p <- phenos[, c(1, cT)]
  g <- geno[, c(1, i)]
  temp_123456789 <- merge(p, g, by = "Taxa")
  colnames(temp_123456789) <- c("Taxa", "y", "snp")
  temp_123456789 <- na.omit(temp_123456789)
  vGWA1 <- leveneTest(y ~ as.factor(snp), center = "median", data = temp_123456789)
  p.value <- vGWA1$`Pr(>F)`[1] #populate the p-value column
  DF.gt <- vGWA1$Df[1]
  DF.ind <- vGWA1$Df[2]
  Fstat <- vGWA1$`F value`[1]
  out <- data.frame(snp = as.character(mymap$snp), chr = as.character(mymap$chr), pos = as.character(mymap$pos), maf = as.numeric(mymap$MAF), DF.genotype = DF.gt, DF.individ = DF.ind, F.value = Fstat, p = p.value)
  return(out)
}

#Test run just to make sure everything is working
my.bft(cT = 4, phenos = residuals, i = 2, geno = jus.geno.list[[1]], map = jus.map.list[[1]])

#for loop
for (cT in 2:ncol(residuals)) {
  TF.list <- vector(mode = "list", length(seq_along(jus.geno.list)))
  for (j in seq_along(jus.geno.list)) {
    geno <- jus.geno.list[[j]]
    map <- jus.map.list[[j]]
    TF <- matrix(NA, nrow = dim(map)[1], ncol = 8)
    TF <- as.data.frame(TF)
    colnames(TF) <- c("snp", "chr", "pos", "maf", "DF.genotype", "DF.individ", "F.value", "p.value")
    tic()
    for (k in 2:dim(geno)[2]) {
      try({
        outm <- my.bft(cT = cT, phenos = residuals, i = k, geno = geno, map = map)
        TF[k - 1, ] <- outm
        if (k %% 1000 == 0) {
             print(k)
        }
    }, silent = T)
    }
    toc()
    TF.list[[j]] <- TF
  }
  bft_results <- rbindlist(TF.list)
  bft_results <- bft_results[order(p.value), ]
  write_csv(bft_results, file = paste0("BFT_results", colnames(residuals[cT]), ".csv"))
}

############################# END #############################