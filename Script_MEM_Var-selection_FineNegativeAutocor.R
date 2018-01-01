# Bauman et al. 2018. Disentangling good from bad practices of spatial and phylogenetic
# variable selection in eigenvector-based methods. - Ecography.

# Appendix A2: R code used to compute the type I error rate, power and R² estimation accuracy 
# for the different eigenvector selection methods.
# *******************************************************************************************

################################################################################################
################################################################################################
################################################################################################
# ******       Test of the power and R² estimation accuracy of the three EV selection # ****** #
# ******                                    approaches                                # ****** #
################################################################################################
################################################################################################
################################################################################################


# Useful packages and functions:
# ******************************

library(vegan)
library(adespatial)
library(spdep)

# Function allowing retrieving the pvalue of a function lm():

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# Construction of a results matrix for each scale and for MEM selection methods:
# ******************************************************************************
# B, M, F stand for Broad, Medium and Fine; AIC, FWD and MIR stand for AIC-based selection,
# forward selection of Blanchet et al. 2008, and Minimisation of Moran's I in Residuals.

resultsF_FWD <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsF_FWD) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                            paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsF_FWD[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsF_FWD[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsF_I <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsF_I) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                          paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsF_I[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsF_I[,2] <- c("1-(D/4t)^2", NA, NA, NA)

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "negative"    ; autocor <- "neg"

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 10000

# Generation of the 117 quadrats:
#################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats
  
  # We define the quadrat coordinates
  X1 <- c()
  Y1 <- c()
  for(i in 1:50){ X1 <- c(X1, rep(i, 25))}
  for(i in 1:50){ Y1 <- c(Y1, c(1:25))}
  
  C[,1] <- X1
  C[,2] <- Y1
  
  # We choose the 117 regularly spaced quadrats of the grid
  
  tx <- seq(from = 1, to = 50, by = 4)
  ty <- seq(from = 1, to = 25, by = 3)
  C <- C[C[,1] %in% tx, ]
  C <- C[C[,2] %in% ty, ]
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
} else {
  
  # We choose 117 irregularly spaced quadrats inside the grid 
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats
  
  set.seed(12)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25)
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
}

xy.d1 <- dist(C)

# We generate the MEM variables with the db-MEM (PCNM) and generate nperm simulated species 
# structured at 1) broad, 2) intermediate, and 3) fine scale.
###########################################################################################

funPCNM <- function (D, t) {1-(D/(4*t))^2}          

# Minimum spanning tree
(thresh <- give.thresh(dist(C)))

list <- dnearneigh(thresh, x = as.matrix(C), d1 = 0)

Y.DB.lw <- nb2listw(list)

Y.DBMEM <- scores.listw(Y.DB.lw, MEM.autocor = MEM_model)
MEM <- as.data.frame(Y.DBMEM)

# MEM is the reference used for building the response variables. 

spesimB <- vector("list", nperm)   # List of the nperm simulated species (Broad scale) 
spesimM <- vector("list", nperm)   # List of the nperm simulated species (Medium scale)
spesimF <- vector("list", nperm)   # List of the nperm simulated species (Fine scale)

n <- nrow(C)

# Generation of nperm realisations of a structured species at broad, medium and fine scale:

if(design == "regular"){
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.5) + (MEM[, 2] * 0.5) + (MEM[, 3] * 0.5) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 25] * 0.6) + (MEM[, 26] * 1) + (MEM[, 27] * 0.8) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 56] * 0.5) + (MEM[, 57] * 1) + (MEM[, 58] * 1) + rnorm(n, mean = 0, sd = 1)
  }  
} else {                  # Irregular design
  for(i in 1:nperm){
    set.seed(i)
    spesimB[[i]] <- (MEM[, 1] * 0.5) + (MEM[, 2] * 0.5) + (MEM[, 3] * 0.5) + rnorm(n, mean = 0, sd = 1)
    spesimM[[i]] <- (MEM[, 19] * 0.6) + (MEM[, 20] * 1) + (MEM[, 21] * 0.8) + rnorm(n, mean = 0, sd = 1)
    spesimF[[i]] <- (MEM[, 37] * 0.5) + (MEM[, 38] * 1) + (MEM[, 39] * 1) + rnorm(n, mean = 0, sd = 1)
  }
}


# We now use the structured response variables we generated above to test the power and 
# accuracy of the three MEM subsets (three EV selection approaches):
# *************************************************************************************

# I. Forward selection with two stopping criteria (Blanchet et al. 2008):
# #######################################################################
# #######################################################################
# #######################################################################

# Fine scale
##############
##############

if(design == "regular") x <- MEM[, 56:58] else x <- MEM[, 37:39]

for(i in 1:nperm){
  
  Y <- spesimF[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsF_FWD[2, 10005+i] <- R2adj
  resultsF_FWD[3, 5+i] <- lmp(lm)
  
  Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  R2adj <- RsquareAdj(rda(Y, Y.MEM$best$MEM))$adj.r.squared
  if(anova.cca(rda(Y, Y.MEM$best$MEM))$Pr[1] <= 0.05){
    class <- class(try(fsel <- forward.sel(Y, Y.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999)
                       , TRUE))
    if(class != "try-error"){
      sign <- sort(fsel$order)
      MEM.FwdSel <- as.data.frame(Y.MEM$best$MEM)[, c(sign)]
      resultsF_FWD[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
      resultsF_FWD[1, i+10005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared - resultsF_FWD[2, 10005+i]
    } } else{ 
      resultsF_FWD[1, i+5] <- 1   # p-val made equal to 1 
      resultsF_FWD[1, i+10005] <- NA
    }
  
}

# Power, median and sd of R2adj:
################################

resultsF_FWD[1, 3] <- length(which(resultsF_FWD[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsF_FWD[1, 4] <- median(na.omit(as.numeric(resultsF_FWD[1, c(10006:(nperm + 10005))])))
resultsF_FWD[1, 5] <- sd(na.omit(as.numeric(resultsF_FWD[1, c(10006:(nperm + 10005))])))
resultsF_FWD[2, 4] <- median(na.omit(as.numeric(resultsF_FWD[2, c(10006:(nperm + 10005))])))
resultsF_FWD[2, 5] <- sd(na.omit(as.numeric(resultsF_FWD[2, c(10006:(nperm + 10005))])))
resultsF_FWD[3, 3] <- length(which(resultsF_FWD[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Fine_FWD", paste(design, "_", MEM_model,
                                                  ".txt", sep = ""), sep = "_")

write.table(resultsF_FWD, file = res_file_name, sep = "\t")

# III. Minimisation of the Moran's I in the Residuals (MIR approach):
# ###################################################################
# ###################################################################
# ###################################################################

# Fine scale
#############
#############

if(design == "regular") x <- MEM[, 56:58] else x <- MEM[, 37:39]

for(i in 1:nperm){
  
  Y <- spesimF[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsF_I[2, 10005+i] <- R2adj
  resultsF_I[3, 5+i] <- lmp(lm)
  
  # Selection based on the minimum number of eigenvectors minimizing the Moran's I of Y's residuals:
  
  moransel <- MEM.moransel(Y, C, Y.DB.lw, MEM.autocor = MEM_model)
  
  if (class(moransel) == "list") {
    resultsF_I[1, i+5] <- 0
    resultsF_I[1, i+10005] <- RsquareAdj(rda(Y, moransel$MEM.select))$adj.r.squared - resultsF_I[2, 10005+i]
    resultsF_I[4, i+5] <- ncol(moransel$MEM.select)
  } else {
    resultsF_I[1, i+5] <- 1
    resultsF_I[1, i+10005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsF_I[1, 3] <- length(which(resultsF_I[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsF_I[1, 4] <- median(na.omit(as.numeric(resultsF_I[1, c(10006:(nperm + 10005))])))
resultsF_I[1, 5] <- sd(na.omit(as.numeric(resultsF_I[1, c(10006:(nperm + 10005))])))
resultsF_I[2, 4] <- median(na.omit(as.numeric(resultsF_I[2, c(10006:(nperm + 10005))])))
resultsF_I[2, 5] <- sd(na.omit(as.numeric(resultsF_I[2, c(10006:(nperm + 10005))])))
resultsF_I[3, 3] <- length(which(resultsF_I[3, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsF_I[4, 4] <- median(na.omit(as.numeric(resultsF_I[4, c(6:(nperm + 5))])))
resultsF_I[4, 5] <- sd(na.omit(as.numeric(resultsF_I[4, c(6:(nperm + 5))])))

# Output of the results:
# **********************

res_file_name <- paste("Power", "Fine_Moran.Res", paste(design, "_", MEM_model,
                                                        ".txt", sep = ""), sep = "_")

write.table(resultsF_I, file = res_file_name, sep = "\t")


