# Bauman et al. 2018. Disentangling good from bad practices of spatial and phylogenetic
# variable selection in eigenvector-based methods. - Ecography.

# Appendix A2: R code used to compute the type I error rate, power and R² estimation accuracy 
# for the different eigenvector selection methods.
# *******************************************************************************************

###############################################################################################
###############################################################################################
###############################################################################################
#*****# Test and comparison of the type I error rate of the FWD, AIC and MIR approaches #*****#
#*****#     for selecting an optimal subset of eigenvectors within a given W matrix.    #*****#
###############################################################################################
###############################################################################################
###############################################################################################


# Useful packages:
# ****************

library(vegan)
library(adespatial)
library(spdep)

# Construction of a result matrix:
# ********************************

# One line = Matrix W created using the db-MEM corresponding to the PCNM criteria (see Methods
# section). For the columns: column 3 = type I error; column 4 = median R2adj; column 5 = sd of 
# the R2adj; 10000 permutations, so that columns 6 to 10005 contain p-values, and columns 
# 10006 to 20005 contain R2adj.

results <- as.data.frame(matrix(nrow = 1, ncol = 2005))

colnames(results) <- c("Matrix B", "Matrix A", "type I error", "mean R2adj", "sd R2adj",
                       paste("p-val", c(1:1000), sep = ""), paste("R2_", c(1:1000), sep = ""))

results[, 1] <- "Thresh MST"   # Two quadrats further away from one another than the largest 
# edge of the minimun spanning tree are not connected.

results[, 2] <- "1-(D/4t)^2"

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 1000   # Max 1000 (otherwise, the result matrix needs to be adapted)

# Generation of the 117 sampled quadrats:
#########################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats
  
  # We define the quadrat coordinates
  X1 <- c()
  Y1 <- c()
  for(i in 1:50){ X1 <- c(X1, rep(i, 25))}
  for(i in 1:50){ Y1 <- c(Y1, c(1:25))}
  
  C[, 1] <- X1
  C[, 2] <- Y1
  
  # We choose the 117 regularly spaced quadrats of the grid
  
  tx <- seq(from = 1, to = 50, by = 4)
  ty <- seq(from = 1, to = 25, by = 3)
  C <- C[C[, 1] %in% tx, ]
  C <- C[C[, 2] %in% ty, ]
  
} else {
  
  # We choose 117 randomly spaced quadrats inside the grid 
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats
  
  set.seed(123)
  C[, 1] <- runif(117, min = 1, max = 50)
  C[, 2] <- runif(117, min = 1, max = 25)
}

xy.d1 <- dist(C)

# Connectivity and weighting matrices based on the PCNM criteria:
# ***************************************************************

funPCNM <- function (D, t) {1-(D/(4*t))^2}          

# Minimum spanning tree
(thresh <- give.thresh(dist(C)))

list <- dnearneigh(thresh, x = as.matrix(C), d1 = 0)

# *******************************************************************************
# The simulation begins here 
# *******************************************************************************

# I. AIC-based approach:
########################
########################
########################

# Simulation procedure:
#######################

for(i in 1:nperm){   
  
  set.seed(i)
  
  Y <- runif(nrow(C), min = 0, max = 20) ; ran <- "runif"             # Random (uniform)
  #   Y <- rnorm(nrow(C), mean = 0, sd = runif(1, 1, 3)) ; ran <- "rnorm" # Random (normal)
  #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                        # Exponential (1)
  #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"                    # Exponential cubed
  
  # Now we can apply the function test.W()
  
  
  Y.thresh.res <- test.W(list, Y = Y, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)

  R2adj <- RsquareAdj(rda(Y, Y.thresh.res$best$MEM))$adj.r.squared
  if(anova.cca(rda(Y, Y.thresh.res$best$MEM))$Pr[1] <= 0.05){
 
    MEMid <- Y.thresh.res$best$AIC$ord[1:which.min(Y.thresh.res$best$AIC$AICc)]
    MEM.aic <- as.data.frame(Y.thresh.res$best$MEM[, sort(c(MEMid))])
    type <- is.vector(MEM.aic)
    # We apply the control of the second stopping criterion of the FWD sel. (R2adj):
    if (type == TRUE) n <- 1 else n <- ncol(MEM.aic)
      for (i in 1:n) {
        if (i == 1) { 
          me <- MEM.aic[, i]
          MEM.select <- MEM.aic[, i]
        } else me <- cbind(me, MEM.aic[, i])
        r2 <- RsquareAdj(rda(Y, me))$adj.r.squared
        if (r2 > R2adj) { break
        } else MEM.select <- me
      }
    results[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
    results[1, i+1005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared
  }
}

# Type I error, median and sd of R2adj:
#######################################

results[1, 3] <- length(which(results[1, c(6:(nperm + 5))] <= 0.05)) / nperm
results[1, 4] <- median(as.numeric(results[1, c(1006:(nperm + 1005))]))
results[1, 5] <- sd(as.numeric(results[1, c(1006:(nperm + 1005))]))

# Output of the results:
# **********************

res_file_name <- paste("Results", framework, ran, paste(design,".mixAicFwd", ".txt", 
                                                        sep = ""), sep = "_")

write.table(results, file = res_file_name, sep = "\t")

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

resultsB_AIC <- as.data.frame(matrix(nrow = 4, ncol = 2005))
colnames(resultsB_AIC) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                            paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsB_AIC[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsB_AIC[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsM_AIC <- as.data.frame(matrix(nrow = 4, ncol = 2005))
colnames(resultsM_AIC) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                            paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsM_AIC[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsM_AIC[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsF_AIC <- as.data.frame(matrix(nrow = 4, ncol = 2005))
colnames(resultsF_AIC) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                            paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsF_AIC[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsF_AIC[,2] <- c("1-(D/4t)^2", NA, NA, NA)

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 1000

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

# Broad scale
#############
#############

x <- MEM[, 1:3]

for(i in 1:nperm){
  
  Y <- spesimB[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsB_AIC[2, 10005+i] <- R2adj
  resultsB_AIC[3, 5+i] <- lmp(lm)
  
  Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  R2adj <- RsquareAdj(rda(Y, Y.MEM$best$MEM))$adj.r.squared
  if(anova.cca(rda(Y, Y.MEM$best$MEM))$Pr[1] <= 0.05){
    
    MEMid <- Y.MEM$best$AIC$ord[1:which.min(Y.MEM$best$AIC$AICc)]
    MEM.aic <- as.data.frame(Y.MEM$best$MEM[, sort(c(MEMid))])
    type <- is.vector(MEM.aic)
    # We apply the control of the second stopping criterion of the FWD sel. (R2adj):
    if (type == TRUE) n <- 1 else n <- ncol(MEM.aic)
    for (i in 1:n) {
      if (i == 1) {
        me <- MEM.aic[, i]
        MEM.select <- MEM.aic[, i]
      } else me <- cbind(me, MEM.aic[, i])
      r2 <- RsquareAdj(rda(Y, me))$adj.r.squared
      if (r2 > R2adj) { break
      } else MEM.select <- me
    }
    resultsB_AIC[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
    resultsB_AIC[1, i+1005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared - resultsB_AIC[2, 1005+i]
  } else {
    resultsB_AIC[1, i+5] <- 1   # p-val made equal to 1 
    resultsB_AIC[1, i+1005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsB_AIC[1, 3] <- length(which(resultsB_AIC[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsB_AIC[1, 4] <- median(na.omit(as.numeric(resultsB_AIC[1, c(1006:(nperm + 1005))])))
resultsB_AIC[1, 5] <- sd(na.omit(as.numeric(resultsB_AIC[1, c(1006:(nperm + 1005))])))
resultsB_AIC[2, 4] <- median(na.omit(as.numeric(resultsB_AIC[2, c(1006:(nperm + 1005))])))
resultsB_AIC[2, 5] <- sd(na.omit(as.numeric(resultsB_AIC[2, c(1006:(nperm + 1005))])))
resultsB_AIC[3, 3] <- length(which(resultsB_AIC[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Broad_mixFWD-AIC", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsB_AIC, file = res_file_name, sep = "\t")


# Medium scale
##############
##############

if(design == "regular") x <- MEM[, 25:27] else x <- MEM[, 19:21]

for(i in 1:nperm){
  
  Y <- spesimM[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsM_AIC[2, 1005+i] <- R2adj
  resultsM_AIC[3, 5+i] <- lmp(lm)
  
  Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  R2adj <- RsquareAdj(rda(Y, Y.MEM$best$MEM))$adj.r.squared
  if(anova.cca(rda(Y, Y.MEM$best$MEM))$Pr[1] <= 0.05){
    
    MEMid <- Y.MEM$best$AIC$ord[1:which.min(Y.MEM$best$AIC$AICc)]
    MEM.aic <- as.data.frame(Y.MEM$best$MEM[, sort(c(MEMid))])
    type <- is.vector(MEM.aic)
    # We apply the control of the second stopping criterion of the FWD sel. (R2adj):
    if (type == TRUE) n <- 1 else n <- ncol(MEM.aic)
    for (i in 1:n) {
      if (i == 1) {
        me <- MEM.aic[, i]
        MEM.select <- MEM.aic[, i]
      } else me <- cbind(me, MEM.aic[, i])
      r2 <- RsquareAdj(rda(Y, me))$adj.r.squared
      if (r2 > R2adj) { break
      } else MEM.select <- me
    }
    resultsM_AIC[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
    resultsM_AIC[1, i+1005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared - resultsM_AIC[2, 1005+i]
  } else {
    resultsM_AIC[1, i+5] <- 1   # p-val made equal to 1 
    resultsM_AIC[1, i+1005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsM_AIC[1, 3] <- length(which(resultsM_AIC[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsM_AIC[1, 4] <- median(na.omit(as.numeric(resultsM_AIC[1, c(1006:(nperm + 1005))])))
resultsM_AIC[1, 5] <- sd(na.omit(as.numeric(resultsM_AIC[1, c(1006:(nperm + 1005))])))
resultsM_AIC[2, 4] <- median(na.omit(as.numeric(resultsM_AIC[2, c(1006:(nperm + 1005))])))
resultsM_AIC[2, 5] <- sd(na.omit(as.numeric(resultsM_AIC[2, c(1006:(nperm + 1005))])))
resultsM_AIC[3, 3] <- length(which(resultsM_AIC[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Medium_mixFWD-AIC", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsM_AIC, file = res_file_name, sep = "\t")


# Fine scale
##############
##############

if(design == "regular") x <- MEM[, 56:58] else x <- MEM[, 37:39]

for(i in 1:nperm){
  
  Y <- spesimF[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsF_AIC[2, 1005+i] <- R2adj
  resultsF_AIC[3, 5+i] <- lmp(lm)
  
  Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  R2adj <- RsquareAdj(rda(Y, Y.MEM$best$MEM))$adj.r.squared
  if(anova.cca(rda(Y, Y.MEM$best$MEM))$Pr[1] <= 0.05){
    
    MEMid <- Y.MEM$best$AIC$ord[1:which.min(Y.MEM$best$AIC$AICc)]
    MEM.aic <- as.data.frame(Y.MEM$best$MEM[, sort(c(MEMid))])
    type <- is.vector(MEM.aic)
    # We apply the control of the second stopping criterion of the FWD sel. (R2adj):
    if (type == TRUE) n <- 1 else n <- ncol(MEM.aic)
    for (i in 1:n) {
      if (i == 1) {
        me <- MEM.aic[, i]
        MEM.select <- MEM.aic[, i]
      } else me <- cbind(me, MEM.aic[, i])
      r2 <- RsquareAdj(rda(Y, me))$adj.r.squared
      if (r2 > R2adj) { break
      } else MEM.select <- me
    }
    resultsF_AIC[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
    resultsF_AIC[1, i+1005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared - resultsF_AIC[2, 1005+i]
  } else {
    resultsF_AIC[1, i+5] <- 1   # p-val made equal to 1 
    resultsF_AIC[1, i+1005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsF_AIC[1, 3] <- length(which(resultsF_AIC[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsF_AIC[1, 4] <- median(na.omit(as.numeric(resultsF_AIC[1, c(1006:(nperm + 1005))])))
resultsF_AIC[1, 5] <- sd(na.omit(as.numeric(resultsF_AIC[1, c(1006:(nperm + 1005))])))
resultsF_AIC[2, 4] <- median(na.omit(as.numeric(resultsF_AIC[2, c(1006:(nperm + 1005))])))
resultsF_AIC[2, 5] <- sd(na.omit(as.numeric(resultsF_AIC[2, c(1006:(nperm + 1005))])))
resultsF_AIC[3, 3] <- length(which(resultsF_AIC[3, c(6:(nperm + 5))] <= 0.05)) / nperm 

# Output of the results:
# **********************

res_file_name <- paste("Power", "Fine_mixFWD-AIC", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsF_AIC, file = res_file_name, sep = "\t")
