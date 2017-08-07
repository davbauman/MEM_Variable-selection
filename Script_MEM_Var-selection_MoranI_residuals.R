# Eigenvector selection based on the minimalization of the response variable Moran's I
# (same procedure as in Griffith and Peres-Neto 2006, but without environmental variables)
# ****************************************************************************************

# I. Type I error rate:
# *********************

# Useful packages and functions:
# ******************************

library(vegan)
library(adespatial)
library(spdep)

# Function for the selection of a subset of MEM variables based on the minimization of the
# Moran's I of the response residuals (no environmental dataset considered here).
# The function is based on the Moran's I index and on a Mantel correlogram for uni- and
# multivariate response data, respectively.

MEM.moransel <- function (y, coord, MEM, nperm = 999, style = "B", alpha = 0.05,
                          response.transform = "hellinger") {
  SPATIAL = "FALSE"
  neigh <- dnearneigh(x = as.matrix(coord), d1 = 0, d2 = give.thresh(dist(coord)))
  listw <- nb2listw(neigh, style = style)
  
  if (is.vector(y) == "TRUE") {  # The response is univariate --> Moran's I test (permutation)
    I <- moran.mc(y, listw, nperm)
    if (I$p.value <= alpha) {
      SPATIAL <- "TRUE"
      MEM.sel <- data.frame(row.names = row.names(MEM))
    }
  
    nbloop <- c()
    while (I$p.value <= alpha) {
      nbloop <- c(nbloop, 1)                   # Loop counter
      I.vector <- vector("numeric", ncol(MEM)) # For the I computed with each MEM variable
      for (i in 1:ncol(MEM)) {
        mod <- lm(y ~ MEM[, i])
        I.vector[i] <- moran(residuals(mod), listw, length(neigh), Szero(listw))$I
      }
      min.moran <- which.min(I.vector)
      # Selection of the MEM variable(s) best minimizing the Moran's I value of the residuals:
      MEM.sel[, sum(nbloop)] <- MEM[, min.moran]
      colnames(MEM.sel)[sum(nbloop)] <- colnames(MEM)[min.moran]
      y <- residuals(lm(y ~ MEM.sel[, sum(nbloop)]))
      I <- moran.mc(y, listw, nperm)
    }
  } else {   # The response is multivariate --> Mantel correlogram
    y <- decostand(y, method = response.transform)
    y.D1 <- dist(y)
    M <- mantel.correlog(y.D1, XY = coord, nperm = nperm)
    sub <- as.numeric(which(M$mantel.res[, 3] > 0))   # Only positive spatial correlation
    signif <- length(which(M$mantel.res[sub, 5] <= alpha))
    if (signif > 0) {
      SPATIAL <- "TRUE"
      MEM.sel <- data.frame(row.names = row.names(MEM))
    }
    
    nbloop <- c()
    while (signif > 0) {
      nbloop <- c(nbloop, 1)                   # Loop counter
      M.vector <- vector("numeric", ncol(MEM))
      for (i in 1:ncol(MEM)) {
        mod <- rda(y, MEM[, i])
        ymod.D1 <- dist(residuals(mod))
        M <- mantel.correlog(ymod.D1, XY = coord, nperm = 1)
        M.vector[i] <- sum(M$mantel.res[as.numeric(which(M$mantel.res[, 3] > 0)), 3])
      }
      min.mantel <- which.min(M.vector)
      # Selection of the MEM variable(s) best minimizing the Moran's I value of the residuals:
      MEM.sel[, sum(nbloop)] <- MEM[, min.mantel]
      colnames(MEM.sel)[sum(nbloop)] <- colnames(MEM)[min.mantel]
      y <- residuals(rda(y, MEM.sel[, sum(nbloop)]))
      
      y.D1 <- dist(y)
      M <- mantel.correlog(y.D1, XY = coord, nperm = nperm)
      sub <- as.numeric(which(M$mantel.res[, 3] > 0))   # Only positive spatial correlation
      signif <- length(which(M$mantel.res[sub, 5] <= alpha))      
    }
  }
  
  if (SPATIAL == "FALSE") return("No significant spatial structure")
  else return(MEM.sel)
  # By David Bauman
}

# Construction of a result matrix:
# ********************************

# One line = Matrix W created using the db-MEM corresponding to the PCNM criteria (see Material and methods).
# For the columns: column 3 = type I error; column 4 = median R2adj; column 5 = sd of the R2adj;
# 1000 permutations, so that columns 6 to 1005 contain p-values, and columns 1006 to 2005 
# contain R2adj.

results <- as.data.frame(matrix(nrow = 1, ncol = 2005))

colnames(results) <- c("Matrix B", "Matrix A", "type I error", "mean R2adj", "sd R2adj",
                       paste("p-val", c(1:1000), sep = ""), paste("R2_", c(1:1000), sep = ""))

results[, 1] <- "Thresh MST"   # Two quadrats further away from one another than the smallest edge of the
# minimun spanning tree are not connected.

results[, 2] <- "1-(D/4t)^2"

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"
#MEM_model = "negative"   ; autocor <- "neg"
#MEM_model = "all"   ; autocor <- "all"  

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 1000

# Generation of the 117 quadrats:
#################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats de 20 m
  
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
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats de 20 m
  
  set.seed(123)
  C[, 1] <- runif(117, min = 1, max = 50)
  C[, 2] <- runif(117, min = 1, max = 25)
}

xy.d1 <- dist(C)

# Connectivity and weighting matrices based on the PCNM method:
# *************************************************************

funPCNM <- function (D, t) {1-(D/(4*t))^2}          

# Minimum spanning tree
(thresh <- give.thresh(dist(C)))

list <- dnearneigh(thresh, x = as.matrix(C), d1 = 0)

# *******************************************************************************
# The simulation begins here 
# *******************************************************************************

# Simulation procedure:
#######################

for(i in 1:nperm){   
  
  if(framework == "univariate"){
    
    set.seed(i)
    
    Y <- runif(nrow(C), min = 0, max = 20) ; ran <- "runif"             # Random (uniform)
    #   Y <- rnorm(nrow(C), mean = 0, sd = runif(1, 1, 3)) ; ran <- "rnorm" # Random (normal)
    #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                        # Exponential (1)
    #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"                    # Exponential cubed
    
  } else {
    
    Y <- matrix(ncol = 5, nrow = nrow(C))
    if(i == 1) CountSeed <- 1
    for(b in 1:ncol(Y)){
      
      set.seed(CountSeed + b)
      
      Y[, b] <- runif(nrow(C), min = 0, max = 20) ; ran <- "runif"    
      #     Y[, b] <- rnorm(nrow(C), mean = 0, sd = runif(1, 1, 3)) ; ran <- "rnorm" 
      #     Y[, b] <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                       
      #     Y[, b] <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3" 
    }  
    CountSeed <- CountSeed + ncol(Y)                 
  }
  
  Y.thresh.res <- test.W(list, Y = Y, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  moransel <- MEM.moransel(Y, C, Y.thresh.res$best$MEM)

  if (class(moransel) == "data.frame") {
    results[1, i+5] <- 0
    results[1, i+1005] <- RsquareAdj(rda(Y, moransel))$adj.r.squared
  } else {
    results[1, i+5] <- 1
    results[1, i+1005] <- NA
  }
}
  
# Type I error, median and sd of R2adj:
#######################################
  
results[1, 3] <- length(which(results[1, c(6:(nperm + 5))] <= 0.05)) / nperm
results[1, 4] <- median(na.omit(as.numeric(results[1, c(1006:(nperm + 1005))])))
results[1, 5] <- sd(na.omit(as.numeric(results[1, c(1006:(nperm + 1005))])))

# Output of the results:
# **********************

res_file_name <- paste("Results", framework, ran, paste(design,".MoranRes", ".txt", sep = ""), sep = "_")
write.table(results, file = res_file_name, sep = "\t")
  
##############################################################################################
  
# II. Power and accuracy:
# ***********************

# Useful packages and functions:
# ******************************

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
# B, M, F stand for Broad, Medium and Fine; AIC et FWD stand for AIC-based selection and forward selection

resultsB_I <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsB_I) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                            paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsB_I[,1] <- c("Thresh MST", "R2adjReal", "pvalReal")
resultsB_I[,2] <- c("1-(D/4t)^2", NA, NA)

resultsM_I <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsM_I) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                            paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsM_I[,1] <- c("Thresh MST", "R2adjReal", "pvalReal")
resultsM_I[,2] <- c("1-(D/4t)^2", NA, NA)

resultsF_I <- as.data.frame(matrix(nrow = 3, ncol = 2005))
colnames(resultsF_I) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                            paste("p-val", c(1:1000), sep = ""), paste("deltaR2_", c(1:1000), sep = ""))
resultsF_I[,1] <- c("Thresh MST", "R2adjReal", "pvalReal")
resultsF_I[,2] <- c("1-(D/4t)^2", NA, NA)

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"
#MEM_model = "negative"   ; autocor <- "neg"
#MEM_model = "positive"   ; autocor <- "all"  

# Do we work with univariate or multivariate random dependent variables?

framework <- "univariate"    # "univariate" or "multivariate" 

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 1000

# Generation of the 117 quadrats:
#################################

if(design == "regular"){
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 1250))   # 1250 = nb quadrats de 20 m
  
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
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats de 20 m
  
  set.seed(12)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25)
  
  # C <- C*20
  # par(mar = c(2, 2, 1, 1))
  # plot(C, pch=15, cex = 2, xlab = "", ylab = "")
  
}

xy.d1 <- dist(C)

# We generate the MEM variables with the db-MEM (PCNM) and generate nperm simulated species structured
# at 1) broad, 2) intermediate, and 3) fine scale.
######################################################################################################

funPCNM <- function (D, t) {1-(D/(4*t))^2}          

# Minimum spanning tree
(thresh <- give.thresh(dist(C)))

list <- dnearneigh(thresh, x = as.matrix(C), d1 = 0)

Y.DB.lw <- nb2listw(list)

Y.DBMEM <- scores.listw(Y.DB.lw, MEM.autocor = MEM_model)
MEM <- as.data.frame(Y.DBMEM)

# MEM is the r?f?rence used for building the response variables. 

spesimB <- vector("list", nperm)   # List of the nperm simulated structured species (Broad scale) 
spesimM <- vector("list", nperm)   # List of the nperm simulated structured species (Medium scale)
spesimF <- vector("list", nperm)   # List of the nperm simulated structured species (Fine scale)

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
    spesimF[[i]] <- (MEM[, 38] * 0.5) + (MEM[, 39] * 1) + (MEM[, 40] * 1) + rnorm(n, mean = 0, sd = 1)
  }
}

# We now use the structured response variables we generated above to test the power and accuracy of both selection
# **************************************************************************************************************** 
# methods:
# ********

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
  resultsB_I[2, 1005+i] <- R2adj
  resultsB_I[3, 5+i] <- lmp(lm)
  
  Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  # Selection based on the minimum number of eigenvectors minimizing the Moran's I of Y's residuals:

  moransel <- MEM.moransel(Y, C, Y.MEM$best$MEM)
  
  if (class(moransel) == "data.frame") {
    resultsB_I[1, i+5] <- 0
    resultsB_I[1, i+1005] <- RsquareAdj(rda(Y, moransel))$adj.r.squared - resultsB_I[2, 1005+i]
  } else {
    resultsB_I[1, i+5] <- 0
    resultsB_I[1, i+1005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsB_I[1, 3] <- length(which(resultsB_I[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsB_I[1, 4] <- median(na.omit(as.numeric(resultsB_I[1, c(1006:(nperm + 1005))])))
resultsB_I[1, 5] <- sd(na.omit(as.numeric(resultsB_I[1, c(1006:(nperm + 1005))])))
resultsB_I[2, 4] <- median(na.omit(as.numeric(resultsB_I[1, c(1006:(nperm + 1005))])))
resultsB_I[2, 5] <- sd(na.omit(as.numeric(resultsB_I[1, c(1006:(nperm + 1005))])))
resultsB_I[3, 3] <- length(which(resultsB_I[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Broad_Moran.Res", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsB_I, file = res_file_name, sep = "\t")

# Medium scale
#############
#############

if(design == "regular") x <- MEM[, 25:27] else x <- MEM[, 19:21]

for(i in 1:nperm){
  
  Y <- spesimM[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsM_I[2, 1005+i] <- R2adj
  resultsM_I[3, 5+i] <- lmp(lm)
  
  Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  # Selection based on the minimum number of eigenvectors minimizing the Moran's I of Y's residuals:
  
  moransel <- MEM.moransel(Y, C, Y.MEM$best$MEM)
  
  if (class(moransel) == "data.frame") {
    resultsM_I[1, i+5] <- 0
    resultsM_I[1, i+1005] <- RsquareAdj(rda(Y, moransel))$adj.r.squared - resultsM_I[2, 1005+i]
  } else {
    resultsM_I[1, i+5] <- 0
    resultsM_I[1, i+1005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsM_I[1, 3] <- length(which(resultsM_I[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsM_I[1, 4] <- median(na.omit(as.numeric(resultsM_I[1, c(1006:(nperm + 1005))])))
resultsM_I[1, 5] <- sd(na.omit(as.numeric(resultsM_I[1, c(1006:(nperm + 1005))])))
resultsM_I[2, 4] <- median(na.omit(as.numeric(resultsM_I[1, c(1006:(nperm + 1005))])))
resultsM_I[2, 5] <- sd(na.omit(as.numeric(resultsM_I[1, c(1006:(nperm + 1005))])))
resultsM_I[3, 3] <- length(which(resultsM_I[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Medium_Moran.Res", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsM_I, file = res_file_name, sep = "\t")

# Fine scale
#############
#############

if(design == "regular") x <- MEM[, 56:58] else x <- MEM[, 38:40]

for(i in 1:nperm){
  
  Y <- spesimF[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsF_I[2, 1005+i] <- R2adj
  resultsF_I[3, 5+i] <- lmp(lm)
  
  Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)
  
  # Selection based on the minimum number of eigenvectors minimizing the Moran's I of Y's residuals:
  
  moransel <- MEM.moransel(Y, C, Y.MEM$best$MEM)
  
  if (class(moransel) == "data.frame") {
    resultsF_I[1, i+5] <- 0
    resultsF_I[1, i+1005] <- RsquareAdj(rda(Y, moransel))$adj.r.squared - resultsF_I[2, 1005+i]
  } else {
    resultsF_I[1, i+5] <- 0
    resultsF_I[1, i+1005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsF_I[1, 3] <- length(which(resultsF_I[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsF_I[1, 4] <- median(na.omit(as.numeric(resultsF_I[1, c(1006:(nperm + 1005))])))
resultsF_I[1, 5] <- sd(na.omit(as.numeric(resultsF_I[1, c(1006:(nperm + 1005))])))
resultsF_I[2, 4] <- median(na.omit(as.numeric(resultsF_I[1, c(1006:(nperm + 1005))])))
resultsF_I[2, 5] <- sd(na.omit(as.numeric(resultsF_I[1, c(1006:(nperm + 1005))])))
resultsF_I[3, 3] <- length(which(resultsF_I[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Fine_Moran.Res", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsF_I, file = res_file_name, sep = "\t")

