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

results <- as.data.frame(matrix(nrow = 1, ncol = 20005))

colnames(results) <- c("Matrix B", "Matrix A", "type I error", "mean R2adj", "sd R2adj",
   paste("p-val", c(1:10000), sep = ""), paste("R2_", c(1:10000), sep = ""))

results[, 1] <- "Thresh MST"   # Two quadrats further away from one another than the largest 
# edge of the minimun spanning tree are not connected.

results[, 2] <- "1-(D/4t)^2"

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 10000   # Max 10000 (otherwise, the result matrix needs to be adapted)

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

# Retrieval of the MEM eigenvectors and significance test by anova.cca:
# *********************************************************************

MEMid <- Y.thresh.res$best$AIC$ord[1:which.min(Y.thresh.res$best$AIC$AICc)]
 MEM.select <- as.data.frame(Y.thresh.res$best$MEM[, sort(c(MEMid))])
  results[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
   results[1, i+10005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared

}

# Type I error, median and sd of R2adj:
#######################################

   results[1, 3] <- length(which(results[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   results[1, 4] <- median(as.numeric(results[1, c(10006:(nperm + 10005))]))
   results[1, 5] <- sd(as.numeric(results[1, c(10006:(nperm + 10005))]))

# Output of the results:
# **********************

res_file_name <- paste("Results", framework, ran, paste(design,".AIC", ".txt", sep = ""), sep = "_")

write.table(results, file = res_file_name, sep = "\t")



# II. Forward selection approach (Blanchet et al. 2008):
#################################
#################################
#################################

# Construction of a result matrix:
# ********************************

results <- as.data.frame(matrix(nrow = 1, ncol = 20005))

colnames(results) <- c("Matrix B", "Matrix A", "type I error", "mean R2adj", "sd R2adj",
   paste("p-val", c(1:10000), sep = ""), paste("R2_", c(1:10000), sep = ""))

results[, 1] <- "Thresh MST"   
results[, 2] <- "1-(D/4t)^2"

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
 class <- class(try(fsel <- forward.sel(Y, Y.thresh.res$best$MEM, adjR2thresh = R2adj, 
                                        nperm = 999), TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- as.data.frame(Y.thresh.res$best$MEM)[, c(sign)]
   results[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   results[1, i+10005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared
 } } else{ 
          results[1, i+5] <- 1   # p-val made equal to 1 
          results[1, i+10005] <- NA
     }

}


# Type I error, median and sd of R2adj:
#######################################

   results[1, 3] <- length(which(results[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   results[1, 4] <- median(na.omit(as.numeric(results[1, c(10006:(nperm + 10005))])))
   results[1, 5] <- sd(na.omit(as.numeric(results[1, c(10006:(nperm + 10005))])))

# Output of the results:
# **********************

res_file_name <- paste("Results", framework, ran, paste(design,".FwdSel", ".txt", sep = ""), 
                       sep = "_")

write.table(results, file = res_file_name, sep = "\t")

# III. Moran's I minimization approach:
#######################################
#######################################
#######################################

MEM.moransel <- function (y, coord, listw, MEM.autocor = c("positive", "negative", "all"), 
                          nperm = 999, alpha = 0.05) {
  
  # The function computes MEM based on any given listw provided by the user, and performs a
  # MEM variable selection based on the minimization of the
  # Moran's I of the response residuals (no environmental dataset considered here).
  # The function is based on the Moran's I index for univariate response data.
  
  SPATIAL = "FALSE"
  # number of regions:
  nb_sites <- length(y)
  
  MEM.autocor <- match.arg(MEM.autocor) 
  MEM <- scores.listw(listw, MEM.autocor = MEM.autocor)
  
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
      I.vector[i] <- moran(residuals(mod), listw, nb_sites, Szero(listw))$I
    }
    min.moran <- which.min(I.vector)
    # Selection of the MEM variable(s) best minimizing the Moran's I value of the residuals:
    MEM.sel[, sum(nbloop)] <- MEM[, min.moran]
    colnames(MEM.sel)[sum(nbloop)] <- colnames(MEM)[min.moran]
    y <- residuals(lm(y ~ MEM.sel[, sum(nbloop)]))
    I <- moran.mc(y, listw, nperm)
  }
  
  if (SPATIAL == "FALSE") return("No significant spatial structure")
  else list(MEM.all = MEM, MEM.select = MEM.sel)
  # By David Bauman
}

# Construction of a result matrix:
# ********************************

results <- as.data.frame(matrix(nrow = 1, ncol = 20005))

colnames(results) <- c("Matrix B", "Matrix A", "type I error", "mean R2adj", "sd R2adj",
                       paste("p-val", c(1:10000), sep = ""), paste("R2_", c(1:10000), sep = ""))

results[, 1] <- "Thresh MST"

results[, 2] <- "1-(D/4t)^2"

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors

MEM_model = "positive"    ; autocor <- "pos"

# Regular or irregular sampling design:

design <- "regular"   # or "irregular"

nperm <- 10000

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

if (design == "regular") {
  list <- dnearneigh(thresh, x = as.matrix(C), d1 = 0)
} else list <- dnearneigh(thresh+0.00001, x = as.matrix(C), d1 = 0)

listw <- nb2listw(list, style = "B")

# *******************************************************************************
# The simulation begins here 
# *******************************************************************************

# Simulation procedure:
#######################

for(i in 1:nperm){   

    set.seed(i)
    
    Y <- runif(nrow(C), min = 0, max = 20) ; ran <- "runif"             # Random (uniform)
    #   Y <- rnorm(nrow(C), mean = 0, sd = runif(1, 1, 3)) ; ran <- "rnorm" # Random (normal)
    #   Y <- rexp(nrow(C), rate = 1) ; ran <- "rexp"                        # Exponential (1)
    #   Y <- rexp(nrow(C), rate = 1)^3  ; ran <- "rexp3"                    # Exponential cubed
  
  moransel <- MEM.moransel(Y, C, listw, MEM.autocor = MEM_model)
  
  if (class(moransel) == "list") {
    results[1, i+5] <- 0
    results[1, i+10005] <- RsquareAdj(rda(Y, moransel$MEM.select))$adj.r.squared
  } else {   # No spatial structure in the response
    results[1, i+5] <- 1
    results[1, i+10005] <- NA
  }
}

# Type I error, median and sd of R2adj:
#######################################

results[1, 3] <- length(which(results[1, c(6:(nperm + 5))] <= 0.05)) / nperm
results[1, 4] <- median(na.omit(as.numeric(results[1, c(10006:(nperm + 10005))])))
results[1, 5] <- sd(na.omit(as.numeric(results[1, c(10006:(nperm + 10005))])))

# Output of the results:
# **********************

res_file_name <- paste("Results", framework, ran, paste(design,".MoranRes", ".txt", sep = ""), sep = "_")
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

resultsB_AIC <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsB_AIC) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsB_AIC[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsB_AIC[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsB_FWD <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsB_FWD) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsB_FWD[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsB_FWD[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsB_I <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsB_I) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                          paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsB_I[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsB_I[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsM_AIC <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsM_AIC) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsM_AIC[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsM_AIC[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsM_FWD <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsM_FWD) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsM_FWD[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsM_FWD[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsM_I <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsM_I) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
                          paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsM_I[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsM_I[,2] <- c("1-(D/4t)^2", NA, NA, NA)

resultsF_AIC <- as.data.frame(matrix(nrow = 4, ncol = 20005))
colnames(resultsF_AIC) <- c("Matrix B", "Matrix A", "Power", "MedianDeltaR2", "sd DeltaR2",
   paste("p-val", c(1:10000), sep = ""), paste("deltaR2_", c(1:10000), sep = ""))
resultsF_AIC[,1] <- c("Thresh MST", "R2adjReal", "pvalReal", "NbVar")
resultsF_AIC[,2] <- c("1-(D/4t)^2", NA, NA, NA)

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

MEM_model = "positive"    ; autocor <- "pos"

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

   # Broad scale
   #############
   #############

x <- MEM[, 1:3]

for(i in 1:nperm){

   Y <- spesimB[[i]]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsB_FWD[2, 10005+i] <- R2adj
   resultsB_FWD[3, 5+i] <- lmp(lm)

Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)

# Forward selection with double stopping criteria:

R2adj <- RsquareAdj(rda(Y, Y.MEM$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.MEM$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999)
   , TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- as.data.frame(Y.MEM$best$MEM)[, c(sign)]
   resultsB_FWD[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   resultsB_FWD[1, i+10005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared - resultsB_FWD[2, 10005+i]
 } } else{ 
          resultsB_FWD[1, i+5] <- 1   # p-val made equal to 1 
          resultsB_FWD[1, i+10005] <- NA
     }

}

# Power, median and sd of DeltaR2adj:
#####################################

   resultsB_FWD[1, 3] <- length(which(resultsB_FWD[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsB_FWD[1, 4] <- median(na.omit(as.numeric(resultsB_FWD[1, c(10006:(nperm + 10005))])))
   resultsB_FWD[1, 5] <- sd(na.omit(as.numeric(resultsB_FWD[1, c(10006:(nperm + 10005))])))
   resultsB_FWD[2, 4] <- median(na.omit(as.numeric(resultsB_FWD[2, c(10006:(nperm + 10005))])))
   resultsB_FWD[2, 5] <- sd(na.omit(as.numeric(resultsB_FWD[2, c(10006:(nperm + 10005))])))
   resultsB_FWD[3, 3] <- length(which(resultsB_FWD[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Broad_FWD", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsB_FWD, file = res_file_name, sep = "\t")


   # Medium scale
   ##############
   ##############

if(design == "regular") x <- MEM[, 25:27] else x <- MEM[, 19:21]

for(i in 1:nperm){

   Y <- spesimM[[i]]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsM_FWD[2, 10005+i] <- R2adj
   resultsM_FWD[3, 5+i] <- lmp(lm)

Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)

R2adj <- RsquareAdj(rda(Y, Y.MEM$best$MEM))$adj.r.squared
 if(anova.cca(rda(Y, Y.MEM$best$MEM))$Pr[1] <= 0.05){
 class <- class(try(fsel <- forward.sel(Y, Y.MEM$best$MEM, adjR2thresh = R2adj, nperm = 999)
   , TRUE))
  if(class != "try-error"){
   sign <- sort(fsel$order)
   MEM.FwdSel <- as.data.frame(Y.MEM$best$MEM)[, c(sign)]
   resultsM_FWD[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.FwdSel)))$Pr[1]
   resultsM_FWD[1, i+10005] <- RsquareAdj(rda(Y, MEM.FwdSel))$adj.r.squared - resultsM_FWD[2, 10005+i]
 } } else{ 
          resultsM_FWD[1, i+5] <- 1   # p-val made equal to 1 
          resultsM_FWD[1, i+10005] <- NA
     }

}

# Power, median and sd of R2adj:
################################

   resultsM_FWD[1, 3] <- length(which(resultsM_FWD[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsM_FWD[1, 4] <- median(na.omit(as.numeric(resultsM_FWD[1, c(10006:(nperm + 10005))])))
   resultsM_FWD[1, 5] <- sd(na.omit(as.numeric(resultsM_FWD[1, c(10006:(nperm + 10005))])))
   resultsM_FWD[2, 4] <- median(na.omit(as.numeric(resultsM_FWD[2, c(10006:(nperm + 10005))])))
   resultsM_FWD[2, 5] <- sd(na.omit(as.numeric(resultsM_FWD[2, c(10006:(nperm + 10005))])))
   resultsM_FWD[3, 3] <- length(which(resultsM_FWD[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Medium_FWD", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsM_FWD, file = res_file_name, sep = "\t")


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

res_file_name <- paste("Power", "Fine_FWD", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsF_FWD, file = res_file_name, sep = "\t")


# II. Generation of the MEM variables with the data-driven procedure (Dray et al. 2006):
########################################################################################
########################################################################################
########################################################################################

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

MEMid <- Y.MEM$best$AIC$ord[1:which.min(Y.MEM$best$AIC$AICc)]
 MEM.select <- as.data.frame(Y.MEM$best$MEM)[, sort(c(MEMid))]
  resultsB_AIC[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
   resultsB_AIC[1, i+10005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared - resultsB_AIC[2, i+10005]

}

# Power, median and sd of R2adj:
################################

   resultsB_AIC[1, 3] <- length(which(resultsB_AIC[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsB_AIC[1, 4] <- median(na.omit(as.numeric(resultsB_AIC[1, c(10006:(nperm + 10005))])))
   resultsB_AIC[1, 5] <- sd(na.omit(as.numeric(resultsB_AIC[1, c(10006:(nperm + 10005))])))
   resultsB_AIC[2, 4] <- median(na.omit(as.numeric(resultsB_AIC[2, c(10006:(nperm + 10005))])))
   resultsB_AIC[2, 5] <- sd(na.omit(as.numeric(resultsB_AIC[2, c(10006:(nperm + 10005))])))
   resultsB_AIC[3, 3] <- length(which(resultsB_AIC[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Broad_AIC", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsB_AIC, file = res_file_name, sep = "\t")


   # Medium scale
   ##############
   ##############

if(design == "regular") x <- MEM[, 25:27] else x <- MEM[, 19:21]

for(i in 1:nperm){

   Y <- spesimM[[i]]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsM_AIC[2, 10005+i] <- R2adj
   resultsM_AIC[3, 5+i] <- lmp(lm)

Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)

MEMid <- Y.MEM$best$AIC$ord[1:which.min(Y.MEM$best$AIC$AICc)]
 MEM.select <- as.data.frame(Y.MEM$best$MEM)[, sort(c(MEMid))]
  resultsM_AIC[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
   resultsM_AIC[1, i+10005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared - resultsM_AIC[2, i+10005]

}

# VII. Power, median and sd of R2adj:
##########################################

   resultsM_AIC[1, 3] <- length(which(resultsM_AIC[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsM_AIC[1, 4] <- median(na.omit(as.numeric(resultsM_AIC[1, c(10006:(nperm + 10005))])))
   resultsM_AIC[1, 5] <- sd(na.omit(as.numeric(resultsM_AIC[1, c(10006:(nperm + 10005))])))
   resultsM_AIC[2, 4] <- median(na.omit(as.numeric(resultsM_AIC[2, c(10006:(nperm + 10005))])))
   resultsM_AIC[2, 5] <- sd(na.omit(as.numeric(resultsM_AIC[2, c(10006:(nperm + 10005))])))
   resultsM_AIC[3, 3] <- length(which(resultsM_AIC[3, c(6:(nperm + 5))] <= 0.05)) / nperm

# Output of the results:
# **********************

res_file_name <- paste("Power", "Medium_AIC", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsM_AIC, file = res_file_name, sep = "\t")


   # Fine scale
   ############
   ############

if(design == "regular") x <- MEM[, 56:58] else x <- MEM[, 37:39]

for(i in 1:nperm){

   Y <- spesimF[[i]]
   lm <- lm(Y ~ ., data = x)
   R2adj <- summary(lm)$adj.r.squared
   resultsF_AIC[2, 10005+i] <- R2adj
   resultsF_AIC[3, 5+i] <- lmp(lm)

Y.MEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = funPCNM, t = thresh)

MEMid <- Y.MEM$best$AIC$ord[1:which.min(Y.MEM$best$AIC$AICc)]
 MEM.select <- as.data.frame(Y.MEM$best$MEM)[, sort(c(MEMid))]
  resultsF_AIC[1, i+5] <- as.data.frame(anova.cca(rda(Y, MEM.select)))$Pr[1]
   resultsF_AIC[1, i+10005] <- RsquareAdj(rda(Y, MEM.select))$adj.r.squared - resultsF_AIC[2, i+10005]

}

# Power, median and sd of R2adj:
################################

   resultsF_AIC[1, 3] <- length(which(resultsF_AIC[1, c(6:(nperm + 5))] <= 0.05)) / nperm
   resultsF_AIC[1, 4] <- median(na.omit(as.numeric(resultsF_AIC[1, c(10006:(nperm + 10005))])))
   resultsF_AIC[1, 5] <- sd(na.omit(as.numeric(resultsF_AIC[1, c(10006:(nperm + 10005))])))
   resultsF_AIC[2, 4] <- median(na.omit(as.numeric(resultsF_AIC[2, c(10006:(nperm + 10005))])))
   resultsF_AIC[2, 5] <- sd(na.omit(as.numeric(resultsF_AIC[2, c(10006:(nperm + 10005))])))
   resultsF_AIC[3, 3] <- length(which(resultsF_AIC[3, c(6:(nperm + 5))] <= 0.05)) / nperm
   
# Output of the results:
# **********************

res_file_name <- paste("Power", "Fine_AIC", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsF_AIC, file = res_file_name, sep = "\t")

# III. Minimisation of the Moran's I in the Residuals (MIR approach):
# ###################################################################
# ###################################################################
# ###################################################################

# Broad scale
#############
#############

x <- MEM[, 1:3]

for(i in 1:nperm){
  
  Y <- spesimB[[i]]
  lm <- lm(Y ~ ., data = x)
  R2adj <- summary(lm)$adj.r.squared
  resultsB_I[2, 10005+i] <- R2adj
  resultsB_I[3, 5+i] <- lmp(lm)
  
  # Selection based on the minimum number of eigenvectors minimizing the Moran's I of Y's residuals:
  
  moransel <- MEM.moransel(Y, C, Y.DB.lw, MEM.autocor = MEM_model)
  
  if (class(moransel) == "list") {
    resultsB_I[1, i+5] <- 0
    resultsB_I[1, i+10005] <- RsquareAdj(rda(Y, moransel$MEM.select))$adj.r.squared - resultsB_I[2, 10005+i]
    resultsB_I[4, i+5] <- ncol(moransel$MEM.select)
  } else {  # No spatial structure in the response
    resultsB_I[1, i+5] <- 1
    resultsB_I[1, i+10005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsB_I[1, 3] <- length(which(resultsB_I[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsB_I[1, 4] <- median(na.omit(as.numeric(resultsB_I[1, c(10006:(nperm + 10005))])))
resultsB_I[1, 5] <- sd(na.omit(as.numeric(resultsB_I[1, c(10006:(nperm + 10005))])))
resultsB_I[2, 4] <- median(na.omit(as.numeric(resultsB_I[2, c(10006:(nperm + 10005))])))
resultsB_I[2, 5] <- sd(na.omit(as.numeric(resultsB_I[2, c(10006:(nperm + 10005))])))
resultsB_I[3, 3] <- length(which(resultsB_I[3, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsB_I[4, 4] <- median(na.omit(as.numeric(resultsB_I[4, c(6:(nperm + 5))])))
resultsB_I[4, 5] <- sd(na.omit(as.numeric(resultsB_I[4, c(6:(nperm + 5))])))

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
  resultsM_I[2, 10005+i] <- R2adj
  resultsM_I[3, 5+i] <- lmp(lm)
  
  # Selection based on the minimum number of eigenvectors minimizing the Moran's I of Y's residuals:
  
  moransel <- MEM.moransel(Y, C, Y.DB.lw, MEM.autocor = MEM_model)
  
  if (class(moransel) == "list") {
    resultsM_I[1, i+5] <- 0
    resultsM_I[1, i+10005] <- RsquareAdj(rda(Y, moransel$MEM.select))$adj.r.squared - resultsM_I[2, 10005+i]
    resultsM_I[4, i+5] <- ncol(moransel$MEM.select)
  } else {
    resultsM_I[1, i+5] <- 1
    resultsM_I[1, i+10005] <- NA
  }
}

# Power, median and sd of DeltaR2adj:
#####################################

resultsM_I[1, 3] <- length(which(resultsM_I[1, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsM_I[1, 4] <- median(na.omit(as.numeric(resultsM_I[1, c(10006:(nperm + 10005))])))
resultsM_I[1, 5] <- sd(na.omit(as.numeric(resultsM_I[1, c(10006:(nperm + 10005))])))
resultsM_I[2, 4] <- median(na.omit(as.numeric(resultsM_I[2, c(10006:(nperm + 10005))])))
resultsM_I[2, 5] <- sd(na.omit(as.numeric(resultsM_I[2, c(10006:(nperm + 10005))])))
resultsM_I[3, 3] <- length(which(resultsM_I[3, c(6:(nperm + 5))] <= 0.05)) / nperm
resultsM_I[4, 4] <- median(na.omit(as.numeric(resultsM_I[4, c(6:(nperm + 5))])))
resultsM_I[4, 5] <- sd(na.omit(as.numeric(resultsM_I[4, c(6:(nperm + 5))])))

# Output of the results:
# **********************

res_file_name <- paste("Power", "Medium_Moran.Res", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsM_I, file = res_file_name, sep = "\t")

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

res_file_name <- paste("Power", "Fine_Moran.Res", paste(design, ".txt", sep = ""), sep = "_")

write.table(resultsF_I, file = res_file_name, sep = "\t")


