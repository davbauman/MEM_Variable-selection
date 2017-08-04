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

MEM.moransel <- function (y, neigh, MEM, nperm = 999, style = "B", alpha = 0.05) {
  SPATIAL = "FALSE"
  listw <- nb2listw(neigh, style = style)
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
  
  if (SPATIAL == "FALSE") return("No significant spatial structure")
  else return(MEM.sel)
  # Written by David Bauman
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

  nb <- dnearneigh(x = as.matrix(C), d1 = 0, d2 = give.thresh(dist(C)))  
  
  moransel <- MEM.moransel(Y, nb, Y.thresh.res$best$MEM)

  