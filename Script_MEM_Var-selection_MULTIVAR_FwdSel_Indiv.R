# Useful packages:
# ****************

library(vegan)
library(adespatial)
library(spdep)

# Result matrix:
# **************
# FwdSel of Blanchet et al. 2008:
results <- as.data.frame(matrix(nrow = 4, ncol = 10004))
colnames(results) <- c("EV selection", "TypeIerror", "Median", "sd",
                            paste("p-val", c(1:10000), sep = ""))
results[,1] <- c("FWD", "FWDind", "NbVar_W_FWD", "NbVar_W_FWDind")

# Definition of the simulation parameters:
##########################################

# Define if we want positive, negative or all eigenvectors
MEM_model = "positive"    ; autocor <- "pos"

# Regular or irregular sampling design:
#design <- "regular"   # or "irregular"

nb_sp <- 10   # Nb of species for the multivariate response
intensity_spe <- 0.6
intensity_env <- 0.5
nperm <- 100

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
  
} else {
  
  # We choose 117 irregularly spaced quadrats inside the grid 
  
  C <- as.matrix(matrix(0, ncol = 2, nrow = 117))   # 1250 = nb quadrats
  
  set.seed(12)
  C[,1] <- runif(117, min = 1, max = 50)
  C[,2] <- runif(117, min = 1, max = 25)
  
}

xy.d1 <- dist(C)

# Generation of the structured Y multivariate responses:
########################################################

# We generate the MEM variables with the db-MEM (PCNM) and generate nperm simulated communities 
# of nb_sp species structured at a medium scale. Each species is structured by three MEM
# variables, some of which can structure more than one species.

funPCNM <- function (D, t) {1-(D/(4*t))^2}          

# Minimum spanning tree
(thresh <- give.thresh(dist(C)))
list <- dnearneigh(thresh, x = as.matrix(C), d1 = 0)
Y.DB.lw <- nb2listw(list)
Y.DBMEM <- scores.listw(Y.DB.lw, MEM.autocor = MEM_model)
MEM <- as.data.frame(Y.DBMEM)
# MEM is the reference used for building the response matrices. 

spesim <- vector("list", nperm)   # List of the nperm simulated species (Medium scale)

# Generation of nperm realisations of a community structured at medium spatial scale:
spe <- matrix(nrow = nrow(C), ncol = nb_sp)

  for(i in 1:nperm){
    set.seed(i)
    for(j in 1:nb_sp) {
      if (j < round(nb_sp/2)) {   # Structured species
        coefs <- sample(c(rep(0, 10), runif(10, 0, 1)))
        mem <- scale(rowSums(sweep(MEM[, c(11:30)], 2, coefs, "*")))
        noise <- scale(rnorm(nrow(C), mean = 0, sd = 1))
        spe[, j] <- (intensity_spe * mem) + ((1-intensity_spe) * noise)
      } else {   # Non structured species
          spe[, j] <- scale(rnorm(nrow(C), mean = 0, sd = 1))
      }
    }
    spesim[[i]] <- spe
  }

# Generation of the structured environmental variables (matrix X):
##################################################################

envsim <- vector("list", nperm)   # List of the nperm simulated species (Medium scale)

# Generation of nperm realisations of a structured (medium scale) environmental matrix X:
# The number of explanatory variables is the species number.

env <- matrix(nrow = nrow(C), ncol = nb_sp)

for(i in 1:nperm){
#  set.seed(i)
  for(j in 1:nb_sp) {
    coefs <- sample(c(rep(0, 10), runif(10, 0, 1)))
    mem <- scale(rowSums(sweep(MEM[, 1:20], 2, coefs, "*")))
    noise <- scale(rnorm(nrow(C), mean = 0, sd = 1))
    env[, j] <- (intensity_env * mem) + ((1-intensity_env) * noise)
  }
  envsim[[i]] <- env
}

# Eigenvector selection based on 1) the normal FWD selection of Blanchet et al. (2008), and
# 2) on the same selection performed for each species individually to keep the union of the
# eigenvectors selected at least once. ####################################################
###########################################################################################

MEM_FS    <- vector("list", nperm)
MEM_FSind <- vector("list", nperm)

for (i in 1:nperm) {
  
  Y <- spesim[[i]]
  
  # I. Forward selection of Blanchet et al. 2008:
  # *********************************************
  R2adj <- RsquareAdj(rda(Y, MEM))$adj.r.squared
  if (anova.cca(rda(Y, MEM))$Pr[1] <= 0.05) {
    #class <- class(try(fsel <- forward.sel(Y, MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
    class <- class(try(fsel <- forward.sel(Y, MEM, nperm = 999), TRUE))
    if (class != "try-error") {
      sign <- sort(fsel$order)
      MEM.FwdSel <- as.data.frame(MEM)[, c(sign)]
      MEM_FS[[i]] <- MEM.FwdSel
      } else MEM_FS[[i]] <- MEM
  } else MEM_FS[[i]] <- MEM
  results[3, (i+4)] <- ncol(MEM_FS[[i]])

  # II. Forward selection of Blanchet et al. 2008 performed for each species separately:
  # ************************************************************************************
  MEM_union <- c()
  for (h in 1:nb_sp) {
    y <- Y[, h]
    R2adj <- RsquareAdj(rda(y, MEM))$adj.r.squared
    if (anova.cca(rda(y, MEM))$Pr[1] <= 0.05) {
      #class <- class(try(fsel <- forward.sel(y, MEM, adjR2thresh = R2adj, nperm = 999), TRUE))
      class <- class(try(fsel <- forward.sel(y, MEM, nperm = 999), TRUE))
      if (class != "try-error") MEM_union <- c(MEM_union, sort(fsel$order))
    }
  }
  sort <- sort(MEM_union[-which(duplicated(MEM_union) == TRUE)])
  MEM_FSind[[i]] <- MEM[, sort]
  results[4, (i+4)] <- ncol(MEM_FSind[[i]])
}

# Variation partitioning and type I error rate of fraction [a] after the two different types 
# of Forward selection (Blanchet et al. 2008): #############################################
############################################################################################

for(i in 1:nperm) {
  
  Y   <- spesim[[i]]
  X   <- envsim[[i]]
  
  # I. Forward selection of Blanchet et al. 2008:
  # *********************************************
  W   <- MEM_FS[[i]]
  
  if (is.vector(W) == FALSE) {
    if (ncol(W) != ncol(MEM)) {
      results[1, i+4] <- anova.cca(rda(Y, X, W))$Pr[1]
      #results[1, i+4] <- anova.cca(rda(Y, X))$Pr[1]
    } else results[1, i+4] <- NA   # We did not manage to get a significantly structured Y
  } else results[1, i+4] <- anova.cca(rda(Y, X, W))$Pr[1]

   # II. Forward selection of Blanchet et al. 2008 performed for each species separately:
   # ************************************************************************************
  W   <- MEM_FSind[[i]]
  
   results[2, i+4] <- anova.cca(rda(Y, X, W))$Pr[1]
    #results[2, i+4] <- anova.cca(rda(Y, X))$Pr[1]
}

# Type I error rate of the fraction [a]:
########################################
results[1, 2] <- length(which(results[1, c(5:nperm+4)] <= 0.05)) / nperm
results[2, 2] <- length(which(results[2, c(5:nperm+4)] <= 0.05)) / nperm
# Median and sd of the number of MEM selected with both FWD selections:
#######################################################################
results[3, 3] <- median(as.numeric(results[3, c(5:nperm+4)][which(results[3, c(5:nperm+4)] != ncol(MEM))]))
results[3, 4] <- sd(as.numeric(results[3, c(5:nperm+4)][which(results[3, c(5:nperm+4)] != ncol(MEM))]))
results[4, 3] <- median(as.numeric(results[4, c(5:nperm+4)][which(results[4, c(5:nperm+4)] != ncol(MEM))]))
results[4, 4] <- sd(as.numeric(results[4, c(5:nperm+4)][which(results[4, c(5:nperm+4)] != ncol(MEM))]))
