#################################################################
#*****# I. Variation partitioning of Mikembo tree species #*****#
#################################################################

# Useful packages:
# ****************

library(vegan)
library(adespatial)
library(spdep)
library(usdm)

##################
# Fix parameters #
##################

# Abondance totale minimum sous laquelle on ne considère pas les espèces : 
ab <- 20

##############
# Data input #
##############

C <- read.table("spa_25kr.txt", h = T, sep = "\t", row.names = 1)

env <- read.table("env_25kr.txt", h = T, sep = "\t", row.names = 1)

spe <- read.table("spe_25kr.txt", h = T, sep = "\t", row.names = 1)
sum <- c() ; for(i in 1:ncol(spe)) sum <- c(sum, sum(spe[,i]))
spe <- spe[, which(sum >= ab)]

# Parcimonious environmental model selection:
# *******************************************

vif(env)
vifcor(env, th = 0.7)
env <- env[, -c(2,4,6,8,10,12,13,15:17,19:22,25,27,29:31)]   # Pour cor < 70 %

###################
# Result matrices #
###################

results_FWD <- as.data.frame(matrix(nrow = ncol(spe), ncol = 14))
row.names(results_FWD) <- colnames(spe)
colnames(results_FWD) <- c("Matrix B", "Matrix A", "Nb_MEM", "X1", "X2", "a", "b", "c", "d", "pval_X1", 
                           "pval_X2", "pval_a", "pval_b", "pval_c")
results_FWD[, 1] <- "DB (PCNM)"
results_FWD[, 2] <- "1-(D/4t)^2"

results_AIC <- as.data.frame(matrix(nrow = ncol(spe), ncol = 14))
row.names(results_AIC) <- colnames(spe)
colnames(results_AIC) <- c("Matrix B", "Matrix A", "Nb_MEM", "X1", "X2", "a", "b", "c", "d", "pval_X1", 
                           "pval_X2", "pval_a", "pval_b", "pval_c")
results_AIC[, 1] <- "DB (PCNM)"
results_AIC[, 2] <- "1-(D/4t)^2"

results_MIR <- as.data.frame(matrix(nrow = ncol(spe), ncol = 14))
row.names(results_MIR) <- colnames(spe)
colnames(results_MIR) <- c("Matrix B", "Matrix A", "Nb_MEM", "X1", "X2", "a", "b", "c", "d", "pval_X1", 
                           "pval_X2", "pval_a", "pval_b", "pval_c")
results_MIR[, 1] <- "DB (PCNM)"
results_MIR[, 2] <- "1-(D/4t)^2"

# Lists to keep all the MEM selections based on the three methods and for all species.
# This will allow to display graphical outputs of the species spatial patterns according to
# different MEM selection approaches:

listMEMsel_FWD <- vector("list", ncol(spe))
names(listMEMsel_FWD) <- colnames(spe)
listMEMsel_AIC <- vector("list", ncol(spe))
names(listMEMsel_AIC) <- colnames(spe)
listMEMsel_MIR <- vector("list", ncol(spe))
names(listMEMsel_MIR) <- colnames(spe)

for (a in 1:ncol(spe)) {
  
  Y <- spe[, a]

   ###################
   # I. MEM Analysis #
   ###################

MEM_model = "positive"

xy.d1 <- dist(C)

# Distance-based B (radius around points):
thresh <- give.thresh(xy.d1)
f4 <- function (D, t) {1-(D/(4*t))^2}           # PCNM criterion

list <- dnearneigh(thresh+1, x = as.matrix(C), d1 = 0)
Y.DB.lw <- nb2listw(list, style = "B")
Y.DBMEM <- scores.listw(Y.DB.lw, MEM.autocor = MEM_model)

   ###################
   # II. FWD approach:
   ###################

# Retrieval of the MEM eigenvectors of all final models (when a matrix A is used)
# and significance test by anova.cca:
# ***********************************

R2adj <- RsquareAdj(rda(Y, Y.DBMEM))$adj.r.squared
if (anova.cca(rda(Y, Y.DBMEM))$Pr[1] <= 0.05) {
  class <- class(try(fsel <- forward.sel(Y, Y.DBMEM, adjR2thresh = R2adj, nperm = 999,
                                         R2more = 0.01), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DBMEM[, c(sign)]
#    if(is.matrix(MEM.FwdSel) == TRUE){ 
      results_FWD[a, 3] <- ncol(MEM.FwdSel)
#    } else { 
#      results_FWD[a, 3] <- 1 
#      }
    MEM <- MEM.FwdSel
  }
} else { 
  results_FWD[a, 3] <- ncol(Y.DBMEM)
  MEM <- Y.DBMEM
}

listMEMsel_FWD[[a]] <- MEM


# Variation Partitioning #
##########################
    
    varpart.real <- varpart(Y, env, MEM)
    
    # Tests of significance (X1, X2, a, c)
    # *********************
    
    # X1 (env)
    results_FWD[a, 10] <- anova.cca(rda(Y, env))$Pr[1]
    results_FWD[a, 4] <- RsquareAdj(rda(Y, env))$adj.r.squared
    
    # X2 (space)
    results_FWD[a, 11] <- anova.cca(rda(Y, MEM))$Pr[1]
    results_FWD[a, 5] <- RsquareAdj(rda(Y, MEM))$adj.r.squared
    
    # Fraction [a], pure environmental
    results_FWD[a, 12] <- anova.cca(rda(Y, env, MEM))$Pr[1]
    results_FWD[a, 6] <- RsquareAdj(rda(Y, env, MEM))$adj.r.squared
    
    # Fraction [c], pure spatial
    results_FWD[a, 14] <- anova.cca(rda(Y, MEM, env))$Pr[1]
    results_FWD[a, 8] <- RsquareAdj(rda(Y, MEM, env))$adj.r.squared
    
    # Fraction [d], residuals
    results_FWD[a, 9] <- 1 - RsquareAdj(rda(Y, cbind(MEM, env)))$adj.r.squared
  
    # Fraction [b], not tested (0 degree of freedom, see Borcard et al. 2011)
    results_FWD[a, 7] <- results_FWD[a, 4] - results_FWD[a, 6]


   # III. AIC approach:
   ####################
   ####################
   ####################

Y.DBMEM.AIC <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = f4, t = thresh)

MEMid <- Y.DBMEM.AIC$best$AIC$ord[1:which.min(Y.DBMEM.AIC$best$AIC$AICc)]
MEM.select <- Y.DBMEM.AIC$best$MEM[, sort(c(MEMid))]
if (length(class(MEM.select)) == 3) { results_AIC[a, 3] <- ncol(MEM.select)
} else { results_AIC[a, 3] <- 1 }
MEM <- MEM.select

listMEMsel_AIC[[a]] <- MEM

# Variation Partitioning #
##########################

    varpart.real <- varpart(Y, env, MEM)
    
    # Tests of significance (X1, X2, a, c)
    # *********************
    
    # X1 (env)
    results_AIC[a, 10] <- anova.cca(rda(Y, env))$Pr[1]
    results_AIC[a, 4] <- RsquareAdj(rda(Y, env))$adj.r.squared
    
    # X2 (space)
    results_AIC[a, 11] <- anova.cca(rda(Y, MEM))$Pr[1]
    results_AIC[a, 5] <- RsquareAdj(rda(Y, MEM))$adj.r.squared
    
    # Fraction [a], pure environmental
    results_AIC[a, 12] <- anova.cca(rda(Y, env, MEM))$Pr[1]
    results_AIC[a, 6] <- RsquareAdj(rda(Y, env, MEM))$adj.r.squared
    
    # Fraction [c], pure spatial
    results_AIC[a, 14] <- anova.cca(rda(Y, MEM, env))$Pr[1]
    results_AIC[a, 8] <- RsquareAdj(rda(Y, MEM, env))$adj.r.squared
    
    # Fraction [d], residuals
    results_AIC[a, 9] <- 1 - RsquareAdj(rda(Y, cbind(MEM, env)))$adj.r.squared
    
    # Fraction [b], not tested (0 degree of freedom, see Borcard et al. 2011)
    results_AIC[a, 7] <- results_AIC[a, 4] - results_AIC[a, 6]    
    
    
   ###################
   # IV. MIR approach:
   ###################
   
    source("MEM.moransel.R") 
    moransel <- MEM.moransel(Y, C, Y.DB.lw, MEM.autocor = MEM_model)
    
    if (class(moransel) == "list") {
      results_MIR[a, 3] <- ncol(moransel$MEM.select)
      MEM <- moransel$MEM.select
    } else {
      results_MIR[a, 3] <- ncol(Y.DBMEM)
      MEM <- Y.DBMEM
    }
    
    listMEMsel_MIR[[a]] <- MEM
 
    # Variation Partitioning #
    ##########################
    
    varpart.real <- varpart(Y, env, MEM)
    
    # Tests of significance (X1, X2, a, c)
    # *********************
    
    # X1 (env)
    results_MIR[a, 10] <- anova.cca(rda(Y, env))$Pr[1]
    results_MIR[a, 4] <- RsquareAdj(rda(Y, env))$adj.r.squared
    
    # X2 (space)
    results_MIR[a, 11] <- anova.cca(rda(Y, MEM))$Pr[1]
    results_MIR[a, 5] <- RsquareAdj(rda(Y, MEM))$adj.r.squared
    
    # Fraction [a], pure environmental
    results_MIR[a, 12] <- anova.cca(rda(Y, env, MEM))$Pr[1]
    results_MIR[a, 6] <- RsquareAdj(rda(Y, env, MEM))$adj.r.squared
    
    # Fraction [c], pure spatial
    results_MIR[a, 14] <- anova.cca(rda(Y, MEM, env))$Pr[1]
    results_MIR[a, 8] <- RsquareAdj(rda(Y, MEM, env))$adj.r.squared
    
    # Fraction [d], residuals
    results_MIR[a, 9] <- 1 - RsquareAdj(rda(Y, cbind(MEM, env)))$adj.r.squared
    
    # Fraction [b], not tested (0 degree of freedom, see Borcard et al. 2011)
    results_MIR[a, 7] <- results_MIR[a, 4] - results_MIR[a, 6]         

}   # Fin du for de l'espèce 'a'

write.table(results_FWD, "results_Mikembo_FWD.txt", sep = "\t")
write.table(results_AIC, "results_Mikembo_AIC.txt", sep = "\t")
write.table(results_MIR, "results_Mikembo_MIR.txt", sep = "\t")

# write.table(detrend, file = "detrending.by.species.txt", sep = "\t")

# ****************************************************************************************
# Visualisation of the spatial patterns highlighted by the different eigenvector selection
# methods:
# ****************************************************************************************

listMEMsel_AIC
listMEMsel_FWD
listMEMsel_MIR

colnames(spe)
a <- 13

Y <- spe[, a]

MEM_AIC <- listMEMsel_AIC[[a]]
MEM_FWD <- listMEMsel_FWD[[a]]
MEM_MIR <- listMEMsel_MIR[[a]]

lm_AIC <- lm(Y ~., data = MEM_AIC)
lm_FWD <- lm(Y ~., data = MEM_FWD)
lm_MIR <- lm(Y ~., data = MEM_MIR)

par(mfrow = c(3, 1))
s.value(C, lm_AIC$fitted.values)
s.value(C, lm_FWD$fitted.values)
s.value(C, lm_MIR$fitted.values)

################
# Simulated
# Broad (i = 1000)
x <- MEM[, c(1:3)]
MEM_AIC <- MEM.select
MEM_FWD <- MEM.FwdSel
MEM_MIR <- moransel$MEM.select

lm_ref <- lm(Y ~., data= x)
lm_AIC <- lm(Y ~., data = MEM_AIC)
lm_FWD <- lm(Y ~., data = MEM_FWD)
lm_MIR <- lm(Y ~., data = MEM_MIR)

par(mfrow = c(4, 1))
s.value(C, lm_ref$fitted.values)
s.value(C, lm_AIC$fitted.values)
s.value(C, lm_FWD$fitted.values)
s.value(C, lm_MIR$fitted.values)
