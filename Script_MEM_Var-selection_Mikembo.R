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

# S?lection mod?le environnemental plus parsimonieux :
# ****************************************************

vif(env)
vifcor(env, th = 0.7)
env <- env[, -c(2,4,6,8,10,12,13,15:17,19:22,25,27,29:31)]   # Pour cor < 70 %

#R2adj <- RsquareAdj(rda(Y, env))$adj.r.squared
#(fwd <- forward.sel(Y, env, adjR2thresh = R2adj, R2more = 0.01))
#order <- fwd$order
#env <- env[, order]

#detrend <- c()

for (a in 1:ncol(spe)) {
  
  Y <- spe[, a] ; spename <- colnames(spe)[a]

#################
# Result matrix #
#################

# Valeurs de p des simulations :
# ******************************
results_FWD <- as.data.frame(matrix(nrow = 1, ncol = 14))   # Valeurs de p des simulations

colnames(results_FWD) <- c("Matrix B", "Matrix A", "Nb_MEM", "X1", "X2", "a", "b", "c", "d", "pval_X1", 
   "pval_X2", "pval_a", "pval_b", "pval_c")
results_FWD[,1] <- "DB (PCNM)"
results_FWD[,2] <- "1-(D/4t)^2"

results_AIC <- as.data.frame(matrix(nrow = 1, ncol = 14))   # Valeurs de p des simulations

colnames(results_AIC) <- c("Matrix B", "Matrix A", "Nb_MEM", "X1", "X2", "a", "b", "c", "d", "pval_X1", 
                           "pval_X2", "pval_a", "pval_b", "pval_c")
results_AIC[,1] <- "DB (PCNM)"
results_AIC[,2] <- "1-(D/4t)^2"

results_MIR <- as.data.frame(matrix(nrow = 1, ncol = 14))   # Valeurs de p des simulations

colnames(results_MIR) <- c("Matrix B", "Matrix A", "Nb_MEM", "X1", "X2", "a", "b", "c", "d", "pval_X1", 
                           "pval_X2", "pval_a", "pval_b", "pval_c")
results_MIR[,1] <- "DB (PCNM)"
results_MIR[,2] <- "1-(D/4t)^2"

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

# Retrieval of the MEM eigenvectors of all final models (when a matrix A is used)
# and significance test by anova.cca:
# ***********************************

R2adj <- RsquareAdj(rda(Y, Y.DBMEM))$adj.r.squared
if(anova.cca(rda(Y, Y.DBMEM))$Pr[1] <= 0.05){
  class <- class(try(fsel <- forward.sel(Y, Y.DBMEM, adjR2thresh = R2adj, nperm = 999,
                                         R2more = 0.01), TRUE))
  if(class != "try-error"){
    sign <- sort(fsel$order)
    MEM.FwdSel <- Y.DBMEM[, c(sign)]
    if(is.matrix(MEM.FwdSel) == TRUE){ results_FWD[1, 3] <- ncol(MEM.FwdSel)
    } else { results_FWD[1, 3] <- 1 }
    MEM <- MEM.FwdSel
  }
} else { 
  results_FWD[1, 3] <- ncol(Y.DBMEM)
  MEM <- Y.DBMEM
}

##############################
# II. Variation Partitioning #
##############################
  
  if (length(MEM) > 0) {
    
    (varpart.real <- varpart(Y, env, MEM))
    
    # Tests of significance (X1, X2, a, c)
    # *********************
    
    # X1 (env)
    results_FWD[1, 10] <- anova.cca(rda(Y, env))$Pr[1]
    results_FWD[1, 4] <- RsquareAdj(rda(Y, env))$adj.r.squared
    
    # X2 (space)
    results_FWD[1, 11] <- anova.cca(rda(Y, MEM))$Pr[1]
    results_FWD[1, 5] <- RsquareAdj(rda(Y, MEM))$adj.r.squared
    
    # Fraction [a], pure environmental
    results_FWD[1, 12] <- anova.cca(rda(Y, env, MEM))$Pr[1]
    results_FWD[1, 6] <- RsquareAdj(rda(Y, env, MEM))$adj.r.squared
    
    # Fraction [c], pure spatial
    results_FWD[1, 14] <- anova.cca(rda(Y, MEM, env))$Pr[1]
    results_FWD[1, 8] <- RsquareAdj(rda(Y, MEM, env))$adj.r.squared
    
    # Fraction [d], residuals
    results_FWD[1, 9] <- 1 - RsquareAdj(rda(Y, cbind(MEM, env)))$adj.r.squared
  
    # Fraction [b], not tested (0 degree of freedom, see Borcard et al. 2011)
    results_FWD[1, 7] <- results_FWD[1, 4] - results_FWD[1, 6]
    
  }   # Fin du if(length(MEM) > 0)

fileFWD <- paste("Results_FWD___", spename, paste(".txt", sep = ""), sep = "")
write.table(results_FWD, file = fileFWD, sep = "\t")

# IV. Generation of the MEM variables with the AIC selection:
#############################################################
#############################################################
#############################################################

Y.DBMEM <- test.W(Y = Y, nb = list, xy = C, MEM.autocor = MEM_model, f = f4, t = thresh)

MEMid <- Y.DBMEM$best$AIC$ord[1:which.min(Y.DBMEM$best$AIC$AICc)]
MEM.select <- Y.DBMEM$best$MEM[, sort(c(MEMid))]
if (length(class(MEM.select)) == 3) { results_AIC[1, 3] <- ncol(MEM.select)
} else { results_AIC[1, 3] <- 1 }
MEM <- MEM.select


##############################
# II. Variation Partitioning #
##############################

for(i in 1:nrow(results_AIC)){
  
  MEM <- list_MEM[[i]]
  
  if(length(MEM) > 0) {
    
    (varpart.real <- varpart(Y, env, MEM))
    
    # Tests of significance (X1, X2, a, c)
    # *********************
    
    # X1 (env)
    results_AIC[i, 10] <- anova.cca(rda(Y, env))$Pr[1]
    results_AIC[i, 4] <- RsquareAdj(rda(Y, env))$adj.r.squared
    
    # X2 (space)
    results_AIC[i, 11] <- anova.cca(rda(Y, MEM))$Pr[1]
    results_AIC[i, 5] <- RsquareAdj(rda(Y, MEM))$adj.r.squared
    
    # Fraction [a], pure environmental
    results_AIC[i, 12] <- anova.cca(rda(Y, env, MEM))$Pr[1]
    results_AIC[i, 6] <- RsquareAdj(rda(Y, env, MEM))$adj.r.squared
    
    # Fraction [c], pure spatial
    results_AIC[i, 14] <- anova.cca(rda(Y, MEM, env))$Pr[1]
    results_AIC[i, 8] <- RsquareAdj(rda(Y, MEM, env))$adj.r.squared
    
    # Fraction [d], residuals
    results_AIC[i, 9] <- 1 - RsquareAdj(rda(Y, cbind(MEM, env)))$adj.r.squared
    
    # Test of the [b] fraction (Vleminckx et al. 2016):
    # *************************************************
    
    # Pour sauver les valeurs simul?es
    E.b     = c() ## vecteur pour stocker le R? de la fraction b
    
    # M = matrice  => colonnes 1 et 2              = coord. x et y
    #              => colonnes 3:1+2         = abond. des sp
    #              => colonnes 1+3 ? ncol(M) = var. en
    
    M <- cbind(C, Y, env)
    
    for(k in 1:ntoro){
      set.seed(k)
      M2 <- M 
      
      rx <- ceiling(runif(1, 0.01, 0.99) * Xmax)
      ry <- ceiling(runif(1, 0.01, 0.99) * Ymax) 
      
      # Pour d?cider si effet miroir ou pas (retournement 180? de la grille)
      sens = runif(1, min = 0, max = 1)  
      if(sens <= 0.5) { M2[,1] <-                (M2 [,1] + rx) %% max(M2 [,1]) 
      M2[,2] <-                (M2 [,2] + ry) %% max(M2 [,2]) }
      if(sens > 0.5)  { M2[,1] <- max(M2 [,1]) - (M2 [,1] + rx) %% max(M2 [,1]) 
      M2[,2] <- max(M2 [,2]) - (M2 [,2] + ry) %% max(M2 [,2]) }
      
      # R?ordonnons pour assembler coordonn?es initiales (dans M3 = M) aux variables
      # environnementales 'toroidalis?es' (dans M2):
      
      M3 <- M ; M2 <- M2[order(M2[,2]), ] ; M2 <- M2[order(M2[,1]), ]
      M3[, (1+3):ncol(M)] <- M2[ ,(1+3):ncol(M)]
      
      # Matrice de variables environnementales translat?es
      
      E.toro <- M3[,(1+3):ncol(M)]   # variables env. translat?es
      
      # Pour sauver r?sultats des partitions de variations des dataset simul?s
      
      E.b <- c(E.b, varpart(M[,c(3:(1+2))], E.toro, MEM)$part$indfract$Adj.R.square[2])
      
    }  # end "for(k in 1:ntoro){"
    
    # Valeurs r?elles
    # ???????????????
    
    R2RDA <- RsquareAdj(rda(Y, env))$adj.r.squared
    a <- RsquareAdj(rda(Y, env, MEM))$adj.r.squared
    b <- R2RDA - a   # R2 de la fraction b
    
    results_AIC[i, 7] <- b
    
    ## P-value R? RDA simul?e (b = valeur observ?e du R?)
    results_AIC[i, 13] <- length ( E.b[E.b > b]) / (ntoro+1) 
    
  }   # Fin du if(length(MEM) > 0)
}   # Fin de la boucle for

fileAIC <- paste("Results_AIC___", spename, paste(".txt", sep = ""), sep = "")
write.table(results_AIC[-c(9:16, 18, 19), ], file = fileAIC, sep = "\t")

}   # Fin du for de l'esp?ce 'a'

write.table(detrend, file = "detrending.by.species.txt", sep = "\t")

