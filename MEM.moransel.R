# Bauman, D. et al. 2017. Disentangling good from bad practices in the selection of spatial or 
# phylogenetic eigenvectors. – Ecography 000: 000–000.      
# Appendix A3: Function MEM.moransel used to the perform the MIR procedures.

# The function computes MEM variables based on any given listw provided by the user, and 
# select the smallest subset of spatial eigenvectors minimising the autocorrelation (Moran's I) 
# of the residuals of a model containing only an intercept term (no environmental variables).

# 'y' is a univariate response vector; 'coord' is a matrix of spatial coordinates (cartesian);
# listw is a spatial weighting matrix of class 'listw'.

MEM.moransel <- function (y, coord, listw, MEM.autocor = c("positive", "negative", "all"), 
                          nperm = 999, alpha = 0.05) {
  
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