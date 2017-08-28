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
