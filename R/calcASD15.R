# Gap Analysis package
# CIAT, 2009-2014
# Author: Julian Ramirez-Villegas
# Contributors: Nora Castaneda, Harold Achicanoy

# Calculate proportion of the dist. range with SD above 0.15 (ASD15) for a list of species
# and writes a csv file

calcASD15 <- function(idir) {
  
  odir <- paste(idir, "/summary-files", sep="")
  if (!file.exists(odir)) {
    dir.create(odir)
  }
  
  spList <- list.files(paste(idir, "/occurrence_files", sep=""))
  
  sppC <- 1
  for (spp in spList) {
    spp <- unlist(strsplit(spp, ".", fixed=T))[1]
    fdName <- spp #paste("sp-", spp, sep="")
    spFolder <- paste(idir, "/maxent_modeling/models/", fdName, sep="")
    
    if (file.exists(spFolder)) {
      
      res <- spASD15(idir, spp)
      
      metFile <- paste(spFolder, "/metrics/ASD15.csv", sep="")
      metrics <- read.csv(metFile)
      
      if (sppC == 1) {
        outSum <- metrics
      } else {
        outSum <- rbind(outSum, metrics)
      }
      sppC <- sppC + 1
    }
  }
  outFile <- paste(idir, "/maxent_modeling/summary-files/ASD15.csv", sep="")
  write.csv(outSum, outFile, quote=F, row.names=F)
  return(outSum)
}
