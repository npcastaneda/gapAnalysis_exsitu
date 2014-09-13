# Gap Analysis package
# CIAT, 2009-2014
# Author: Julian Ramirez-Villegas
# Contributors: Nora Castaneda, Harold Achicanoy

# Calculates proportion of the dist. range with SD above 0.15 (ASD15) for a given species

spASD15 <- function(idir, spID) {
  cat("Taxon", spID, "\n")
  spFolder <- paste(idir, "/maxent_modeling/models/", spID, sep="")
  projFolder <- paste(spFolder, "/projections", sep="")
  
  esdCpt <- paste(spID, "_worldclim2_5_ESD.asc.gz", sep="") #s.d. of prob (of k-folds)
  esdThr <- paste(spID, "_worldclim2_5_ESD_PR.asc.gz", sep="") #prob. s.d. (of k-folds) but only in areas of presence
  dumm <- paste(spID, "_worldclim2_5_EMN.asc.gz", sep="") #mean probability (of k-folds)
  
  if (file.exists(paste(projFolder,"/",dumm,sep=""))) {
    cat("..Reading raster files \n")
    dumm <- zipRead(projFolder, dumm) #load mean prob.
    esdCpt <- zipRead(projFolder, esdCpt) #load s.d. prob
    esdThr <- zipRead(projFolder, esdThr) #load thresholded s.d.
    
    esdCpt[which(dumm[] < 0.001)] <- NA #from 'all area' remove areas with ~0 mean prob
    rm(dumm); g=gc(); rm(g) #remove mean prob. raster
    
    esdThr[which(esdThr[] == 0)] <- NA #set to NA areas of zero (absence)
    
    cat("..Calculating \n")
    szCpt <- length(which(esdCpt[] >= 0)) #total no. locations with an s.d. value
    szCptUncertain <- length(which(esdCpt[] >= 0.15)) #locations above threshold
    rateCpt <- szCptUncertain / szCpt * 100 #percentage
    
    szThr <- length(which(esdThr[] >= 0)) #locations with an s.d. value in thresholded s.d. prob
    szThrUncertain <- length(which(esdThr[] >= 0.15)) #locations above threshold
    rateThr <- szThrUncertain / szThr * 100 #percentage
    
    #putting into a data.frame for writing an output file
    cat("..Writing results \n")
    dfOut <- data.frame(taxon=spID, sizeComplete=szCpt, sizeCompleteUncertain=szCptUncertain, rateComplete=rateCpt, sizeThresholded=szThr, sizeThresholdedUncertain=szThrUncertain, rateThresholded=rateThr)
  } else {
    cat("..Writing results \n")
    dfOut <- data.frame(taxon=spID, sizeComplete=NA, sizeCompleteUncertain=NA, rateComplete=NA, sizeThresholded=NA, sizeThresholdedUncertain=NA, rateThresholded=NA)
  }
  
  oFile <- paste(spFolder, "/metrics/ASD15.csv", sep="")
  write.csv(dfOut, oFile, quote=F, row.names=F)
  
  return(dfOut)
}
