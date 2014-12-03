# An R code to calculate the areas of potential distribution range, herbarium buffered samples, germplasm
# buffered samples and convex hulls for crop wild relatives.
#
# "Adapting crops to climate change: collecting, protecting and preparing crop wild relatives"
# www.cwrdiversity.org
#
# J. Ramirez, N. Castaneda,  H. Achicanoy - 2013

require(rgdal)
require(raster)
require(sp)
require(maptools)
require(rgeos)
#Calculate the size of the DR, of the convexhull in km2, of the native area, and of the herbarium samples
#based on the area of the cells

sizeDR <- function(bdir, crop, spID) {
	idir <- paste(bdir, "/maxent_modeling", sep="")
	ddir <- paste(bdir, "/samples_calculations", sep="")
	naDir = paste(bdir, "/biomod_modeling/native-areas/polyshps", sep="") #Now restrict data to known native areas
	
	#Creating the directories
	if (!file.exists(ddir)) {
		dir.create(ddir)
	}
	
	spOutFolder <- paste(ddir, "/", spID, sep="")
	if (!file.exists(spOutFolder)) {
		dir.create(spOutFolder)
	}
  
  #Start estimating areas
	cat("Taxon", spID, "\n")
	
	mskArea <- paste(bdir, "/masks/cellArea.tif", sep="")
	mskArea <- raster(mskArea, values=T)
	msk <- paste(bdir, "/masks/mask.tif", sep="")
	msk <- raster(msk)
	
	#Size of the convex-hull
	spList <- read.csv(paste(bdir, "/summary-files/taxaForRichness.csv", sep=""))
  isValid <- spList$IS_VALID[which(spList$TAXON == paste(spID))]
	occ <- paste(bdir, "/occurrence_files_narea/", spID, ".csv", sep="")
	
	if (isValid == 0){
	  cat("Reading occurrences \n")

	  if(file.exists(occ)){
	    occ <- read.csv(occ) #Only uses records within native area
	    
	   if (!file.exists(paste(spOutFolder, "/convex-hull.tif",sep=""))) {
	     cat("Creating the convex hull \n")
	     #     ch <- occ[chull(occ$lon,occ$lat)]
	     ch <- occ[chull(cbind(occ$lon,occ$lat)),1:2]
	     
	     #   	ch <- occ[chull(cbind(occ$lon, occ$lat)),2:3]
	     ch <- rbind(ch, ch[1,])
	     
	     cat("Transforming to polygons \n")
	     pol <- SpatialPolygons(list(Polygons(list(Polygon(ch)), 1)))
	     grd <- rasterize(pol, msk)
	     
	     cat("Final fixes \n")
	     grd[which(!is.na(grd[]))] <- 1
	     grd[which(is.na(grd[]) & msk[] == 1)] <- 0
	     grd[which(is.na(msk[]))] <- NA
	     
	     cat("Writing convex hull \n")
	     chName <- writeRaster(grd,paste(spOutFolder,"/convex-hull.tif", sep=""))
	     
	     cat("Size of the convex hull \n")
	     grd <- grd * mskArea
	     areaCH <- sum(grd[which(grd[] != 0)])
	     rm(grd)
	   } else {
       grd <- raster(paste(spOutFolder,"/convex-hull.tif", sep=""))
       grd <- grd * mskArea
       areaCH <- sum(grd[which(grd[] != 0)])
       rm(grd)
	   }
	  } else {
	    areaCH <- NA
	  }
	} else {
	  areaCH <- NA
	}
	
	# Size of the native area
	  areaNA <- NA
	
	#Load all occurrences
	allOcc <- read.csv(paste(bdir, "/occurrences/",crop,".csv", sep=""))
	allOcc <- allOcc[which(allOcc$Taxon == spID),]
	
	#Size of the herbarium samples CA50
	cat("Size of the h-samples buffer \n")
	hOcc <- allOcc[which(allOcc$H == 1),]
	if (nrow(hOcc) != 0) {
	  xy <- cbind(hOcc$lon, hOcc$lat)
    xy <- data.frame(xy)
    xy <- unique(xy)
	  occ <- SpatialPoints(xy)
	  
	  cat("Loading", spID, "native areas \n")
	  narea = paste(naDir, "/", spID, "/narea.shp", sep="")
    if(!file.exists(narea)){
      hOcc <- as.data.frame(cbind(as.character(hOcc$Taxon), hOcc$lon, hOcc$lat))
      names(hOcc) <- c("taxon", "lon", "lat")
      write.csv(hOcc, paste(spOutFolder, "/hsamples.csv", sep=""), quote=F, row.names=F)
      buff <- gBuffer(occ, width=0.5,byid=T)
      occBuff <- rasterize(buff,msk)
      grd <- reclassify(occBuff,c(minValue(occBuff),maxValue(occBuff),1))
      grd <- writeRaster(grd,paste(spOutFolder,"/hsamples-buffer.tif",sep=""))      
   }else{
      narea = readShapeSpatial(narea)
      cat ("Projecting files \n")
      proj4string(occ) = CRS("+proj=longlat +datum=WGS84")
      proj4string(narea) = CRS("+proj=longlat +datum=WGS84")
      cat("Selecting occurrences within native area \n")
      x <- over(narea, occ)
      x <- sum(x, na.rm=T)
      if(x==0){
        cat("No points within native area \n")
      } else {
        occ = occ[narea]
        occ = as.data.frame(occ)
        occBuff <- occ
        names(occ) = c("lon","lat")
        occ["taxon"] <- spID
        write.csv(occ, paste(spOutFolder, "/hsamples.csv", sep=""), quote=F, row.names=F)
        occBuff <- SpatialPoints(occBuff)
        buff <- gBuffer(occBuff, width=0.5,byid=T)
        occBuff <- rasterize(buff,msk)
        grd <- reclassify(occBuff,c(minValue(occBuff),maxValue(occBuff),1))
        grd <- writeRaster(grd,paste(spOutFolder,"/hsamples-buffer.tif",sep=""))
      }
    }
	}
    rm(hOcc)
		rm(occ)
    
    if (file.exists(paste(spOutFolder, "/hsamples-buffer.tif",sep=""))) {
      grd <- raster(paste(spOutFolder,"/hsamples-buffer.tif",sep=""))
      grd <- grd * mskArea
      areaHB <- sum(grd[which(grd[] != 0)])
    } else {
      areaHB <- 0
    }
	  
	#Size of the germplasm samples CA50
	cat("Size of the g-samples buffer \n")
	gOcc <- allOcc[which(allOcc$G == 1),]
	if (nrow(gOcc) != 0) {
		xy <- cbind(gOcc$lon, gOcc$lat)
		occ <- SpatialPoints(xy)
		
		cat("Loading", spID, "native areas \n")
		narea = paste(naDir, "/", spID, "/narea.shp", sep="")
		if(!file.exists(narea)){
		  gOcc <- as.data.frame(cbind(as.character(gOcc$Taxon), gOcc$lon, gOcc$lat))
		  names(gOcc) <- c("taxon", "lon", "lat")
		  write.csv(gOcc, paste(spOutFolder, "/gsamples.csv", sep=""), quote=F, row.names=F)
      
		  buff <- gBuffer(occ, width=0.5,byid=T)
		  occBuff <- rasterize(buff,msk)
		  grd <- reclassify(occBuff,c(minValue(occBuff),maxValue(occBuff),1))
		  grd <- writeRaster(grd,paste(spOutFolder,"/gsamples-buffer.tif",sep=""))  
      
		}else{
		  narea = readShapeSpatial(narea)
		  cat ("Projecting files \n")
		  proj4string(occ) = CRS("+proj=longlat +datum=WGS84")
		  proj4string(narea) = CRS("+proj=longlat +datum=WGS84")
		  cat("Selecting occurrences within native area \n")
		  x <- over(narea, occ)
		  x <- sum(x, na.rm=T)
		  if(x==0){
		    cat("No points within native area \n")
		  } else {
		    occBuff <- occ
		    occ = occ[narea]
		    occ = as.data.frame(occ)
		    names(occ) = c("lon","lat")
		    occ["taxon"] <- spID
		    write.csv(occ, paste(spOutFolder, "/gsamples.csv", sep=""), quote=F, row.names=F)

		    occBuff <- SpatialPoints(occBuff)
		    buff <- gBuffer(occBuff, width=0.5,byid=T)
		    occBuff <- rasterize(buff,msk)
		    grd <- reclassify(occBuff,c(minValue(occBuff),maxValue(occBuff),1))
		    grd <- writeRaster(grd,paste(spOutFolder,"/gsamples-buffer.tif",sep=""))
        
		  }
		}
	}
		rm(gOcc)
		rm(occ)
    
		if (file.exists(paste(spOutFolder,"/gsamples-buffer.tif",sep=""))) {
		  grd <- raster(paste(spOutFolder,"/gsamples-buffer.tif",sep=""))
		  grd <- grd * mskArea
		  areaGB <- sum(grd[which(grd[] != 0)])
		} else {
		  areaGB <- 0
		}
  
	#Size of the DR
	spFolder <- paste(bdir, "/maxent_modeling/models/", spID, sep="")
	projFolder <- spFolder
	spList <- read.csv(paste(bdir, "/summary-files/taxaForRichness.csv", sep=""))
	isValid <- spList$IS_VALID[which(spList$TAXON == paste(spID))]
	
	if (isValid == 1){
	  cat("Reading raster files \n")
	  grd <- paste(projFolder, "_worldclim2_5_EMN_PA.tif", sep="")
	  grd <- raster(grd)
    
	  cat("Size of the DR \n")
	  grd <- grd * mskArea
	  areaDR <- sum(grd[which(grd[] != 0)])
	  rm(grd) 
	  
	} else {
	  areaDR <- NA
	}
	
	outDF <- data.frame(DRSize=areaDR, CHSize=areaCH, NASize=areaNA, HBSize=areaHB, GBSize=areaGB)
	write.csv(outDF, paste(spOutFolder, "/areas.csv", sep=""), quote=F, row.names=F)
# 	return(outDF)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~
# LOOPS
# ~~~~~~~~~~~~~~~~~~~~~~~~
inputDir <- crop_dir

spList <- list.files(paste(inputDir, "/occurrence_files", sep=""),pattern=".csv")
for(sp in spList){
  sp <- unlist(strsplit(sp, ".", fixed=T))[1]
  out <- sizeDR(inputDir, crop, sp)
}

