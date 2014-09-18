# N. Castaneda, J. Ramirez-Villegas, H. Achicanoy
# 2014

createMasks <- function(shpMask,env_dir){
  
  if (!file.exists("./masks")) {dir.create("./masks")}
  outdir <- "./masks"
  
  cat("Preparing the mask! \n")
  
  shpName <- shpMask
  
  #Reading the polygon shapefile
  cat("Reading and converting \n")
  pol <- readShapePoly(shpName)
  ls <- list.files(env_dir, pattern="bio")    
  rs <- raster(paste(env_dir, "/", ls[1],sep=""))
  
  pa <- rasterize(pol,rs)
  pa <- trim(pa)
  
  #     pa[which(!is.na(pa[]))] <- 1
  pa <- reclassify(pa,c(minValue(pa), maxValue(pa), 1))
  #pa[which(is.na(pa[]) & rs[] == 1)] <- 0
  
  # Prepare cellArea.asc
  rs_a <- area(pa)
  rs_a <- mask(rs_a, pa)
  
  cat("Writing outputs \n")
  writeRaster(rs_a,paste(outdir,"/cellArea.tif",sep=""), overwrite=T)
  writeRaster(pa,paste(outdir,"/mask.tif",sep=""), overwrite=T)
}