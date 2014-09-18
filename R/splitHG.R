# N. Castaneda, J. Ramirez, H. Achicanoy
# 2014

# This code create separate files for germplasm accessions and reference (or herbarium) sightings

splitHG <- function(occ, crop, lon, lat){

  data(wrld_simpl) #nice world map
  
  if (!file.exists("./occurrences")) {dir.create("./occurrences")}
  
  occ <- occ[which(!is.na(occ$lat)),] #select only records with coordinates
  occ <- occ[which(!is.na(occ$lon)),] #select only records with coordinates
  
  h <- occ[which(occ$H==1),]
  g <- occ[which(occ$G==1),]
  
  write.csv(h,paste("./occurrences/",crop,"_h.csv",sep=""),quote=F,row.names=F) #georeferenced herbarium records
  write.csv(g,paste("./occurrences/",crop,"_g.csv",sep=""),quote=F,row.names=F) #georeferenced germplasm records
  write.csv(occ,paste("./occurrences/",crop,".csv",sep=""),quote=F,row.names=F) #all georeferenced records

  #==map H/G densities ==#
  r <- raster()
  res(r) <- 1
  
  g_ras <- g[,c(lon,lat)]
  g_ras <- rasterize(g_ras,r,fun="count")
  
  h_ras <- h[,c(lon,lat)]
  h_ras <- rasterize(h_ras,r,fun="count")
  
  h_ras[which(h_ras[]==0)] <- NA; g_ras[which(g_ras[]==0)] <- NA
  
  brks <- unique(quantile(c(h_ras[],g_ras[]),na.rm=T,probs=c(seq(0,1,by=0.05))))
  cols <- colorRampPalette(c("dark green","yellow","orange","red"))(length(brks)-1)
  brks.lab <- round(brks,0)
  
  if (!file.exists("./figures")) {dir.create("./figures")}
  
  z <- extent(h_ras)
  aspect <- (z@ymax-z@ymin)*1.4/(z@xmax-z@xmin)
  
  #herbarium map
  tiff("./figures/h_samples_count.tif",
       res=300,pointsize=5,width=1500,height=1500*aspect,units="px",compression="lzw")
  par(mar=c(2.5,2.5,1,1),cex=0.8,lwd=0.8)
  plot(h_ras,col=cols,zlim=c(min(brks),max(brks)), main = "Reference samples (H)",
       breaks=brks,lab.breaks=brks.lab,useRaster=F,
       horizontal=T,
       legend.width=1,
       legend.shrink=0.99)
  plot(wrld_simpl,add=T,lwd=0.5, border="azure4")
  grid()
  dev.off()
  
  #germplasm map
  tiff("./figures/g_samples_count.tif",
       res=300,pointsize=5,width=1500,height=1500*aspect,units="px",compression="lzw")
  par(mar=c(2.5,2.5,1,1),cex=0.8, lwd=0.8)
  plot(g_ras,col=cols,zlim=c(min(brks),max(brks)),useRaster=F, main="Genebank accessions (G)",
       breaks=brks,lab.breaks=brks.lab,
       horizontal=T,
       legend.width=1,
       legend.shrink=0.99)
  plot(wrld_simpl,add=T,lwd=0.5, border="azure4")
  grid()
  dev.off()
}
