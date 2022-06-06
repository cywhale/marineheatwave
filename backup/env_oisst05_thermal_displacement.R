# Thermal displacement by marine heatwaves By Jacox et al. 2020 https://www.nature.com/articles/s41586-020-2534-z
# ref: https://github.com/mjacox/Thermal_Displacement/blob/master/thermal_displacement.m
library(stars)
library(data.table)
library(magrittr)
library(abind)
library(dplyr)
library(geosphere)
library(sp)
library(sf)
library(rgdal)
library(maptools)
library(spatstat)
library(MASS)
library(rgeos)
library(raster)
library(ggplot2)
library(ggrepel)
library(viridis)
library(rnaturalearth)
library(grid)
library(gridExtra)
library(geojsonsf)
library(geojsonio)
library(odbapi) #use intersect_polyx

sp_coast <- ne_coastline() ## for ppp analysis, to check result
# for 0-360 degree, not break whole contour, so require a continuous world...
# https://gis.stackexchange.com/questions/266535/change-a-raster-from-longitude-display-180-180-to-0-360
sx <- rasterize(sp_coast, raster())
x1 <- crop(sx, extent(-180, 0, -90, 90))
x2 <- crop(sx, extent(0, 180, -90, 90))   
x3 <- copy(x1)
extent(x1) <- c(180, 360, -90, 90)
mcoa <- merge(x3, x1, x2)


# == function to output required geojson layers by analysing point patterns of thermal displacement from heatwavs ==
# code ref: https://www.r-bloggers.com/aggregating-spatial-points-by-clusters/
ppp_estimator <- function(yr, mon, tdf=NULL, heatwave=NULL, coast=NULL, xPlot=FALSE, xGjson=TRUE) {
  jt <- as.integer(mon)
  mon <- fifelse(jt<10, paste0("0",jt), paste0(jt))
  if (is.null(tdf)) {
    tdf <- paste0("../data_src/oisst/monthly_td/", yr, mon, "_td.csv")
  }
  if (is.character(tdf)) {
    if (!file.exists(tdf)) {
      print(paste0("Fail to output Year-month: ", yr," - ", mon, " because input file: ", tdf, " not exist!!"))
      return(NULL)    
    }
     zt <- fread(tdf) %>% .[td_flag==TRUE, .(lng, lat, sst, anom, xd, yd, sst_dis, ddis)] #"../data_src/oisst/monthly_td/201908_td.csv"
  } else if (is.data.frame(tdf)) {
    zt <- setDT(tdf)[td_flag==TRUE, .(lng, lat, sst, anom, xd, yd, sst_dis, ddis)]
  } else {
    print(paste0("Fail to output Year-month: ", yr," - ", mon, " because input file or data not exist!!"))
    return(NULL)    
  }
  
  zt[,`:=`(rowid=.I)]
  zcm <- zt[lng<0, .(lng, lat, sst, rowid)]
  zc1 <- zt[,.(lng, lat, sst, anom, rowid)] 
  zc <- st_as_sf(rbindlist(list(zc1[,.(lng,lat,sst,rowid)] %>% .[,lng:=fifelse(lng<0, lng+360, lng)], zcm), use.names = T), coords=c("lng", "lat"))
  #st_crs(zc) <- 4326
  ptc<- as.ppp(zc)
  Window(ptc) <- owin(xrange=c(-180,360), yrange=c(-90,90))
  kzc <- density(ptc, sigma=15,  kernel = "quartic", adjust = 0.2, cut=6) #, ext=extent(x))
  kzct <- copy(kzc)
  kzct[kzct<5] <- NA
  
  zhm <- zt[xd<0, .(xd, yd, sst_dis, rowid)] %>% setnames(1:3, c("lng","lat","sst"))
  zh1 <- zt[,.(xd, yd, sst_dis, ddis,rowid)] %>% setnames(1:3, c("lng","lat","sst"))
  zh <- st_as_sf(rbindlist(list(zh1[,.(lng,lat,sst,rowid)] %>% .[,lng:=fifelse(lng<0, lng+360, lng)], zhm), use.names = T), coords=c("lng", "lat"))
  pth<- as.ppp(zh)
  Window(pth) <- owin(xrange=c(-180,360), yrange=c(-90,90))
  kzh <- density(pth, sigma=15,  kernel = "quartic", adjust = 0.2, cut=6)
  kzht <- copy(kzh)
  kzht[kzht<5] <- NA

  if (xPlot) {
    plot(kzc, main=NULL, las=1, box=FALSE,axes=FALSE, ribbon =FALSE, useRaster=FALSE, xlim=c(-180,360), ylim=c(-90,90))
    contour(kzct, nlevels = 6, col="red", add=T) 
    contour(kzht, nlevels = 6, col="green",add=T) 
    if (!is.null(coast)) plot(coast, col="grey", add=T, legend=F, axes=F, box=F)
    abline(h=0, col="lightgrey", lty=2, lwd=1)
  }
  #--------------------------------------------------------------------
  ## thermal displacement, need extend -180 - 180 - 360 to do raster analysis

  gzc <- as(kzct, "SpatialGridDataFrame")  # convert to spatial grid class
  igzc<- as.image.SpatialGridDataFrame(gzc)  # convert again to an image
  cgzc <- contourLines(igzc, nlevels = 6)  # create contour object - change 8 for more/fewer levels
  sldc <- ContourLines2SLDF(cgzc, #proj4string= CRS("+init=epsg:4326")) #got the same as CRS(+proj=...)
                            CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  # convert to SpatialLinesDataFrame
  
  gzh <- as(kzht, "SpatialGridDataFrame")  # convert to spatial grid class
  igzh<- as.image.SpatialGridDataFrame(gzh)  # convert again to an image
  cgzh <- contourLines(igzh, nlevels = 6)  # create contour object - change 8 for more/fewer levels
  sldh <- ContourLines2SLDF(cgzh, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  # convert to SpatialLinesDataFrame

  Polyc <- gPolygonize(sldc[-1, ])
  gac <- gArea(Polyc, byid = T)/10000
  Polyc <- SpatialPolygonsDataFrame(Polyc, data = data.frame(gac), match.ID = F)
  sf1 <- st_as_sf(Polyc) %>% st_union()
  
  Polyc2 <- gPolygonize(sldh)
  gah <- gArea(Polyc2, byid = T)/10000
  Polyc2 <- SpatialPolygonsDataFrame(Polyc2, data = data.frame(gah), match.ID = F)
  sf2 <- st_as_sf(Polyc2) %>% st_union()

  if (xPlot) {
    plot(sldc, col = heat.colors(8), main=NULL, las=1, axes=FALSE, xlim=c(-180,360), ylim=c(-90,90)) #add=T
    plot(sldh, col = terrain.colors(8), add=T)
    plot(sf1, col="#F8766D80", main=NULL, las=1, axes=FALSE, xlim=c(-180,360), ylim=c(-90,90), add=T)
    plot(sf2, col="#00BFC480", add=T) #Polyc2
    #plot(cen1, col="red", pch=10, add=T)
    #plot(cen2, col="green", pch=7, add=T)
    #abline(h=0, col="lightgrey", lty=2, lwd=1)
  }

  sf1x <- st_cast(sf1, "POLYGON")
  sf2x <- st_cast(sf2, "POLYGON")
  
  ## Use sf, geojson output
  xrng <- c(-90, 270) #c(0.360)
  sflc <- st_as_sf(sldc, crs=4326) %>% st_crop(xmin=xrng[1], ymin=-90, xmax=xrng[2], ymax=90) %>% ## clip [0, 360] range
    # sflc$geometry <- (st_geometry(sflc) + c(360,90)) %% c(360) - c(0,90) 
    sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  st_crs(sflc) <- 4326
  
  sflh <- st_as_sf(sldh, crs=4326) %>% st_crop(xmin=xrng[1], ymin=-90, xmax=xrng[2], ymax=90) %>%
    sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  st_crs(sflh) <- 4326
  
  sf1y <- sf1x %>% #st_crop(xmin=xrng[1], ymin=-90, xmax=xrng[2], ymax=90) %>%
    sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  sf2y <- sf2x %>% #st_crop(xmin=xrng[1], ymin=-90, xmax=xrng[2], ymax=90) %>%
    sf::st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))

  cen1 <- st_centroid(sf1x)
  cen2 <- st_centroid(sf2x)
  # nrx <- st_nearest_feature(cen1, cen2) ## cannot solve bug GEOS has differ version in Linking (3.7.1) & runtime (3.4.2)
  nrx <- st_distance(cen1, cen2) 
  nmin<- apply(nrx, 1, which.min)
  mmin0 <- apply(nrx, 2, which.min)
  mmin<- mmin0
  mmin[nmin]<- seq_along(nmin) # if cen1(n) -> cen2(m) can multiple-to-1 but should not 1-to-multiple
  if (any(duplicated(mmin))) {
    munq <- unique(sort(c(which(!duplicated(mmin)), nmin[unique(mmin[which(duplicated(mmin))])]))) # select unique mmin, Note: mmin content is n, but munq is index of mmin
  } else {
    munq <- seq_along(mmin)
  } # still will have duplicated(mmin[munq])

  #zht <- copy(zh1)[,lng:=fifelse(lng<0, lng+360, lng)] #Nope, sf2x extend from -180 to 360 not only (0,360)
  mply <- odbapi::intersect_polyx(zh1, sf2y[munq,], xy=c("lng", "lat"), polyid=NA, onlyInPoly=TRUE, crs=4326)
  mply[,gid:=munq[gid]][,overlap:=NULL]
  mmed <- mply[!is.na(gid), .(dmed=median(ddis)), by=.(gid)][,`:=`(nid=mmin[gid])]
  jy <- mmin[munq]
  if (any(duplicated(jy))) {
    jm <- c(which(duplicated(jy)), which(duplicated(jy, fromLast = TRUE))) #index to mmin[munq] which is duplicated in its content (n)
    selm <- mmed[nid %in% jy[jm], .(selm=gid[which.min(dmed)]), by=.(nid)]$selm
## if need furthur check some poly in sf1x is missed after delete duplicated mminp[] (prevent sf1x -> sf2x become 1-to-multiple)  
    #currm <- mmed[gid %in% sort(c(munq[seq_along(jy)[-jm]], selm)),]$gid
    #currn <- mmed[gid %in% sort(c(munq[seq_along(jy)[-jm]], selm)),]$nid
    #missn <- seq_along(nmin)[-currn] %in% mmin0[-currm]
    #if (any(missn)) {
    #  nt <- seq_along(nmin)[-currn][missn==TRUE]
    #  mt <- seq_along(mmin0)[-currm]
    #  idxmt <- which(mmin0[mt] %in% nt)
    #  munq2<- mt[idxmt]
    #  mply2 <- odbapi::intersect_polyx(zh1, sf2y[munq2,], xy=c("lng", "lat"), polyid=NA, onlyInPoly=TRUE, crs=4326)
    #... not yet } 
    mmed <- mmed[gid %in% sort(c(munq[seq_along(jy)[-jm]], selm)),]
  }
  nunq <- sort(unique(mmed$nid))
  nply <- odbapi::intersect_polyx(zc1, sf1y[nunq,], xy=c("lng", "lat"), polyid=NA, onlyInPoly=TRUE, crs=4326)
  nply[,`:=`(nid=nunq[gid])][,`:=`(overlap=NULL, gid=NULL)] 
  
  setnames(mply, 1:3, c("xd","yd","sst_dis"))
  mnx <- merge(mmed, merge(mply, nply, by=c("rowid"), all.x=T) %>% .[!is.na(nid),], # & sst>sst_dis,] <- MUST BE
               by=c("gid","nid"), all=T) %>% .[!is.na(dmed),] %>% .[ ## just find an example in each (gid, nid)
    ,.(lng=lng[which.max(ddis)], lat=lat[which.max(ddis)],
       sst=sst[which.max(ddis)], anom=anom[which.max(ddis)],
       xd=xd[which.max(ddis)], yd=yd[which.max(ddis)],
       sst_dis=sst_dis[which.max(ddis)], ddis=ddis[which.max(ddis)]), by=.(gid, nid, dmed)] 
  mnx[,`:=`(xm1=fifelse(sign(lng)*sign(xd) == -1 & (abs(lng)>150 | abs(xd)>150), fifelse(sign(lng) == -1, -180.0, 180.0), NA_real_), ym1=0.5*(lat+yd))][
      ,`:=`(xm2=fifelse(sign(lng)*sign(xd) == -1 & (abs(lng)>150 | abs(xd)>150), fifelse(sign(lng) == -1, 180.0, -180.0), NA_real_), ym2=ym1)]  
  
  # length(unique(mply$gid)) == length(munq) ## check it 
  # if need to plot displacement of centroid of two poly sf1x, sf2x
  #nidx <- mned$nid
  #sf1n <- copy(sf1x)
  #kn <- lapply(nidx, function(n, sfx) {
  #  jy <- mmin[munq]
  #  jn <- munq[which(jy==n)]
  #  return(st_difference(sfx[n], st_union(sf2x[jn])))
  #}, sfx=sf1n)
  #for (i in seq_along(nidx)) { sf1n[nidx[i]] <- kn[[i]]}
  #crd1 <- st_coordinates(st_centroid(sf1n))
  #crd2 <- st_coordinates(cen2)
  #mmed[
  #  ,`:=`(x1=fifelse(crd1[nid,1]>180, crd1[nid,1]-360, crd1[nid,1]), y1=crd1[nid,2], 
  #        x2=fifelse(crd2[gid,1]>180, crd2[gid,1]-360, crd2[gid,1]), y2=crd2[gid,2])][
  #  ,`:=`(xm1=fifelse(sign(x1)*sign(x2) == -1 & (abs(x1)>150 | abs(x2)>150), fifelse(sign(x1) == -1, -180.0, 180.0), NA_real_), ym1=0.5*(y1+y2))][
  #  ,`:=`(xm2=fifelse(sign(x1)*sign(x2) == -1 & (abs(x1)>150 | abs(x2)>150), fifelse(sign(x1) == -1, 180.0, -180.0), NA_real_), ym2=ym1)]  

  sfl <- copy(sflc)
  sfl$level <- rev(as.integer(as.character(sflc$level))) - min(as.integer(as.character(sflh$level)))
  sfl2 <- copy(sflh)
  sfl2$level <- as.integer(as.character(sfl2$level))
  sfl <- rbind(sfl, sfl2)

  if (xPlot) {
    if (is.null(heatwave)) {
      heatwave <- read_stars(paste0("../data_src/oisst/monthly_heatwave/", yr, mon, "_heatwave.nc"))
      names(heatwave)[1] <- "heatwave"      
    }
    ht <- as.data.table(heatwave) %>% 
      .[,`:=`(longitude=fifelse(x>180, x-360, x), latitude=y)] %>% .[,.(longitude, latitude, heatwave)]

    g1 <- ggplot() +  
      #geom_tile(data=tt, aes_string(x="longitude", y="latitude", fill="heatwave"), alpha=0.8) + 
      #scale_fill_viridis() +
      geom_point(data=ht[heatwave==1,], aes_string(x="longitude", y="latitude"), alpha=0.65, col="gold", size = 0.01, stroke = 0, shape = 16) +
      geom_point(data=zc1, aes(x=lng, y=lat), alpha=0.75, col="orange", size = 0.01, stroke = 0, shape = 16) + 
      geom_point(data=zh1, aes(x=lng, y=lat), alpha=0.75, col="purple", size = 0.01, stroke = 0, shape = 16) + 
      geom_sf(data = ne_coastline(scale = "large", returnclass = "sf"), color = 'darkgray', size = .3) +
      geom_sf(data = sfl, aes(color=factor(level))) +
      scale_colour_brewer(palette = "RdYlGn") +
      geom_sf(data = sf1y, fill="#F8766D80") + 
      geom_sf(data = sf2y, fill="#00BFC480") +
      guides(fill=FALSE, color=FALSE) + coord_sf() + xlim(c(-180, 180)) + ylim(c(-90, 90)) +
      #coord_sf() + xlim(c(-160, -140)) + ylim(c(46, 62))  
      labs(x="Longitude",y="Latitude") +
      theme(
        panel.background = element_rect(fill="white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.75),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0),#,face = "italic"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),#, #"white"),
        legend.position = NULL
      )
    
    #if (any(!is.na(mmed$xm1))) {
    #  g1 <- g1 +
    #    geom_segment(data=mmed[!is.na(xm1),], aes(x=x1, xend=xm1, y=y1, yend=ym1), size = 0.75, color="purple") + #, color=factor(seamask)
    #    geom_segment(data=mmed[!is.na(xm1),], aes(x=xm2, xend=x2, y=ym2, yend=y2), size = 0.75, arrow = arrow(length = unit(0.08, "cm")), color="purple") + 
    #    geom_segment(data=mmed[is.na(xm1),], aes(x=x1, xend=x2, y=y1, yend=y2), size = 0.75, arrow = arrow(length = unit(0.08, "cm")), color="purple") + 
    #    geom_text_repel(data=mmed[!is.na(xm1),], aes(x=x2, y=y2, label=paste0(round(dmed,0), " km")), color = "black", size=3, hjust=1, vjust=1) + #geom_text
    #    geom_text_repel(data=mmed[is.na(xm1),], aes(x=0.5*(x1+x2), y=0.5*(y1+y2), label=paste0(round(dmed,0), " km")), color = "black", size=3, hjust=1, vjust=1)
    #} else {
    #  g1 <- g1 +
    #    geom_segment(data=mmed, aes(x=x1, xend=x2, y=y1, yend=y2), size = 0.75, arrow = arrow(length = unit(0.08, "cm")), color="purple") + 
    #    geom_text_repel(data=mmed, aes(x=0.5*(x1+x2), y=0.5*(y1+y2), label=paste0(round(dmed,0), " km")), color = "black", size=3, hjust=1, vjust=1) 
    #}
    if (any(!is.na(mnx$xm1))) {
      g1 <- g1 +
        geom_segment(data=mnx[!is.na(xm1),], aes(x=lng, xend=xm1, y=lat, yend=ym1), size = 0.75, color="purple") + #, color=factor(seamask)
        geom_segment(data=mnx[!is.na(xm1),], aes(x=xm2, xend=xd, y=ym2, yend=yd), size = 0.75, arrow = arrow(length = unit(0.08, "cm")), color="purple") + 
        geom_segment(data=mnx[is.na(xm1),], aes(x=lng, xend=xd, y=lat, yend=yd), size = 0.75, arrow = arrow(length = unit(0.08, "cm")), color="purple") + 
        geom_text_repel(data=mnx[!is.na(xm1),], aes(x=xd, y=yd, label=paste0(round(ddis,0), " km")), color = "black", size=3, hjust=1, vjust=1) + #geom_text
        geom_text_repel(data=mnx[is.na(xm1),], aes(x=0.5*(lng+xd), y=0.5*(lat+yd), label=paste0(round(ddis,0), " km")), color = "black", size=3, hjust=1, vjust=1)
    } else {
      g1 <- g1 +
        geom_segment(data=mnx, aes(x=lng, xend=xd, y=lat, yend=yd), size = 0.75, arrow = arrow(length = unit(0.08, "cm")), color="purple") + 
        geom_text_repel(data=mnx, aes(x=0.5*(lng+xd), y=0.5*(lat+yd), label=paste0(round(ddis,0), " km")), color = "black", size=3, hjust=1, vjust=1) 
    }
    #dev.off()
    #g1
  }  
  # ==== Output as geojson
  if (xGjson) {
    outdir <- "../data_src/oisst/thermal_displace/"
    #pthw <- merge(tt[heatwave==1,], copy(zc1) %>% setnames(1:2, c("longitude", "latitude")), by= c("longitude", "latitude"), all.x=TRUE) %>%
    #  merge(copy(zh1) %>% setnames(1:3, c("longitude", "latitude", "sst_dis")), by= c("longitude", "latitude"), all.x=TRUE)
    #pthw[, heatwave:=NULL] ## MUST 1
    sfc1 <- st_as_sf(zc1[,.(sst, anom, lng, lat)], coords=c("lng", "lat"), crs=4326) ## add sst anom too big for geojson 15Mb
    jfc1 <- sf_geojson(sfc1)
    geojson_write(jfc1, file = paste0(outdir, yr, mon, "_hwv_occ.geojson"))
    
    sfc2 <- st_as_sf(zh1[,.(sst, ddis, lng, lat)], coords=c("lng", "lat"), crs=4326) ## add sst anom too big for geojson 15Mb
    jfc2 <- sf_geojson(sfc2)
    geojson_write(jfc2, file = paste0(outdir, yr, mon, "_hwv_dis.geojson"))

    geojson_write(sf_geojson(sflc), file = paste0(outdir, yr, mon, "_hwv_occ_contour.geojson"))
    geojson_write(sf_geojson(sflh), file = paste0(outdir, yr, mon, "_hwv_dis_contour.geojson"))
    
    sf1d <- st_as_sf(data.frame(ID=seq_along(sf1y), geometry=st_geometry(sf1y)), crs=4326)
    sf2d <- st_as_sf(data.frame(ID=seq_along(sf2y), geometry=st_geometry(sf2y)), crs=4326)
    geojson_write(sf1d, file = paste0(outdir, yr, mon, "_hwv_occ_poly.geojson"))
    geojson_write(sf2d, file = paste0(outdir, yr, mon, "_hwv_dis_poly.geojson"))
    
    mnxf <- st_as_sf(mnx, coords=c("lng", "lat"), crs=4326)
    geojson_write(mnxf, file = paste0(outdir, yr, mon, "_hwv_dis_correspond.geojson"))
  }
  print(paste0("Successfully output Year-month: ", yr," - ", mon, " geojson for thermal displacement"))
  if (xPlot) {
    return(g1)
  }
}
# g1 <- ppp_estimator(yr=2019, mon=9, coast=mcoa, xPlot=TRUE, xGjson=TRUE)
# test: tdf, heatwave can be assigned automatically, just leave NULL, or manually assign as
# heatwave <- read_stars("../data_src/oisst/monthly_heatwave/201912_heatwave.nc")
# names(heatwave)[1] <- "heatwave" 
# ppp_estimator(yr=2019, mon=12, tdf="../data_src/oisst/monthly_td/201912_td.csv", heatwave=heatwave, coast=mcoa, xPlot=TRUE, xGjson=TRUE)

# arbitary get one file to get lon-lat matrix
dt <- read_stars("../data_src/oisst/monthly_heatwave/198201_heatwave.nc") %>% as.data.table() %>%
  .[,`:=`(lng=fifelse(x>180, x-360, x), lat=y)] %>% .[,.(lng, lat)]

latx0 <- sort(unique(dt$lat))
lngx0 <- sort(unique(dt$lng))

Detrend <- TRUE ## selece result in monthly_heatwave_detrend
RE <- 6371 #km of Earth radius
res <- lngx0[2]-lngx0[1]

bbox= c(115.0, 20.0, 135.0, 35.0) #c(-180.0, -90.0, 180.0, 90.0) #first try a smaller bounding box
grd_x <- res
grd_y <- res
# hridx <- data.table(lngidx= c(which(round(lngx0,2)<=(bbox[1]+grd_x/2)) %>% .[length(.)], which(round(lngx0,2)>=(bbox[3]-grd_x/2))[1]), 
#                     latidx= c(which(round(latx0,2)<=(bbox[2]+grd_y/2)) %>% .[length(.)], which(round(latx0,2)>(bbox[4]-grd_y/2))[1]))
# dt <- CJ(hridx$lngidx[1]:hridx$lngidx[2], 
#           hridx$latidx[1]:hridx$latidx[2]) %>% setnames(1:2, c("lngx", "latx")) %>%
#   .[,`:=`(lng=lngx0[lngx], lat=latx0[latx])]
## Just test
# lngx1 <- c(-359.75, -359.5, -270.25, -270, -269.75, -180.25, -180, -179.75, 
#            -90.25, -90, -89.75, -45.25, -45, -44.75, -30.25, -30, -29.75, -1, -0.75, -0.5, -0.25, 
#            0, 0,25, 0.5, 0.75, 1, 29.75, 30, 30.25, 44.75, 45, 45.25, 89.75, 90, 90.25, 
#            179.75, 180, 180.25, 269.75, 270, 270.25, 359.5, 359.75)
# latx1 <- c(-89.875, -89.625, -45.125, -44.875, -30.125, -29.875, -1.125, -0.875, -0.625, -0.375, -0.125, 
#            0.125, 0.375, 0.625, 0.875, 1.125, 29.875, 30.125, 44.875, 45.125, 89.625, 89.875)

lngx1 <- seq(as.integer(-360*1000+grd_x*1000), 
             as.integer(360*1000-grd_x*1000), by=as.integer(grd_x*1000))/1000.0
latx1 <- rev(seq(as.integer(-90*1000+(grd_x/2)*1000), 
                 as.integer(90*1000), by=as.integer(grd_x*1000))/1000.0)

lngxt <- sort(fifelse(lngx0<0, lngx0+360, lngx0))
latxt <- rev(latx0) #in NetCDF, lat is from 89.875 to -89.875

## NOTE: very large dm
## print(object.size(dm), units="Gb") ## 11.1 Gb ## consider generate later in parallel code

# dm <- array(rep(NA_real_, length(latx1)*length(latx1)*length(lngx1)), dim=c(length(latx1),length(latx1),length(lngx1)))
# for (i in seq_along(latx1)) {
#  dm[i,,] <- Re(RE*acos(sin(latx1[i]*pi/180)*sin(latx1*pi/180) + 
#                        cos(latx1[i]*pi/180)*cos(latx1*pi/180) %*% t(cos(lngx1*pi/180))))
#}
## Just check ##
# using marked example: dim(dm) == c(22,44); tt <- dm[1,,] #i.e i=1
# tt[,1] == tt[,21] #i.e -359.75 (two longitude difference, ldif) equal ldif=-0.25, and
# tt[,2] == tt[,20]; #tt[]
# distGeo(c(-0.125, -89.875), c(-179.875, -0.125))/1000 # 10002.11
# tt[11, 37] #10007.54
# distGeo(c(-0.125, -89.875), c(-179.875, 89.875))/1000 # 20003.87
# tt[22, 37] #20015.03

# apply_oisst_masks code ref: https://github.com/mjacox/Thermal_Displacement/blob/master/apply_oisst_masks.m
# used data from env_oisst02_landsea.R
dt <- fread("../data_src/oisst/sea_icemask_025d.csv")
setnames(dt, 1:2, c("lng", "y"))
dt[,x:=fifelse(lng<0, lng+360, lng)]
setorder(dt, -y, x)
dt[,ry:=rowid(x)] ## NOTE that in OISST.nc file, y is from Northest (89.875) to Southest (-89.875)
#################### Actually, x is 0-360 in nc so is 0-180 -180-0 arranged 
mask <- dcast(dt, x ~ ry, value.var = "seamask")[,-1] %>% as.matrix()
mlon <- dcast(dt, x ~ ry, value.var = "lng")[,-1] %>% as.matrix()
mlat <- dcast(dt, x ~ ry, value.var = "y")[,-1] %>% as.matrix()

apply_oisst_masks <- function(ii, jj, d_mask) { #, mask, mlat, mlon) {
# ====================================================
# Inputs:
#       ii: longitude index on OISST grid
#       jj: latitude index on OISST grid
#       d:  matrix of distances to all other OISST grid cells
#           from point [ii, jj]. Dimensions are [lon lat]
#       mask: mask defining regions, created by make_oisst_masks.m
#             Dimensions are [lon lat]
#       lat:  Matrix of OISST latitudes
#       lon:  Matrix of OISST longitudes
# Output:
#       d_mask: matrix of distances with unavailable locations masked out
  # d_mask <- d
  # Exclude ice-surrounded areas
  d_mask[mask==3] <- NA_real_
  if (!is.na(mask[ii, jj])) {
    maskij <- letters[as.integer(mask[ii, jj])]
    # Handle regional cases
    switch(maskij,
           d = { d_mask[mask != 4] <- NA_real_ }, # maskij ==  4
           e = { d_mask[mask != 5] <- NA_real_ }, # maskij ==  5
           f = {
             d_mask[mask<=5 | mask==7 | mask==8 | mlat>48 | (mlat>43 & mlon>351)] <- NA_real_
             if (mlon[ii,jj]>12 & mlon[ii,jj]<20 & mlat[ii,jj]>42.3 & mlat[ii,jj]<46) {
               d_mask[mask!=6] <- NA_real_
             }
           }, # maskij ==  6
           g = { d_mask[!(mask==7 | mask==11)] <- NA_real_ },  # maskij ==  7
           h = { d_mask[!(mask==8 | mask==9)] <- NA_real_ },   # maskij ==  8
           i = { d_mask[!(mask==7 | mask==8 | mask==9 | mask==11)] <- NA_real_ }, # maskij ==  9
           j = { d_mask[!(mask==10 | mask==11)] <- NA_real_ }, # maskij ==  10
           k = { d_mask[(mask>=4 & mask<=6) | mask==12] <- NA_real_ }, # maskij ==  11
           l = { d_mask[mask>=9 & mask<=11] <- NA_real_ }, # maskij ==  12
           m = { d_mask[!(mask==13 | mask==14)] <- NA_real_ }, # maskij ==  13
           n = { d_mask[(mask>=15 & mask<=17) | mlat>283] <- NA_real_ }, # maskij ==  14
           o = { d_mask[mask==2 | mask==13 | mask==14 | mask==17 | mlon<260 | mlon>280] <- NA_real_ }, # maskij ==  15
           p = { d_mask[mask==13 | mask==14 | mlon<260] <- NA_real_ }, #maskij ==  16
           q = { d_mask[mask==15] <- NA_real_ } # maskij ==  17
    )
  }
  return(d_mask)
}

# code ref: https://github.com/mjacox/Thermal_Displacement/blob/master/thermal_displacement.m
if (Detrend) {
  indir <- "../data_src/oisst/monthly_heatwave_detrend/"
  andir <- "../data_src/oisst/monthly_anom_detrend/"
} else {
  indir <- "../data_src/oisst/monthly_heatwave/"
  andir <- "../data_src/oisst/monthly_anom/"
}

td_max = 3000 # km
latn <- length(latx1) #dim(dm)[2] #720
dlon <- length(lngx1) #dim(dm)[3] #2879 (difference of longitude)
lonn <- dim(mask)[1] #1440

## Update 20200814: only fetch till 20200727 ([1] "All resolved BUT NOT exist: 28,29,30,31") 
yrng <- seq(1982,2020)
trackdate <- seq.Date(as.IDate("2020-07-27")-6, as.IDate("2020-07-27"), by="day")
curryr <- year(as.IDate("2020-07-27"))
currmo <- month(as.IDate("2020-07-27"))
#monstr <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

setkey(dt, x, y)
rowGrp <- 3L

library(future.apply)
plan(multisession)
options(future.globals.maxSize= 1048576000*2)

for (i in yrng) {
  mmx <- fifelse(i==curryr, currmo-1L, 12L)
  for (j in 1:mmx) {
    monj <- fifelse(j<10, paste0("0",j), paste0(j))
    #jstr <- monstr[j]
    print(paste0("Now in Year-month: ", i," - ", monj, " to calculate thermal displacement"))
    filet <- paste0("../data_src/oisst/monthly_td/", i, monj, "_td.csv")
    if (!file.exists(filet)) {
      z <- read_stars(paste0(indir, i, monj, "_heatwave.nc"))
      names(z)[1] <- "heatwave" 
      sst <- read_stars(paste0("../data_src/oisst/monthly_sst/", i, monj, "_sst.nc"))
      names(sst)[1] <- "sst" 
      anom <- read_stars(paste0(andir, i, monj, "_anom.nc"))
      names(anom)[1] <- "anom"
      zd <- merge(as.data.table(z) %>% setkey(x,y), dt[,.(x,y,lng,seamask)], by=c("x","y"), all=T)[
        heatwave==1 & (is.na(seamask) | seamask >3), .(lng, y, heatwave, seamask)][
          ,`:=`(xd=NA_real_, yd=NA_real_, ddis=NA_real_, 
                sst=NA_real_, anom=NA_real_, sst_dis=NA_real_, td_flag=NA, 
                rowgrp=cut(seq_len(.N), rowGrp, labels = FALSE))] #%>%
        #.[sample.int(nrow(.), 3000),] ## Just test performance bottelneck
      setnames(zd, 2, "lat")
      
      #print("Start get_tdx...")
      get_tdx <- function(lng, lat) {
        lngt <- fifelse(lng<0, lng+360, lng)
        ii <- which(lngxt==lngt)
        it <- which(lngx0==lng)
        jt <- which(latx0==lat)
        jj <- latn-jt+1 ## NetCDF lat in rev order
        #if (is.na(mask[ii, jj]) | mask[ii, jj] > 3) { #z[[1]][ii, jj]==1 & () ## already check in zd[]
          ## (old version) Note that dm (aboving), diff of lng order is reversed to match NetCDF, not the same of ref matlab
          ## (old version use whole dm 11GB too large)
          ## d_lat <- apply_oisst_masks(ii, jj, t(dm[jj, ,seq(it, it+lonn-1)])) #, mask, latm, lngm) #t(dm[jj, ,seq(dlon-lonn-it+2, dlon-it+1)])
          
          rlng <- seq(as.integer((0.125-lngt)*1000), 
                      as.integer((359.875-lngt)*1000), by=as.integer(grd_x*1000))/1000.0
          #rlng <- lngx1[seq(dlon-lonn-it+2, dlon-it+1)] #seq(it, it+lonn-1)
          #dm <- array(rep(NA_real_, length(latx1)*length(rlng)), dim=c(length(latx1),length(rlng)))
          dm <- t(Re(RE*acos(sin(latx1[jj]*pi/180)*sin(latx1*pi/180) + 
                             cos(latx1[jj]*pi/180)*cos(latx1*pi/180) %*% t(cos(rlng*pi/180)))))
          # print(object.size(dm), units="Mb") ## 7.9 Mb << 11Gb !!
          # print(paste0("After matrix multiplication: ", rowid))
          
          d_lat <- apply_oisst_masks(ii, jj, dm) #, mask, latm, lngm) #t(dm[jj, ,seq(dlon-lonn-it+2, dlon-it+1)])
          sst_norm <- sst[[1]][ii,jj] - anom[[1]][ii,jj]; # i.e. climatology, because sst - climatology = anomaly
          d_lat[is.na(sst[[1]]) | sst[[1]]>sst_norm | d_lat>td_max] <- NA_real_;
          idx <- which(!is.na(d_lat))
          if (any(idx)) {
            minidx <- which.min(d_lat[idx])
            xt <- idx[minidx] %% lonn
            yt <- as.integer(idx[minidx]/lonn)
            yt <- fifelse(xt==0, yt, yt+1)
            xt <- fifelse(xt==0, lonn, xt)
            xx <- lngxt[xt]
            xx <- fifelse(xx>180, xx-360, xx)
            yy <- latxt[yt]
            #zd[lng==lngx0[ii] & lat==latx0[jt], `:=`
            return(list(xd=xx, yd=yy, ddis=d_lat[idx[minidx]], 
                        sst=sst[[1]][ii,jj], anom=anom[[1]][ii,jj], sst_dis=sst[[1]][xt,yt], td_flag=TRUE)) #]
          } else {
            #zd[lng==lngx0[ii] & lat==latx0[jt], `:=`
            return(list(xd=NA_real_, yd=NA_real_, ddis=NA_real_, 
                        sst=sst[[1]][ii,jj], anom=anom[[1]][ii,jj], sst_dis=NA_real_, td_flag=FALSE)) #]
          }
        #}
      }
      
      #for (ii in seq_along(lngx0)) {
      #  for (jj in seq_along(latx0)) {
      zdx <- rbindlist(future_lapply(seq_len(rowGrp), function(grp) {
        x <- zd[rowgrp==grp, ]
        return(
          x[,c("xd", "yd", "ddis", "sst", "anom", "sst_dis", "td_flag"):=get_tdx(lng, lat), by = 1:nrow(x)] 
        )
      }), use.names = TRUE, fill = TRUE)
      #} #use rowGrp =3 with future_lapply for nrow=3000 user  system elapsed 0.742   0.033  70.476 
         #use rowGrp =4 for nrow=3000  user  system elapsed 0.867   0.048  60.149 
         #For one month data 1400x720 rows need  user   system  elapsed 21.244    1.086 3309.989 
      #zd[#heatwave==1 & (is.na(seamask) | seamask >3)
      #  , c("xd", "yd", "ddis", "sst", "anom", "sst_dis", "td_flag"):=get_tdx(lng, lat, rowid), by = 1:nrow(zd)] #by=.(lng, lat)]
      setkey(zdx, lng, lat)
      fwrite(zdx[,rowgrp:=NULL], file= filet)
      #} #for-loop is very slow
    }
    
    ppp_estimator(yr=i, mon=j, xPlot=FALSE, xGjson=TRUE)
  }
}


## Just check with plot --------------------------------------------------------------------------------------
if (FALSE) {
  zt <- fread("../data_src/oisst/monthly_td/201908_td.csv") %>%
    .[td_flag==TRUE,] 
  # .[lat>=46 & lat<=62 & lng >= -160 & lng <= -140,]
  
  zd <- zt[,.(lng,lat,sst)][,label:="A"] %>%
    list(zt[,.(xd,yd,sst_dis)][,label:="B"] %>% setnames(1:3,c("lng","lat","sst"))) %>%
    rbindlist(use.names = T, fill=T)
  
  # https://stackoverflow.com/questions/48282989/show-only-high-density-areas-with-ggplot2s-stat-density-2d
  # nzd <- zd %>% group_by(label) %>% do(Dens=kde2d(.$lng, .$lat, n=100, lims=c(c(-180,180),c(-90,90))))
  # nzd %<>% do(label=.$label, V=expand.grid(.$Dens$x,.$Dens$y), Value=c(.$Dens$z)) %>% 
  #   do(data.frame(label=.$label,x=.$V$Var1, y=.$V$Var2, Value=.$Value))
  
  ggplot(nzd, aes(x,y, z=Value, fill=label, alpha=..level..)) + stat_contour(geom="polygon")
  
  tt <- read_stars(paste0("../data_src/oisst/monthly_heatwave_detrend/", 2019, "08", "_heatwave.nc"))
  names(tt)[1] <- "heatwave"
  tt <- as.data.table(tt) %>% 
    .[,`:=`(longitude=fifelse(x>180, x-360, x), latitude=y)] %>% .[,.(longitude, latitude, heatwave)]
  
  ggplot() +  
    geom_tile(data=tt, aes_string(x="longitude", y="latitude", fill="heatwave"), alpha=0.8) + 
    scale_fill_viridis() +
    #stat_contour(data=nzd, aes(x=x, y=y, z=Value, color=label), geom="polygon", alpha=0.5) +
    #stat_contour(data=zd, aes(x=lng, y=lat, z=sst, fill=after_stat(level), color=label), geom="polygon", alpha=0.5) +
    #stat_density2d(data=zd, aes(x=lng, y=lat, fill=after_stat(level), color=label), alpha=0.6, geom="polygon") + 
    #geom_point(data=zt, aes(x=lng, y=lat), size=2) + #, color=factor(seamask)
    #geom_segment(data=zt,aes(x=lng, xend=xd, y=lat, yend=yd), size = 1, arrow = arrow(length = unit(0.08, "cm"))) + #, color=factor(seamask)
    geom_sf(data = ne_coastline(scale = "large", returnclass = "sf"), color = 'darkgray', size = .3) +
    guides(color=FALSE, fill=FALSE) + 
    # coord_sf() + xlim(c(-160, -140)) + ylim(c(46, 62))  
    coord_sf() + xlim(c(-180, 180)) + ylim(c(-90, 90))  
  
  # ------------------------------------ Other trials: export as a tiff, But x, y seems biased and rescaled???
  # gdal_translate -of GTiff -a_srs WGS84 -a_ullr -180 90 360 -90 201910_td.tiff 201910_td.gtiff
  rx <- raster("../data_src/raster_layer/201910_td.gtiff")
  rdt <-  setDT(as.data.frame(rx, xy = TRUE)) %>% setnames(3, "val") %>% .[x>=0 & val!=255,] %>% .[,`:=`(x=fifelse(x>180, x-360, x))]
  
  ggplot() +
    geom_tile(data = rdt , aes(x = x, y = y, fill = val)) + #geom_raster
    scale_fill_viridis_c() +
    geom_sf(data = ne_coastline(scale = "large", returnclass = "sf"), color = 'darkgray', size = .3) +
    guides(color=FALSE, fill=FALSE) + 
    coord_sf() + xlim(c(-180, 180)) + ylim(c(-90, 90))  
  
  #-------------------------------------------------------------------- Other trials
  # mdist <- distm(as_Spatial(zf$geometry)) ## too big
  zf <- st_as_sf(setDF(zt[,.(lng, lat, sst, anom)]), coords=c("lng", "lat"))
  pts<- as.ppp(zf)
  Window(pts) <- owin(xrange=c(-180,180), yrange=c(-90,90))
  
  #zf.km <- rescale(zf, 1000, "km")  
  #x <- list(x=c(-180, 180), y=c(-90,90))
  #extent(x)
  kzf <- density(pts, sigma=15,  kernel = "quartic", adjust = 0.2, cut=6) #, ext=extent(x))
  #kzf[kzf<1e-3] <- NA_real_ # then can plot transparent color
  plot(kzf, main=NULL, las=1)
  #contour(density(pts, adjust = 0.2), nlevels = 6) 
  #pts.km <- rescale(pts, 40075, "km")
  #kzf.km <- density(pts.km, sigma=50, kernel = "quartic", edge=T)
  kzft <- copy(kzf)
  kzft[kzft<5] <- NA
  contour(kzft, nlevels = 6, col="red", add=T) 
  
  # https://maczokni.github.io/crimemapping_textbook_bookdown/studying-spatial-point-patterns.html
  bw.diggle(pts) # sigma 0.286144 # too small
  ## bw.ppl(pts) ## not run bw.scott as well, too slow...
  # kzf2 <- density.ppp(pts, sigma = 1.5, kernel = "quartic", edge=T) #sigma = bw.diggle(pts)
  # plot(kzf2, main=NULL, las=1, box=FALSE,axes=FALSE, ribbon =FALSE, useRaster=FALSE, xlim=c(-180,180), ylim=c(-90,90))
  # plot(sp_coast, col="grey", add=T)
  # good but contour cannot compeltely be a circle, in a mess
  
  zg <- st_as_sf(setDF(zt[,.(xd, yd, sst_dis)] %>% setnames(1:3, c("lng","lat","sst"))), coords=c("lng", "lat"))
  ptg<- as.ppp(zg)
  Window(ptg) <- owin(xrange=c(-180,180), yrange=c(-90,90))
  kzg <- density(ptg, sigma=15,  kernel = "quartic", adjust = 0.2, cut=6)
  #kzg2<-density.ppp(ptg, sigma = 2, kernel = "quartic", edge=T) #sigma = bw.diggle(pts)
  kzgt <- copy(kzg)
  kzgt[kzgt<5] <- NA
  contour(kzgt, nlevels = 4, col="green",add=T) 
  abline(h=0, col="lightgrey", lty=2, lwd=1)
  plot(sp_coast, col="grey", add=T)
  
  contour(density(ptg, adjust = 0.2), nlevels = 3) 
  contour(kzf2, nlevels = 3, col="red", add=T) 
  #abline(v=0, col="lightgrey", lty=2, lwd=1)
  abline(h=0, col="lightgrey", lty=2, lwd=1)
  #spdf_world <- ne_countries()
  #plot(spdf_world, )
  
  gzf <- as(kzft, "SpatialGridDataFrame")  # convert to spatial grid class
  igzf<- as.image.SpatialGridDataFrame(gzf)  # convert again to an image
  cgzf <- contourLines(igzf, nlevels = 6)  # create contour object - change 8 for more/fewer levels
  sldf <- ContourLines2SLDF(cgzf, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  # convert to SpatialLinesDataFrame "+proj=longlat +ellps=WGS84 +datum=WGS84"
  plot(sldf, col = heat.colors(8))
  
  gzg <- as(kzgt, "SpatialGridDataFrame")  # convert to spatial grid class
  igzg<- as.image.SpatialGridDataFrame(gzg)  # convert again to an image
  cgzg <- contourLines(igzg, nlevels = 6)  # create contour object - change 8 for more/fewer levels
  sldg <- ContourLines2SLDF(cgzg, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  # convert to SpatialLinesDataFrame
  plot(sldg, col = terrain.colors(8), add=T)
  
  Polyclust <- gPolygonize(sldf[3, ])
  gas <- gArea(Polyclust, byid = T)/10000
  Polyclust <- SpatialPolygonsDataFrame(Polyclust, data = data.frame(gas), match.ID = F)
  plot(Polyclust, border="red", add=T)
  
  Polyclust2 <- gPolygonize(sldg[1, ])
  gas2 <- gArea(Polyclust2, byid = T)/10000
  Polyclust2 <- SpatialPolygonsDataFrame(Polyclust2, data = data.frame(gas2), match.ID = F)
  plot(Polyclust2, border="blue", add=T)
  
  #e <- drawExtent()
  
  zsf <- st_sf(zf, crs=4326)
  #zst <- as(zsf$geometry,"Spatial")
  #crs(zst) <- "+proj=utm +zone=19 +datum=NAD83 +units=m +no_defs"
  zst <- SpatialPointsDataFrame(st_coordinates(zf), data.frame(ID=1:nrow(zsf)), 
                                proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  cag <- aggregate(zst, by = Polyclust, FUN = length)
  plot(kzf, main=NULL, las=1, box=FALSE,axes=FALSE, ribbon =FALSE, useRaster=FALSE, xlim=c(-180,180), ylim=c(-90,90))
  plot(sldf, col = terrain.colors(8), add = T)
  plot(cag, col = "red", border = "white", add = T)
  #graphics::text(coordinates(cag) + 1000, labels = cAg$CODE)
  
  sin <- zst[cag, ]  # select the stations inside the clusters
  sout<- zst[!row.names(zst) %in% row.names(sin), ]  # stations outside the clusters
  #plot(sout) #, add=T)  # the more sparsely distributed points - notice the 'holes' of low density
  #plot(cag, border = "red", lwd = 3, add = T)
  

  ## Just test (old version)
  # lat = latx0; lon = lngx0
  # which(lon>=26 & lon<=37) #825 ..... 868
  # which(lat>=30.5 & lat<=39.5) #483.. 518
  # ii=826; jj=484
  # lngx0[ii]; latx0[jj]  #26.375 #30.875
  # lngx1[range(seq(dlon-lonn-ii+2, dlon-ii+1))]  #lngx1[c(615, 2054)] => c(-206.25  153.50) => -179.875-26.375, 179.875-26.375
  # BUT it's land, so seamask==0
  # change to c(35.875 35.375) #ii=which(lon==35.875); jj=which(lat==35.375) #864, 502 
  # d_mask <- t(dm[jj, ,seq(dlon-lonn-ii+2, dlon-ii+1)])
}

