if (!require("curl")) install.packages("curl"); library(curl)
if (!require("magrittr")) install.packages("magrittr"); library(magrittr)
if (!require("data.table")) install.packages("data.table"); library(data.table)
if (!require("stars")) install.packages("stars"); library(stars)
if (!require("abind")) install.packages("abind"); library(abind)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("future.apply")) install.packages("future.apply"); library(future.apply)

if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)
if (!require("ggExtra")) install.packages("ggExtra"); library(ggExtra)
if (!require("viridis")) install.packages("viridis"); library(viridis)
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library(RColorBrewer)
if (!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if (!require("rnaturalearth")) install.packages("rnaturalearth"); library(rnaturalearth)

# plot functions the same as R/01-1_OISST_monthly.Rmd (but modify plot and panel margins)
gplotx <- function(z, var="sst", returnx=FALSE, minz=NA_real_, maxz=NA_real_, no_coord=FALSE, fillopts="D", xlab=NULL, ylab=NULL, maxlim=360, legend_pos="bottom", legend_label=var, legend_direction="horizontal") {
  df <- as.data.frame(z)
  if (is.na(var) | var=="") { var = colnames(df)[3] }
  zcol <- chmatch(var, colnames(df))
  if (!length(zcol) | is.na(zcol)) { #ensym(var)
    print("Warning: change z column name")
    colnames(df)[3] <- var
    zcol <- grep(var, colnames(df))
  }  
  
  df[,zcol] <- unclass(df[,zcol])
  if (any(is.na(c(minz, maxz)))) {
    rng <- range(na.omit(df[,zcol]))
    if (is.na(minz)) {
      minz <- floor(rng[1]) - 0.5
    }
    if (is.na(maxz)) {
      maxz <- floor(rng[2]) + 1.5
    }
  }
  
  if (!is.na(maxlim) & maxlim!=360) {
    df[,"x"] <- fifelse(df[,"x"]>180, df[,"x"]-360, df[,"x"])
  }
  
  if ("topright" %in% legend_pos) {
    lt <- c(0.9,0.75)
  } else if ("bottomright" %in% legend_pos) {
    lt <- c(0.9,0.35)
  } else if ("topleft" %in% legend_pos) {
    lt <- c(0.1,0.75)
  } else if ("bottomleft" %in% legend_pos) {
    lt <- c(0.1,0.35)
  } else {
    lt <- legend_pos
  }
  
  gx <- ggplot() +  
    geom_tile(data=df, aes_string(x="x", y="y", fill=var), alpha=0.8) + 
    xlab(xlab) + ylab(ylab) 
  
  if (!is.na(maxlim) & maxlim!=360) {
    if (!no_coord) {
      gx <- gx + coord_sf(xlim=c(-180, 180), ylim=c(-90, 90))
    }
    gx <- gx + 
      scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) + 
      scale_y_continuous(limits = c(-90, 90), expand = c(0, 0))
  } else {
    if (!no_coord) {
      gx <- gx + coord_sf(xlim=c(0, 360), ylim=c(-90, 90))
    }
    gx <- gx + 
      scale_x_continuous(limits = c(0, 360), expand = c(0, 0)) + 
      scale_y_continuous(limits = c(-90, 90), expand = c(0, 0))
  } 
  if (!is.na(fillopts) & fillopts %in% LETTERS[1:5]) {
    if (any(is.na(c(minz, maxz)))) {
      gx <- gx + scale_fill_viridis(legend_label, option=fillopts) 
    } else {
      gx <- gx + scale_fill_viridis(legend_label, limits=c(minz, maxz),option=fillopts,
                                    breaks=seq(100*minz, 100*maxz, by = as.integer((maxz-minz)*100/3))/100)
    }
  } else {
    jet.colors <-
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (any(is.na(c(minz, maxz)))) {
      gx <- gx + scale_fill_gradientn(legend_label, colors = jet.colors(16))
    } else {
      gx <- gx + scale_fill_gradientn(legend_label, colors = jet.colors(16), limits=c(minz, maxz),
                                      breaks=seq(100*minz, 100*maxz, by = as.integer((maxz-minz)*100/3))/100)
    }
  }
  
  gx <- gx +
    ggExtra::removeGrid(x=TRUE, y=TRUE) +
    theme(
      panel.background = element_blank(), 
      panel.border = element_rect(colour = "black", fill=NA, size=0.75),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0,0,0,0), "cm"),
      panel.spacing = unit(0,"cm"),
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0),
      axis.text.x  = element_text(family = "sans"),
      axis.title.x = element_text(family = "sans"),
      axis.title.y = element_text(family = "sans"),
      axis.text.y = element_text(family = "sans"),
      legend.key = element_blank(),
      legend.key.size = unit(0.8,"line"),
      legend.box.background = element_blank(),
      legend.title = element_text(size=10),
      legend.text =  element_text(family = "sans", size=7),
      legend.background = element_rect(fill = "transparent", colour = "transparent"),
      legend.direction=legend_direction,
      legend.key.height = unit(fifelse(legend_direction=="horizontal",0.2,2),'cm'),
      legend.position = lt
    )
  
  if (returnx) {return(gx)}
  gx
}
