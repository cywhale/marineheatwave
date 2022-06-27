library(curl)
baseurl <- "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr"

year_month <- "198801"
date <- paste0(year_month, "01")
fileprex <- "oisst-avhrr-v02r01"
filename <- paste(fileprex, date, "nc", sep=".")
oisstfile <- paste(baseurl, year_month, filename, sep="/")
outputfile <- paste0("../data_src/oisst/daily/", filename)

print(oisstfile)

if (!file.exists(outputfile)) {
  tryCatch({
    curl::curl_download(oisstfile, destfile = outputfile)
  }, error = function(err) print(err))
} else {
  print(paste0(oisstfile, " existed!"))
}

library(stars)
sts <- read_stars(outputfile)
str(sts)

library(ggplot2)
library(viridis)
library(rnaturalearth)
library(sf)
library(data.table)
library(magrittr)

yr=2018
mon="06"
outdir <- "../data_src/oisst/thermal_displace/"
sfc1 <- st_read(paste0(outdir, yr, mon, "_hwv_occ.geojson"))

sf_use_s2(FALSE)
ggplot() +
  geom_sf(data = ne_coastline(scale = "large", returnclass = "sf"), color = 'darkgray', size = .3) +
  geom_sf(data=sfc1, aes(colour = sst), alpha=0.5, size = 0.01, shape = 16) + 
  scale_colour_viridis(option="H") + 
  #guides(fill="none", color="none") + 
  coord_sf() + 
  scale_x_continuous(limits = c(115, 155), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(10, 50), expand = c(0, 0)) +
  labs(x=NULL, y=NULL) +
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


crd <- st_coordinates(sfc1) %>% as.data.table
sfc1[which(crd$X>=115 & crd$X<=155 & crd$Y>=10 & crd$Y<=50), c("sst","geometry")] %>% 
  setDT %>% .[1:1000,]
