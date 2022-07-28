if (!require("magrittr")) install.packages("magrittr"); library(magrittr)
if (!require("data.table")) install.packages("data.table"); library(data.table)
if (!require("stars")) install.packages("stars"); library(stars)
if (!require("abind")) install.packages("abind"); library(abind)
if (!require("dplyr")) install.packages("dplyr"); library(dplyr)
if (!require("starsExtra")) install.packages("starsExtra"); library(starsExtra)

load("../data_src/stats/test_stars_predict00.RData")
seq1 <- seqx[c(3,2,1), rev(seq(1,dim(seqx)[2]))]
seq1[1, c(4,10)] <- rep(NA_real_,2)

for (j in 1:dim(seqx)[2]) {
  x <- matrix(c(seqx[,j], seq1[,j]), ncol=2) %>% starsExtra::matrix_to_stars()
  names(x)[1] <- "anom"
  if (j==1) {
    stx <- x
  } else {
    stx <- c(stx, x)
  }
}

mrt <- merge(stx) %>% st_set_dimensions(3, values = as.POSIXct(dimnames(seqx)[[2]], format="%Y_%m_%d"), names = "time") 
timex <- st_get_dimension_values(mrt, 3)

predt <- function(x) {
  if (all(is.na(x))) return(list(as.numeric(x)))
  x[!is.na(x)] <- as.numeric(stats::predict(lm(x ~ timex))) 
  return(list(as.numeric(x))) 
} 

pmrt <- st_apply(adrop(mrt), c(1,2), predt)
tk <- array(unlist(pmrt[[1]]), dim = c(length(timex), 2, dim(seqx)[1])) 

haty <- apply(seqx, MARGIN=1, FUN=function(y){ 
  if (all(is.na(y))) return(y)
  yh <- as.numeric(stats::predict(lm(y ~ seq_along(timex)))) 
  y[!is.na(y)] <- yh
  return(c(as.numeric(y)))
}) %>% as.data.frame()  

all.equal(round(haty$test,5), round(tk[,1,3],5)) #almost the same(slightly difference)
