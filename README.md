# Marine Heatwaves (MHW)

* Original reference by Jacox 2020 

  - Jocox (2020): https://www.nature.com/articles/s41586-020-2534-z
    Jacox, Michael G., et al. "Thermal displacement by marine heatwaves." Nature 584.7819 (2020): 82-86.

  - Note that in original paper by Jacox 2020: "For each grid cell we calculated time series of SST anomalies relative to 
    the 1982–2011 climatology and classified MHWs as periods with SST anomalies above a seasonally varying 90th-percentile 
    threshold (Extended Data Fig. 5). Our analysis differs from those used in some other studies in that 
    we used monthly averaged SST rather than daily data..." AND Here we do not implement "detrend" global trend in paper 

  - original code ref (in Matlab): https://github.com/mjacox/Thermal_Displacement/

  - code in R, and data use NOAA OISST v2.1: https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/ 


* Usage

  - Install R, (and Rstudio recommended) from their official website first

  - Find souce codes in R/, \*.Rmd are for Rmarkdown files with codes enclosed in \``` section \```; \*.html is the results generated from .Rmd

  - To preview .html (generated from Rmarkdown) in this site:

    > In chrome (or other browser) URL input https://htmlpreview.github.io/? + HTML URL in this repository. 
    > For example: https://htmlpreview.github.io/?https://github.com/cywhale/marineheatwave/blob/main/R/01_OISST_data.html
