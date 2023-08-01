# Open API of Marine Heatwaves (MHWs), ODB

Marine heatwaves (MHWs) evaluated from 0.25-degree gridded NOAA OISST v2.1

#### Reference

Jacox, Alexander, Bograd, and Scott (2020), Thermal Displacement by Marine Heatwaves, Nature, 584, 82–86, doi:10.1038/s41586-020-2534-z",

#### Swagger API doc

https://eco.odb.ntu.edu.tw/api/swagger


#### Endpoint: /api/mhw

###### One-point MHWs without time-span limitation:

- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=135&lat0=15

###### Bonding-box <= 90 * 90 in degrees: 1-year time-span limitation for fetching MHWs

- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=135&lon1&=140&lat0=15&lat1=30&start=2021-01-01

###### Bonding-box > 90 * 90 in degrees: 1-month time-span limitation for fetching MHWs

- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=-180&lon1&=180&lat0=-90&lat1=90&start=2021-01-01


###### The parameter 'append' can append more data-variables (level, sst_anomaly, td) in data output

1. level: Severity of MHWs (ref: Hobday et al., 2018) 0 - 4: None, Moderate, Strong, Severe, Extreme

2. sst: Monthly SST

3. sst_anomaly: Monthly SST Anomalies

4. td: Thermal Displacement (ref: Jacox et al., 2020)
    
- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=135&lat0=15&append=sst_anomaly,level (default is level)

- The parameters lon1, lat1, start, end, and append are optional which can be none.

#### /api/mhw/csv

- The same usage as /api/mhw, but in CSV file response

#### /api/mhw/plot

- parameters see: https://eco.odb.ntu.edu.tw/api/swagger

- Note that wider time/spatial span specified by {period} and {bbox} lead to much slower response.

- {mode} = series: Time-series e.g. /api/mhw/plot?bbox=(-90,-10,-80,0),(-150,-5,-90,5)&period=2015-202306&sstvar=sst_anomaly

- {mode} = month: Month climatology. /api/mhw/plot?bbox=(-90,-60,0,60)&sstvar=sst&period=1982-2011,2012-2021,2022,2023&mode=month


#### Other resources

###### MHW evaluation source code

- https://github.com/yeh-tc/Python_project-Thermal_Displacement

###### QGIS tutorial for MHW WMS layer

- https://github.com/oceandatabank/MHW_QGIS

###### Introduction of MHW database by ODB (e-paper in pdf)

- https://www.odb.ntu.edu.tw/wp-content/uploads/2023/05/ODB_Newsletter_21.pdf
