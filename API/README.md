# Open API of Marine Heatwaves (MHWs), ODB

Marine heatwaves (MHWs) evaluated from 0.25-degree gridded NOAA OISST v2.1.\n" +

#### Reference

Jacox, Alexander, Bograd, and Scott (2020), Thermal Displacement by Marine Heatwaves, Nature, 584, 82–86, doi:10.1038/s41586-020-2534-z",

#### Swagger API doc

https://eco.odb.ntu.edu.tw/api/swagger


#### Usage

###### One-point MHWs without time-span limitation:

- e.g. https://eco.odb.ntu.edu.tw/api/swagger?lon0=135&lat0=15

###### Bonding-box MHWs with 1-year time-span limitation:

- e.g. https://eco.odb.ntu.edu.tw/api/swagger?lon0=135&lon1&=140&lat0=15&lat1=30&start=2021-01-01

###### The parameter 'append' can append more data-variables (level, sst_anomaly, td) in data output

1. level: Severity of MHWs (ref: Hobday et al., 2018) 0 - 4: None, Moderate, Strong, Severe, Extreme

2. sst_anomaly: Monthly SST Anomalies

3. td: Thermal Displacement (ref: Jacox et al., 2020)
    
- e.g. https://eco.odb.ntu.edu.tw/api/swagger?lon0=135&lat0=15&append=sst_anomaly,level (default is level)

- The parameters lon1, lat1, start, end, and append are optional which can be none.

#### Other resource

###### MHW evaluation source code

- https://github.com/yeh-tc/Python_project-Thermal_Displacement

###### QGIS tutorial for MHW WMS layer

- https://github.com/oceandatabank/MHW_QGIS

##### Introduction of MHW database by ODB (e-paper in pdf)

- https://www.odb.ntu.edu.tw/wp-content/uploads/2023/05/ODB_Newsletter_21.pdf
