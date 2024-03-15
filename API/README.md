# Open API of Marine Heatwaves (MHWs), ODB

Marine heatwaves (MHWs) evaluated from 0.25-degree gridded NOAA OISST v2.1

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.8250030.svg)](https://doi.org/10.5281/zenodo.8250030)


#### ODB MHW open data official website

- https://eco.odb.ntu.edu.tw/pub/MHW

#### Reference

- Jacox, Alexander, Bograd, and Scott (2020), Thermal Displacement by Marine Heatwaves, Nature, 584, 82–86, doi:10.1038/s41586-020-2534-z",

#### Swagger API doc

- [ODB Marine Heatwaves (MHW) API manual/online try-out](https://api.odb.ntu.edu.tw/hub/swagger?node=odb_mhw_v1)

#### Usage: Endpoint /api/mhw

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

#### Demo on <a href="https://eco.odb.ntu.edu.tw/pub/MHW/" target="_blank">ODB Marine Heatwaves</a>

<p align="center"><img src="https://eco.odb.ntu.edu.tw/pub/MHW/assets/202401_level.jpg" width=640 alt="Marine Heatwaves level 202401" /></p>

<p align="center"><em>Global marine heatwave levels at 0.25 degree grid for January 2024.</em></p>

#### Other resources

###### MHW evaluation source code

- https://github.com/yeh-tc/Python_project-Thermal_Displacement

###### QGIS tutorial for MHW WMS layer

- https://github.com/oceandatabank/MHW_QGIS

###### Introduction of MHW database by ODB (e-paper in pdf)

- https://www.odb.ntu.edu.tw/wp-content/uploads/2023/05/ODB_Newsletter_21.pdf
