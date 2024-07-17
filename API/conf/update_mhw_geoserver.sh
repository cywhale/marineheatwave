#!/bin/bash
GEOSERVER_URL="localhost:8080/geoserver"
GEOSERVER_ADMIN='MyAdmin'
GEOSERVER_PASS='MyPass'
logg="update_geoserver.log"
echo "#------------------------------- Update MHW geoserver start" |tee -a "$logg"
echo "#------------------------------- $(date '+%Y%m%d %H:%M:%S')" |tee -a "$logg"

# API endpoint and credentials
WORKSPACE="marineheatwave"
STORE="mhw"
LAYER="mhw"
ENDPOINT="$GEOSERVER_URL/rest/workspaces/$WORKSPACE/datastores/$STORE/featuretypes/$LAYER.xml"
Modify_XML="TRUE"
Put_XML="TRUE"

# TO get the first and last days of the previous month
current_date=$(date +"%Y-%m-%d")

# Extract the year and month from the current date
year=$(date -d "$current_date" +"%Y")
month=$(date -d "$current_date" +"%-m")  # Use %-m to remove leading zeros

# Calculate the first day of the previous month
if [ "$month" -eq 1 ]; then
    prev_month_year=$((year - 1))
    prev_month=12
else
    prev_month_year="$year"
    prev_month=$((month - 1))
fi

# Add leading zero if month is a single digit
if [ ${#prev_month} -eq 1 ]; then
    prev_month="0$prev_month"
fi

PREVIOUS_MONTH_FIRST_DAY="${prev_month_year}-${prev_month}-01"

# Calculate the last day of the previous month
PREVIOUS_MONTH_LAST_DAY=$(date -d "${prev_month_year}-${prev_month}-01 + 1 month - 1 day" +"%Y-%m-%d")
echo "Update dimension to first_day $PREVIOUS_MONTH_FIRST_DAY, and end_date $PREVIOUS_MONTH_LAST_DAY by end_point $ENDPOINT" |tee -a "$logg"

# Fetch the existing XML from GeoServer
curl -u "$GEOSERVER_ADMIN:$GEOSERVER_PASS" -X GET "$ENDPOINT" -o "feature.xml"
cp feature.xml bak_xml/"feature_$current_date.xml"

if [[ "$Modify_XML" == "TRUE" ]]; then
  # Modify the XML using xmlstarlet. For instance, to update the referenceValue:
  xmlstarlet ed -L -u "/featureType/metadata/entry[@key='time']/dimensionInfo/defaultValue/referenceValue" -v "NEW_DATE_VALUE" feature.xml

  # Modify the XML using xmlstarlet. Update the referenceValue and endValue
  xmlstarlet ed -L -u "/featureType/metadata/entry[@key='time']/dimensionInfo/defaultValue/referenceValue" -v "$PREVIOUS_MONTH_FIRST_DAY" feature.xml
  xmlstarlet ed -L -u "/featureType/metadata/entry[@key='time']/dimensionInfo/endValue" -v "$PREVIOUS_MONTH_LAST_DAY" feature.xml
fi

echo "Starting curl PUT feature.xml to $ENDPOINT" |tee -a "$logg"
if [[ "$Put_XML" == "TRUE" ]]; then
  # POST the modified XML back to GeoServer
  curl --http1.1 -v -i -u "$GEOSERVER_ADMIN:$GEOSERVER_PASS" -X PUT -H "Content-Type: application/xml" -d "@feature.xml" "$ENDPOINT" >> $logg 2>&1
fi

# Clean up
rm feature.xml
echo "#------------------------------- Update End" |tee -a "$logg"
echo "#------------------------------- $(date '+%Y%m%d %H:%M:%S')" |tee -a "$logg"
