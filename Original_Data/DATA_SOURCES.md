# Data Sources
Data Derived from ZIP file transmitted from GZA to LCWMD and Curtis Bohlen via
on-line transfer.

Original form of the data in excel is poorly organized, and took up a lot of 
disk space, memory, and computation to access.

1.  It contains a high fraction of replicate data, and
2.  It uses inefficient lookup code to create consolidated data sheets. 

This is a substantial amount of data, which contains a fair amount of
duplication. We omitted some portions of the original data transmission from GZA
that we did not use from the GIT archive to save disk space.


## Weather Data
CBEP uses a custom Python program to download data from NOAA's online data
repositories.  Specifically, data were accessed through NOAA's National Centers
for Environmental Information.

Here, we have downloaded daily (GHCND) weather summaries via API v2. Information
on this API is available here: https://www.ncdc.noaa.gov/cdo-web/webservices/v2

Documentation on specific datasets is available at
https://www.ncdc.noaa.gov/cdo-web/datasets

Portland Jetport weather data was  downloaded by Curtis C. Bohlen using a custom
python script, titled "noaaweatherdataGUI.py" on April 23, 2020.

### Units
Data is in SI units, except that NOAA provides some data in tenths of the
nominal units.  This is not well documented through the API, but obvious in 
context. Temperatures are reported in tenths of degrees C, and precipitation in
tenths of a millimeter.  For this analysis, we disregard trace rainfall
events.


### Software Access
CBEP's python program is available in this repository. This is not a polished
program, and may or may not work "out of the box." Some functionality implied by
the User Interface is incomplete. This software is in ongoing development by
CBEP staff, but the version here will not be regularly updated.  For the most
up-to-date version, contact CBEP.

To use the program, you will need to get an access token from NOAA (see the
website describing the API, above), and modify the code slightly on line 43 to
set MYTOKEN equal to the value of the token you receive from NOAA.

