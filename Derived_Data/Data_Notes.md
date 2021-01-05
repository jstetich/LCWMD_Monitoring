# Data Notes
==============

<img
  src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
  style="position:absolute;top:10px;right:50px;" />
  
## Data Files in This Folder
1.  **Sonde_Data.csv**  --  (Omitted from Github Repository because of its size)
    Contains raw data assembled from original Excel files.  See code in 
    *Import_Data.Rmd*  for details.  
2.  **Daily_Data.csv**   -- Daily summaries (min, max, median, mean, sd, IQr, n)
    derived from *Sonde_Data.csv*. See related code in 
    *Make_Daily_Summaries.Rmd* for details.  
3.  **Exceeds_Data.csv**  -- Data derived from *Daily_Data.csv* containing flags
    that indicate whether conditions that day exceeded acute or chronic exposure
    thresholds.   Thresholds for Chlorides, DO, and Percent Saturation are
    derived from Maine water quality criteria.  Temperature thresholds are
    derived from a study of brook trout habitat use in the upper midwest.  See
    the code in *Make_Daily_Summaries.Rmd* for threshold values and other
    details.  
4.  **Full_Data.csv** -- Data derived from *Daily_data.csv* containing lags and
    weighted sums for time series analysis.  The data contains missing values
    where there are gaps in the data, so that autocorrelation-based regression
    models rely only on data where lag data are available.  See the code in
    *Make_Complete_Data.Rmd* for details.  
5.  **Site_IC_Data.csv**  and  **Site_IC_Data.xlsx** --  Simplified data of
    direct and cumulative subwatershed area and imperviousness for each Long
    Creek subwatershed / Monitoring Location.  Entered by hand from a table in a
    report by GZA to LCWMD.  
6.  **Weather Data.csv**  --  Selected Weather data (including average wind
    speeds, precipitation, snow, and temperature), along with lagged and
    weighted precipitation data (mostly lagged and weighted log of
    precipitation).  See the code in *Make_Complete_Data.Rmd* for details.  
7.  **grab_data** a folder containing several csv files.  Files contain subsets
    by subject area of data derived from the original excel data sheets. Data
    includes grab sample data only (no data logger or other "continuous" data).
    See the code in  *Import_Grab_Sample_Data.Rmd* for details.  

## Data omitted from this folder
The LCWMD sonde data, with data collected generally every fifteen minutes to
half an hour across multiple sites, amounts to a substantial amount of data. In
order to reduce the size of this Archive, we have omitted some intermediate
data.

In particular, we  omitted:

-  **Sonde_Data.csv** which is about 54 MB in size.  This is an aggregated and
   simplified version of the raw sonde data received from GZA, and included in
   the folder "Original_Data". If you need this file, you can generate it by
   downloading the archive and running "Import_Data.Rmd".

-  **extracted_grab_data.xlsm** --  file contains copies of the original
   spreadsheets, without the graphics.  it is ~ 235 MB but much of that space is
   wasted.  Much of that can be saved by converting to "long" data format.
