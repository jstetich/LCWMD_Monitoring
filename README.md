# LCWMD_Monitoring

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />


Analysis of ten years of water quality monitoring data from the Long Creek
Watershed Management District in Maine.

# Introduction
This archive includes R scripts and related products analyzing approximately ten
years of water quality monitoring data from the Long Creek Watershed Management
District (LCWMD).The District is a Maine quasi-municipal non-profit corporation
that manages water quality in the Long Creek watershed on behalf of about 130
land owners.

The Long Creek Watershed, a bit over two square miles in area, is dominated by
commercial land use. The Maine Mall is one of the largest landowners in the
watershed, and it is surrounded by a range of commercial businesses, from
medical offices, to car washes.  About a third of the watershed is impervious
surfaces like roads, parking lots, and rooftops.  Landowners with an acre or
more of impervious area are required to get a Clean Water Act permit for
stormwater discharges from their property.  The LCWMD provides an alternative
for landowners to attain an individual permit.  Landowners who elect
to participate in the The Long Creek Watershed Management District receive a
General Permit, in return for providing funding to the District, and
facilitating the work of the district by permitting access to their property for
certain activities.

For more information on LCWMD, see [their web site](restorelongcreek.org).

# Statement of Purpose
CBEP is committed to the ideal of open science.  Our State of the Bay data
archives ensure the science underlying the 2020 State of the Bay report is
documented and reproducible by others. The purpose of these archives is to
release raw data and data analysis code whenever possible to allow others to
review, critique, learn from, and build upon CBEP science.

# Archive Structure
CBEP 2020 State of the Bay data analysis repositories are divided into from two
to four sub-folders.  All archives contain at least an "Original_Data" and a
"Graphics" folder.  The other two folders are only included if strictly
necessary.

- Original Data.  Original data, with a "DATA_SOURCES.md" or "READ ME.txt" file 
that documents data sources.
**DATA IN THIS FOLDER IS AS ORIGINALLY PROVIDED OR ACCESSED.** 

- Derived Data.  Data derived from the original raw data.  Includes
documentation of data reorganization steps, either in the form of files (R
notebooks, Excel files, etc.) that embody data transformations, or via README.md
or DATA_NOTES.md files.

- Analysis.  Contains one or more R Notebooks proceeding through the data
analysis steps. This often includes both preliminary data analysis --
principally graphical, and detailed analysis where necessary.

- Graphics.  Contains R Notebooks stepping through development of graphics, and
also copies of resulting graphics, usually in \*.png and \*.pdf formats.  These
graphics may differ from graphics as they appear in final State of the Bay
graphical layouts.

# Summary of Data Sources
All data analyzed here is derived from data provided by GZA to LCWMD.  GZA has
been the primary monitoring, data management, and data analysis contractor for
the Long Creek Watershed Management District for the last several years. GZA
manages a large water quality data database on behalf of LCWMD. Tee data
provided constitutes a "snapshot" of the database as of mid 2019.

GZA produced a series of detailed data analysis reports for LCWMD based on these 
data.

Data was received from GZA in multiple formats, including large Excel Files and
a very large tab-delimited text file.  Data access in R is complicated by the
size of the raw data files, which often contain large volumes of data collected
with automated water quality sensor systems.
