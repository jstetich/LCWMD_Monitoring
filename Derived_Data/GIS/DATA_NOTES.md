# GIS Data Notes
## Shapefile "LC_Mon_Loc_Nov_2020" 
was received by Curtis C. Bohlen via e-mail from staff at the Cumberland County
Soil and Water Conservation District, which manages the LCWMD GIS data, 
in November of 2020.

## Shapefile "LC_Main_Mon_Loc"
was derived from the previous shapefile by selecting "Tier 1" monitoring 
locations only.

## Shapefile "LongCreekCatchments" 
was derived from NHD2+ WBD data by selecting
the relevant catchments and exporting a separate data file.

## The Shapefile "LongCreekOutline_UTM_smooth"
was derived from "LongCreekCatchments" by:
1. merging the Long Creek subwatersheds,
2. converting to UTM coordinates,
3. running a polygon smooth on the resulting "jagged" polygon; and 
4. hand editing the polygon just south of the Maine mall and adjacent to 
   Home Depot to address obvious inconsistencies where stream channels crossed
   nominal watershed boundaries.
   
** "LongCreekCatchments" is not an accurate representation of the Long
Creek watershed, and should not be used if precision is 
needed. It is intended only to show the general location of the watershed.   
In particular, the shapefile is derived from remotely sensed data, and 
does not include ground-truthed information on drainage.  Below-ground drainage
infrastructure significantly alters the boundaries of the effective watershed 
from what would be inferred from surface features alone.  Furthermore, obvious
errors regarding location of watershed boundaries suggest users treat the
boundary with a certain degree fo skepticism.**


## The Shapefile "Long_Creek" 
was derived from the National Hydrography Database
NHDFlowlines, by selecting features that intersect "LongCreekCatchments."  The
below ground portions of a flowline that crosses the Main Mall parking area,
south of the Maine Mall were deleted, leaving only the portions of a
related channel that are exposed above ground.
