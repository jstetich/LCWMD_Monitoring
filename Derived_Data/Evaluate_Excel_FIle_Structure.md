Examining the Structure of the Excel Files from GZA
================
Curtis C. Bohlen, Casco Bay Estuary Partnership
7/22/2020

-   [Load Libraries](#load-libraries)
-   [Introduction](#introduction)
-   [List Excel Files COntaining Grab Sample Data](#list-excel-files-containing-grab-sample-data)
-   [Check File Structure](#check-file-structure)
    -   [Check Sheet Names](#check-sheet-names)
    -   [Check Names of Data Columns](#check-names-of-data-columns)
-   [Read Test Data](#read-test-data)
-   [Identifying Sonde and Pressure transducer Data](#identifying-sonde-and-pressure-transducer-data)
-   [What Data is Consistently Missing After we Remove Sonde Data?](#what-data-is-consistently-missing-after-we-remove-sonde-data)
-   [How to remove selected data](#how-to-remove-selected-data)

<img
  src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
  style="position:absolute;top:10px;right:50px;" />

Load Libraries
==============

``` r
library(readxl)
library(readr)
#library(lubridate)
library(tidyverse)
```

    ## -- Attaching packages ------------------------------------------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v dplyr   1.0.0
    ## v tibble  3.0.1     v stringr 1.4.0
    ## v tidyr   1.1.0     v forcats 0.5.0
    ## v purrr   0.3.4

    ## -- Conflicts ---------------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

Introduction
============

The Excel files provided by GZA are hard to work with. They are very large and slow to load. Our first attempt was to simplify these files by pooling them into separate tabs in an Excel Spreadsheet using an Excel Macro. it worked, but generated a huge file. What became clear is that the Excel Files contain BOTH sonde and grab sample data, so for loading the grab sample data we are going to if we read things in one at a time and filter them immediately.

This R Notebook explores the structure of the Excel files from GZA, in preparation for assembling a single consistent data set of grab sample data.

List Excel Files COntaining Grab Sample Data
============================================

``` r
sibfldnm    <- 'Original_Data'
subfldnm    <- 'Data_Package_to_LCWMD'
subsubfldnm <- 'Data_by_Sites_Types'

parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

targetfldr <- file.path(sibling, subfldnm, subsubfldnm)
fns <- list.files(targetfldr)
(fns <- fns[nchar(fns)<=9])
```

    ## [1] "S01.xlsx"  "S03.xlsx"  "S04.xlsx"  "S05.xlsx"  "S06B.xlsx" "S07.xlsx" 
    ## [7] "S11.xlsx"  "S12.xlsx"  "S17.xlsx"

Check File Structure
====================

Now, lets check to see if the Excel Spreadsheets contain the same data columns (and in the same order....)

Check Sheet Names
-----------------

``` r
for (fn in fns) {
  cat(fn)
  cat('\n')
  fpath <- file.path(sibling, subfldnm, subsubfldnm, fn)
  sheets <- excel_sheets(fpath)
  print(sheets)
  cat('\n')
}
```

    ## S01.xlsx
    ## [1] "S01_halfRL"
    ## 
    ## S03.xlsx
    ## [1] "S03_halfRL"
    ## 
    ## S04.xlsx
    ## [1] "S04"
    ## 
    ## S05.xlsx
    ## [1] "S05"
    ## 
    ## S06B.xlsx
    ## [1] "S06B"
    ## 
    ## S07.xlsx
    ## [1] "S07"    "Sheet2"
    ## 
    ## S11.xlsx
    ## [1] "S11"
    ## 
    ## S12.xlsx
    ## [1] "S12"
    ## 
    ## S17.xlsx
    ## [1] "S17"

In all cases, the first sheet contains the data we want.

Check Names of Data Columns
---------------------------

We learned from directly examining the files that each has two rows that should be skipped before the real data starts. In each case, the third row contains the column names.

``` r
# A simple file reader that encapuslated reading that third line
myreader <- function(nm) {
  cat(nm)
  cat('\n')
  fpath <- file.path(sibling, subfldnm, subsubfldnm, nm)
  line <- read_excel(fpath, sheet=1, col_names = FALSE, range = cell_rows(3), col_types = 'text')
  line$site = strsplit(nm, '.')[1]
  line
}

#Using lapply, we generate a list of one-row tibbles, and bind_rows them togehter into one tible.
nms <- lapply(fns, myreader)
```

    ## S01.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S03.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S04.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S05.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S06B.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S07.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S11.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S12.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

    ## S17.xlsx

    ## New names:
    ## * `` -> ...1
    ## * `` -> ...2
    ## * `` -> ...3
    ## * `` -> ...4
    ## * `` -> ...5
    ## * ...

``` r
names(nms) <- fns
bind_rows(nms, .id= 'source')
```

    ## # A tibble: 9 x 81
    ##   source ...1  ...2  ...3  ...4  ...5  ...6  ...7  ...8  ...9  ...10 ...11 ...12
    ##   <chr>  <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr> <chr>
    ## 1 S01.x~ samp~ Samp~ Samp~ Upda~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~
    ## 2 S03.x~ samp~ Samp~ Samp~ Upda~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~
    ## 3 S04.x~ samp~ Samp~ Samp~ Task~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~
    ## 4 S05.x~ samp~ Samp~ Samp~ Upda~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~
    ## 5 S06B.~ samp~ Samp~ Samp~ upda~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~
    ## 6 S07.x~ samp~ Samp~ Samp~ Upda~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~
    ## 7 S11.x~ samp~ Samp~ Samp~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~ Anti~
    ## 8 S12.x~ samp~ Samp~ Samp~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~ Anti~
    ## 9 S17.x~ samp~ Samp~ Samp~ Upda~ Task~ 2-Ch~ 2-Me~ Acen~ Acen~ Alum~ Ammo~ Anth~
    ## # ... with 68 more variables: ...13 <chr>, ...14 <chr>, ...15 <chr>,
    ## #   ...16 <chr>, ...17 <chr>, ...18 <chr>, ...19 <chr>, ...20 <chr>,
    ## #   ...21 <chr>, ...22 <chr>, ...23 <chr>, ...24 <chr>, ...25 <chr>,
    ## #   ...26 <chr>, ...27 <chr>, ...28 <chr>, ...29 <chr>, ...30 <chr>,
    ## #   ...31 <chr>, ...32 <chr>, ...33 <chr>, ...34 <chr>, ...35 <chr>,
    ## #   ...36 <chr>, ...37 <chr>, ...38 <chr>, ...39 <chr>, ...40 <chr>,
    ## #   ...41 <chr>, ...42 <chr>, ...43 <chr>, ...44 <chr>, ...45 <chr>,
    ## #   ...46 <chr>, ...47 <chr>, ...48 <chr>, ...49 <chr>, ...50 <chr>,
    ## #   ...51 <chr>, ...52 <chr>, ...53 <chr>, ...54 <chr>, ...55 <chr>,
    ## #   ...56 <chr>, ...57 <chr>, ...58 <chr>, ...59 <chr>, ...60 <chr>,
    ## #   ...61 <chr>, ...62 <chr>, ...63 <chr>, ...64 <chr>, ...65 <chr>,
    ## #   ...66 <chr>, ...67 <chr>, ...68 <chr>, ...69 <chr>, ...70 <chr>,
    ## #   ...71 <chr>, ...72 <chr>, ...73 <chr>, ...74 <chr>, ...75 <chr>,
    ## #   ...76 <chr>, ...77 <chr>, ...78 <chr>, ...79 <chr>, site <list>

So, what we see is that the structure appears similar -- two blank lines, then column names, then data.

The data column names are similar but not identical from one spreadsheet to the next. In particular: 1. S11 and S12 lack the "Updated Task Code" column, 2. "Updated Task Code" is called "Task Code Updated" in S04. 3. "Precipitation" is "Precipitation (inches) in S04, S11, and S12.

Other Notes: 1. Non-standard names (with spaces, parentheses, etc.) 2. A wide range of parameters, alphabetical, not sorted by categories 3. What was the 'Updated Task Code'? 4. Several Parameters are entered in multiple ways \* "CHLORIDE (AS Cl)" \* "CHLORIDE (AS CL)" \* "Chloride (Chlorine)" \* "Chloride, Calculated" \* "Dissolved Oxygen" \* "DISSOLVED OXYGEN" \* "DISSOLVED OXYGEN SAT PERCENT" \* "Dissolved oxygen saturation" \* "Dissolved Oxygen Saturation" 5. No data specifying the sample site. We'll have to add that.

Read Test Data
==============

``` r
fn = fns[1]
fpath <- file.path(sibling, subfldnm, subsubfldnm, fn)
test_data <- read_excel(fpath, sheet = 1, skip=2)
```

That's a BIG dataframe, accounting for \`r prettyNum(object.size(test\_data),big.mark = ",")' bytes. Ouch.

``` r
cat(prettyNum(object.size(test_data), big.mark = ','))
```

    ## 114,725,264

Identifying Sonde and Pressure transducer Data
==============================================

``` r
xtabs(~`Task Code` + `Updated Task Code`, data = test_data , addNA = TRUE)
```

    ##              Updated Task Code
    ## Task Code     Baseflow Pressure Sonde Data Spring Melt Stormwater
    ##   Baseflow          84        0          0           0          0
    ##   Pressure           0     7459          0           0          0
    ##   Sonde Data         0   100419      57661           0          0
    ##   Spring Melt        0        0          0          60          0
    ##   Stormwater         0        0          0           0         87

So, the Updated Task Code converted a large number of "Sonde Data" to "Pressure" codes. But that's all, so we can filter on either one safely.

What Data is Consistently Missing After we Remove Sonde Data?
=============================================================

``` r
test_data_2<- test_data %>%
  filter(`Task Code` != 'Sonde Data' & `Task Code` != 'Pressure' )

(nms <- names(test_data_2[,sapply(test_data_2, function(x) all(is.na(x)))]))
```

    ##  [1] "Barometric Pressure"          "CHLORIDE (AS CL)"            
    ##  [3] "Chloride (Chlorine)"          "CHLORIDE, CALCULATED"        
    ##  [5] "Cyanide"                      "Depth"                       
    ##  [7] "DISSOLVED OXYGEN"             "DISSOLVED OXYGEN SAT PERCENT"
    ##  [9] "Dissolved oxygen saturation"  "Pressure"                    
    ## [11] "Rainfall"

Most of those are parameters principally associated with sondes, pressure transducers, or weather records. These are appropriately removed from the data for grab samples. But it does not include ALL chloride or DO measures, presumably because of data QA/QC data collected.

How to remove selected data
===========================

``` r
(nms <- nms[c(1,4,6,7,8,9,10,11)])
```

    ## [1] "Barometric Pressure"          "CHLORIDE, CALCULATED"        
    ## [3] "Depth"                        "DISSOLVED OXYGEN"            
    ## [5] "DISSOLVED OXYGEN SAT PERCENT" "Dissolved oxygen saturation" 
    ## [7] "Pressure"                     "Rainfall"

``` r
test_data_2 <- test_data_2 %>% select( - any_of(nms))
```

That's better. Only \`r prettyNum(object.size(test\_data\_2),big.mark = ",")' bytes.

``` r
cat(prettyNum(object.size(test_data_2), big.mark = ','))
```

    ## 161,696
