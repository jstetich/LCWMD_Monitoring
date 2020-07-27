LCWMD Data Import
================
Curtis C. Bohlen, Casco Bay Estuary Partnership
(Revised 7/21/2020)

  - [Import Libraries](#import-libraries)
  - [Data Review](#data-review)
  - [Load Final Data](#load-final-data)
      - [First Site](#first-site)
      - [Iterate over other sites](#iterate-over-other-sites)
      - [Rename and Reorder](#rename-and-reorder)
  - [Export the Data](#export-the-data)

<img
  src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
  style="position:absolute;top:10px;right:50px;" />

# Import Libraries

``` r
library(readxl)
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------------------------------------------------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.1     v dplyr   1.0.0
    ## v tidyr   1.1.0     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ------------------------------------------------------------------------------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

# Data Review

We begin by reviewing content of one of the Excel spreadsheets provided
to us., for site S)1.

We have some trouble reading in these files, because of their size. If
we mis-specify data types, it slows down loading. But if we get that
right, we can read these in fairly quickly. It pays to specify data
columns to speed the data load process.

``` r
sibfldnm    <- 'Original_Data'
subfldnm    <- 'Data_Package_to_LCWMD'
subsubfldnm <- 'Data_by_Sites_Types'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

fn <- "S01 Merged Pressure and Sonde.xlsx"

fpath <- file.path(sibling, subfldnm, subsubfldnm, fn)

test.data <- read_excel(fpath, 
        sheet = "Consolidated 15 min",
        col_types = c("date", 
                "numeric", "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric"),
        skip = 2,
        na='NA#')
#View(test.data)
summary(test.data)
```

    ##  Sample Date (nearest 15 minutes) Chloride, Calculated CHLORIDE, CALCULATED
    ##  Min.   :2010-06-04 16:30:00      Min.   :   0.117     Min.   : NA         
    ##  1st Qu.:2011-09-15 06:07:30      1st Qu.: 273.417     1st Qu.: NA         
    ##  Median :2012-12-24 03:15:00      Median : 381.346     Median : NA         
    ##  Mean   :2013-12-07 14:46:26      Mean   : 422.439     Mean   :NaN         
    ##  3rd Qu.:2016-05-08 01:15:00      3rd Qu.: 532.810     3rd Qu.: NA         
    ##  Max.   :2019-06-10 08:30:00      Max.   :2796.682     Max.   : NA         
    ##                                   NA's   :14795        NA's   :148819      
    ##      Depth       Dissolved Oxygen DISSOLVED OXYGEN Dissolved oxygen saturation
    ##  Min.   :0.00    Min.   : 0.00    Min.   : NA      Min.   : NA                
    ##  1st Qu.:0.11    1st Qu.: 7.14    1st Qu.: NA      1st Qu.: NA                
    ##  Median :0.14    Median : 8.49    Median : NA      Median : NA                
    ##  Mean   :0.17    Mean   : 8.81    Mean   :NaN      Mean   :NaN                
    ##  3rd Qu.:0.18    3rd Qu.:10.93    3rd Qu.: NA      3rd Qu.: NA                
    ##  Max.   :2.28    Max.   :17.41    Max.   : NA      Max.   : NA                
    ##  NA's   :40941   NA's   :18117    NA's   :148819   NA's   :148819             
    ##  Dissolved Oxygen Saturation       pH            Pressure     
    ##  Min.   : -1.60              Min.   : NA      Min.   : 99.43  
    ##  1st Qu.: 75.50              1st Qu.: NA      1st Qu.:102.13  
    ##  Median : 84.00              Median : NA      Median :102.69  
    ##  Mean   : 80.06              Mean   :NaN      Mean   :102.85  
    ##  3rd Qu.: 91.20              3rd Qu.: NA      3rd Qu.:103.39  
    ##  Max.   :149.80              Max.   : NA      Max.   :122.18  
    ##  NA's   :16614               NA's   :148819   NA's   :40886   
    ##  Specific Conductivity  Temperature     Precipitation  
    ##  Min.   :   0          Min.   :-4.600   Min.   :0.00   
    ##  1st Qu.: 990          1st Qu.: 7.782   1st Qu.:0.00   
    ##  Median :1373          Median :14.325   Median :0.00   
    ##  Mean   :1513          Mean   :13.074   Mean   :0.19   
    ##  3rd Qu.:1910          3rd Qu.:18.770   3rd Qu.:0.15   
    ##  Max.   :9888          Max.   :35.435   Max.   :6.43   
    ##  NA's   :13034         NA's   :1299     NA's   :47554

``` r
# rm(test.data)
```

The structure of the source Excel files is a bit confusing. The
consolidated data in the first tab of each data sheet are derived from
lookup formulas that extract data from the source data, present in the
“Sonde” and “Pressure” tabs. The consolidated data sheet FIRST tries
to pull matching data from the Sonde data, and if that is not available,
it then pulls data from the pressure transducer.

That order of precedence is arguably incorrect for water depth, which is
based on pressure, and thus should be more accurate if derived from the
pressure transducers, not the sondes. However, it turns out that no
pressure or depth data was included in the sonde data (even though some
sondes DO have pressure sensors). We are correctly pulling depth data
from the pressure transducers only, despite the ambiguous formulas.

For some parameters, notably calculated chlorides, dissolved oxygen, and
percent saturation, each observation is split into two columns, the
first being the observation, and the second being labeled as “ND” in the
original data.

Lets figure out what type of values those include.

``` r
fn <- "S01 Merged Pressure and Sonde.xlsx"
fpath <- file.path(sibling, subfldnm, subsubfldnm, fn)
test.data <- read_excel(fpath, 
        sheet = "Consolidated 15 min",
        col_types = c("skip", "skip", "text", "skip", 
                "skip", "text", "text", "skip", 
                "skip", "skip", "skip", "skip","skip"),
        skip = 2,
        na='NA#') %>%
    mutate_all(~factor(.))
test.data %>% summarise()
```

    ## # A tibble: 1 x 0

So, there is No data in those columns in the first spreadsheet. We need
to check all available data, so I’ll go through each spreadsheet in turn
and see whether these columns contain any non-missing values.

``` r
results <- matrix(c(length(levels(test.data[,1])),
             length(levels(test.data[,2])),
             length(levels(test.data[,3]))), nrow=1)

othersites = c("S03", "S05", "S06B", "S07",  "S17")

for (site in othersites) {
    cat(site)
    cat("...")
    fn <- paste(site, "Merged Pressure and Sonde.xlsx")
    fpath <- file.path(sibling, subfldnm, subsubfldnm, fn)
    test.data <- read_excel(fpath, 
        sheet = "Consolidated 15 min",
        col_types = c("skip", "skip", "text", "skip", 
                "skip", "text", "text", "skip", 
                "skip", "skip", "skip", "skip","skip"),
        skip = 2,
        na='NA#') %>%
      mutate_all(~factor(.))
      current <- c(length(levels(test.data[,1])),
                   length(levels(test.data[,2])),
                   length(levels(test.data[,3])))
    results<- rbind(results, current)
}
```

    ## S03...S05...S06B...S07...S17...

``` r
rownames(results) <- c("S01", othersites)
results
```

    ##      [,1] [,2] [,3]
    ## S01     0    0    0
    ## S03     0    0    0
    ## S05     0    0    0
    ## S06B    0    0    0
    ## S07     0    0    0
    ## S17     0    0    0

So, there is no data in the the “ND” columns anywhere. We can safely
drop them.

``` r
rm(test.data, results, current, site)
```

# Load Final Data

## First Site

Using column designations speeds up file import – unless there are
errors. Preparing and printing error messages drastically slows this
code.

We load one data set first, to set up data format, before loading the
other files. This is not strictly necessary, but it’s convenient to have
a first dataframe to add others to using bind\_rows.

``` r
fn <- "S01 Merged Pressure and Sonde.xlsx"
fpath <- file.path(sibling, subfldnm, subsubfldnm, fn)

the_data <- read_excel(fpath, 
        sheet = "Consolidated 15 min",
        col_types = c("date", 
                "numeric", "skip", "numeric", "numeric", 
                "skip", "skip", "numeric", "numeric", 
                "numeric", "numeric", "numeric", 
                "numeric"),
        skip = 2,
        na='NA#')
```

## Iterate over other sites

``` r
othersites = c("S03", "S05", "S06B", "S07",  "S17")
the_data <- the_data %>%
    mutate(Site = "S01" )
for (site in othersites) {
    cat(site)
    cat("...\n")
    fn <- paste(site, "Merged Pressure and Sonde.xlsx")
    fpath <- file.path(sibling, subfldnm, subsubfldnm, fn)
    
    tmp.data <- read_excel(fpath, 
            sheet = "Consolidated 15 min",
            col_types = c("date", 
                    "numeric", "skip", "numeric", "numeric", 
                    "skip", "skip", "numeric", "numeric", 
                    "numeric", "numeric", "numeric", 
                    "numeric"),
            skip = 2,
            na='NA#') %>%
        mutate(Site = site)
    the_data <- the_data %>% bind_rows(tmp.data)
}
```

    ## S03...
    ## S05...
    ## S06B...
    ## S07...
    ## S17...

## Rename and Reorder

It is worth noting that the data are NOT sorted by date and time, but
apparently by Month, day, and time. That is, at least some of the time,
Data with the same calendar dates and times from different years are
stored together. I address that in this next block of code.

``` r
newnames <- c('DT', 'Chl', 'D', 'DO', 'PctSat', 'pH', 'Press', 'SpCond', 'T','Precip', 'Site')
names(the_data) <- newnames

the_data <- the_data %>%
    mutate(Site = factor(Site)) %>%
    arrange(Site, DT)

summary(the_data)
```

    ##        DT                           Chl                D         
    ##  Min.   :2010-06-04 13:00:00   Min.   :   0.12   Min.   : 0.00   
    ##  1st Qu.:2011-12-17 06:30:00   1st Qu.: 150.64   1st Qu.: 0.11   
    ##  Median :2014-05-30 17:00:00   Median : 257.24   Median : 0.17   
    ##  Mean   :2014-07-06 16:30:16   Mean   : 298.56   Mean   : 0.22   
    ##  3rd Qu.:2016-10-27 06:30:00   3rd Qu.: 397.15   3rd Qu.: 0.28   
    ##  Max.   :2019-06-10 13:00:00   Max.   :9110.04   Max.   :18.18   
    ##                                NA's   :160603    NA's   :155695  
    ##        DO             PctSat             pH             Press       
    ##  Min.   : 0.00    Min.   : -1.60   Min.   :3.6      Min.   :  6.01  
    ##  1st Qu.: 6.82    1st Qu.: 70.80   1st Qu.:6.9      1st Qu.:102.15  
    ##  Median : 8.30    Median : 82.30   Median :7.2      Median :102.94  
    ##  Mean   : 8.59    Mean   : 78.85   Mean   :7.1      Mean   :102.89  
    ##  3rd Qu.:10.54    3rd Qu.: 90.40   3rd Qu.:7.3      3rd Qu.:103.86  
    ##  Max.   :24.93    Max.   :217.90   Max.   :8.0      Max.   :181.60  
    ##  NA's   :176393   NA's   :198729   NA's   :688697   NA's   :146565  
    ##      SpCond             T              Precip         Site       
    ##  Min.   :    0    Min.   :-4.600   Min.   :0.00     S01 :148819  
    ##  1st Qu.:  560    1st Qu.: 7.983   1st Qu.:0.00     S03 :149047  
    ##  Median :  936    Median :14.038   Median :0.01     S05 :134512  
    ##  Mean   : 1080    Mean   :12.936   Mean   :0.19     S06B: 68552  
    ##  3rd Qu.: 1430    3rd Qu.:18.426   3rd Qu.:0.17     S07 :144486  
    ##  Max.   :32134    Max.   :41.820   Max.   :6.43     S17 : 48974  
    ##  NA's   :151997   NA's   :6920     NA's   :251436

``` r
rm(tmp.data, site, newnames)
```

# Export the Data

``` r
write.csv(the_data, 'Sonde_Data.csv')
```
