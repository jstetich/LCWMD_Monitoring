Analysis of LCWMD ‘Chloride’ Data based on Conductivity Measurement
================
Curtis C. Bohlen, Casco Bay Estuary Partnership.
01/06/2021

  - [Introduction](#introduction)
  - [Import Libraries](#import-libraries)
  - [Data Preparation](#data-preparation)
      - [Initial Folder References](#initial-folder-references)
      - [Load Weather Data](#load-weather-data)
      - [Update Folder References](#update-folder-references)
      - [Load Data on Sites and Impervious
        Cover](#load-data-on-sites-and-impervious-cover)
      - [Load Main Data](#load-main-data)
          - [Cleanup](#cleanup)
      - [Data Correction](#data-correction)
          - [Anomolous Depth Values](#anomolous-depth-values)
          - [Single S06B Chloride Observation from
            2017](#single-s06b-chloride-observation-from-2017)
      - [Add Stream Flow Index](#add-stream-flow-index)
      - [Create Working Data](#create-working-data)
      - [Cleanup](#cleanup-1)
  - [Linear Models](#linear-models)
      - [Two Starting Models](#two-starting-models)
      - [Initial Practical Model](#initial-practical-model)
          - [Visualizing Estimated Marginal
            Means](#visualizing-estimated-marginal-means)
          - [Visualizing Trends](#visualizing-trends)
      - [Model Selection Using `step()`](#model-selection-using-step)
          - [Visualizing Estimated Marginal
            Means](#visualizing-estimated-marginal-means-1)
          - [Visualizing Trends](#visualizing-trends-1)
      - [Cleanup](#cleanup-2)
  - [GLS Analyses](#gls-analyses)
      - [Initial Model](#initial-model)
          - [Refit the Simple Model with
            REML](#refit-the-simple-model-with-reml)
          - [Visualizing Estimated Marginal
            Means](#visualizing-estimated-marginal-means-2)
          - [Visualizing Trends](#visualizing-trends-2)
      - [Model Selection Using `step()`](#model-selection-using-step-1)
      - [Manual Identification of Simpler
        Models](#manual-identification-of-simpler-models)
          - [Refit the Selected Model with
            REML](#refit-the-selected-model-with-reml)
          - [Visualizing Estimated Marginal
            Means](#visualizing-estimated-marginal-means-3)
          - [Visualizing Trends](#visualizing-trends-3)
      - [Cleanup](#cleanup-3)
  - [GAMM Analysis](#gamm-analysis)
      - [Initial Model](#initial-model-1)
          - [Diagnostic Plots](#diagnostic-plots)
          - [Visualizing Estimated Marginal
            Means](#visualizing-estimated-marginal-means-4)
          - [Visualizing Trends](#visualizing-trends-4)
      - [Full Model](#full-model)
          - [Visualizing Estimated Marginal
            Means](#visualizing-estimated-marginal-means-5)
          - [Visualizing Trends](#visualizing-trends-5)
      - [Simplified Models](#simplified-models)
          - [Visualizing Estimated Marginal
            Means](#visualizing-estimated-marginal-means-6)
          - [Visualizing Trends](#visualizing-trends-6)
  - [Seeking a Simplifies GLS Model](#seeking-a-simplifies-gls-model)
  - [Comparison Of Estimated Marginal
    Means](#comparison-of-estimated-marginal-means)
      - [Calculate Monthly Marginal
        Means](#calculate-monthly-marginal-means)
          - [the\_gls](#the_gls)
          - [selected\_gls](#selected_gls)
          - [highflow\_gls\_1](#highflow_gls_1)
          - [the\_gamm](#the_gamm)
          - [alt\_gamm\_1](#alt_gamm_1)
          - [full\_gamm](#full_gamm)
      - [Compare Results](#compare-results)
          - [Look at Correlations](#look-at-correlations)
          - [Look at Maximum Percent
            Differences](#look-at-maximum-percent-differences)
          - [Look at Graphic Comparison of Point
            Estimates](#look-at-graphic-comparison-of-point-estimates)
  - [Conclusion](#conclusion)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

This R Notebook explores various models for analyzing chloride levels in
Long Creek. Our principal goal is NOT to predict future values, but to
assess the contribution of certain (time varying) predictors to
explaining the pattern in the time series.

This notebook looks principally at chlorides, which are frequently
elevated in Maine urban streams because of use of salt for deicing of
roads, parking areas, sidewalks, and paths in winter. While water
quality standards are often written in terms of chlorides, it may be
better to think about chlorides as a measure of the salinity of the
water in the stream. Physiologically, it is probably salinity or
osmolarity that principally affects organisms in the stream, not
chlorides *per se*. The data we examine here is based on measurement of
conductivity, which is converted to an estimate of in-stream chlorides
based on a robust regression relationship developed over several years.

We are especially interested in looking at whether there have been
long-term trends in water quality over the duration of the LCWMD
monitoring program.

We develop various models that predict chlorides based on site,
precipitation, and stream flow,to reduce unexplained variability in the
data, and better characterize seasonal patterns and long term trends.

Simple linear models of the LCWMD data are based on the assumption that
observations are independent, however, we know both on principal and
from working with the data, that the data are both auto-correlated and
cross-correlated in complex ways. More robust analysis needs to address
that reality.

We sequentially work through 1. Linear models.  
2\. Generalized Least Squares Models with autocorrelated (ar1) error.  
3\. Generalized Additive Models with autocorrelated (ar1) error.

In each case we explore several models that estimate 1. Effect of time
of year (Month, or Day of Year) 2. Effect of Site 3. Differences (if
any) in trends from site to site.

To approach this question, we consider models of the form:

\[ Chlorides = f(Covariates) + Predictors + Error\]

Where: \* precipitation and flow covariates enter into the model via
linear or smooth functions. Precipitation terms are uniform across
sites, but the flow correction can vary site to site.

  - The predictors include linear functions of:  
    – site,  
    – year,  
    – a site by year interaction, and  
    – time of year.

  - The error includes either a simple normal error or a combination of
    a normal error and an AR(1) autocorrelated error.

We abuse the autocorrelation models slightly, since we use sequential
autocorrelations (not time-based) and we don’t fit separate
autocorrelations for each site and season. That should have little
impact on results, as transitions are relatively rare in a dense data
set, and missing values at the beginning of each season at each site
prevent estimation near season and site transitions in the sequential
data anyway.

Note that rather than work with time series methods, we use linear
models that incorporate lag terms, especially autocorrelated errors and
lagged precipitation variables instead.

The analysis may appear somewhat cavalier, but that is because much of
the usual statistical legwork was carried out previously. We explored
dozens of models before settling on the ones we present here. We have
checked model diagnostics.

On the whole, these models are OK, but not great. They tend to have
heavy tailed, skewed residuals. We should not trust the asymptotic p
values calculated by default. But since sample sizes are large and
results tend to have high statistical significance, p values are not
much use anyway. Even traditional methods using information criteria are
relatively little help for identifying minimal models that make
efficient use of the flow and precipitation covariates.

# Import Libraries

``` r
library(tidyverse)
#> -- Attaching packages --------------------------------------- tidyverse 1.3.0 --
#> v ggplot2 3.3.2     v purrr   0.3.4
#> v tibble  3.0.4     v dplyr   1.0.2
#> v tidyr   1.1.2     v stringr 1.4.0
#> v readr   1.4.0     v forcats 0.5.0
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
library(readr)

library(emmeans) # Provides tools for calculating marginal means
library(nlme)
#> 
#> Attaching package: 'nlme'
#> The following object is masked from 'package:dplyr':
#> 
#>     collapse

#library(zoo)     # here, for the `rollapply()` function

library(mgcv)    # generalized additive models. Function gamm() allows
#> This is mgcv 1.8-33. For overview type 'help("mgcv-package")'.
                 # autocorrelation.

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())
```

# Data Preparation

## Initial Folder References

``` r
sibfldnm    <- 'Original_Data'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)

dir.create(file.path(getwd(), 'figures'), showWarnings = FALSE)
dir.create(file.path(getwd(), 'models'), showWarnings = FALSE)
```

## Load Weather Data

``` r
fn <- "Portland_Jetport_2009-2019.csv"
fpath <- file.path(sibling, fn)

weather_data <- read_csv(fpath, 
 col_types = cols(.default = col_skip(),
        date = col_date(),
        PRCP = col_number(), PRCPattr = col_character() #,
        #SNOW = col_number(), SNOWattr = col_character(), 
        #TMIN = col_number(), TMINattr = col_character(), 
        #TAVG = col_number(), TAVGattr = col_character(), 
        #TMAX = col_number(), TMAXattr = col_character(), 
        )) %>%
  rename(sdate = date) %>%
  mutate(pPRCP = dplyr::lag(PRCP))
```

## Update Folder References

``` r
sibfldnm    <- 'Derived_Data'
parent      <- dirname(getwd())
sibling     <- file.path(parent,sibfldnm)
```

## Load Data on Sites and Impervious Cover

These data were derived from Table 2 from a GZA report to the Long Creek
Watershed Management District, titled “Re: Long Creek Watershed Data
Analysis; Task 2: Preparation of Explanatory and Other Variables.” The
Memo is dated November 13, 2019 File No. 09.0025977.02.

Cumulative Area and IC calculations are our own, based on the GZA data
and the geometry of the stream channel.

``` r
# Read in data and drop the East Branch, where we have no data
fn <- "Site_IC_Data.csv"
fpath <- file.path(sibling, fn)

Site_IC_Data <- read_csv(fpath) %>%
  filter(Site != "--") 
#> 
#> -- Column specification --------------------------------------------------------
#> cols(
#>   Site = col_character(),
#>   Subwatershed = col_character(),
#>   Area_ac = col_double(),
#>   IC_ac = col_double(),
#>   CumArea_ac = col_double(),
#>   CumIC_ac = col_double(),
#>   PctIC = col_character(),
#>   CumPctIC = col_character()
#> )

# Now, create a factor that preserves the order of rows (roughly upstream to downstream). 
Site_IC_Data <- Site_IC_Data %>%
  mutate(Site = factor(Site, levels = Site_IC_Data$Site))

# Finally, convert percent covers to numeric values
Site_IC_Data <- Site_IC_Data %>%
  mutate(CumPctIC = as.numeric(substr(CumPctIC, 1, nchar(CumPctIC)-1))) %>%
  mutate(PctIC = as.numeric(substr(PctIC, 1, nchar(PctIC)-1)))
Site_IC_Data
#> # A tibble: 6 x 8
#>   Site  Subwatershed      Area_ac IC_ac CumArea_ac CumIC_ac PctIC CumPctIC
#>   <fct> <chr>               <dbl> <dbl>      <dbl>    <dbl> <dbl>    <dbl>
#> 1 S07   Blanchette Brook     434.  87.7       434.     87.7  20.2     20.2
#> 2 S06B  Upper Main Stem      623.  80.2       623.     80.2  12.9     12.9
#> 3 S05   Middle Main Stem     279.  53.6      1336     222.   19.2     16.6
#> 4 S17   Lower Main Stem      105   65.1      1441     287.   62       19.9
#> 5 S03   North Branch Trib    298. 123         298.    123    41.2     41.2
#> 6 S01   South Branch Trib    427. 240.        427.    240.   56.1     56.1
```

## Load Main Data

Read in the data from the Derived Data folder.

Note that I filter out data from 2019 because that is only a partial
year, which might affect estimation of things like seasonal trends. We
could add it back in, but with care….

*Full\_Data.csv* does not include a field for precipitation from the
previous day. In earlier work, we learned that a weighted sum of recent
precipitation provided better explanatory power. But we also want to
check a simpler model, so we construct a “PPrecip” data field. This is
based on a modification of code in the “Make\_Daily\_Summaries.Rmd”
notebook.

``` r
fn <- "Full_Data.csv"
fpath <- file.path(sibling, fn)

full_data <- read_csv(fpath, 
    col_types = cols(DOY = col_integer(), 
        D_Median = col_double(), Precip = col_number(), 
        X1 = col_skip(), Year = col_integer(), 
        FlowIndex = col_double())) %>%

  mutate(Site = factor(Site, levels=levels(Site_IC_Data$Site))) %>%
  mutate(Month = factor(Month, levels = month.abb)) %>%
  mutate(IC=as.numeric(Site_IC_Data$CumPctIC[match(Site, Site_IC_Data$Site)])) %>%
  mutate(Yearf = factor(Year)) %>%

# We combine data using "match" because we have data for multiple sites and 
# therefore dates are not unique.  `match()` correctly assigns weather
# data by date.
mutate(PPrecip = weather_data$pPRCP[match(sdate, weather_data$sdate)])
#> Warning: Missing column names filled in: 'X1' [1]
#> Warning: The following named parsers don't match the column names: FlowIndex
```

### Cleanup

``` r
rm(Site_IC_Data, weather_data)
rm(fn, fpath, parent, sibling, sibfldnm)
```

## Data Correction

### Anomolous Depth Values

Several depth observations in the record appear highly unlikely. In
particular, several observations show daily median water depths over 15
meters. And those observations were recorded in May or June, at site
S05, with no associated record of significant precipitation, and no
elevated depths at other sites on the stream.

We can trace these observations back to the raw QA/QC’d pressure and
sonde data submitted to LCWMD by GZA, so they are not an artifact of our
data preparation.

A few more observations show daily median depths over 4 meters, which
also looks unlikely in a stream of this size. All these events also
occurred in May or June of 2015 at site S05. Some sort of malfunction of
the pressure transducer appears likely.

We remove these extreme values. The other daily medians in May and June
of 2015 appear reasonable, and we leave them in place, although given
possible instability of the pressure sensors, it might make sense to
remove them all.

``` r
full_data <- full_data %>%
  mutate(D_Median = if_else(D_Median > 4, NA_real_, D_Median),
         lD_Median = if_else(D_Median > 4, NA_real_, lD_Median))
```

### Single S06B Chloride Observation from 2017

The data includes just a single chloride observation from site S06B from
any year other than 2013. While we do not know if the data point is
legitimate or not, it has very high leverage in several models, and we
suspect a transcription error of some sort.

``` r
full_data %>%
  filter(Site == 'S06B') %>%
  select(sdate, Chl_Median) %>%
  ggplot(aes(x = sdate, y = Chl_Median)) + geom_point()
#> Warning: Removed 1214 rows containing missing values (geom_point).
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

We remove the Chloride value from the data.

``` r
full_data <- full_data %>%
  mutate(Chl_Median = if_else(Site == 'S06B' & Year > 2014,
                              NA_real_, Chl_Median))
```

## Add Stream Flow Index

We worked through many models on a site by site basis in which we
included data on water depth, but since the depth coordinate is
site-specific, a 10 cm depth at one site may be exceptional, while at
another it is commonplace. We generally want not a local measure of
stream depth, but a watershed-wide metric of high, medium, or low stream
flow.

Middle and Lower Maine Stem sites would be suitable for a general flow
indicator across the watershed. The monitoring sites in that stretch of
Long Creek include include S05 and S17, however only site S05 has been
in continuous operation throughout the period of record, so we use depth
data from S05 to construct our general stream flow indicator.

Stream flow at S05 is correlated with flow at other sites, although not
all that closely correlated to flow in the downstream tributaries.

``` r
full_data %>%
  select(sdate, Site, lD_Median) %>%
  pivot_wider(names_from = Site, values_from = lD_Median) %>%
  select( -sdate) %>%
  cor(use = 'pairwise', method = 'pearson')
#>            S07      S06B       S05       S17       S03       S01
#> S07  1.0000000 0.5882527 0.7177120 0.7327432 0.4434108 0.5859663
#> S06B 0.5882527 1.0000000 0.8043943 0.8778188 0.7152403 0.6310361
#> S05  0.7177120 0.8043943 1.0000000 0.7906571 0.4250331 0.6668570
#> S17  0.7327432 0.8778188 0.7906571 1.0000000 0.6666414 0.7290077
#> S03  0.4434108 0.7152403 0.4250331 0.6666414 1.0000000 0.4441852
#> S01  0.5859663 0.6310361 0.6668570 0.7290077 0.4441852 1.0000000
```

We use the log of the daily median flow at S05 as a general
watershed-wide stream flow indicator, which we call `FlowIndex`. We use
the log of the raw median, to lessen the effect of the highly skewed
distribution of stream depths on the metric.

``` r
depth_data <- full_data %>%
  filter (Site == 'S05') %>%
  select(sdate, lD_Median)

full_data <- full_data %>%
  mutate(FlowIndex = depth_data$lD_Median[match(sdate, depth_data$sdate)])
  rm(depth_data)
```

Note that because the flow record at S05 has some gaps, any model using
this predictor is likely to have a smaller sample size.

## Create Working Data

Including Site = S06B in the GLS models causes an error, because models
that includes a Site:Year interaction are rank deficient. We only have
one year’s worth of data from that site. (`lm()` handles that case
gracefully, `gls()` does not.)

``` r
xtabs(~ Site + Year, data = full_data)
#>       Year
#> Site   2010 2011 2012 2013 2014 2015 2016 2017 2018
#>   S07   176  311  288  230  262  191  246  265  240
#>   S06B    0    0    0  223  217  193  247  265  240
#>   S05   126  283  182  190  217  231  248  265  241
#>   S17     0    0    0    0    0  223  235  265  241
#>   S03   192  311  298  256  262  223  246  265  241
#>   S01   192  311  298  271  252  217  248  265  241
```

We proceed with analyses that omits Site S06B.

``` r
reduced_data <- full_data %>%
  select (Site, Year, Month, DOY,
          Precip, lPrecip, PPrecip, wlPrecip,
          D_Median, lD_Median,
          Chl_Median, 
          IC, FlowIndex) %>%
  filter (Site != 'S06B' ) %>%
    mutate(Site = droplevels(Site)) %>%
  mutate(Year_f = factor(Year))
```

## Cleanup

``` r
rm(full_data)
```

# Linear Models

We begin by developing several simple linear models. We do not have much
confidence in these models, but they are a starting point for
understanding the analytic choices we face in developing more
sophisticated models.

## Two Starting Models

We use these principally to demonstrate the sizable advantage of our
weighted recent rainfall predictor over a simple lagged rainfall metric.

``` r
simple_lm_tmp <- lm(log(Chl_Median) ~ Site + 
               lPrecip + 
               log(PPrecip+1) +
               FlowIndex +
               Month + 
               Year, 
               data = reduced_data,
               subset = ! is.na(wlPrecip))   # so both models use the same data
AIC(simple_lm_tmp)
#> [1] 9079.1

the_lm <- lm(log(Chl_Median) ~ Site + 
               lPrecip + 
               wlPrecip +
               FlowIndex +
               Month + 
               Year,
             data = reduced_data)
AIC(the_lm)
#> [1] 8710.248
anova(the_lm)
#> Analysis of Variance Table
#> 
#> Response: log(Chl_Median)
#>             Df  Sum Sq Mean Sq  F value    Pr(>F)    
#> Site         4  526.86  131.72  574.697 < 2.2e-16 ***
#> lPrecip      1   70.73   70.73  308.606 < 2.2e-16 ***
#> wlPrecip     1  420.55  420.55 1834.923 < 2.2e-16 ***
#> FlowIndex    1   22.73   22.73   99.171 < 2.2e-16 ***
#> Month        9  238.66   26.52  115.702 < 2.2e-16 ***
#> Year         1  117.93  117.93  514.540 < 2.2e-16 ***
#> Residuals 6350 1455.36    0.23                       
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

So the model using a weighted sum of recent rainfall does much better
than one based only on the prior day’s rainfall. We continue, looking
only at models  
using weighted recent rainfall.

## Initial Practical Model

This model is the linear model prototype for similar GLS and GAM models
we develop later.

``` r
the_lm <- lm(log(Chl_Median) ~ Site + 
               lPrecip + 
               wlPrecip +
               FlowIndex +
               Month + 
               Year +
               Site: Year,
             data = reduced_data)
AIC(the_lm)
#> [1] 8115.919
anova(the_lm)
#> Analysis of Variance Table
#> 
#> Response: log(Chl_Median)
#>             Df  Sum Sq Mean Sq F value    Pr(>F)    
#> Site         4  526.86  131.72  631.31 < 2.2e-16 ***
#> lPrecip      1   70.73   70.73  339.01 < 2.2e-16 ***
#> wlPrecip     1  420.55  420.55 2015.69 < 2.2e-16 ***
#> FlowIndex    1   22.73   22.73  108.94 < 2.2e-16 ***
#> Month        9  238.66   26.52  127.10 < 2.2e-16 ***
#> Year         1  117.93  117.93  565.23 < 2.2e-16 ***
#> Site:Year    4  131.35   32.84  157.39 < 2.2e-16 ***
#> Residuals 6346 1324.02    0.21                      
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Visualizing Estimated Marginal Means

``` r
month.emm <- emmeans(the_lm, ~Month, 
                     cov.reduce = median, 
                     type = 'response')

labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Sites Combined'

plot(month.emm) +
  annotate('text', 450, 7, label = labl, size = 3) +
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab('') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip()
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

``` r
site.emm <- emmeans(the_lm, ~Site, 
                    cov.reduce = median, 
                    type='response')
#> NOTE: Results may be misleading due to involvement in interactions

labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation'

plot(site.emm) +
  annotate('text', 450, 2.5, label = labl, size = 3) +
  ylab("Upstream        Main Stem                 Lower Tribs         ") +
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip() 
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

### Visualizing Trends

We use `emmeans()` and `emtrends()` to extract

``` r
a <- summary(emmeans(the_lm, 'Site', 
                     cov.reduce = median))
#> NOTE: Results may be misleading due to involvement in interactions
b <- summary(emtrends(the_lm, 'Site', 'Year', 
                      cov.reduce = median))

lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])
rm(a,b)
```

And create a temporary dataframe for the results.

``` r
df <- tibble(Site = rep(levels(reduced_data$Site), each = 9), 
              Year = rep(2010:2018, 5)) %>%
  mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
         iintercept = lookup$Intercept[match(Site, lookup$Site)],
         pred = exp((Year - 2014) * sslope + iintercept)) %>%
  select(-sslope, -iintercept)
```

``` r
ggplot(df, aes(x = Year, y = pred, color = Site)) +
         geom_step(direction = 'mid') +
  ylab('Flow-Adjusted Annual Mean Chloride')
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

``` r
rm(df)
```

## Model Selection Using `step()`

``` r
full_lm <- lm(log(Chl_Median) ~ Site + 
               lPrecip * wlPrecip * FlowIndex +
               FlowIndex : Site +   # We expect flow to play out differently by site
               Month + 
               Year + 
               Year:Site,
              data = reduced_data)

step_lm <- step(full_lm,
                upper = )
#> Start:  AIC=-10634.69
#> log(Chl_Median) ~ Site + lPrecip * wlPrecip * FlowIndex + FlowIndex:Site + 
#>     Month + Year + Year:Site
#> 
#>                              Df Sum of Sq    RSS      AIC
#> - lPrecip:wlPrecip:FlowIndex  1      0.02 1187.5 -10636.6
#> <none>                                    1187.5 -10634.7
#> - Site:Year                   4    108.24 1295.7 -10087.2
#> - Site:FlowIndex              4    130.20 1317.7  -9980.2
#> - Month                       9    204.58 1392.1  -9640.5
#> 
#> Step:  AIC=-10636.58
#> log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + Month + 
#>     Year + lPrecip:wlPrecip + lPrecip:FlowIndex + wlPrecip:FlowIndex + 
#>     Site:FlowIndex + Site:Year
#> 
#>                      Df Sum of Sq    RSS      AIC
#> - lPrecip:wlPrecip    1     0.003 1187.5 -10638.6
#> <none>                            1187.5 -10636.6
#> - wlPrecip:FlowIndex  1     0.434 1188.0 -10636.3
#> - lPrecip:FlowIndex   1     4.083 1191.6 -10616.7
#> - Site:Year           4   108.224 1295.7 -10089.2
#> - Site:FlowIndex      4   130.202 1317.7  -9982.1
#> - Month               9   205.193 1392.7  -9639.6
#> 
#> Step:  AIC=-10638.57
#> log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + Month + 
#>     Year + lPrecip:FlowIndex + wlPrecip:FlowIndex + Site:FlowIndex + 
#>     Site:Year
#> 
#>                      Df Sum of Sq    RSS      AIC
#> <none>                            1187.5 -10638.6
#> - wlPrecip:FlowIndex  1     0.452 1188.0 -10638.1
#> - lPrecip:FlowIndex   1     6.209 1193.7 -10607.4
#> - Site:Year           4   108.222 1295.7 -10091.2
#> - Site:FlowIndex      4   130.202 1317.7  -9984.1
#> - Month               9   206.337 1393.8  -9636.4
anova(step_lm)
#> Analysis of Variance Table
#> 
#> Response: log(Chl_Median)
#>                      Df  Sum Sq Mean Sq  F value    Pr(>F)    
#> Site                  4  526.86  131.72  703.214 < 2.2e-16 ***
#> lPrecip               1   70.73   70.73  377.618 < 2.2e-16 ***
#> wlPrecip              1  420.55  420.55 2245.259 < 2.2e-16 ***
#> FlowIndex             1   22.73   22.73  121.348 < 2.2e-16 ***
#> Month                 9  238.66   26.52  141.576 < 2.2e-16 ***
#> Year                  1  117.93  117.93  629.605 < 2.2e-16 ***
#> lPrecip:FlowIndex     1    6.51    6.51   34.742 3.958e-09 ***
#> wlPrecip:FlowIndex    1    0.15    0.15    0.787     0.375    
#> Site:FlowIndex        4  152.97   38.24  204.178 < 2.2e-16 ***
#> Site:Year             4  108.22   27.06  144.446 < 2.2e-16 ***
#> Residuals          6340 1187.51    0.19                       
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The best fit linear model includes Site terms related to the Flow Index,
and interaction terms relating precipitation and recent precipitation to
flow conditions. The sums of squares associated with the interaction
terms, however, are quite low. We could probably drop the interaction
terms for the sake of simplicity, but we retain it for now, as the GLS
analyses are likely to be slightly more parsimonious.

``` r
plot(step_lm)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-24-1.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-24-2.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-24-3.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-24-4.png" style="display: block; margin: auto;" />

The distribution of errors is skewed, and there are a couple of
outliers, luckily all with low leverage. None of the points have
especially high leverage. If anything the model is slightly light
tailed, which suggests an underfit model.

``` r
old <- options() 
options(digits=3,scipen=2)
summary(step_lm)
#> 
#> Call:
#> lm(formula = log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + 
#>     Month + Year + lPrecip:FlowIndex + wlPrecip:FlowIndex + Site:FlowIndex + 
#>     Site:Year, data = reduced_data)
#> 
#> Residuals:
#>    Min     1Q Median     3Q    Max 
#> -5.630 -0.232  0.036  0.277  2.563 
#> 
#> Coefficients:
#>                       Estimate  Std. Error t value Pr(>|t|)    
#> (Intercept)        -264.294104    9.292337  -28.44  < 2e-16 ***
#> SiteS05             128.854545   15.272278    8.44  < 2e-16 ***
#> SiteS17             548.775223   32.519425   16.88  < 2e-16 ***
#> SiteS03             248.874162   12.394102   20.08  < 2e-16 ***
#> SiteS01             164.362131   12.464722   13.19  < 2e-16 ***
#> lPrecip              -0.000645    0.009871   -0.07     0.95    
#> wlPrecip             -0.121601    0.007628  -15.94  < 2e-16 ***
#> FlowIndex            -2.516665    0.134801  -18.67  < 2e-16 ***
#> MonthApr             -0.169833    0.038764   -4.38  1.2e-05 ***
#> MonthMay             -0.268812    0.038764   -6.93  4.5e-12 ***
#> MonthJun             -0.461909    0.038970  -11.85  < 2e-16 ***
#> MonthJul             -0.557996    0.040456  -13.79  < 2e-16 ***
#> MonthAug             -0.726729    0.039598  -18.35  < 2e-16 ***
#> MonthSep             -0.749232    0.039629  -18.91  < 2e-16 ***
#> MonthOct             -0.644141    0.038899  -16.56  < 2e-16 ***
#> MonthNov             -0.558096    0.038767  -14.40  < 2e-16 ***
#> MonthDec             -0.418712    0.049835   -8.40  < 2e-16 ***
#> Year                  0.134488    0.004609   29.18  < 2e-16 ***
#> lPrecip:FlowIndex    -0.178739    0.031045   -5.76  8.9e-09 ***
#> wlPrecip:FlowIndex    0.037568    0.024196    1.55     0.12    
#> SiteS05:FlowIndex     0.870368    0.129134    6.74  1.7e-11 ***
#> SiteS17:FlowIndex     0.971933    0.154628    6.29  3.5e-10 ***
#> SiteS03:FlowIndex     3.048621    0.120257   25.35  < 2e-16 ***
#> SiteS01:FlowIndex     1.320442    0.119840   11.02  < 2e-16 ***
#> SiteS05:Year         -0.064179    0.007577   -8.47  < 2e-16 ***
#> SiteS17:Year         -0.272266    0.016126  -16.88  < 2e-16 ***
#> SiteS03:Year         -0.123689    0.006151  -20.11  < 2e-16 ***
#> SiteS01:Year         -0.081396    0.006186  -13.16  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 0.433 on 6340 degrees of freedom
#>   (3377 observations deleted due to missingness)
#> Multiple R-squared:  0.584,  Adjusted R-squared:  0.582 
#> F-statistic:  329 on 27 and 6340 DF,  p-value: <2e-16
options(old)
rm(old)
```

So: 1. Present-day precipitation suppresses chlorides, but effect is
tiny compared to the effect of recent precipitation, and depends on flow
conditions. The Effect of precipitation is lower under high flow
conditions. It’s not yet clear how to include those terms in more
complex models. 3. Chlorides have been INCREASING over the period of
record, but the rate of change varies by Site. 6. The highest (flow
adjusted) chlorides are in the winter months, especially February and
March. Flow adjusted Chlorides are lowest in summer\!

### Visualizing Estimated Marginal Means

``` r
month.emm <- emmeans(step_lm, ~Month, 
                     cov.reduce = median, 
                     type = 'response')

labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Sites Combined'

plot(month.emm) +
  annotate('text', 400, 6, label = labl, size = 3) +
  ylab("Month")+
  xlab('Flow-Adjusted Chlorides (mg/l)') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip()
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-26-1.png" style="display: block; margin: auto;" />

``` r
site.emm <- emmeans(step_lm, ~Site, 
                    cov.reduce = median, 
                    type='response')
#> NOTE: Results may be misleading due to involvement in interactions

labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation'


plot(site.emm) +
  annotate('text', 450, 2.5, label = labl, size = 3) +
  ylab("Upstream        Main Stem                 Lower Tribs         ") +
  xlab('Chlorides (mg/l)') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip() 
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

### Visualizing Trends

We extract marginal means and marginal slopes using `emmeans()` and
`emtrends()`. These calculate marginal means at “typical” conditions -
usually the means of all other predictors, but we specify medians, which
handle skewed predictors better.

``` r
a <- summary(emmeans(step_lm, 'Site', 
                     cov.reduce = median))
#> NOTE: Results may be misleading due to involvement in interactions
b <- summary(emtrends(step_lm, 'Site', 'Year', 
                      cov.reduce = median))

lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])
rm(a,b)
```

And create a temporary dataframe for the results.

``` r
df <- tibble(Site = rep(levels(reduced_data$Site), each = 9), 
              Year = rep(2010:2018, 5)) %>%
  mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
         iintercept = lookup$Intercept[match(Site, lookup$Site)],
         pred = exp((Year - 2014) * sslope + iintercept)) %>%
  select(-sslope, -iintercept)
```

``` r
ggplot(df, aes(x = Year, y = pred, color = Site)) +
  geom_step(direction = 'mid') +
  ggtitle('Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  ylab('Flow-Adjusted Mean Chloride')
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

``` r
rm(df)
```

The big outlier here is S17, which is a partial record, and thus can not
say much about the longer-term trends. S03 is a full record, but that is
the South Branch, which has separate issues.

``` r
plot(emtrends(step_lm, 'Site', 'Year', cov.reduce = median)) +
  geom_vline(xintercept = 0) +
 xlab('Slope of Predicted Log of Chloride (mg/l)' )
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

## Cleanup

``` r
rm(full_lm, simple_lm_tmp, step_lm)
```

# GLS Analyses

We can improve on those models two ways: 1. By taking into account
temporal autocorrelation and 2. By allowing relationships with
predictors to be non-linear

Here we use generalize least squares to accommodate autocorrelated
errors.

We abuse the autocorrelation model slightly, since we don’t fit separate
autocorrelations for each site and season. That should have little
impact on results, as transitions are rare, and missing values at
beginning of most seasonal time series prevent estimation near season
and site transitions in the sequential data anyway.

## Initial Model

This model mirrors the simple linear model. This model takes about a
minute to run.

``` r
the_gls_ml <- gls(log(Chl_Median) ~ Site + 
               lPrecip + 
               wlPrecip +
               FlowIndex +
               Month + 
               Year +
               Site: Year,
               correlation = corAR1(0.8),
               na.action = na.omit, 
               method = 'ML',
               data = reduced_data)
anova(the_gls_ml)
#> Denom. DF: 6346 
#>             numDF  F-value p-value
#> (Intercept)     1 68754.56  <.0001
#> Site            4    43.53  <.0001
#> lPrecip         1     6.03  0.0141
#> wlPrecip        1  2213.57  <.0001
#> FlowIndex       1  1403.19  <.0001
#> Month           9    26.81  <.0001
#> Year            1   111.21  <.0001
#> Site:Year       4     4.74  0.0008
```

This produces coefficients that are rather different from the simple
linear model, although differences in predictions will likely be small.
The big differences relate to what is absorbed into the intercepts and
what turns up elsewhere. Site S17 continues to cause trouble because of
its relative short record. For now, we leave it in place.

``` r
cbind(coef(the_lm),coef(the_gls_ml))
#>                       [,1]          [,2]
#> (Intercept)  -279.79616602 -273.21966576
#> SiteS05       138.77732736   91.79621489
#> SiteS17       568.46646674  239.09769513
#> SiteS03       282.54722541  185.56863602
#> SiteS01       178.15142087  116.72930096
#> lPrecip        -0.04838229   -0.03728404
#> wlPrecip       -0.11001639   -0.09869611
#> FlowIndex      -1.30315792   -1.76630810
#> MonthApr       -0.17583431   -0.09873760
#> MonthMay       -0.28141077   -0.25320624
#> MonthJun       -0.47277133   -0.43643139
#> MonthJul       -0.57543237   -0.50116278
#> MonthAug       -0.74350639   -0.62697035
#> MonthSep       -0.77019715   -0.68812521
#> MonthOct       -0.65188930   -0.58287250
#> MonthNov       -0.55628552   -0.41922031
#> MonthDec       -0.38350020   -0.27366991
#> Year            0.14205147    0.13879278
#> SiteS05:Year   -0.06900814   -0.04569495
#> SiteS17:Year   -0.28192375   -0.11862222
#> SiteS03:Year   -0.14006572   -0.09192408
#> SiteS01:Year   -0.08809458   -0.05763025
```

### Refit the Simple Model with REML

The function `emmeans()` is a bit happier working with models fit with
REML.

``` r
the_gls <- gls(log(Chl_Median) ~ Site + 
               lPrecip + 
               wlPrecip +
               FlowIndex +
               Month + 
               Year +
               Site : Year,
               correlation = corAR1(0.8),
               na.action = na.omit, 
               method = 'REML',
               data = reduced_data)
summary(the_gls)
#> Generalized least squares fit by REML
#>   Model: log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + Month +      Year + Site:Year 
#>   Data: reduced_data 
#>        AIC      BIC    logLik
#>   188.2198 350.3538 -70.10992
#> 
#> Correlation Structure: AR(1)
#>  Formula: ~1 
#>  Parameter estimate(s):
#>       Phi 
#> 0.8602781 
#> 
#> Coefficients:
#>                   Value Std.Error   t-value p-value
#> (Intercept)  -273.07758  35.30401  -7.73503  0.0000
#> SiteS05        89.84925  56.79114   1.58210  0.1137
#> SiteS17       226.13403 111.03171   2.03666  0.0417
#> SiteS03       181.01034  46.83263   3.86505  0.0001
#> SiteS01       113.99565  47.38149   2.40591  0.0162
#> lPrecip        -0.03729   0.00312 -11.96408  0.0000
#> wlPrecip       -0.09871   0.00501 -19.70745  0.0000
#> FlowIndex      -1.76734   0.04540 -38.92576  0.0000
#> MonthApr       -0.09830   0.05094  -1.92977  0.0537
#> MonthMay       -0.25223   0.05748  -4.38785  0.0000
#> MonthJun       -0.43456   0.05943  -7.31267  0.0000
#> MonthJul       -0.49773   0.06201  -8.02673  0.0000
#> MonthAug       -0.62254   0.06183 -10.06895  0.0000
#> MonthSep       -0.68404   0.05995 -11.41099  0.0000
#> MonthOct       -0.58013   0.05739 -10.10822  0.0000
#> MonthNov       -0.41744   0.05390  -7.74441  0.0000
#> MonthDec       -0.27199   0.07050  -3.85780  0.0001
#> Year            0.13872   0.01752   7.91776  0.0000
#> SiteS05:Year   -0.04473   0.02818  -1.58719  0.1125
#> SiteS17:Year   -0.11219   0.05506  -2.03767  0.0416
#> SiteS03:Year   -0.08966   0.02325  -3.85712  0.0001
#> SiteS01:Year   -0.05627   0.02352  -2.39282  0.0167
#> 
#>  Correlation: 
#>              (Intr) SitS05 SitS17 SitS03 SitS01 lPrecp wlPrcp FlwInd MnthAp
#> SiteS05      -0.682                                                        
#> SiteS17      -0.312  0.173                                                 
#> SiteS03      -0.747  0.515  0.196                                          
#> SiteS01      -0.738  0.506  0.233  0.533                                   
#> lPrecip      -0.005  0.005  0.024  0.003  0.002                            
#> wlPrecip     -0.007  0.010  0.039  0.003  0.002  0.640                     
#> FlowIndex     0.004 -0.003 -0.002 -0.001 -0.001 -0.371 -0.519              
#> MonthApr     -0.023 -0.010 -0.010  0.005  0.014 -0.040 -0.052  0.064       
#> MonthMay     -0.041 -0.017 -0.036  0.010  0.019 -0.054 -0.060  0.084  0.728
#> MonthJun     -0.047 -0.013 -0.017  0.009  0.018 -0.068 -0.095  0.111  0.659
#> MonthJul     -0.043 -0.018 -0.023  0.006  0.015 -0.055 -0.078  0.129  0.608
#> MonthAug     -0.036 -0.033 -0.042 -0.004  0.007 -0.066 -0.074  0.128  0.593
#> MonthSep     -0.057 -0.028 -0.046 -0.008  0.005 -0.037 -0.053  0.103  0.583
#> MonthOct     -0.055 -0.008 -0.036 -0.001  0.007 -0.078 -0.122  0.122  0.574
#> MonthNov     -0.060  0.010 -0.009  0.011  0.008 -0.036 -0.064  0.069  0.554
#> MonthDec     -0.084  0.018  0.012 -0.003 -0.006 -0.064 -0.086  0.066  0.453
#> Year         -1.000  0.682  0.312  0.747  0.738  0.005  0.007 -0.005  0.022
#> SiteS05:Year  0.682 -1.000 -0.173 -0.515 -0.507 -0.005 -0.010  0.003  0.010
#> SiteS17:Year  0.313 -0.174 -1.000 -0.196 -0.233 -0.024 -0.039  0.002  0.010
#> SiteS03:Year  0.747 -0.515 -0.196 -1.000 -0.532 -0.003 -0.003  0.001 -0.005
#> SiteS01:Year  0.738 -0.506 -0.233 -0.533 -1.000 -0.002 -0.002  0.001 -0.014
#>              MnthMy MnthJn MnthJl MnthAg MnthSp MnthOc MnthNv MnthDc Year  
#> SiteS05                                                                    
#> SiteS17                                                                    
#> SiteS03                                                                    
#> SiteS01                                                                    
#> lPrecip                                                                    
#> wlPrecip                                                                   
#> FlowIndex                                                                  
#> MonthApr                                                                   
#> MonthMay                                                                   
#> MonthJun      0.797                                                        
#> MonthJul      0.706  0.815                                                 
#> MonthAug      0.671  0.733  0.814                                          
#> MonthSep      0.656  0.693  0.726  0.808                                   
#> MonthOct      0.637  0.664  0.672  0.712  0.806                            
#> MonthNov      0.626  0.636  0.627  0.646  0.699  0.761                     
#> MonthDec      0.465  0.482  0.457  0.456  0.473  0.486  0.548              
#> Year          0.039  0.046  0.042  0.035  0.056  0.054  0.058  0.083       
#> SiteS05:Year  0.017  0.013  0.018  0.033  0.028  0.008 -0.010 -0.018 -0.682
#> SiteS17:Year  0.036  0.017  0.023  0.042  0.046  0.036  0.009 -0.012 -0.313
#> SiteS03:Year -0.010 -0.009 -0.006  0.004  0.008  0.001 -0.011  0.003 -0.747
#> SiteS01:Year -0.019 -0.018 -0.015 -0.007 -0.005 -0.007 -0.008  0.006 -0.738
#>              SS05:Y SS17:Y SS03:Y
#> SiteS05                          
#> SiteS17                          
#> SiteS03                          
#> SiteS01                          
#> lPrecip                          
#> wlPrecip                         
#> FlowIndex                        
#> MonthApr                         
#> MonthMay                         
#> MonthJun                         
#> MonthJul                         
#> MonthAug                         
#> MonthSep                         
#> MonthOct                         
#> MonthNov                         
#> MonthDec                         
#> Year                             
#> SiteS05:Year                     
#> SiteS17:Year  0.174              
#> SiteS03:Year  0.515  0.196       
#> SiteS01:Year  0.507  0.233  0.532
#> 
#> Standardized residuals:
#>         Min          Q1         Med          Q3         Max 
#> -11.6648396  -0.5199376   0.1206259   0.6390631   3.7213273 
#> 
#> Residual standard error: 0.4762795 
#> Degrees of freedom: 6368 total; 6346 residual
```

### Visualizing Estimated Marginal Means

`emmeans()` fails to estimate Satterthwaite estimates of the effective
degrees of freedom from GLS models, so we provide the `mode =
"df.error"` parameter. This may slightly over estimate the degrees of
freedom, and thus cause tests to be overly optimistic. See the ‘models’
vignette of `emmeans` for details.

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Sites Combined'

month.emm <- emmeans(the_gls, ~ Month,
                     cov.reduce = median,
                     type = 'response',
                     mode = "df.error")

plot(month.emm) +
  annotate('text', 400, 7, label = labl, size = 3) +
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab ('') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip()
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-37-1.png" style="display: block; margin: auto;" />

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation'

site.emm <- emmeans(the_gls, ~Site,
                    cov.reduce = median,
                    type ='response',
                    mode = "df.error")
#> NOTE: Results may be misleading due to involvement in interactions
plot(site.emm) +
  annotate('text', 450, 2.5, label = labl, size = 3) +
  ylab("Upstream        Main Stem                 Lower Tribs         ") +
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip() 
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-38-1.png" style="display: block; margin: auto;" />

### Visualizing Trends

``` r
a <- summary(emmeans(the_gls, 'Site',
             cov.reduce = median, 
             mode = "df.error"))
#> NOTE: Results may be misleading due to involvement in interactions

b <- summary(emtrends(the_gls, 'Site', 'Year',
                      cov.reduce = median,
                      mode = "df.error"))

lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])


df <- tibble(Site = rep(levels(reduced_data$Site), each = 9), 
             Year = rep(2010:2018, 5)) %>%
      mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
             iintercept = lookup$Intercept[match(Site, lookup$Site)],
             pred = exp((Year - 2014) * sslope + iintercept)) %>%
      select(-sslope, -iintercept)

ggplot(df, aes(x = Year, y = pred, color = Site)) +
         geom_step(direction = 'mid') +
  ylab('Predicted Mean Chloride') +
  theme_cbep(base_size = 12)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-39-1.png" style="display: block; margin: auto;" />
Which shows fairly consistent behavior across sites. Again, the slope
for S17 reflect only a few recent years.

``` r
rm(a,b)
```

## Model Selection Using `step()`

Fitting a larger GLS model takes a few minutes.

``` r
full_gls <- gls(log(Chl_Median) ~ Site + 
                  lPrecip * wlPrecip * FlowIndex +
                  FlowIndex : Site +   # We expect flow to play out differently by site
                  Month + 
                  Year + 
                  Site:Year,
               correlation = corAR1(0.8),
               na.action = na.omit,
               method = 'ML',
               data = reduced_data)
```

The following does not work if we use the function `step()`, so we use
the similar function from the `MASS` package. Since we don’t need that
package otherwise, we load the function by a direct call to the package
namespace.

As this fits several GLS models, it takes several minutes to run to
completion.

``` r
step_gls <- MASS::stepAIC(full_gls)
#> Start:  AIC=-463.46
#> log(Chl_Median) ~ Site + lPrecip * wlPrecip * FlowIndex + FlowIndex:Site + 
#>     Month + Year + Site:Year
#> 
#>                              Df     AIC
#> <none>                          -463.46
#> - Site:Year                   4 -458.51
#> - lPrecip:wlPrecip:FlowIndex  1 -457.02
#> - Month                       9 -339.24
#> - Site:FlowIndex              4  -47.41
anova(step_gls)
#> Denom. DF: 6338 
#>                            numDF  F-value p-value
#> (Intercept)                    1 66061.99  <.0001
#> Site                           4    41.51  <.0001
#> lPrecip                        1     6.25  0.0124
#> wlPrecip                       1  2412.47  <.0001
#> FlowIndex                      1  1546.72  <.0001
#> Month                          9    27.77  <.0001
#> Year                           1   117.98  <.0001
#> lPrecip:wlPrecip               1    50.68  <.0001
#> lPrecip:FlowIndex              1    43.89  <.0001
#> wlPrecip:FlowIndex             1    34.48  <.0001
#> Site:FlowIndex                 4   109.48  <.0001
#> Site:Year                      4     3.79  0.0044
#> lPrecip:wlPrecip:FlowIndex     1     8.41  0.0038
```

The stepwise procedure retains all interactions among the covariates.
Those sharply increase model complexity, but hugely increase the log
likelihood as well.

``` r
anova(the_gls_ml, step_gls)
#>            Model df       AIC       BIC    logLik   Test  L.Ratio p-value
#> the_gls_ml     1 24   71.6142  233.8312 -11.80711                        
#> step_gls       2 32 -463.4563 -247.1670 263.72813 1 vs 2 551.0705  <.0001
```

Although stepwise selection uses an automated procedure to identify the
“best” model for predictive purposes, in our setting, a simpler model
is easier to understand, explain, and communicate. We want to find a
model that has a better balance between complexity and adequacy. It
would be nice if most of the improvement in model fit is associated with
only a couple of the covariates.

We need to check.

## Manual Identification of Simpler Models

We look at the impact of dropping single terms from the “best” model.
This fits several alternate models, and so takes some time.

``` r
drop1(step_gls, test = 'Chisq')
#> Single term deletions
#> 
#> Model:
#> log(Chl_Median) ~ Site + lPrecip * wlPrecip * FlowIndex + FlowIndex:Site + 
#>     Month + Year + Site:Year
#>                            Df     AIC    LRT  Pr(>Chi)    
#> <none>                        -463.46                     
#> Month                       9 -339.24 142.22 < 2.2e-16 ***
#> Site:FlowIndex              4  -47.41 424.05 < 2.2e-16 ***
#> Site:Year                   4 -458.51  12.95  0.011537 *  
#> lPrecip:wlPrecip:FlowIndex  1 -457.02   8.44  0.003671 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

So dropping the three-way interaction has the smallest effect on
likelihood. Interestingly, dropping the Site : Year interaction term
also has relatively low effect of likelihood, but as that is a focus of
the analysis, rather than a covariate, we retain it.

``` r
gls_1 <- update(step_gls, . ~ . - lPrecip:wlPrecip:FlowIndex )
drop1(gls_1, test = 'Chisq')
#> Single term deletions
#> 
#> Model:
#> log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + Month + 
#>     Year + lPrecip:wlPrecip + lPrecip:FlowIndex + wlPrecip:FlowIndex + 
#>     Site:FlowIndex + Site:Year
#>                    Df     AIC    LRT  Pr(>Chi)    
#> <none>                -457.02                     
#> Month               9 -336.99 138.03 < 2.2e-16 ***
#> lPrecip:wlPrecip    1 -457.39   1.63   0.20156    
#> lPrecip:FlowIndex   1 -396.76  62.26 3.011e-15 ***
#> wlPrecip:FlowIndex  1 -432.91  26.10 3.236e-07 ***
#> Site:FlowIndex      4  -41.77 423.25 < 2.2e-16 ***
#> Site:Year           4 -452.11  12.91   0.01175 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
gls_2 <- update(gls_1, . ~ . - lPrecip:wlPrecip )
drop1(gls_2, test = 'Chisq')
#> Single term deletions
#> 
#> Model:
#> log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + Month + 
#>     Year + lPrecip:FlowIndex + wlPrecip:FlowIndex + Site:FlowIndex + 
#>     Site:Year
#>                    Df     AIC    LRT  Pr(>Chi)    
#> <none>                -457.39                     
#> Month               9 -337.21 138.17 < 2.2e-16 ***
#> lPrecip:FlowIndex   1 -348.47 110.92 < 2.2e-16 ***
#> wlPrecip:FlowIndex  1 -428.01  31.37 2.131e-08 ***
#> Site:FlowIndex      4  -43.02 422.37 < 2.2e-16 ***
#> Site:Year           4 -452.44  12.94   0.01156 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

``` r
gls_3 <- update(gls_2, . ~ . - wlPrecip:FlowIndex )
drop1(gls_3, test = 'Chisq')
#> Single term deletions
#> 
#> Model:
#> log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + Month + 
#>     Year + lPrecip:FlowIndex + Site:FlowIndex + Site:Year
#>                   Df     AIC    LRT Pr(>Chi)    
#> <none>               -428.01                    
#> wlPrecip           1   15.95 445.96   <2e-16 ***
#> Month              9 -296.49 149.52   <2e-16 ***
#> lPrecip:FlowIndex  1 -349.69  80.32   <2e-16 ***
#> Site:FlowIndex     4   -8.32 427.69   <2e-16 ***
#> Site:Year          4 -422.78  13.23   0.0102 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

At this point, any further simplification of the covariate terms in the
model has relatively large impact on likelihoods.

``` r
anova(step_gls, gls_3)
#>          Model df       AIC       BIC   logLik   Test  L.Ratio p-value
#> step_gls     1 32 -463.4563 -247.1670 263.7281                        
#> gls_3        2 29 -428.0140 -232.0018 243.0070 1 vs 2 41.44231  <.0001
```

### Refit the Selected Model with REML

The function `emmeans()` is a bit happier working with models fit with
REML.

``` r
selected_gls <- update(gls_3,  method = 'REML')
anova(selected_gls)
#> Denom. DF: 6341 
#>                   numDF  F-value p-value
#> (Intercept)           1 63218.51  <.0001
#> Site                  4    39.62  <.0001
#> lPrecip               1     6.13  0.0133
#> wlPrecip              1  2398.64  <.0001
#> FlowIndex             1  1542.70  <.0001
#> Month                 9    27.21  <.0001
#> Year                  1   116.46  <.0001
#> lPrecip:FlowIndex     1    88.22  <.0001
#> Site:FlowIndex        4   110.82  <.0001
#> Site:Year             4     3.47  0.0077
```

### Visualizing Estimated Marginal Means

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Sites Combined'

month.emm <- emmeans(selected_gls, ~ Month,
                     cov.reduce = median,
                     type = 'response',
                     mode = "df.error")

plot(month.emm) +
  annotate('text', 400, 7, label = labl, size = 3) +
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab ('') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip()
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-50-1.png" style="display: block; margin: auto;" />

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation'


site.emm <- emmeans(selected_gls, ~Site,
                    cov.reduce = median,
                    type ='response',
                    mode = "df.error")
#> NOTE: Results may be misleading due to involvement in interactions
plot(site.emm) +
  annotate('text', 450, 2.5, label = labl, size = 3) +
  ylab("Upstream        Main Stem                 Lower Tribs         ") +
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ggtitle('Marginal Geometric Means of Daily Medians') +
  theme_cbep(base_size = 12) +
  coord_flip() 
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-51-1.png" style="display: block; margin: auto;" />

### Visualizing Trends

``` r
a <- summary(emmeans(selected_gls, 'Site',
             cov.reduce = median, 
             mode = "df.error"))
#> NOTE: Results may be misleading due to involvement in interactions

b <- summary(emtrends(selected_gls, 'Site', 'Year',
                      cov.reduce = median,
                      mode = "df.error"))

lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])


df <- tibble(Site = rep(levels(reduced_data$Site), each = 9), 
             Year = rep(2010:2018, 5)) %>%
      mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
             iintercept = lookup$Intercept[match(Site, lookup$Site)],
             pred = exp((Year - 2014) * sslope + iintercept)) %>%
      select(-sslope, -iintercept)

ggplot(df, aes(x = Year, y = pred, color = Site)) +
         geom_step(direction = 'mid') +
  ylab('Predicted Mean Chloride') +
  theme_cbep(base_size = 12)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-52-1.png" style="display: block; margin: auto;" />
Which shows fairly consistent behavior across sites. Again, the slope
for S17 reflect only a few recent years.

``` r
rm(a,b)
```

## Cleanup

``` r
rm(the_gls_ml, step_gls, gls_3, gls_2, gls_1)
```

# GAMM Analysis

Here we use more sophisticated “General Additive Models” that allow
non-linear (smoother) fits for some parameters. Our emphasis is on using
smoothers to better account for non-linearities in relationships between
weather or flow-related predictors and chlorides.

We use the function `gamm()` because it has a relatively simple
interface for incorporating autocorrelated errors.

We abuse the autocorrelation model slightly, since we don’t fit separate
autocorrelations for each site and season. That should have little
impact on results, as missing values at beginning and end of most time
series prevent estimation anyway.

## Initial Model

Our first GAMM simply fits smoothers for each of the major
weather-related covariates. Arguably, we should fit separate smoothers
by `FlowIndex` for each site, but we did not include interaction terms
in our earlier base models, so we leave that out here as well.

This model takes several minutes to run (more than 5, less than 15)

``` r
if (! file.exists("models/the_gamm.rds")) {
  the_gamm <- gamm(log(Chl_Median) ~ Site + 
                     s(lPrecip) + 
                     s(wlPrecip) +
                     s(FlowIndex) +
                     Month +
                     Year +
                     Site : Year,
                   correlation = corAR1(0.8),
                   na.action = na.omit, 
                   method = 'REML',
                   data = reduced_data)
  saveRDS(the_gamm, file="models/the_gamm.rds")
} else {
  the_gamm <- readRDS("models/the_gamm.rds")
}
```

``` r
anova(the_gamm$gam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> log(Chl_Median) ~ Site + s(lPrecip) + s(wlPrecip) + s(FlowIndex) + 
#>     Month + Year + Site:Year
#> 
#> Parametric Terms:
#>           df      F  p-value
#> Site       4  2.162   0.0706
#> Month      9 20.325  < 2e-16
#> Year       1 49.773 1.91e-12
#> Site:Year  4  2.148   0.0723
#> 
#> Approximate significance of smooth terms:
#>                edf Ref.df      F p-value
#> s(lPrecip)   7.035  7.035  26.80  <2e-16
#> s(wlPrecip)  4.142  4.142  57.59  <2e-16
#> s(FlowIndex) 8.459  8.459 231.41  <2e-16
```

Interestingly, differences between sites and differences in slopes are
marginally not significant in this simplified model.

``` r
plot(the_gamm$gam)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-56-1.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-56-2.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-56-3.png" style="display: block; margin: auto;" />
Note that the function for recent weighted precipitation is nearly
linear, while the effect of present-day precipitation is near zero for
low to moderate rainfall, but drops quickly for rainfall over about 4 cm
or 1.5 inches (rare events). Chlorides drop with increasing water depth,
up to a point, but then climb again at the highest (very rare) flow
levels.

What these smoothers show is that sticking with linear terms for many of
our covariates should work fairly well, except at the highest flow
conditions. We might also consider adding a “high rainfall” term, rather
than fitting a a linear or smoothed predictor term for today’s rain. The
cost of such model simplification would be a drop in ability to
accurately predict chloride levels under the highest flow, highest
rainfall conditions.

### Diagnostic Plots

The help files for `gam.check()` suggest using care when interpreting
results for GAMM models, since the function does not correctly
incorporate the error correlations structure. However, for our purposes,
this is probably sufficient, since our focus is not on statistical
significance, but on estimation.

``` r
gam.check(the_gamm$gam)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-57-1.png" style="display: block; margin: auto;" />

    #> 
    #> 'gamm' based fit - care required with interpretation.
    #> Checks based on working residuals may be misleading.
    #> Basis dimension (k) checking results. Low p-value (k-index<1) may
    #> indicate that k is too low, especially if edf is close to k'.
    #> 
    #>                k'  edf k-index p-value    
    #> s(lPrecip)   9.00 7.03    0.99     0.2    
    #> s(wlPrecip)  9.00 4.14    0.93  <2e-16 ***
    #> s(FlowIndex) 9.00 8.46    0.90  <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

As with the linear model, we have a skewed, slightly heavy tailed
distribution of residuals, with a couple of very large outliers. There
is perhaps slight evidence for lack of complete independence between
residuals and predictors. T his model is adequate, but not great. For
careful work, we should probably use bootstrapped confidence intervals
or something similar, but for our purposes, that is probably overkill.

### Visualizing Estimated Marginal Means

Reliably calling `emmeans()` for these large `gamm()` models appears to
require creating a call object and associating it with the model (e.g.,
as `the_gamm$gam$call`). (See the `emmeans` models vignette for more
info, although not all strategies recommended there worked for us).

We first create the call object, then associate it with the model, and
finally manually construct a reference grid before calling `emmeans()`
to extract marginal means. This workflow has the advantage that it
requires us to think carefully about the structure of the reference
grid.

Note also that we explicitly specify that we want the marginal means
estimated at Year = 2014. This is largely to be explicit, and avoid
possible confusion from here on out. The default method creates a
reference grid where marginal means are keyed to mean values of all
predictors, which would be some value slightly larger than 2014.
However, we specified `cov.reduce = median`, and the median Year
predictor is precisely 2014. Although this setting is probably
unnecessary, we chose to be explicit from here on out.

``` r
the_call <-  quote(gamm(log(Chl_Median) ~ Site + 
                          s(lPrecip) + 
                          s(wlPrecip) +
                          s(FlowIndex) +
                          Month +
                          Year +
                          Site : Year,
                        correlation = corAR1(0.8),
                        na.action = na.omit, 
                        method = 'REML',
                        data = reduced_data))
the_gamm$gam$call <- the_call

my_ref_grid <- ref_grid(the_gamm, at = list(Year = 2014), cov.reduce = median) 
a <- emmeans(my_ref_grid, ~ Month, type = 'response')

labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Sites Combined'

plot(a) + 
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab ('') +
  annotate('text', 400, 6, label = labl, size = 3) +
  xlim(0,500) +
  geom_vline(xintercept =  230, color = 'orange') +
  geom_vline(xintercept =  860, color = 'red') +
  coord_flip() +
  theme_cbep(base_size = 12)
#> Warning: Removed 1 rows containing missing values (geom_vline).
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-58-1.png" style="display: block; margin: auto;" />

Note that all estimated monthly means are somewhat lower than from the
GLS model.

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Dates Combined'

a <- emmeans(my_ref_grid, ~ Site, type = 'response')
#> NOTE: Results may be misleading due to involvement in interactions
plot(a) + 
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab("Upstream                  Main Stem                                 Lower Tribs                   ") +
  annotate('text', 400, 2.5, label = labl, size = 3) +
  xlim(0,500) +
  geom_vline(xintercept =  230, color = 'orange') +
  geom_vline(xintercept =  860, color = 'red') +
  coord_flip() +
  theme_cbep(base_size = 12)
#> Warning: Removed 1 rows containing missing values (geom_vline).
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-59-1.png" style="display: block; margin: auto;" />

### Visualizing Trends

We extract results on the log scale, so we can calculate the linear
predictor by hand, then back transform.

``` r
my_ref_grid <- ref_grid(the_gamm, at = list(Year = 2014, Month = 'Jul'),
                        cov.reduce = median)

a <- summary(emmeans(my_ref_grid, 'Site'))
#> NOTE: Results may be misleading due to involvement in interactions
b <- summary(emtrends(the_gamm, 'Site', 'Year'))

lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])
rm(a,b)

df <- tibble(Site = rep(levels(reduced_data$Site), each = 9), 
              Year = rep(2010:2018, 5)) %>%
  mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
         iintercept = lookup$Intercept[match(Site, lookup$Site)],
         pred = exp((Year - 2014) * sslope + iintercept)) %>%
  select(-sslope, -iintercept)
```

``` r
ggplot(df, aes(x = Year, y = pred, color = Site)) +
         geom_step(direction = 'mid') +
  ylab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  xlab('') +
  ylim(0,600) +
  geom_hline(yintercept =  230, color = 'black') +
  #geom_hline(yintercept =  860, color = 'red') +

  theme_cbep(base_size = 12)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-61-1.png" style="display: block; margin: auto;" />

## Full Model

The Full model explores a series of interaction smoothers using tensor
smooth terms. This allows testing of which terms are important int eh
model, much as one might do by looking at interactions in linear models.
A similar (but not identical) final model can be refit more compactly
and with fewer parameters by fitting equivalent joint smooths that only
include important components.

Unfortunately, using `stepAIC()` for GAMM models is tricky, as we need
to use tensor terms to explore model structure. And neither `step()` nor
`stepAIC()` appear to work reliably with a GAMM model.

The full GAMM takes over an hour to run. To speed up analysis, we save a
copy of the model for reuse. But this may pose a problem if the model or
the underlying data changes, in which case you need to delete the old
version of the model manually.

``` r
if (! file.exists("models/full_gamm.rds")) {
  full_gamm <- gamm(log(Chl_Median) ~ Site + 
                      ti(lPrecip) +
                      ti(wlPrecip) + 
                      ti(FlowIndex, by = Site) +
                      ti(lPrecip,  wlPrecip) +
                      ti(lPrecip, FlowIndex) +
                      ti(wlPrecip, FlowIndex) +
                      ti(lPrecip, wlPrecip, FlowIndex) +
                      Month +
                      Year + 
                      Site:Year,
                    correlation = corAR1(0.8),
                    na.action = na.omit,
                    method = 'REML',
                    data = reduced_data)
saveRDS(full_gamm, file="models/full_gamm.rds")
} else {
  full_gamm <- readRDS("models/full_gamm.rds")
}
```

``` r
anova(full_gamm$gam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> log(Chl_Median) ~ Site + ti(lPrecip) + ti(wlPrecip) + ti(FlowIndex, 
#>     by = Site) + ti(lPrecip, wlPrecip) + ti(lPrecip, FlowIndex) + 
#>     ti(wlPrecip, FlowIndex) + ti(lPrecip, wlPrecip, FlowIndex) + 
#>     Month + Year + Site:Year
#> 
#> Parametric Terms:
#>           df      F  p-value
#> Site       4  0.569    0.685
#> Month      9 15.940  < 2e-16
#> Year       1 39.214 4.05e-10
#> Site:Year  4  0.567    0.687
#> 
#> Approximate significance of smooth terms:
#>                                   edf Ref.df       F  p-value
#> ti(lPrecip)                     3.059  3.059  14.280  < 2e-16
#> ti(wlPrecip)                    1.816  1.816 142.756  < 2e-16
#> ti(FlowIndex):SiteS07           3.671  3.671 161.143  < 2e-16
#> ti(FlowIndex):SiteS05           1.927  1.927  18.239 4.76e-06
#> ti(FlowIndex):SiteS17           3.536  3.536  73.228  < 2e-16
#> ti(FlowIndex):SiteS03           3.731  3.731   9.747 2.02e-07
#> ti(FlowIndex):SiteS01           3.845  3.845 167.522  < 2e-16
#> ti(lPrecip,wlPrecip)           13.301 13.301  13.763  < 2e-16
#> ti(lPrecip,FlowIndex)           8.760  8.760   6.106  < 2e-16
#> ti(wlPrecip,FlowIndex)          6.424  6.424   6.468  < 2e-16
#> ti(lPrecip,wlPrecip,FlowIndex) 27.869 27.869   5.670  < 2e-16
```

All terms are highly significant by ANOVA, except differences among
sites. but several interaction terms terms have high standard errors /
low T values, despite appearing statistically significant by ANOVA.

``` r
summary(full_gamm$lme)
#> Linear mixed-effects model fit by REML
#>  Data: strip.offset(mf) 
#>         AIC       BIC   logLik
#>   -1115.436 -791.2289 605.7181
#> 
#> Random effects:
#>  Formula: ~Xr - 1 | g
#>  Structure: pdTens
#>              Xr1     Xr2      Xr3
#> StdDev: 21.10682 6.74354 1.590896
#> 
#>  Formula: ~Xr.0 - 1 | g.0 %in% g
#>  Structure: pdTens
#>            Xr.01    Xr.02    Xr.03
#> StdDev: 107.2626 42.53465 8.338999
#> 
#>  Formula: ~Xr.1 - 1 | g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.11    Xr.12     Xr.13
#> StdDev: 6.495293 2.332564 0.5266659
#> 
#>  Formula: ~Xr.2 - 1 | g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.21    Xr.22    Xr.23
#> StdDev: 35.86388 12.87929 2.907995
#> 
#>  Formula: ~Xr.3 - 1 | g.3 %in% g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.31    Xr.32     Xr.33
#> StdDev: 6.026581 2.164242 0.4886607
#> 
#>  Formula: ~Xr.4 - 1 | g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.41    Xr.42     Xr.43
#> StdDev: 6.120255 2.197882 0.4962562
#> 
#>  Formula: ~Xr.5 - 1 | g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.51    Xr.52     Xr.53
#> StdDev: 4.261418 1.530344 0.3455338
#> 
#>  Formula: ~Xr.6 - 1 | g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.61    Xr.62    Xr.63    Xr.64    Xr.65    Xr.66    Xr.67
#> StdDev: 2.609149 2.058366 1.990198 1.986107 1.866367 1.699081 1.692051
#>             Xr.68    Xr.69    Xr.610    Xr.611   Xr.612    Xr.613    Xr.614
#> StdDev: 0.9552702 0.797844 0.7875847 0.5622217 0.540603 0.2002675 0.1544075
#>            Xr.615
#> StdDev: 0.1275359
#> 
#>  Formula: ~Xr.7 - 1 | g.7 %in% g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.71    Xr.72    Xr.73    Xr.74    Xr.75    Xr.76    Xr.77    Xr.78
#> StdDev: 18.81276 6.088098 1.746146 1.021901 18.78857 18.78516 18.78498 6.012931
#>            Xr.79   Xr.710   Xr.711   Xr.712   Xr.713     Xr.714   Xr.715
#> StdDev: 1.462676 0.366981 6.002293 6.001721 1.418313 0.08286005 1.415891
#> 
#>  Formula: ~Xr.8 - 1 | g.8 %in% g.7 %in% g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.81    Xr.82    Xr.83    Xr.84    Xr.85   Xr.86    Xr.87    Xr.88
#> StdDev: 33.25122 32.55351 14.77441 32.45428 32.44894 7.68595 7.260152 13.12901
#>            Xr.89   Xr.810   Xr.811   Xr.812   Xr.813    Xr.814   Xr.815
#> StdDev: 12.88099 12.86753 3.627908 2.607237 2.590479 0.5886839 2.522703
#> 
#>  Formula: ~Xr.9 - 1 | g.9 %in% g.8 %in% g.7 %in% g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g
#>  Structure: pdTens
#>            Xr.91    Xr.92    Xr.93    Xr.94    Xr.95    Xr.96    Xr.97    Xr.98
#> StdDev: 4.986063 2.814897 2.470902 2.449121 4.674353 4.766477 4.628873 4.626417
#>         Xr.99   Xr.910   Xr.911   Xr.912   Xr.913   Xr.914   Xr.915   Xr.916
#> StdDev: 4.726 4.724375 2.216393 2.404605 1.759008 1.728279 1.990939 1.963842
#>           Xr.917   Xr.918   Xr.919   Xr.920  Xr.921   Xr.922   Xr.923  Xr.924
#> StdDev: 2.118794 2.113423 2.323337 4.439376 2.32003 1.634319 1.601199 1.62735
#>           Xr.925   Xr.926  Xr.927   Xr.928   Xr.929   Xr.930   Xr.931   Xr.932
#> StdDev: 1.594085 1.891985 1.86345 1.887922 1.859324 4.391463 4.388874 4.395888
#>          Xr.933   Xr.934   Xr.935   Xr.936   Xr.937   Xr.938    Xr.939   Xr.940
#> StdDev: 4.39414 4.347495 4.345729 4.344881 4.343113 1.664596 0.9760075 0.919472
#>           Xr.941   Xr.942   Xr.943 Xr.944   Xr.945    Xr.946    Xr.947
#> StdDev: 1.532243 1.524808 1.544879 1.5399 0.727653 0.6498592 0.7118635
#>            Xr.948    Xr.949    Xr.950    Xr.951    Xr.952   Xr.953   Xr.954
#> StdDev: 0.6321296 0.7538975 0.6791168 0.7436415 0.6677132 1.401263 1.395772
#>           Xr.955   Xr.956    Xr.957    Xr.958    Xr.959    Xr.960    Xr.961
#> StdDev: 1.393129 1.387606 0.3811163 0.1951611 0.3604039 0.3500293 0.1507619
#>            Xr.962    Xr.963  Residual
#> StdDev: 0.1239302 0.3273558 0.5078005
#> 
#> Correlation Structure: AR(1)
#>  Formula: ~1 | g/g.0/g.1/g.2/g.3/g.4/g.5/g.6/g.7/g.8/g.9 
#>  Parameter estimate(s):
#>       Phi 
#> 0.9065866 
#> Fixed effects: y ~ X - 1 
#>                                         Value Std.Error   DF   t-value p-value
#> X(Intercept)                       -274.17194  44.88140 6338 -6.108809  0.0000
#> XSiteS05                             72.34453  70.67474 6338  1.023627  0.3061
#> XSiteS17                            133.75644 129.49581 6338  1.032902  0.3017
#> XSiteS03                             70.78263  58.80149 6338  1.203756  0.2287
#> XSiteS01                             61.00246  59.69571 6338  1.021890  0.3069
#> XMonthApr                            -0.06068   0.04786 6338 -1.267781  0.2049
#> XMonthMay                            -0.19018   0.05574 6338 -3.412091  0.0006
#> XMonthJun                            -0.34220   0.05827 6338 -5.872328  0.0000
#> XMonthJul                            -0.41651   0.06143 6338 -6.779815  0.0000
#> XMonthAug                            -0.52427   0.06135 6338 -8.545774  0.0000
#> XMonthSep                            -0.57820   0.05959 6338 -9.703033  0.0000
#> XMonthOct                            -0.54536   0.05569 6338 -9.793228  0.0000
#> XMonthNov                            -0.37561   0.05150 6338 -7.292978  0.0000
#> XMonthDec                            -0.22531   0.06561 6338 -3.434094  0.0006
#> XYear                                 0.13894   0.02227 6338  6.238204  0.0000
#> XSiteS05:Year                        -0.03605   0.03507 6338 -1.027910  0.3040
#> XSiteS17:Year                        -0.06641   0.06422 6338 -1.034136  0.3011
#> XSiteS03:Year                        -0.03494   0.02919 6338 -1.197031  0.2313
#> XSiteS01:Year                        -0.03000   0.02963 6338 -1.012393  0.3114
#> Xti(lPrecip)Fx1                      -0.01459   0.08490 6338 -0.171793  0.8636
#> Xti(wlPrecip)Fx1                     -0.57446   0.05922 6338 -9.699971  0.0000
#> Xti(FlowIndex):SiteS07Fx1            -1.11828   0.16238 6338 -6.886951  0.0000
#> Xti(FlowIndex):SiteS05Fx1            -0.46972   0.14239 6338 -3.298859  0.0010
#> Xti(FlowIndex):SiteS17Fx1            -1.10061   0.17382 6338 -6.331935  0.0000
#> Xti(FlowIndex):SiteS03Fx1            -0.77150   0.16125 6338 -4.784359  0.0000
#> Xti(FlowIndex):SiteS01Fx1            -1.46971   0.17686 6338 -8.310135  0.0000
#> Xti(lPrecip,wlPrecip)Fx1             -0.93837   1.26870 6338 -0.739626  0.4596
#> Xti(lPrecip,FlowIndex)Fx1             0.26739   0.49109 6338  0.544475  0.5861
#> Xti(wlPrecip,FlowIndex)Fx1            0.90814   0.29689 6338  3.058816  0.0022
#> Xti(lPrecip,wlPrecip,FlowIndex)Fx1   -0.33432   1.72996 6338 -0.193252  0.8468
#>  Correlation: 
#>                                    X(Int) XStS05 XStS17 XStS03 XStS01 XMnthAp
#> XSiteS05                           -0.735                                    
#> XSiteS17                           -0.335  0.181                             
#> XSiteS03                           -0.753  0.564  0.185                      
#> XSiteS01                           -0.742  0.551  0.256  0.520               
#> XMonthApr                          -0.031 -0.004  0.004  0.003  0.021        
#> XMonthMay                          -0.047 -0.015 -0.022  0.007  0.022  0.742 
#> XMonthJun                          -0.058 -0.011 -0.001  0.007  0.022  0.660 
#> XMonthJul                          -0.060 -0.014 -0.009  0.007  0.020  0.595 
#> XMonthAug                          -0.056 -0.029 -0.025 -0.004  0.011  0.568 
#> XMonthSep                          -0.078 -0.022 -0.032 -0.007  0.009  0.548 
#> XMonthOct                          -0.077  0.000 -0.030  0.006  0.013  0.537 
#> XMonthNov                          -0.077  0.024  0.002  0.023  0.017  0.527 
#> XMonthDec                          -0.096  0.025  0.027  0.006  0.004  0.450 
#> XYear                              -1.000  0.735  0.335  0.753  0.742  0.030 
#> XSiteS05:Year                       0.735 -1.000 -0.181 -0.564 -0.551  0.004 
#> XSiteS17:Year                       0.335 -0.182 -1.000 -0.185 -0.256 -0.004 
#> XSiteS03:Year                       0.753 -0.564 -0.185 -1.000 -0.520 -0.003 
#> XSiteS01:Year                       0.742 -0.551 -0.256 -0.520 -1.000 -0.021 
#> Xti(lPrecip)Fx1                     0.001  0.000  0.002  0.000  0.001  0.009 
#> Xti(wlPrecip)Fx1                    0.002  0.000  0.014 -0.004 -0.002 -0.020 
#> Xti(FlowIndex):SiteS07Fx1          -0.001  0.004  0.007 -0.001  0.000  0.016 
#> Xti(FlowIndex):SiteS05Fx1           0.003  0.008  0.007 -0.004 -0.003  0.036 
#> Xti(FlowIndex):SiteS17Fx1           0.000  0.000  0.004  0.000 -0.001  0.017 
#> Xti(FlowIndex):SiteS03Fx1           0.000  0.001  0.013  0.001  0.000  0.016 
#> Xti(FlowIndex):SiteS01Fx1          -0.001  0.001  0.006  0.001  0.003  0.016 
#> Xti(lPrecip,wlPrecip)Fx1           -0.001 -0.001 -0.002 -0.001  0.000  0.010 
#> Xti(lPrecip,FlowIndex)Fx1           0.000  0.000  0.009 -0.003 -0.002  0.007 
#> Xti(wlPrecip,FlowIndex)Fx1         -0.001  0.001  0.007 -0.004 -0.002  0.014 
#> Xti(lPrecip,wlPrecip,FlowIndex)Fx1  0.001  0.000  0.006  0.000 -0.001 -0.001 
#>                                    XMnthM XMnthJn XMnthJl XMnthAg XMnthS XMnthO
#> XSiteS05                                                                       
#> XSiteS17                                                                       
#> XSiteS03                                                                       
#> XSiteS01                                                                       
#> XMonthApr                                                                      
#> XMonthMay                                                                      
#> XMonthJun                           0.804                                      
#> XMonthJul                           0.704  0.830                               
#> XMonthAug                           0.654  0.737   0.831                       
#> XMonthSep                           0.623  0.680   0.730   0.824               
#> XMonthOct                           0.599  0.637   0.656   0.712   0.813       
#> XMonthNov                           0.589  0.609   0.604   0.634   0.696  0.775
#> XMonthDec                           0.458  0.478   0.451   0.453   0.473  0.493
#> XYear                               0.046  0.057   0.059   0.056   0.077  0.076
#> XSiteS05:Year                       0.015  0.011   0.014   0.029   0.022  0.000
#> XSiteS17:Year                       0.022  0.001   0.009   0.025   0.032  0.030
#> XSiteS03:Year                      -0.007 -0.007  -0.007   0.004   0.007 -0.006
#> XSiteS01:Year                      -0.022 -0.022  -0.020  -0.011  -0.009 -0.013
#> Xti(lPrecip)Fx1                     0.021 -0.035  -0.024  -0.022  -0.014 -0.016
#> Xti(wlPrecip)Fx1                   -0.021 -0.069  -0.076  -0.057  -0.048 -0.069
#> Xti(FlowIndex):SiteS07Fx1           0.022  0.022   0.032   0.030   0.030  0.019
#> Xti(FlowIndex):SiteS05Fx1           0.047  0.046   0.052   0.052   0.052  0.041
#> Xti(FlowIndex):SiteS17Fx1           0.022  0.021   0.031   0.030   0.030  0.021
#> Xti(FlowIndex):SiteS03Fx1           0.022  0.023   0.031   0.030   0.031  0.018
#> Xti(FlowIndex):SiteS01Fx1           0.021  0.025   0.033   0.032   0.033  0.021
#> Xti(lPrecip,wlPrecip)Fx1            0.010 -0.004  -0.003  -0.003   0.001  0.004
#> Xti(lPrecip,FlowIndex)Fx1           0.020  0.007   0.014   0.014   0.017 -0.002
#> Xti(wlPrecip,FlowIndex)Fx1          0.023  0.014   0.012   0.015   0.019  0.006
#> Xti(lPrecip,wlPrecip,FlowIndex)Fx1  0.000  0.012   0.012   0.011   0.008  0.000
#>                                    XMnthN XMnthD XYear  XSS05: XSS17: XSS03:
#> XSiteS05                                                                    
#> XSiteS17                                                                    
#> XSiteS03                                                                    
#> XSiteS01                                                                    
#> XMonthApr                                                                   
#> XMonthMay                                                                   
#> XMonthJun                                                                   
#> XMonthJul                                                                   
#> XMonthAug                                                                   
#> XMonthSep                                                                   
#> XMonthOct                                                                   
#> XMonthNov                                                                   
#> XMonthDec                           0.559                                   
#> XYear                               0.076  0.095                            
#> XSiteS05:Year                      -0.024 -0.025 -0.735                     
#> XSiteS17:Year                      -0.002 -0.027 -0.335  0.182              
#> XSiteS03:Year                      -0.023 -0.006 -0.753  0.564  0.185       
#> XSiteS01:Year                      -0.017 -0.004 -0.742  0.551  0.256  0.520
#> Xti(lPrecip)Fx1                    -0.007 -0.010 -0.001  0.000 -0.002  0.000
#> Xti(wlPrecip)Fx1                   -0.029 -0.036 -0.002  0.000 -0.014  0.004
#> Xti(FlowIndex):SiteS07Fx1           0.019  0.015  0.001 -0.004 -0.007  0.001
#> Xti(FlowIndex):SiteS05Fx1           0.039  0.032 -0.003 -0.008 -0.007  0.004
#> Xti(FlowIndex):SiteS17Fx1           0.020  0.016  0.000  0.000 -0.004  0.000
#> Xti(FlowIndex):SiteS03Fx1           0.018  0.017  0.000 -0.001 -0.013 -0.001
#> Xti(FlowIndex):SiteS01Fx1           0.020  0.019  0.001 -0.001 -0.006 -0.001
#> Xti(lPrecip,wlPrecip)Fx1            0.004  0.007  0.001  0.001  0.002  0.001
#> Xti(lPrecip,FlowIndex)Fx1           0.004  0.010  0.000  0.000 -0.009  0.003
#> Xti(wlPrecip,FlowIndex)Fx1          0.017  0.012  0.001 -0.001 -0.007  0.004
#> Xti(lPrecip,wlPrecip,FlowIndex)Fx1  0.006  0.007 -0.001  0.000 -0.006  0.000
#>                                    XSS01: Xt(lP)F1 Xt(wP)F1 X(FI):SS07
#> XSiteS05                                                              
#> XSiteS17                                                              
#> XSiteS03                                                              
#> XSiteS01                                                              
#> XMonthApr                                                             
#> XMonthMay                                                             
#> XMonthJun                                                             
#> XMonthJul                                                             
#> XMonthAug                                                             
#> XMonthSep                                                             
#> XMonthOct                                                             
#> XMonthNov                                                             
#> XMonthDec                                                             
#> XYear                                                                 
#> XSiteS05:Year                                                         
#> XSiteS17:Year                                                         
#> XSiteS03:Year                                                         
#> XSiteS01:Year                                                         
#> Xti(lPrecip)Fx1                    -0.001                             
#> Xti(wlPrecip)Fx1                    0.002  0.191                      
#> Xti(FlowIndex):SiteS07Fx1           0.000 -0.026   -0.060             
#> Xti(FlowIndex):SiteS05Fx1           0.003 -0.016   -0.071    0.628    
#> Xti(FlowIndex):SiteS17Fx1           0.001 -0.014   -0.054    0.603    
#> Xti(FlowIndex):SiteS03Fx1           0.000 -0.026   -0.061    0.669    
#> Xti(FlowIndex):SiteS01Fx1          -0.003 -0.038   -0.061    0.593    
#> Xti(lPrecip,wlPrecip)Fx1            0.000  0.400    0.123   -0.013    
#> Xti(lPrecip,FlowIndex)Fx1           0.002  0.084    0.003    0.605    
#> Xti(wlPrecip,FlowIndex)Fx1          0.002 -0.001    0.158    0.637    
#> Xti(lPrecip,wlPrecip,FlowIndex)Fx1  0.001 -0.241   -0.096    0.203    
#>                                    X(FI):SS05 X(FI):SS1 X(FI):SS03 X(FI):SS01
#> XSiteS05                                                                     
#> XSiteS17                                                                     
#> XSiteS03                                                                     
#> XSiteS01                                                                     
#> XMonthApr                                                                    
#> XMonthMay                                                                    
#> XMonthJun                                                                    
#> XMonthJul                                                                    
#> XMonthAug                                                                    
#> XMonthSep                                                                    
#> XMonthOct                                                                    
#> XMonthNov                                                                    
#> XMonthDec                                                                    
#> XYear                                                                        
#> XSiteS05:Year                                                                
#> XSiteS17:Year                                                                
#> XSiteS03:Year                                                                
#> XSiteS01:Year                                                                
#> Xti(lPrecip)Fx1                                                              
#> Xti(wlPrecip)Fx1                                                             
#> Xti(FlowIndex):SiteS07Fx1                                                    
#> Xti(FlowIndex):SiteS05Fx1                                                    
#> Xti(FlowIndex):SiteS17Fx1           0.581                                    
#> Xti(FlowIndex):SiteS03Fx1           0.632      0.606                         
#> Xti(FlowIndex):SiteS01Fx1           0.552      0.514     0.598               
#> Xti(lPrecip,wlPrecip)Fx1           -0.013     -0.022    -0.014     -0.017    
#> Xti(lPrecip,FlowIndex)Fx1           0.577      0.556     0.610      0.532    
#> Xti(wlPrecip,FlowIndex)Fx1          0.605      0.565     0.643      0.618    
#> Xti(lPrecip,wlPrecip,FlowIndex)Fx1  0.202      0.200     0.208      0.191    
#>                                    X(P,P) Xt(lP,FI)F1 Xt(wP,FI)F1
#> XSiteS05                                                         
#> XSiteS17                                                         
#> XSiteS03                                                         
#> XSiteS01                                                         
#> XMonthApr                                                        
#> XMonthMay                                                        
#> XMonthJun                                                        
#> XMonthJul                                                        
#> XMonthAug                                                        
#> XMonthSep                                                        
#> XMonthOct                                                        
#> XMonthNov                                                        
#> XMonthDec                                                        
#> XYear                                                            
#> XSiteS05:Year                                                    
#> XSiteS17:Year                                                    
#> XSiteS03:Year                                                    
#> XSiteS01:Year                                                    
#> Xti(lPrecip)Fx1                                                  
#> Xti(wlPrecip)Fx1                                                 
#> Xti(FlowIndex):SiteS07Fx1                                        
#> Xti(FlowIndex):SiteS05Fx1                                        
#> Xti(FlowIndex):SiteS17Fx1                                        
#> Xti(FlowIndex):SiteS03Fx1                                        
#> Xti(FlowIndex):SiteS01Fx1                                        
#> Xti(lPrecip,wlPrecip)Fx1                                         
#> Xti(lPrecip,FlowIndex)Fx1          -0.012                        
#> Xti(wlPrecip,FlowIndex)Fx1         -0.023  0.603                 
#> Xti(lPrecip,wlPrecip,FlowIndex)Fx1 -0.406  0.301       0.260     
#> 
#> Standardized Within-Group Residuals:
#>        Min         Q1        Med         Q3        Max 
#> -9.7084084 -0.5519914  0.1366988  0.6879176  3.4128525 
#> 
#> Number of Observations: 6368
#> Number of Groups: 
#>                                                                                           g 
#>                                                                                           1 
#>                                                                                  g.0 %in% g 
#>                                                                                           1 
#>                                                                         g.1 %in% g.0 %in% g 
#>                                                                                           1 
#>                                                                g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1 
#>                                                       g.3 %in% g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1 
#>                                              g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1 
#>                                     g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1 
#>                            g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1 
#>                   g.7 %in% g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1 
#>          g.8 %in% g.7 %in% g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1 
#> g.9 %in% g.8 %in% g.7 %in% g.6 %in% g.5 %in% g.4 %in% g.3 %in% g.2 %in% g.1 %in% g.0 %in% g 
#>                                                                                           1
```

### Visualizing Estimated Marginal Means

``` r
the_call <-  quote(gamm(log(Chl_Median) ~ Site + 
                      ti(lPrecip) +
                      ti(wlPrecip) + 
                      ti(FlowIndex, by = Site) +
                      ti(lPrecip,  wlPrecip) +
                      ti(lPrecip, FlowIndex) +
                      ti(wlPrecip, FlowIndex) +
                      ti(lPrecip, wlPrecip, FlowIndex) +
                      Month +
                      Year + 
                      Site:Year,
                    correlation = corAR1(0.8),
                    na.action = na.omit,
                    method = 'REML',
                    data = reduced_data))

full_gamm$gam$call <- the_call

my_ref_grid <- ref_grid(full_gamm, at = list(Year = 2014), cov.reduce = median) 
a <- emmeans(my_ref_grid, ~ Month, type = 'response')
```

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Sites Combined'

plot(a) + 
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab ('') +
  annotate('text', 400, 6, label = labl, size = 3) +
  xlim(0,500) +
  geom_vline(xintercept =  230, color = 'orange') +
  geom_vline(xintercept =  860, color = 'red') +
  
  coord_flip() +
  
  theme_cbep(base_size = 12)
#> Warning: Removed 1 rows containing missing values (geom_vline).
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-65-1.png" style="display: block; margin: auto;" />

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Dates Combined'

a <- emmeans(my_ref_grid, ~ Site, type = 'response')
#> NOTE: Results may be misleading due to involvement in interactions
plot(a) + 
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab("Upstream                  Main Stem                                 Lower Tribs                   ") +
  annotate('text', 400, 2.5, label = labl, size = 3) +
  xlim(0,500) +
  geom_vline(xintercept =  230, color = 'orange') +
  geom_vline(xintercept =  860, color = 'red') +
  coord_flip() +
  theme_cbep(base_size = 12)
#> Warning: Removed 1 rows containing missing values (geom_vline).
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-66-1.png" style="display: block; margin: auto;" />

### Visualizing Trends

We extract results on the log scale, so we can calculate the linear
predictor by hand, then back transform.

``` r
my_ref_grid <- ref_grid(full_gamm, at = list(Year = 2014, Month = 'Jul'),
                        cov.reduce = median)
a <- summary(emmeans(my_ref_grid, 'Site'))
#> NOTE: Results may be misleading due to involvement in interactions
b <- summary(emtrends(full_gamm, 'Site', 'Year'))

lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])
rm(a,b)

df <- tibble(Site = rep(levels(reduced_data$Site), each = 9), 
              Year = rep(2010:2018, 5)) %>%
  mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
         iintercept = lookup$Intercept[match(Site, lookup$Site)],
         pred = exp((Year - 2014) * sslope + iintercept)) %>%
  select(-sslope, -iintercept)

ggplot(df, aes(x = Year, y = pred, color = Site)) +
         geom_step(direction = 'mid') +
  ylab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  xlab('') +
  ylim(0,600) +
  geom_hline(yintercept =  230, color = 'black') +
  #geom_hline(yintercept =  860, color = 'red') +

  theme_cbep(base_size = 12)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-67-1.png" style="display: block; margin: auto;" />

## Simplified Models

We examined numerous alternative model structures, on a site by site
basis, only to find predictions of most reasonable models were highly
correlated, with a correlation coefficient over 0.98, and often over
0.99. Thus for most purposes, the performance differences between “best”
models and “good” models will be small. That matters principally because
these complex models are so slow to run that we can not practically
explore multiple model alternatives.

Here we develop models that omits all interaction terms. This is the
simplest model that incorporates different site by FlowIndex smoothers,
making it more directly comparable to the Full model

This took about 20 minutes to run.

``` r
if (! file.exists("Models/alt_gamm_1.rds")) {
  alt_gamm_1 <- gamm(log(Chl_Median) ~ Site + 
                       s(lPrecip) +
                       s(wlPrecip) + 
                       s(FlowIndex, by = Site) +
                       Month +
                       Year + 
                       Site:Year,
                     correlation = corAR1(0.8),
                     na.action = na.omit,
                     method = 'REML',
                     data = reduced_data)
  saveRDS(alt_gamm_1, file="models/alt_gamm_1.rds")
} else {
  alt_gamm_1 <- readRDS("models/alt_gamm_1.rds")
}
```

``` r
anova(alt_gamm_1$gam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> log(Chl_Median) ~ Site + s(lPrecip) + s(wlPrecip) + s(FlowIndex, 
#>     by = Site) + Month + Year + Site:Year
#> 
#> Parametric Terms:
#>           df      F  p-value
#> Site       4  1.367    0.243
#> Month      9 19.636  < 2e-16
#> Year       1 44.449 2.83e-11
#> Site:Year  4  1.358    0.246
#> 
#> Approximate significance of smooth terms:
#>                        edf Ref.df      F p-value
#> s(lPrecip)           7.019  7.019  28.52  <2e-16
#> s(wlPrecip)          4.342  4.342  60.91  <2e-16
#> s(FlowIndex):SiteS07 8.021  8.021 155.25  <2e-16
#> s(FlowIndex):SiteS05 6.897  6.897  24.67  <2e-16
#> s(FlowIndex):SiteS17 7.068  7.068  51.60  <2e-16
#> s(FlowIndex):SiteS03 6.947  6.947  38.53  <2e-16
#> s(FlowIndex):SiteS01 7.320  7.320 187.83  <2e-16
```

``` r
plot(alt_gamm_1$gam)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-69-1.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-69-2.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-69-3.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-69-4.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-69-5.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-69-6.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-69-7.png" style="display: block; margin: auto;" />
Those curves show fairly substantial differences in slope for the
relationship between FlowIndex and Chlorides from site to site. However,
over the majority of the range of the data, the relationships are not
too far from linear. This again suggests that and some inconsistency of
response to increasing FlowIndex at each site.

Some of the apparent effect is because these curves are fit to data with
different depth ranges at each site, but are plotted across the entire
observed dept range, thus leading the curves to include spurious values
outside the range of observed depths. In other words, the relationships
that matter within each Site are close to linear, with just a bit of
curvature tossed in, which is being overfit near the ends of the
splines, where data is sparse.

``` r
reduced_data %>%
  select(Site, D_Median, lD_Median) %>%
  group_by(Site) %>%
  summarize(D_min = min(D_Median, na.rm = TRUE),
            D_max = max(D_Median, na.rm = TRUE),
            lD_min = min(lD_Median, na.rm = TRUE),
            lD_max = max(lD_Median, na.rm = TRUE))
#> `summarise()` ungrouping output (override with `.groups` argument)
#> # A tibble: 5 x 5
#>   Site   D_min D_max lD_min lD_max
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 S07   0      1.02  0       0.704
#> 2 S05   0.024  1.62  0.0237  0.961
#> 3 S17   0.0195 1.27  0.0193  0.820
#> 4 S03   0.092  1.20  0.0880  0.788
#> 5 S01   0.039  0.924 0.0383  0.654
```

This is similar to the previous model, but using a Day of Year (‘DOY’)
smoothed predictor instead of the month by month factor terms. This
model is also slow to run….

``` r
if (! file.exists("Models/alt_gamm_2.rds")) {
  alt_gamm_2 <- gamm(log(Chl_Median) ~ Site + 
                       s(lPrecip) +
                       s(wlPrecip) + 
                       s(FlowIndex, by = Site) +
                       s(DOY) +
                       Year +
                       Site:Year,
                     correlation = corAR1(0.8),
                     na.action = na.omit,
                     method = 'REML',
                     data = reduced_data)
  saveRDS(alt_gamm_2, file="models/alt_gamm_2.rds")
} else {
  alt_gamm_2 <- readRDS("models/alt_gamm_2.rds")
}
```

``` r
plot(alt_gamm_2$gam)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-1.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-2.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-3.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-4.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-5.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-6.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-7.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-71-8.png" style="display: block; margin: auto;" />

The effect of time of year remains whether expressed as monthly
parameters or via a smoothing term. The other smoothing terms have not
changed appreciably. The two models are functionally similar. The Month
by Month version provides slightly better graphics.

### Visualizing Estimated Marginal Means

Reliably calling `emmeans()` for these large `gamm()` models appears to
require creating a call object and associating it with the model (e.g.,
as `the_gamm$gam$call`). (See the `emmeans` models vignette for more
info, although not all strategies recommended there worked for us).

We first create the call object, then associate it with the model, and
finally manually construct a reference grid before calling `emmeans()`
to extract marginal means. This workflow has the advantage that it
requires us to think carefully about the structure of the reference
grid.

``` r
the_call <-  quote(gamm(log(Chl_Median) ~ Site + 
                       s(lPrecip) +
                       s(wlPrecip) + 
                       s(FlowIndex, by = Site) +
                       Month +
                       Year + 
                       Site:Year,
                     correlation = corAR1(0.8),
                     na.action = na.omit,
                     method = 'REML',
                     data = reduced_data))

alt_gamm_1$gam$call <- the_call

my_ref_grid <- ref_grid(alt_gamm_1, at = list(Year = 2014), cov.reduce = median) 
a <- emmeans(my_ref_grid, ~ Month, type = 'response')

labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Sites Combined'

plot(a) + 
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab ('') +
  annotate('text', 400, 6, label = labl, size = 3) +
  xlim(0,500) +
  geom_vline(xintercept =  230, color = 'orange') +
  geom_vline(xintercept =  860, color = 'red') +
  
  coord_flip() +
  
  theme_cbep(base_size = 12)
#> Warning: Removed 1 rows containing missing values (geom_vline).
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-72-1.png" style="display: block; margin: auto;" />

Note that all estimated monthly means are somewhat lower than from the
GLS model.

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Dates Combined'

a <- emmeans(my_ref_grid, ~ Site, type = 'response')
#> NOTE: Results may be misleading due to involvement in interactions
plot(a) + 
  xlab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  ylab("Upstream                  Main Stem                                 Lower Tribs                   ") +
  annotate('text', 400, 2.5, label = labl, size = 3) +
  xlim(0,500) +
  geom_vline(xintercept =  230, color = 'orange') +
  geom_vline(xintercept =  860, color = 'red') +
  coord_flip() +
  theme_cbep(base_size = 12)
#> Warning: Removed 1 rows containing missing values (geom_vline).
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-73-1.png" style="display: block; margin: auto;" />

### Visualizing Trends

We extract results on the log scale, so we can calculate the linear
predictor by hand, then back transform.

``` r
my_ref_grid <- ref_grid(alt_gamm_1, at = list(Year = 2014, Month = 'Jul'),
                        cov.reduce = median)

a <- summary(emmeans(my_ref_grid, 'Site'))
#> NOTE: Results may be misleading due to involvement in interactions

b <- summary(emtrends(alt_gamm_1, 'Site', 'Year'))

lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])
rm(a,b)

df <- tibble(Site = rep(levels(reduced_data$Site), each = 9), 
              Year = rep(2010:2018, 5)) %>%
  mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
         iintercept = lookup$Intercept[match(Site, lookup$Site)],
         pred = exp((Year - 2014) * sslope + iintercept)) %>%
  select(-sslope, -iintercept)
```

``` r
ggplot(df, aes(x = Year, y = pred, color = Site)) +
         geom_step(direction = 'mid') +
  ylab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  xlab('') +
  ylim(0,600) +
  geom_hline(yintercept =  230, color = 'black') +
  #geom_hline(yintercept =  860, color = 'red') +

  theme_cbep(base_size = 12)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-75-1.png" style="display: block; margin: auto;" />

# Seeking a Simplifies GLS Model

The GLS models run much more quickly, which is considerably more
convenient that dealing with full GAMM models. Also, there is value in
using a model that does NOT include any interaction terms between
covariates and Site, if only to avoid warning labels on output.

Most to the important nonlinearities in the covariate smoothers from the
GAMM models occur at high rainfall and high flow conditions. We can
generate a simpler linear model that accounts for different behavior at
high flows by including appropriate indicator variables.

``` r
reduced_data <- reduced_data %>%
  mutate(high_precip =  Precip > 4,
         high_flow = FlowIndex > 0.65)

highflow_gls_1 <- gls(log(Chl_Median) ~ Site + 
                        lPrecip +
                        high_precip +
                        wlPrecip +
                        FlowIndex +
                        high_flow +
                        high_flow : lPrecip +
                        high_flow : wlPrecip +
                        Month +
                        Year +
                        Site:Year,
                    correlation = corAR1(0.8),
                    na.action = na.omit, 
                    method = 'ML',
                    data = reduced_data)
```

``` r
anova(highflow_gls_1)
#> Denom. DF: 6342 
#>                    numDF  F-value p-value
#> (Intercept)            1 64935.93  <.0001
#> Site                   4    40.89  <.0001
#> lPrecip                1     6.02  0.0142
#> high_precip            1     8.65  0.0033
#> wlPrecip               1  2293.28  <.0001
#> FlowIndex              1  1456.75  <.0001
#> high_flow              1     8.49  0.0036
#> Month                  9    27.55  <.0001
#> Year                   1   112.38  <.0001
#> lPrecip:high_flow      1    72.09  <.0001
#> wlPrecip:high_flow     1   133.86  <.0001
#> Site:Year              4     4.09  0.0026
```

``` r
summary(highflow_gls_1)
#> Generalized least squares fit by maximum likelihood
#>   Model: log(Chl_Median) ~ Site + lPrecip + high_precip + wlPrecip + FlowIndex +      high_flow + high_flow:lPrecip + high_flow:wlPrecip + Month +      Year + Site:Year 
#>   Data: reduced_data 
#>         AIC      BIC   logLik
#>   -140.7574 48.49575 98.37869
#> 
#> Correlation Structure: AR(1)
#>  Formula: ~1 
#>  Parameter estimate(s):
#>      Phi 
#> 0.862947 
#> 
#> Coefficients:
#>                             Value Std.Error   t-value p-value
#> (Intercept)            -271.56323  35.33393  -7.68562  0.0000
#> SiteS05                  89.63723  56.79022   1.57839  0.1145
#> SiteS17                 219.00562 110.68982   1.97855  0.0479
#> SiteS03                 177.49647  46.85048   3.78857  0.0002
#> SiteS01                 112.16071  47.40778   2.36587  0.0180
#> lPrecip                  -0.02987   0.00570  -5.24054  0.0000
#> high_precipTRUE           0.00547   0.01384   0.39489  0.6929
#> wlPrecip                 -0.08429   0.00510 -16.52633  0.0000
#> FlowIndex                -2.02657   0.05266 -38.48368  0.0000
#> high_flowTRUE             1.37800   0.09586  14.37517  0.0000
#> MonthApr                 -0.11724   0.05019  -2.33603  0.0195
#> MonthMay                 -0.28182   0.05674  -4.96662  0.0000
#> MonthJun                 -0.47251   0.05877  -8.03938  0.0000
#> MonthJul                 -0.54242   0.06148  -8.82223  0.0000
#> MonthAug                 -0.67040   0.06135 -10.92708  0.0000
#> MonthSep                 -0.72600   0.05954 -12.19358  0.0000
#> MonthOct                 -0.61051   0.05665 -10.77770  0.0000
#> MonthNov                 -0.43997   0.05313  -8.28059  0.0000
#> MonthDec                 -0.29960   0.06947  -4.31263  0.0000
#> Year                      0.13800   0.01754   7.86974  0.0000
#> lPrecip:high_flowTRUE    -0.17902   0.01430 -12.51815  0.0000
#> wlPrecip:high_flowTRUE   -0.16892   0.01460 -11.57320  0.0000
#> SiteS05:Year             -0.04462   0.02818  -1.58350  0.1134
#> SiteS17:Year             -0.10866   0.05489  -1.97960  0.0478
#> SiteS03:Year             -0.08792   0.02325  -3.78065  0.0002
#> SiteS01:Year             -0.05536   0.02353  -2.35283  0.0187
#> 
#>  Correlation: 
#>                        (Intr) SitS05 SitS17 SitS03 SitS01 lPrecp hgh_pTRUE
#> SiteS05                -0.684                                             
#> SiteS17                -0.313  0.174                                      
#> SiteS03                -0.747  0.517  0.195                               
#> SiteS01                -0.738  0.508  0.234  0.532                        
#> lPrecip                -0.001  0.005  0.012  0.001  0.000                 
#> high_precipTRUE        -0.002 -0.003  0.001  0.000  0.001 -0.836          
#> wlPrecip               -0.006  0.010  0.038  0.003  0.002  0.360 -0.013   
#> FlowIndex               0.003 -0.005 -0.002 -0.001 -0.002 -0.221  0.028   
#> high_flowTRUE           0.002  0.003  0.002  0.000  0.001  0.131 -0.039   
#> MonthApr               -0.024 -0.010 -0.010  0.004  0.014 -0.031  0.009   
#> MonthMay               -0.041 -0.018 -0.036  0.010  0.019 -0.037  0.006   
#> MonthJun               -0.047 -0.013 -0.017  0.009  0.018 -0.036 -0.004   
#> MonthJul               -0.044 -0.019 -0.023  0.005  0.015 -0.056  0.027   
#> MonthAug               -0.037 -0.033 -0.042 -0.005  0.007 -0.057  0.021   
#> MonthSep               -0.058 -0.028 -0.045 -0.008  0.004 -0.039  0.019   
#> MonthOct               -0.056 -0.008 -0.036 -0.001  0.007 -0.037 -0.010   
#> MonthNov               -0.060  0.010 -0.009  0.011  0.008 -0.023  0.001   
#> MonthDec               -0.084  0.018  0.012 -0.003 -0.006 -0.048  0.014   
#> Year                   -1.000  0.685  0.313  0.747  0.738  0.001  0.002   
#> lPrecip:high_flowTRUE  -0.003  0.001 -0.001  0.001  0.000 -0.214  0.132   
#> wlPrecip:high_flowTRUE -0.001 -0.004 -0.003 -0.001  0.000 -0.059 -0.015   
#> SiteS05:Year            0.684 -1.000 -0.174 -0.517 -0.508 -0.005  0.003   
#> SiteS17:Year            0.314 -0.174 -1.000 -0.196 -0.234 -0.012 -0.001   
#> SiteS03:Year            0.747 -0.517 -0.195 -1.000 -0.532 -0.001  0.000   
#> SiteS01:Year            0.738 -0.508 -0.234 -0.532 -1.000  0.000 -0.001   
#>                        wlPrcp FlwInd hgh_fTRUE MnthAp MnthMy MnthJn MnthJl
#> SiteS05                                                                   
#> SiteS17                                                                   
#> SiteS03                                                                   
#> SiteS01                                                                   
#> lPrecip                                                                   
#> high_precipTRUE                                                           
#> wlPrecip                                                                  
#> FlowIndex              -0.549                                             
#> high_flowTRUE           0.229 -0.395                                      
#> MonthApr               -0.061  0.077 -0.033                               
#> MonthMay               -0.072  0.100 -0.044     0.729                     
#> MonthJun               -0.110  0.132 -0.056     0.659  0.797              
#> MonthJul               -0.097  0.157 -0.068     0.607  0.706  0.816       
#> MonthAug               -0.095  0.160 -0.073     0.591  0.671  0.733  0.816
#> MonthSep               -0.077  0.142 -0.073     0.580  0.654  0.693  0.728
#> MonthOct               -0.128  0.124 -0.039     0.571  0.635  0.662  0.670
#> MonthNov               -0.071  0.081 -0.034     0.553  0.624  0.634  0.626
#> MonthDec               -0.093  0.077 -0.034     0.453  0.465  0.481  0.457
#> Year                    0.006 -0.003 -0.002     0.023  0.040  0.046  0.043
#> lPrecip:high_flowTRUE  -0.068  0.097 -0.673     0.006  0.010  0.008  0.011
#> wlPrecip:high_flowTRUE -0.222  0.303 -0.887     0.029  0.038  0.050  0.056
#> SiteS05:Year           -0.010  0.005 -0.003     0.010  0.018  0.013  0.019
#> SiteS17:Year           -0.038  0.002 -0.002     0.010  0.036  0.017  0.023
#> SiteS03:Year           -0.003  0.001  0.000    -0.004 -0.010 -0.009 -0.005
#> SiteS01:Year           -0.002  0.002 -0.001    -0.014 -0.019 -0.018 -0.015
#>                        MnthAg MnthSp MnthOc MnthNv MnthDc Year   lP:_TR wP:_TR
#> SiteS05                                                                       
#> SiteS17                                                                       
#> SiteS03                                                                       
#> SiteS01                                                                       
#> lPrecip                                                                       
#> high_precipTRUE                                                               
#> wlPrecip                                                                      
#> FlowIndex                                                                     
#> high_flowTRUE                                                                 
#> MonthApr                                                                      
#> MonthMay                                                                      
#> MonthJun                                                                      
#> MonthJul                                                                      
#> MonthAug                                                                      
#> MonthSep                0.810                                                 
#> MonthOct                0.711  0.804                                          
#> MonthNov                0.644  0.697  0.761                                   
#> MonthDec                0.456  0.473  0.486  0.549                            
#> Year                    0.036  0.057  0.055  0.059  0.084                     
#> lPrecip:high_flowTRUE   0.010 -0.004  0.023  0.013  0.011  0.003              
#> wlPrecip:high_flowTRUE  0.061  0.066  0.027  0.024  0.028  0.001  0.414       
#> SiteS05:Year            0.033  0.028  0.008 -0.010 -0.018 -0.685 -0.001  0.004
#> SiteS17:Year            0.042  0.045  0.036  0.009 -0.012 -0.314  0.001  0.003
#> SiteS03:Year            0.005  0.008  0.001 -0.011  0.003 -0.747 -0.001  0.001
#> SiteS01:Year           -0.007 -0.004 -0.007 -0.008  0.006 -0.738  0.000  0.000
#>                        SS05:Y SS17:Y SS03:Y
#> SiteS05                                    
#> SiteS17                                    
#> SiteS03                                    
#> SiteS01                                    
#> lPrecip                                    
#> high_precipTRUE                            
#> wlPrecip                                   
#> FlowIndex                                  
#> high_flowTRUE                              
#> MonthApr                                   
#> MonthMay                                   
#> MonthJun                                   
#> MonthJul                                   
#> MonthAug                                   
#> MonthSep                                   
#> MonthOct                                   
#> MonthNov                                   
#> MonthDec                                   
#> Year                                       
#> lPrecip:high_flowTRUE                      
#> wlPrecip:high_flowTRUE                     
#> SiteS05:Year                               
#> SiteS17:Year            0.174              
#> SiteS03:Year            0.517  0.196       
#> SiteS01:Year            0.508  0.234  0.532
#> 
#> Standardized residuals:
#>         Min          Q1         Med          Q3         Max 
#> -11.8717012  -0.5247764   0.1307798   0.6532462   3.5247628 
#> 
#> Residual standard error: 0.471479 
#> Degrees of freedom: 6368 total; 6342 residual
```

``` r
plot(highflow_gls_1)
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-79-1.png" style="display: block; margin: auto;" />

# Comparison Of Estimated Marginal Means

A quick review of the estimated marginal means from all of these models
show very similar results for all of the GLS and GAMM models. We gather
marginal means for two GLS models and two GAM models to compare.

## Calculate Monthly Marginal Means

### the\_gls

``` r
month.emm_gls_1 <- summary(emmeans(the_gls, ~ Month,
                           at = list(Year = 2014),
                           cov.reduce = median,
                           type = 'response',
                           mode = "df.error"))
```

### selected\_gls

``` r
month.emm_gls_2 <- summary(emmeans(selected_gls, ~ Month,
                           at = list(Year = 2014),
                           cov.reduce = median,
                           type = 'response',
                           mode = "df.error"))
```

### highflow\_gls\_1

``` r
month.emm_gls_3 <- summary(emmeans(highflow_gls_1, ~ Month,
                           at = list(Year = 2014),
                           cov.reduce = median,
                           type = 'response',
                           mode = "df.error"))
```

### the\_gamm

``` r
the_call <-  quote(gamm(log(Chl_Median) ~ Site + 
                          s(lPrecip) + 
                          s(wlPrecip) +
                          s(FlowIndex) +
                          Month +
                          Year +
                          Site : Year,
                        correlation = corAR1(0.8),
                        na.action = na.omit, 
                        method = 'REML',
                        data = reduced_data))

the_gamm$gam$call <- the_call

my_ref_grid <- ref_grid(the_gamm, at = list(Year = 2014), cov.reduce = median)

month.emm_gamm_1 <- summary(emmeans(my_ref_grid, ~ Month, type = 'response'))
```

### alt\_gamm\_1

``` r
the_call <-  quote(gamm(log(Chl_Median) ~ Site + 
                       s(lPrecip) +
                       s(wlPrecip) + 
                       s(FlowIndex, by = Site) +
                       Month +
                       Year + 
                       Site:Year,
                     correlation = corAR1(0.8),
                     na.action = na.omit,
                     method = 'REML',
                     data = reduced_data))

alt_gamm_1$gam$call <- the_call

my_ref_grid <- ref_grid(alt_gamm_1, at = list(Year = 2014), cov.reduce = median)

month.emm_gamm_2 <- summary(emmeans(my_ref_grid, ~ Month, type = 'response'))
```

### full\_gamm

``` r
the_call <-  quote(gamm(log(Chl_Median) ~ Site + 
                      ti(lPrecip) +
                      ti(wlPrecip) + 
                      ti(FlowIndex, by = Site) +
                      ti(lPrecip,  wlPrecip) +
                      ti(lPrecip, FlowIndex) +
                      ti(wlPrecip, FlowIndex) +
                      ti(lPrecip, wlPrecip, FlowIndex) +
                      Month +
                      Year + 
                      Site:Year,
                    correlation = corAR1(0.8),
                    na.action = na.omit,
                    method = 'REML',
                    data = reduced_data))

full_gamm$gam$call <- the_call

my_ref_grid <- ref_grid(full_gamm,
                        at = list(Year = 2014),
                        cov.reduce = median)

month.emm_gamm_3 <- summary(emmeans(my_ref_grid, ~ Month, type = 'response'))
```

## Compare Results

### Look at Correlations

``` r
a <- cbind(month.emm_gls_1$response,
      month.emm_gls_2$response,
      month.emm_gls_2$response,
      month.emm_gamm_1$response,
      month.emm_gamm_2$response,
      month.emm_gamm_3$response)
colnames(a) <- c('GLS', 'Int_GLS', 'Hi_Flow_GLS',
              'Simple_Gam', 'Int_GAM', 'Full_Gam')
cor(a)
#>                   GLS   Int_GLS Hi_Flow_GLS Simple_Gam   Int_GAM  Full_Gam
#> GLS         1.0000000 0.9979218   0.9979218  0.9987213 0.9981917 0.9947305
#> Int_GLS     0.9979218 1.0000000   1.0000000  0.9972925 0.9988194 0.9982829
#> Hi_Flow_GLS 0.9979218 1.0000000   1.0000000  0.9972925 0.9988194 0.9982829
#> Simple_Gam  0.9987213 0.9972925   0.9972925  1.0000000 0.9992152 0.9953991
#> Int_GAM     0.9981917 0.9988194   0.9988194  0.9992152 1.0000000 0.9973084
#> Full_Gam    0.9947305 0.9982829   0.9982829  0.9953991 0.9973084 1.0000000
```

Results are HIGHLY correlated.

### Look at Maximum Percent Differences

``` r
as_tibble(a) %>%
  mutate(Month = month.abb[3:12]) %>%
  rowwise() %>%
  mutate(m =  mean(c_across(GLS:Full_Gam)),
         d =  max(m - (c_across(GLS:Full_Gam)))) %>%
  summarize (Month = first(Month),
             pct_dif = round(d/m,3) * 100,
             .groups = 'drop')
#> # A tibble: 10 x 2
#>    Month pct_dif
#>    <chr>   <dbl>
#>  1 Mar      9.7 
#>  2 Apr      7.8 
#>  3 May      6.6 
#>  4 Jun      4.2 
#>  5 Jul      5   
#>  6 Aug      4.40
#>  7 Sep      4.6 
#>  8 Oct      7.1 
#>  9 Nov      7.6 
#> 10 Dec      6.5
```

Maximum differences never exceed 10% of the estimates. While we don’t
show it here, the error bars of all estimates overlap.

### Look at Graphic Comparison of Point Estimates

``` r
as_tibble(a) %>%
  mutate(Month = month.abb[3:12],
         Month = factor(Month, levels = month.abb[3:12])) %>%
  pivot_longer(-Month, names_to = 'Which', values_to = 'Value') %>%
  ggplot(aes(x = Month, y = Value, color = Which)) +
  geom_point() +
  geom_line(aes(x = as.numeric(Month))) +
  theme(axis.text.x = element_text(angle = 90, vjust = .25))
```

<img src="Chloride_Analysis_files/figure-gfm/unnamed-chunk-88-1.png" style="display: block; margin: auto;" />

GLS models consistently estimate slightly higher marginal means than the
GAMM models. In general, adding more parameters to the model appears to
flatten out differences between months.

# Conclusion

Any of these models does an adequate job of summarize our results. There
is value is sticking with simpler, faster to compute, easier to explain
models.

The simplest GAM model we explored appears to be the best compromise
between sophistication and simplicity.

``` r
the_gls$call
#> gls(model = log(Chl_Median) ~ Site + lPrecip + wlPrecip + FlowIndex + 
#>     Month + Year + Site:Year, data = reduced_data, correlation = corAR1(0.8), 
#>     method = "REML", na.action = na.omit)
the_gamm$gam$call
#> gamm(log(Chl_Median) ~ Site + s(lPrecip) + s(wlPrecip) + s(FlowIndex) + 
#>     Month + Year + Site:Year, correlation = corAR1(0.8), na.action = na.omit, 
#>     method = "REML", data = reduced_data)
```
