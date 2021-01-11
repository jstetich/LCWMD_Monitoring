Summary of Analysis of LCWMD ‘Chloride’ Data
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
  - [GAMM Analysis](#gamm-analysis)
      - [Initial Model](#initial-model)
      - [ANOVA](#anova)
      - [Summary](#summary)
      - [Structure of the GAM](#structure-of-the-gam)
      - [Diagnostic Plots](#diagnostic-plots)
      - [Checking Estimated Marginal
        Means](#checking-estimated-marginal-means)
      - [Visualizing Trends](#visualizing-trends)
  - [Model without the interactions.](#model-without-the-interactions.)
      - [ANOVA](#anova-1)
      - [Summary](#summary-1)
      - [Structure of the GAM](#structure-of-the-gam-1)
      - [Diagnostic Plots](#diagnostic-plots-1)
  - [Model with Separate Years](#model-with-separate-years)
      - [ANOVA](#anova-2)
      - [Diagnostic Plots](#diagnostic-plots-2)
  - [Compare AIC](#compare-aic)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

Chlorides are frequently elevated in Maine urban streams because of use
of salt for deicing of roads, parking areas, sidewalks, and paths in
winter. While water quality standards are often written in terms of
chlorides, it may be better to think about chlorides as a measure of the
salinity. Physiologically, it is probably salinity or osmolarity that
principally affects organisms in the stream, not chlorides *per se*. The
data we examine here is based on measurement of conductivity, which is
converted to an estimate of in-stream chlorides based on a robust
regression relationship developed over several years.

This R Notebook reviews the model we end up using to analyze chloride
levels in Long Creek. We examined numerous models before settling on
this one. Details of some of those models is available in the
“Chloride\_Analysis.Rmd” notebook.

Several more complex models are “better” using conventional measures of
statistical significance or information criteria. We selected a slightly
simpler model, largely as it makes explaining the model more direct.

Our interest focuses on answering three questions:  
1\. What is the effect of time of year (Month, or Day of Year) on
chlorides?  
2\. Do chloride levels differ from site to site?  
3\. Is there a long-term trend in chlorides? 3. Are there differences in
the magnitude of the trend from site to site?

We use a Generalized Additive Model, with autocorrelated errors to
explore these questions. The model has the following form:

\[ 
\begin{align}
log(Chlorides) &= f(Covariates) + \\
&\qquad \beta_{1,i} Site_i + 
\beta_{2,j} Month_j + \beta_3 Year + \beta_{4,i} Site_i * Year + \epsilon
\end{align}
\]

Where: \* covariates include three terms:  
– Daily precipitation  
– Weighted precipitation from the prior nine days  
– Stream flow in the middle of the watershed  
\* The core predictors enter the model as standard linear terms  
\* The error i an AR(1) correlated error.

We abuse the autocorrelation models slightly, since we use sequential
autocorrelations (not time-based) and we don’t fit separate
autocorrelations for each site and season. That should have little
impact on results, as transitions are relatively rare in a dense data
set, and missing values at the beginning of each season at each site
prevent estimation near season and site transitions in the sequential
data anyway.

On the whole, this models is OK, but not great. It has heavy tailed,
skewed residuals. We should not trust the asymptotic p values. But since
sample sizes are large and results tend to have high statistical
significance, p values are not much use anyway.

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

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

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

## ANOVA

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

## Summary

``` r
summary(the_gamm$gam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> log(Chl_Median) ~ Site + s(lPrecip) + s(wlPrecip) + s(FlowIndex) + 
#>     Month + Year + Site:Year
#> 
#> Parametric coefficients:
#>                Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)  -274.07889   39.67904  -6.907 5.42e-12 ***
#> SiteS05        66.50572   63.29515   1.051  0.29343    
#> SiteS17       170.42246  120.32720   1.416  0.15673    
#> SiteS03       141.72115   52.36986   2.706  0.00682 ** 
#> SiteS01        92.93459   53.01566   1.753  0.07966 .  
#> MonthApr       -0.10089    0.04981  -2.026  0.04285 *  
#> MonthMay       -0.25274    0.05681  -4.449 8.77e-06 ***
#> MonthJun       -0.43644    0.05938  -7.349 2.24e-13 ***
#> MonthJul       -0.51555    0.06255  -8.242  < 2e-16 ***
#> MonthAug       -0.63148    0.06251 -10.102  < 2e-16 ***
#> MonthSep       -0.67389    0.06036 -11.164  < 2e-16 ***
#> MonthOct       -0.62198    0.05703 -10.905  < 2e-16 ***
#> MonthNov       -0.40560    0.05294  -7.662 2.10e-14 ***
#> MonthDec       -0.27603    0.06888  -4.007 6.21e-05 ***
#> Year            0.13892    0.01969   7.055 1.91e-12 ***
#> SiteS05:Year   -0.03316    0.03141  -1.056  0.29118    
#> SiteS17:Year   -0.08457    0.05967  -1.417  0.15643    
#> SiteS03:Year   -0.07016    0.02599  -2.699  0.00697 ** 
#> SiteS01:Year   -0.04583    0.02631  -1.742  0.08160 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                edf Ref.df      F p-value    
#> s(lPrecip)   7.035  7.035  26.80  <2e-16 ***
#> s(wlPrecip)  4.142  4.142  57.59  <2e-16 ***
#> s(FlowIndex) 8.459  8.459 231.41  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.464   
#>   Scale est. = 0.24523   n = 6368
```

## Structure of the GAM

``` r
plot(the_gamm$gam)
```

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-17-1.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-17-2.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-17-3.png" style="display: block; margin: auto;" />
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

## Diagnostic Plots

The help files for `gam.check()` suggest using care when interpreting
results for GAMM models, since the function does not correctly
incorporate the error correlations structure. However, for our purposes,
this is probably sufficient, since our focus is not on statistical
significance, but on estimation.

``` r
gam.check(the_gamm$gam)
```

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

    #> 
    #> 'gamm' based fit - care required with interpretation.
    #> Checks based on working residuals may be misleading.
    #> Basis dimension (k) checking results. Low p-value (k-index<1) may
    #> indicate that k is too low, especially if edf is close to k'.
    #> 
    #>                k'  edf k-index p-value    
    #> s(lPrecip)   9.00 7.03    0.97   0.045 *  
    #> s(wlPrecip)  9.00 4.14    0.91  <2e-16 ***
    #> s(FlowIndex) 9.00 8.46    0.88  <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

As with the linear model, we have a skewed, slightly heavy tailed
distribution of residuals, with a couple of very large outliers. There
is perhaps slight evidence for lack of complete independence between
residuals and predictors. T his model is adequate, but not great. For
careful work, we should probably use bootstrapped confidence intervals
or something similar, but for our purposes, that is probably overkill.

## Checking Estimated Marginal Means

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
(a <- emmeans(my_ref_grid, ~ Month, type = 'response'))
#>  Month response    SE   df lower.CL upper.CL
#>  Mar        384 25.15 6329      338      436
#>  Apr        347 19.38 6329      311      387
#>  May        298 15.30 6329      270      330
#>  Jun        248 12.51 6329      225      274
#>  Jul        229 11.67 6329      207      253
#>  Aug        204 10.15 6329      185      225
#>  Sep        196  9.55 6329      178      215
#>  Oct        206 10.41 6329      187      227
#>  Nov        256 13.63 6329      230      284
#>  Dec        291 20.75 6329      253      335
#> 
#> Results are averaged over the levels of: Site 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale
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

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

``` r
labl <- 'Values Adjusted to Median Flow and\nMedian 10 Day Precipitation\nAll Dates Combined'

(a <- emmeans(my_ref_grid, ~ Site, type = 'response'))
#> NOTE: Results may be misleading due to involvement in interactions
#>  Site response   SE   df lower.CL upper.CL
#>  S07       218 12.1 6329      196      243
#>  S05       166 11.1 6329      146      190
#>  S17       239 40.2 6329      172      332
#>  S03       333 16.6 6329      302      367
#>  S01       408 21.0 6329      368      451
#> 
#> Results are averaged over the levels of: Month 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale
```

``` r
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

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-22-1.png" style="display: block; margin: auto;" />

## Visualizing Trends

We extract results on the log scale, so we can calculate the linear
predictor by hand, then back transform.

``` r
my_ref_grid <- ref_grid(the_gamm, at = list(Year = 2014, Month = 'Jul'),
                        cov.reduce = median)

a <- summary(emmeans(my_ref_grid, 'Site'))
#> NOTE: Results may be misleading due to involvement in interactions
(b <- summary(emtrends(the_gamm, 'Site', 'Year')))
#>  Site Year.trend     SE   df lower.CL upper.CL
#>  S07      0.1389 0.0197 6329   0.1003    0.178
#>  S05      0.1058 0.0224 6329   0.0618    0.150
#>  S17      0.0543 0.0565 6329  -0.0563    0.165
#>  S03      0.0688 0.0172 6329   0.0351    0.102
#>  S01      0.0931 0.0177 6329   0.0584    0.128
#> 
#> Results are averaged over the levels of: Month 
#> Confidence level used: 0.95
```

The key insight here is that the trends are significant for all sites
EXCEPT S17, where we have fewer years of data.

``` r
plot(b)
```

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />
And those trends are NOT statistically different.

``` r
lookup <- tibble(Site = a[[1]], Intercept = a[[2]], Slope = b[[2]])
#rm(a,b)
```

``` r
df <- tibble(Site = rep(levels(reduced_data$Site), each = 10), 
              Year = rep(2010:2019, 5)) %>%
  mutate(sslope =     lookup$Slope[match(Site, lookup$Site)],
         iintercept = lookup$Intercept[match(Site, lookup$Site)],
         pred = exp((Year - 2014) * sslope + iintercept)) %>%
  select(-sslope, -iintercept)
```

``` r
ggplot(df, aes(x = Year, y = pred, color = Site)) +
         geom_step() +
  ylab('Chloride (mg/l)\n(Flow and Precipitation Adjusted)') +
  xlab('') +
  ylim(0,600) +
  geom_hline(yintercept =  230, color = 'black') +
  #geom_hline(yintercept =  860, color = 'red') +

  theme_cbep(base_size = 12)
```

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

# Model without the interactions.

This model takes several minutes to run (more than 5, less than 15)

``` r
if (! file.exists("models/revised_gamm.rds")) {
  revised_gamm <- gamm(log(Chl_Median) ~ Site + 
                     s(lPrecip) + 
                     s(wlPrecip) +
                     s(FlowIndex) +
                     Month +
                     Year,
                   correlation = corAR1(0.8),
                   na.action = na.omit, 
                   method = 'REML',
                   data = reduced_data)
  saveRDS(revised_gamm, file="models/revised_gamm.rds")
} else {
  revised_gamm <- readRDS("models/revised_gamm.rds")
}
```

## ANOVA

``` r
anova(revised_gamm$gam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> log(Chl_Median) ~ Site + s(lPrecip) + s(wlPrecip) + s(FlowIndex) + 
#>     Month + Year
#> 
#> Parametric Terms:
#>       df      F p-value
#> Site   4  37.85  <2e-16
#> Month  9  19.01  <2e-16
#> Year   1 117.10  <2e-16
#> 
#> Approximate significance of smooth terms:
#>                edf Ref.df      F p-value
#> s(lPrecip)   7.049  7.049  26.90  <2e-16
#> s(wlPrecip)  4.156  4.156  57.19  <2e-16
#> s(FlowIndex) 8.468  8.468 233.01  <2e-16
```

## Summary

``` r
summary(revised_gamm$gam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> log(Chl_Median) ~ Site + s(lPrecip) + s(wlPrecip) + s(FlowIndex) + 
#>     Month + Year
#> 
#> Parametric coefficients:
#>               Estimate Std. Error t value Pr(>|t|)    
#> (Intercept) -1.995e+02  1.898e+01 -10.508  < 2e-16 ***
#> SiteS05     -3.155e-01  7.859e-02  -4.015 6.02e-05 ***
#> SiteS17     -8.932e-02  9.121e-02  -0.979   0.3275    
#> SiteS03      3.717e-01  7.483e-02   4.967 6.97e-07 ***
#> SiteS01      5.784e-01  7.659e-02   7.553 4.86e-14 ***
#> MonthApr    -1.003e-01  4.994e-02  -2.009   0.0446 *  
#> MonthMay    -2.503e-01  5.707e-02  -4.386 1.18e-05 ***
#> MonthJun    -4.322e-01  5.984e-02  -7.223 5.70e-13 ***
#> MonthJul    -5.067e-01  6.313e-02  -8.027 1.18e-15 ***
#> MonthAug    -6.186e-01  6.308e-02  -9.806  < 2e-16 ***
#> MonthSep    -6.601e-01  6.090e-02 -10.840  < 2e-16 ***
#> MonthOct    -6.144e-01  5.742e-02 -10.700  < 2e-16 ***
#> MonthNov    -4.028e-01  5.315e-02  -7.579 4.00e-14 ***
#> MonthDec    -2.713e-01  6.906e-02  -3.929 8.63e-05 ***
#> Year         1.019e-01  9.416e-03  10.821  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Approximate significance of smooth terms:
#>                edf Ref.df      F p-value    
#> s(lPrecip)   7.049  7.049  26.90  <2e-16 ***
#> s(wlPrecip)  4.156  4.156  57.19  <2e-16 ***
#> s(FlowIndex) 8.468  8.468 233.01  <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> R-sq.(adj) =  0.426   
#>   Scale est. = 0.26072   n = 6368
```

## Structure of the GAM

Interestingly, differences between sites and differences in slopes are
marginally not significant in this simplified model.

``` r
plot(revised_gamm$gam)
```

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-30-1.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-30-2.png" style="display: block; margin: auto;" /><img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-30-3.png" style="display: block; margin: auto;" />

## Diagnostic Plots

The help files for `gam.check()` suggest using care when interpreting
results for GAMM models, since the function does not correctly
incorporate the error correlations structure. However, for our purposes,
this is probably sufficient, since our focus is not on statistical
significance, but on estimation.

``` r
gam.check(revised_gamm$gam)
```

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-31-1.png" style="display: block; margin: auto;" />

    #> 
    #> 'gamm' based fit - care required with interpretation.
    #> Checks based on working residuals may be misleading.
    #> Basis dimension (k) checking results. Low p-value (k-index<1) may
    #> indicate that k is too low, especially if edf is close to k'.
    #> 
    #>                k'  edf k-index p-value    
    #> s(lPrecip)   9.00 7.05    0.99    0.32    
    #> s(wlPrecip)  9.00 4.16    0.90  <2e-16 ***
    #> s(FlowIndex) 9.00 8.47    0.87  <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

No appreciable changes in moidel adequacy.

# Model with Separate Years

``` r
if (! file.exists("models/years_gamm.rds")) {
  years_gamm <- gamm(log(Chl_Median) ~ Site + 
                     s(lPrecip) + 
                     s(wlPrecip) +
                     s(FlowIndex) +
                     Month +
                     Year_f,
                   correlation = corAR1(0.8),
                   na.action = na.omit, 
                   method = 'REML',
                   data = reduced_data)
  saveRDS(years_gamm, file="models/years_gamm.rds")
} else {
  years_gamm <- readRDS("models/years_gamm.rds")
}
```

## ANOVA

``` r
anova(years_gamm$gam)
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> Formula:
#> log(Chl_Median) ~ Site + s(lPrecip) + s(wlPrecip) + s(FlowIndex) + 
#>     Month + Year_f
#> 
#> Parametric Terms:
#>        df     F p-value
#> Site    4 39.28  <2e-16
#> Month   9 13.59  <2e-16
#> Year_f  8 30.23  <2e-16
#> 
#> Approximate significance of smooth terms:
#>                edf Ref.df      F p-value
#> s(lPrecip)   7.056  7.056  27.53  <2e-16
#> s(wlPrecip)  4.085  4.085  59.41  <2e-16
#> s(FlowIndex) 8.494  8.494 232.88  <2e-16
```

## Diagnostic Plots

The help files for `gam.check()` suggest using care when interpreting
results for GAMM models, since the function does not correctly
incorporate the error correlations structure. However, for our purposes,
this is probably sufficient, since our focus is not on statistical
significance, but on estimation.

``` r
gam.check(years_gamm$gam)
```

<img src="Chloride_Analysis_Summary_files/figure-gfm/unnamed-chunk-33-1.png" style="display: block; margin: auto;" />

    #> 
    #> 'gamm' based fit - care required with interpretation.
    #> Checks based on working residuals may be misleading.
    #> Basis dimension (k) checking results. Low p-value (k-index<1) may
    #> indicate that k is too low, especially if edf is close to k'.
    #> 
    #>                k'  edf k-index p-value    
    #> s(lPrecip)   9.00 7.06    0.99    0.18    
    #> s(wlPrecip)  9.00 4.08    0.95  <2e-16 ***
    #> s(FlowIndex) 9.00 8.49    0.91  <2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Compare AIC

``` r
AIC(revised_gamm$lme)
#> [1] -253.4018
AIC(the_gamm$lme)
#> [1] -232.1566
AIC(years_gamm$lme)
#> [1] -339.9495
```

So the a model that fits multiple slopes is not justified, while a model
that includes separate terms for each year fits substantially better. A
linear fit for the year here is probably an oversimplification.
