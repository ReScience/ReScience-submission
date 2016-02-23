# Reproduce anaylses of Beninca et al (2008)
Owen Petchey et al  
29 Jan 2015  


# Introduction

This document contains various reproductions of the anaylses presented in the paper *Chaos in a long-term experiment with a plankton community*, by Elisa Benincà and others ([the paper on the Nature website](http://www.nature.com/nature/journal/v451/n7180/abs/nature06512.html)). Details of the methods are in [the Supplement to the Nature paper](http://www.nature.com/nature/journal/v451/n7180/extref/nature06512-s1.pdf).

This reproduction was made as part of the Reproducible Research in Ecology, Evolution, Behaviour, and Environmental Studies (RREEBES) Course, lead by Owen Petchey at the University of Zurich. More information about the course [here](https://github.com/opetchey/RREEBES/blob/master/README.md) on github.

The code and data for the reproduction are [here](https://github.com/opetchey/RREEBES/tree/Beninca_development) on github.


# Authors

Lead author was Owen Petchey. Frank Pennekamp and Marco Plebani made sizeable contributions. All contributors are detailed on github in the commit history etc.


# Known issues

The known issues / problems with this reproduction are with the prefix "Beninca" [here](https://github.com/opetchey/RREEBES/issues). Please add and or solve issues there.


# The data

The data are available as an Excel file supplement to [an Ecology Letters publication](http://onlinelibrary.wiley.com/doi/10.1111/j.1461-0248.2009.01391.x/abstract). The Excel file contains several datasheets. Two are particularly important, as they are the source of the raw data (one contains original species abundances, the one with the nutrient concentrations). Another datasheet in the ELE supplement contains transformed variables. We also got some data direct from Steve Ellner, see below for details.

In the code below, the data and any other files are read from github, which means there must be a connection to github.


# R session info


```r
rm(list=ls())
devtools::session_info()
```

```
## Session info --------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.2.3 (2015-12-10)
##  system   x86_64, darwin14.5.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       America/New_York            
##  date     2016-02-23
```

```
## Packages ------------------------------------------------------------------
```

```
##  package   * version date       source        
##  devtools    1.9.1   2015-09-11 CRAN (R 3.2.0)
##  digest      0.6.8   2014-12-31 CRAN (R 3.2.0)
##  evaluate    0.8     2015-09-18 CRAN (R 3.2.0)
##  formatR     1.2.1   2015-09-18 CRAN (R 3.2.0)
##  htmltools   0.3     2015-12-29 CRAN (R 3.2.3)
##  knitr       1.11    2015-08-14 CRAN (R 3.2.2)
##  magrittr    1.5     2014-11-22 CRAN (R 3.2.0)
##  memoise     0.2.1   2014-04-22 CRAN (R 3.2.0)
##  rmarkdown   0.9.1   2015-12-31 CRAN (R 3.2.3)
##  stringi     1.0-1   2015-10-22 CRAN (R 3.2.0)
##  stringr     1.0.0   2015-04-30 CRAN (R 3.2.0)
##  yaml        2.1.13  2014-06-12 CRAN (R 3.2.0)
```


# Repository to read data from

The datasets are read from this github repository:

```r
repo <- "https://raw.githubusercontent.com/opetchey/ReScience-submission/petchey-plebani-pennekamp-2016"
```

# First get the raw data into R, and clean and tidy it.

All required libraries:

```r
library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)
library(RCurl)
library(pracma)
library(oce)
library(tseriesChaos)
library(reshape2)
library(mgcv)
library(repmis)
library(magrittr)
```

Read in the data and remove some blank rows and columns.

```r
spp.abund <- read.csv(text=getURL(paste0(repo,"/data/original/species_abundances_original.csv")), skip=7, header=T)

spp.abund <- select(spp.abund, -X, -X.1)
spp.abund <- spp.abund[-804:-920,]
str(spp.abund)
```

```
## 'data.frame':	803 obs. of  12 variables:
##  $ Date               : Factor w/ 803 levels "","01/02/91",..: 306 440 465 498 520 601 628 673 699 778 ...
##  $ Day.number         : int  1 6 7 8 9 12 13 15 16 19 ...
##  $ Cyclopoids         : num  0 0 0.0353 0 0.0353 ...
##  $ Calanoid.copepods  : num  1.04 2.03 1.72 2.41 1.71 ...
##  $ Rotifers           : num  7.7 10.19 8.08 6.06 5.94 ...
##  $ Protozoa           : Factor w/ 330 levels "","0","0,000001",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ Nanophytoplankton  : num  0.106 0.212 0.212 0.212 0.212 ...
##  $ Picophytoplankton  : num  1 2 1.52 1.52 1.98 ...
##  $ Filamentous.diatoms: num  0 0 0 0 0 0 0 0 0 0 ...
##  $ Ostracods          : num  0 0 0 0.0187 0 ...
##  $ Harpacticoids      : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ Bacteria           : num  2.15 1.97 1.79 1.61 1.43 ...
```

The Protozoa variable contains some numbers with comman as the decimal separator. This creates a question about what dataset was used for the original analyses, as it could not have been this one.

```r
spp.abund$Protozoa <- as.numeric(str_replace(spp.abund$Protozoa, ",", "."))
```

Format the dates as dates.

```r
spp.abund$Date <- dmy(spp.abund$Date)
# y123 <- year(spp.abund$Date) %% 100
# year(spp.abund$Date) <- 1900 + y123
# rm(y123)
```

Check dates match the Day.number (should give true):

```r
sum(spp.abund$Day.number == 1+as.numeric((spp.abund$Date - spp.abund$Date[1]) / 24 / 60 / 60)) == length(spp.abund$Date)
```

```
## [1] TRUE
```

Check for duplicate dates:

```r
spp.abund$Date[duplicated(spp.abund$Date)]
```

```
## [1] "1996-10-28 UTC"
```

```r
which(duplicated(spp.abund$Date))
```

```
## [1] 702
```

Original dataset contains a duplicated date: 28/10/1996 (row 709 and 710 in excel sheet). Lets change the date in row 709 to 26/10/1996, which will put it half way between the two surrounding dates:

```r
which(spp.abund$Date==ymd("1996-10-28 UTC"))
```

```
## [1] 701 702
```

```r
spp.abund$Date[701] <- ymd("1996-10-26 UTC")
```

Check dates match the Day.number (should give false, as we changed a date):

```r
sum(spp.abund$Day.number == 1+as.numeric((spp.abund$Date - spp.abund$Date[1]) / 24 / 60 / 60)) == length(spp.abund$Date)
```

```
## [1] FALSE
```

Recalculate the Day.number variable from the dates:

```r
spp.abund$Day.number <- 1+as.numeric((spp.abund$Date - spp.abund$Date[1]) / 24 / 60 / 60)
```

Data is in wide format, so change it to long:

```r
spp.abund <- gather(spp.abund, "variable", "value", 3:12)
str(spp.abund)
```

```
## 'data.frame':	8030 obs. of  4 variables:
##  $ Date      : POSIXct, format: "1990-07-12" "1990-07-17" ...
##  $ Day.number: num  1 6 7 8 9 12 13 15 16 19 ...
##  $ variable  : Factor w/ 10 levels "Cyclopoids","Calanoid.copepods",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ value     : num  0 0 0.0353 0 0.0353 ...
```

Bring in the nutrient data:

```r
nuts <- read.csv(text=getURL(paste0(repo,"/data/original/nutrients_original.csv")), skip=7, header=T)

nuts <- select(nuts, -X, -X.1)
nuts <- nuts[-349:-8163,]
nuts$Date <- dmy(nuts$Date)
nuts <- select(nuts, -NO2, -NO3, -NH4)
#nuts$Date[duplicated(nuts$Date)]
#which(duplicated(nuts$Date))
nuts <- gather(nuts, "variable", "value", 3:4)
str(nuts)
```

```
## 'data.frame':	696 obs. of  4 variables:
##  $ Date      : POSIXct, format: "1990-09-18" "1990-09-24" ...
##  $ Day.number: int  69 75 82 90 96 103 110 117 124 131 ...
##  $ variable  : Factor w/ 2 levels "Total.dissolved.inorganic.nitrogen",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ value     : num  28.32 20.84 11.15 15.5 5.92 ...
```

Now put the two datasets together

```r
all.data <- rbind(spp.abund, nuts)
```

Now select only the date range used in the Nature paper. From the supplement *The analysis in Benincà et al. (Nature 2008) covered all data from 16/06/1991 until 20/10/1997*. 

```r
all.data <- filter(all.data, Date>=dmy("15/06/1991") & Date<=dmy("21/10/1997"))
```




# Reproducing figure 1b through 1g

(No attempt to reproduce Figure 1a, as its a food web diagram.)

First add a column that gives the variable types, same as in figure 1b through 1g. Make a table giving species type:

```r
tt <- data.frame(variable=unique(all.data$variable),
                 Type=c("Cyclopoids", "Herbivore", "Herbivore", "Herbivore",
                        "Phytoplankton",  "Phytoplankton", "Phytoplankton",
                        "Detritivore", "Detritivore", "Bacteria", "Nutrient", "Nutrient"))
```

And add the Type variable to the new dataset:

```r
all.data <- merge(all.data, tt)
```

First lets set the colours as in the original:

```r
species.colour.mapping <- c("Cyclopoids"="pink",
                            "Calanoid.copepods"="red",
                            "Rotifers"="blue",
                            "Protozoa"="green",
                            "Nanophytoplankton"="red",
                            "Picophytoplankton"="black",
                            "Filamentous.diatoms"="green",
                            "Ostracods"="lightblue",
                            "Harpacticoids"="purple",
                            "Bacteria"="black",
                            "Total.dissolved.inorganic.nitrogen"="red",
                            "Soluble.reactive.phosphorus"="black")                            
```

Next change the order of the levels in the Type variable, so plots appear in the same order as in the original figure:

```r
all.data$Type <- factor(all.data$Type, levels=c("Cyclopoids", "Herbivore", "Phytoplankton", "Nutrient",
                                    "Bacteria", "Detritivore"))
```

<!--
Here are some versions of the graphs that are not included in the report, hence the html comment:

Now a version that doesn't try to recreate the "gap" in the y axes of the original figures:

```r
g1 <- qplot(as.numeric(Day.number), value, col=variable, data=all.data) +
  facet_wrap(~Type, ncol=2, scales="free_y") +
  geom_point() + geom_line() +
  scale_colour_manual(values = species.colour.mapping)
g1
```

![](report_files/figure-html/unnamed-chunk-19-1.png) 
Looks reasonably good.

Now a version that approximates the "gap", by removing data above it:

```r
an2 <- filter(all.data, Type=="Cyclopoids" & value<0.6 |
                Type=="Herbivore" & value<13 |
                Type=="Phytoplankton" & value<1400 |
                Type=="Nutrient" & value<50 |
                Type=="Bacteria" & value<10 |
                Type=="Detritivore" & value<0.7) 
g1 <- qplot(as.numeric(Day.number), value, col=variable, data=an2) +
  facet_wrap(~Type, ncol=2, scales="free_y") +
  geom_point() + geom_line() +
  scale_colour_manual(values = species.colour.mapping)
g1
```

![](report_files/figure-html/unnamed-chunk-20-1.png) 
Difficult it look like the data go off the top of the graph in ggplot.

Try logarithmic y-axes:

```r
g1 <- qplot(as.numeric(Day.number), log10(value+0.00001), col=variable, data=all.data) +
  facet_wrap(~Type, ncol=2, scales="free_y") +
  geom_point() + geom_line() +
  scale_colour_manual(values = species.colour.mapping)
g1
```

![](report_files/figure-html/unnamed-chunk-21-1.png) 

End of html comment
-->

The graph with abundances fourth root transformed, as this is the transformation used in analyses in the original publication.

```r
g1 <- ggplot(aes(x=as.numeric(Day.number), y=value^0.25, col=variable), data=all.data) +
  xlab("Time (days)") + ylab("Abundance value ^ 0.25") +
  facet_wrap(~Type, ncol=2, scales="free_y") +
  geom_line() +
  scale_colour_manual(values = species.colour.mapping)
g1
```

![](report_files/figure-html/unnamed-chunk-22-1.png) 




# Data transformation

Now we need to work with transformed data. Details of the transformation, copied from the Supplmentary information are in indented quote style in the following sections... looks like this:

> 3. Transformation of the time series. We transformed the original time series, shown in Fig. 1b-g of the main text, to obtain stationary time series with equidistant data and homogeneous units of measurement. The transformation steps are illustrated for the bacteria (Fig. S1).

The ELE supplement contains the raw data and the transformed data, in separate data sheets.


## Removal of long sequences of zeros

> The time series of cyclopoid copepods, protozoa, filamentous diatoms, harpacticoid copepods and ostracods contained long sequences of zero values. This does not imply that these species were absent from the food web during these periods, but that their concentrations were below the detection limit. Time series dominated by many zeros can bias the statistical analysis. Therefore, these time series were shortened to remove long sequences of zero values, before the data transformation.

The exact method for identifying the "long sequences of zeros" was not given. However, the transformed data with zeros removed was provided, and from this we can get the day range of the retained data:

```r
tr <- read.csv(text=getURL(paste0(repo,"/data/original/transformed_data_Nature2008.csv")),
               skip=7, na.string="")
tr <- tr[,-14:-24] ## remove bad columns
tr <- tr[-693:-694,] ## remove last two rows (contain summary stats)
tr <- gather(tr, key="Species", value="Abundance", 2:13) ## make long format
levels(tr$Species)[levels(tr$Species)=="Calanoids"] <- "Calanoid.copepods" # rename a level (and next line)
levels(tr$Species)[levels(tr$Species)=="Total.Dissolved.Inorganic.Nitrogen"] <- "Total.dissolved.inorganic.nitrogen"
names(tr)[2] <- "variable" # rename a couple of variables...
names(tr)[1] <- "Day.number"
tr$Day.number <- as.numeric(as.character(tr$Day.number)) # make Day.number numeric
## get day range, in order to use this to constrain later data
min.max.day <- group_by(tr, variable) %>%
  summarise(min.day=min(Day.number[!is.na(Abundance)]), max.day=max(Day.number[!is.na(Abundance)]))
```

Now we create a trimmed version of the raw data:

```r
trimmed_raw <- full_join(all.data, min.max.day)
```

```
## Joining by: "variable"
```

```r
trimmed_raw <- filter(trimmed_raw, Day.number>min.day & Day.number<max.day)
```



## Interpolation (zeros trimmed)

> First, the time series were interpolated using cubic hermite interpolation, to obtain data with equidistant time intervals of 3.35 days (Fig. S1a).

Make a sequence of times at which to interpolate. 

```r
#aggregate(Day.number ~ variable, all.data, min)
#aggregate(Day.number ~ variable, all.data, max)
#xout <- seq(343.35, 2657.2, by=3.35)
xout <- seq(343.35, 2658, by=3.35)
#range(xout)
```

Perform the interpolation:

```r
# make a function to do the interpolation
Interpolate <- function(x, y) {
  xout.inner <- xout[xout>=min(x) & xout<=max(x)]
  res <- pracma::interp1(x=as.numeric(x),
                  y=as.numeric(y),
                  xi=xout.inner,
                  method="cubic")
  res.df <- data.frame(Day.number=xout.inner,
                       value=res)
}

# run the function and tidy the data
# (somewhat ashamed of this code...)
# remove NAs before interpolation, so they don't cause problems
trimmed_raw <- na.omit(trimmed_raw)
res <- group_by(trimmed_raw, variable) %>%
  do(le=Interpolate(.$Day.number, .$value))

trimmed_interp <- res %>% do(data.frame(Day.number = .$le[[1]], value = .$le[[2]], variable = .[[1]]))


# for(p in 1:nrow(res)){
#     temp <- data.frame(variable=res[p,1],
#                        res[p,2][[1]])
#     if(p==1)
#       trimmed_interp <- temp
#     if(p>1)
#       trimmed_interp <- rbind(trimmed_interp, temp)
#   }
```


## Fourth root transform  (zeros trimmed)

> Next, because the original time series showed many sharp spikes, the time series were rescaled using a fourth-root power transformation (Fig. S1b). The sharp spikes bias "direct method" estimates of the Lyapunov exponent, because nearby pairs of reconstructed state vectors mostly occurred in the troughs between spikes. The average rate of subsequent trajectory divergence from these pairs is therefore an estimate of the local Lyapunov exponent in the troughs, which may be very different from the global Lyapunov exponent. By making spikes and troughs more nearly symmetric, the power transformation resulted in a much more even spread of nearby state vector pairs across the full range of the data for all functional groups in the food web. The transformation is also useful for fitting nonlinear models of the deterministic skeleton (used for nonlinear predictability and indirect method estimates of the Lyapunov exponent), which was done by least squares and therefore is most efficient when error variances are stabilized. Fourth-root transformation is intermediate between the square-root transformation that would approximately stabilize the measurement error variance in count data from random subsamples, and the log transformation that is usually recommended for stabilizing process noise variance due to stochastic variation in birth and death rates.


```r
trimmed_interp$fr.value <- trimmed_interp$value^0.25
```


## Detrend  (zeros trimmed)

> The time series were then detrended using a Gaussian kernel with a bandwidth of 300 days (red line in Fig. S1b), to obtain stationary time series. Most species did not show long-term trends, except for the bacteria, detritivores (ostracods and harpacticoid copepods), dissolved inorganic nitrogen and soluble reactive phosphorus. One possible explanation for these trends in the microbial loop could be the slow accumulation of refractory organic material in the mesocosm, but we have not measured this component.



```r
## data to detrend
ww.td <- filter(trimmed_interp, variable=="Total.dissolved.inorganic.nitrogen" |
                  variable=="Soluble.reactive.phosphorus" |
                  variable=="Bacteria" |
                  variable=="Ostracods" |
                  variable=="Harpacticoids")
## and to not detrend
ww.ntd <- filter(trimmed_interp, variable!="Total.dissolved.inorganic.nitrogen" &
                  variable!="Soluble.reactive.phosphorus" &
                  variable!="Bacteria" &
                  variable!="Ostracods" &
                  variable!="Harpacticoids")
## detrend:
ww1 <- group_by(ww.td, variable) %>%
  mutate(trend=ksmooth(Day.number,fr.value,bandwidth=300,kernel="normal")$y)
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```r
ww1$dt.value <- ww1$fr.value-ww1$trend

## don't detrend
ww2 <- ww.ntd
ww2$trend <- 0
ww2$dt.value <- ww2$fr.value

## rejoin
detrended <- rbind(ww1, ww2)
```


## Rescale (zeros trimmed)

> Finally, the time series were linearly rescaled to have zero mean and a standard deviation of 1 (Fig. S1c).



```r
## standardise
scaled_trimmed_final <- group_by(detrended, variable) %>%
  mutate(stand.y=as.numeric(scale(dt.value)))
#summarise(scaled_trimmed_final, mean=mean(stand.y), sd=sd(stand.y))
```

## Compare the zeros removed transformed


```
## Joining by: c("Day.number", "variable")
```

![](report_files/figure-html/unnamed-chunk-30-1.png) 

Looks pretty good, unclear why N and P are less well reproduced than others.


```r
final <- scaled_trimmed_final
final$y <- final$stand.y
```


# Spectral analyses


```r
# Calculate and plot raw spectra
spectra <- final %>% group_by(variable) %>%
  do(spectra = spectrum(ts(data=.$y, end=2650.15, deltat=3.35), log='no', method="pgram", detrend=F, plot=F))
spec <- spectra %>%
  do(data.frame(spec = .$spec[[2]], freq = .$spec[[1]], group = .[[1]]))
ggplot(spec, aes(y=spec, x=1/freq, group=group)) +
  geom_line() +
  facet_wrap(~group, scales="free_y") +
  scale_x_log10() + 
  coord_cartesian(ylim=c(0,50), xlim=c(10,240))
```

![](report_files/figure-html/unnamed-chunk-32-1.png) 

```r
# freq.est <- spec %>% group_by(group) %>% mutate(max_spec = max(spec), freq = freq)
# freq.est <- subset(freq.est, max_spec==spec, select=c(freq,group))
# freq.est$freq <- 1/freq.est$freq

# Calculate and plot Welch's periodogram
wspectra <- final %>% group_by(variable) %>%
  do(spectra = pwelch(ts(data=.$y, end=2650.15, deltat=3.35), window=5, method="pgram", plot=F))
wspec <- wspectra %>% do(data.frame(spec = .$spec[[2]], freq = .$spec[[1]], group = .[[1]]))
ggplot(wspec, aes(y=spec, x=1/freq, group=group)) +
  geom_line() +
  facet_wrap(~group, scales="free_y") +
  scale_x_log10() +
  coord_cartesian(ylim=c(0.01,100), xlim=c(10,240)) +
  scale_y_log10()
```

![](report_files/figure-html/unnamed-chunk-32-2.png) 

```r
# freq.est <- wspec %>%
#   group_by(group) %>%
#   mutate(max_spec = max(spec), freq = freq)
# freq.est <- subset(freq.est, max_spec==spec, select=c(freq,group))
# freq.est$freq <- 1/freq.est$freq
# #frequency(final$y)
# ts <- as.ts(final$y, frequency = 0.3)
# #time(ts)
```

# Reproducing Table 1 using ELE supplement data.

Create dataset with zeros removed for this table (note that this is probably not how Beninca et al did this):

```r
# some data formatting
final_nozeros <- final
final_wide <- spread(as.data.frame(final_nozeros[,c("Day.number","variable","stand.y")]), variable, stand.y)
final_wide <- select(final_wide, - Day.number)
```

Calculate correlation coefficients:

```r
cor.coefs <- cor(final_wide, use="pairwise.complete.obs")
```

Only keep the upper triangle of the cor.pvals matrix:

```r
#for(i in 1:10){
#  for(j in 1:10){
#  cor.coefs[i,j] <- ifelse(i<j, cor.coefs[i,j], NA)
#}}
```

Get p-vals:

```r
cor.pvals <- matrix(NA, length(final_wide[1,]), length(final_wide[1,]))
for(i in 1:length(final_wide[1,]))
  for(j in 1:length(final_wide[1,])) {
    if(sum(!is.na(final_wide[,i]) & !is.na(final_wide[,j]))>10)
    cor.pvals[i,j] <- cor.test(final_wide[,i], final_wide[,j])$p.value
}
    #plot(final_wide[,2], final_wide[,3])
```

Only keep the upper triangle of the cor.pvals matrix:

```r
#for(i in 1:10){
#  for(j in 1:10){
#  cor.pvals[i,j] <- ifelse(i<j, cor.pvals[i,j], NA)
#}}
```


Add significance "stars" to cor.coefs from cor.pvals:

```r
cor.stars <- cor.pvals
cor.stars <- ifelse(cor.pvals<0.0001, "***",
                    ifelse(cor.pvals<0.001, "**",
                           ifelse(cor.pvals<0.05, "*", "")))
```

Make the table entries, by pasting the values and stars:

```r
cor.cp <- cor.coefs
for(i in 1:length(final_wide[1,])){
  for(j in 1:length(final_wide[1,])){
  cor.cp[i,j] <- paste(round(cor.coefs[i,j],2), cor.stars[i,j])
}}
```

Cosmetic work on the table:

```r
oo <- c(10, 9, 8, 11, 12, 6, 5, 3, 4, 2)
cor.coefs <- cor.coefs[oo,oo]
cor.cp <- cor.cp[oo,oo]
for(i in 1:length(cor.cp[1,]))
  for(j in 1:i)
  cor.cp[i,j] <- ""
#sn <- c("Bacteria", "Harps", "Ostr", "N", "P", "Picophyt", "Nanophyt", "Rotifers", "Protozoa", "Calanoids")
sn <- substr(rownames(cor.cp), 1, 4)
sn[4:5] <- c("N", "P")
dimnames(cor.cp) <- list(sn, sn)
```

Make it a table:

```r
library(knitr)
table1b <- kable(cor.cp, format="markdown", col.names = colnames(cor.cp), align="c",
                caption="Table 1.'Correlations between the species in the food web. Table entries show the product–moment correlation coefficients, after transformation of the data to stationary time series (see Methods). Significance tests were corrected for multiple hypothesis testing by calculation of adjusted P values using the false discovery rate.' Significant correlations are indicated as follows: *: P<0.05; **: P<0.01; ***: P<0.001. 'The correlation between calanoid copepods and protozoa could not be calculated, because their time series did not overlap. Filamentous diatoms and cyclopoid copepods were not included in the correlation analysis, because their time series contained too many zeros.' (Beninca et al. 2008)")
table1b
```



|     | Bact | Harp |   Ostr    |    N     |    P     | Pico  |   Nano    |   Roti    |  Prot   |   Cala   |
|:----|:----:|:----:|:---------:|:--------:|:--------:|:-----:|:---------:|:---------:|:-------:|:--------:|
|Bact |      | 0.03 | -0.24 *** |   0.05   | 0.18 *** | 0.03  | -0.17 *** |  0.3 ***  | -0.17 * | 0.22 *** |
|Harp |      |      |  0.18 **  | -0.13 *  |   0.09   | -0.04 |   -0.03   |   0.01    |  0.14   |  -0.02   |
|Ostr |      |      |           | -0.16 ** |  -0.05   | -0.04 |   0.01    |   -0.03   | 0.18 *  |  -0.03   |
|N    |      |      |           |          |  0.08 *  |   0   |   -0.02   |   -0.04   |  0.03   |   0.04   |
|P    |      |      |           |          |          | -0.04 |   0.03    |   0.1 *   |  -0.02  |  -0.08   |
|Pico |      |      |           |          |          |       | -0.17 *** |   -0.03   | -0.22 * | 0.28 *** |
|Nano |      |      |           |          |          |       |           | -0.19 *** |  0.13   | -0.14 *  |
|Roti |      |      |           |          |          |       |           |           |  -0.02  |  -0.1 *  |
|Prot |      |      |           |          |          |       |           |           |         |  NA NA   |
|Cala |      |      |           |          |          |       |           |           |         |          |

Compare correlations in table 1 of the original article with those calculated here:

```r
# read in the correlations given in the original publication
original.cors <- read.csv(text=getURL(paste0(repo,"/data/original/table1_original_article.csv")), skip=0, header=T, row.names = 1)
original.cors <- as.matrix(original.cors)

qplot(x=as.vector(original.cors), y=as.vector(cor.coefs), ylim = c(-0.4,0.4), xlim=c(-0.4, 0.4),
      xlab="Correlations in original article", ylab="Correlations calculated\nin this reproduction") + geom_abline(intercept=0,slope=1, colour="red")
```

```
## Warning: Removed 56 rows containing missing values (geom_point).
```

![](report_files/figure-html/unnamed-chunk-42-1.png) 



# Predictability (Figure 2)

This will be done after getting global Lyapunov exponents by the indirect method.

# Lyapunov exponents by direct method (Figure 3)

Estimate the Lyapunov exponents of the time series, via time-delayed embedding. The Nature report used the [Tisean software](http://www.mpipks-dresden.mpg.de/~tisean/), which was available from CRAN [until mid 2014](http://cran.r-project.org/web/packages/RTisean/index.html). Based on this, and being a bit less well integrated with R, we'll instead use the [tseriesChaos](http://cran.r-project.org/web/packages/tseriesChaos/index.html) package, which was *largely inspired by the TISEAN project*. 

Time delay (1), embedding dimension (6), and Theiler window (50) were used in the Nature report. Other parameters are chosen rather randomly!


```r
time.delay <- 1
embedding.dimension <- 6
Theiler.window <- 50
```

Note that a time step is 3.35 days in the transformed data. So to get a graph with 80 days on the x-axis (as in Figure 3 in the Nature report), we need 80/3.35 = 24 time steps for the calculation of Lyapunov exponents.

```r
time.steps <- 24
```

Remove the species that were not analysed in the Nature report, due to too many zeros in the time series:

```r
led <- filter(tr, variable!="Filamentous.diatoms",
                variable!="Protozoa",
                variable!="Cyclopoids")
```

Do the time-delayed embedding, and save the data for the graphs:

```r
all.species <- unique(as.character(led$variable))
diverg <- matrix(NA, time.steps, length(all.species))
colnames(diverg) <- all.species
for(i in 1:length(all.species)) {
  print(all.species[i])
  tr.fs <- filter(final, variable==all.species[i])$y
  diverg[,i] <- as.numeric(try(lyap_k(tr.fs,
                                      m=embedding.dimension,
                                      d=time.delay,
                                      k=10, # number of considered neighbours 20
                                      ref=40, # number of points to take into account 100
                                      t=Theiler.window,
                                      s=time.steps,
                                      eps=10 # radius where to find nearest neighbours 10
                                      )))
}
```

```
## [1] "Calanoid.copepods"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Rotifers"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Nanophytoplankton"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Picophytoplankton"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Ostracods"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Harpacticoids"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Bacteria"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Total.dissolved.inorganic.nitrogen"
## Finding nearests
## Keeping  40  reference points
## Following points
## [1] "Soluble.reactive.phosphorus"
## Finding nearests
## Keeping  40  reference points
## Following points
```

```r
## a bit of a fudge with the translation to days
diverg <- as.data.frame(cbind(days=1:time.steps, diverg))
diverg <- gather(diverg, Species, Difference, 2:10)
diverg$days <- diverg$days*3.35
#str(diverg)
```

Next calculate the Lyapunov exponents, noting that 6 or 7 points were used in the regressions in the Nature report

```r
#diverg$Difference[is.na(diverg$Difference)] <- 0
diverg$Difference[is.infinite(diverg$Difference)] <- NA
diverg.short <- filter(diverg, days<24) ## 24 is about 6 steps, after initial gap
LEs <- group_by(diverg.short, Species) %>%
  summarise(le=coef(lm(Difference[1:6] ~ days[1:6]))[2])
#pval=summary(lm(Difference[1:6] ~ days[1:6]))$coefficients[2,4])
```


Then plot the graphs with LE:

```r
LEs <- mutate(LEs, days=20, Difference=-0.5)
diverg$Species_short <- plyr::revalue(diverg$Species,
                                 c(Soluble.reactive.phosphorus="P",
                                   Total.dissolved.inorganic.nitrogen="N"))
LEs$Species_short <- plyr::revalue(LEs$Species,
                                 c(Soluble.reactive.phosphorus="P",
                                   Total.dissolved.inorganic.nitrogen="N"))
g1 <- ggplot(diverg, aes(x=days, y=Difference)) + geom_point() + facet_wrap(~Species_short) +
  geom_text(data=LEs, aes(x=60, y=0, label=round(le,3)))
g1
```

```
## Warning: Removed 1 rows containing missing values (geom_point).
```

![](report_files/figure-html/unnamed-chunk-48-1.png) 

Not exactly the same at Figure 3 in the Nature report. Qualitatively the same, except for where the time-delayed embedding failed.

Compare graphically with original LEs

```r
## read in LE values given in original publication
original.LEs <- read.csv(text=getURL(paste0(repo,"/data/original/original_direct_LEs.csv")), skip=0, header=T)
both_LEs <- full_join(original.LEs, LEs) # merge the reproduced and original values
```

```
## Joining by: "Species"
```

```
## Warning in outer_join_impl(x, y, by$x, by$y): joining factors with
## different levels, coercing to character vector
```

```r
ggplot(both_LEs, aes(x=LE, y=le)) +
  geom_point() +
  xlim(0.02, 0.09) +
  ylim(0.02, 0.09) +
  geom_abline(intercept=0, slope=1, colour="red") +
  xlab("Lyapunov exponent from original article") + 
  ylab("Reproduced Lyapunov exponent")
```

![](report_files/figure-html/unnamed-chunk-49-1.png) 

Get the mean and standard deviation of the original and reproduced Lyapunov exponents

```r
mean(original.LEs$LE)
```

```
## [1] 0.05755556
```

```r
sd(original.LEs$LE)
```

```
## [1] 0.004798727
```

```r
mean(LEs$le)
```

```
## [1] 0.05971112
```

```r
sd(LEs$le)
```

```
## [1] 0.01411634
```



# Lyapunov exponents by indirect method

The process is to transform the data, to fit GAMs, and then to calculated the global Lyapunov exponent.

It seems that this must have been done on data without zeros removed, because removal of the zeros would leave too many missing values in the GAMs. So we do the transformation steps again, but without removal of zeros.

## Interpolation (zeros retained)

> First, the time series were interpolated using cubic hermite interpolation, to obtain data with equidistant time intervals of 3.35 days (Fig. S1a).

Make a sequence of times at which to interpolate. 

```r
#aggregate(Day.number ~ variable, all.data, min)
#aggregate(Day.number ~ variable, all.data, max)
#xout <- seq(343.35, 2657.2, by=3.35)
xout <- seq(343.35, 2658, by=3.35)
#range(xout)
```

Perform the interpolation:

```r
all.data <- na.omit(all.data)
mt <- plyr::dlply(all.data,
                  "variable",
                  function(xx) pracma::interp1(x=xx$Day.number,
                                               y=xx$value,
                                               xi=xout,
                                               method="cubic"))
## Aside: the duplicated date that was previously fixed was only discovered by a warning message
## given by the pracma::interp1 function!
mt <- as.data.frame(mt)
mt <- cbind(Day.number=xout, mt)
mt <- gather(mt, variable, value, 2:13)
```

Check this against the data direct from Steve:

```r
from.steve <- read.csv(text=getURL(paste0(repo,"/data/original/interp_short_allsystem_original.csv")), header=T)
name_table <- read.csv(text=getURL(paste0(repo,"/data/reproduction/repro_steve_name_table.csv")), header=T)
names(from.steve) <- name_table$repro_names[match(names(from.steve), name_table$steve_names)]


from.steve <- gather(from.steve, Species, Abundance, 2:13)
names(from.steve) <- c("Day.number", "variable", "value")
from.steve$value_steve <- from.steve$value 
from.steve <- select(from.steve, -value)

## correcting an offset in the two datasets Day numbers
from.steve$Day.number <- from.steve$Day.number + 1.65 + 3.35

## join the data, again first making days character, then back to numeric after
mt$Day.number <- as.character(mt$Day.number)
from.steve$Day.number <- as.character(from.steve$Day.number)
ff <- inner_join(mt, from.steve)
```

```
## Joining by: c("Day.number", "variable")
```

```
## Warning in inner_join_impl(x, y, by$x, by$y): joining factors with
## different levels, coercing to character vector
```

```r
mt$Day.number <- as.numeric(mt$Day.number)
from.steve$Day.number <- as.numeric(from.steve$Day.number)
ff$Day.number <- as.numeric(ff$Day.number)

## and plot time series
g1 <- ggplot(mt, aes(x=as.numeric(Day.number), y=sqrt(value))) +
  facet_wrap(~variable, ncol=2, scales="free_y") +
  geom_line(size=0.5, col="black") +
  scale_colour_manual(values = species.colour.mapping)
g2 <- geom_line(data=from.steve, aes(x=Day.number, y=sqrt(value_steve)), colour="red")
g1 + g2
```

![](report_files/figure-html/unnamed-chunk-53-1.png) 

```r
## same but to allow x-y scatter
mt$Day.number <- as.character(mt$Day.number)
from.steve$Day.number <- as.character(from.steve$Day.number)
ff <- inner_join(mt, from.steve)
```

```
## Joining by: c("Day.number", "variable")
```

```
## Warning in inner_join_impl(x, y, by$x, by$y): joining factors with
## different levels, coercing to character vector
```

```r
mt$Day.number <- as.numeric(mt$Day.number)
from.steve$Day.number <- as.numeric(from.steve$Day.number)
ff$Day.number <- as.numeric(ff$Day.number)
#str(ff)
#ff$difference <- ff$value - ff$value_steve
#qplot(log10(ff$difference))
qplot(data=ff, x=log10(value_steve), y=log10(value)) +
  xlab("Values from Ellner\n(log 10 transformed)") +
  ylab("Reproduced values\n(log 10 transformed)") +
  geom_abline(intercept=0, slope=1, colour="red")
```

![](report_files/figure-html/unnamed-chunk-53-2.png) 

Looks good.







## Fourth root transform (zeros retained)

> Next, because the original time series showed many sharp spikes, the time series were rescaled using a fourth-root power transformation (Fig. S1b). The sharp spikes bias "direct method" estimates of the Lyapunov exponent, because nearby pairs of reconstructed state vectors mostly occurred in the troughs between spikes. The average rate of subsequent trajectory divergence from these pairs is therefore an estimate of the local Lyapunov exponent in the troughs, which may be very different from the global Lyapunov exponent. By making spikes and troughs more nearly symmetric, the power transformation resulted in a much more even spread of nearby state vector pairs across the full range of the data for all functional groups in the food web. The transformation is also useful for fitting nonlinear models of the deterministic skeleton (used for nonlinear predictability and indirect method estimates of the Lyapunov exponent), which was done by least squares and therefore is most efficient when error variances are stabilized. Fourth-root transformation is intermediate between the square-root transformation that would approximately stabilize the measurement error variance in count data from random subsamples, and the log transformation that is usually recommended for stabilizing process noise variance due to stochastic variation in birth and death rates.


```r
mt$fr.value <- mt$value^0.25
```


## Detrend (zeros retained)

> The time series were then detrended using a Gaussian kernel with a bandwidth of 300 days (red line in Fig. S1b), to obtain stationary time series. Most species did not show long-term trends, except for the bacteria, detritivores (ostracods and harpacticoid copepods), dissolved inorganic nitrogen and soluble reactive phosphorus. One possible explanation for these trends in the microbial loop could be the slow accumulation of refractory organic material in the mesocosm, but we have not measured this component.


```r
ww.td <- filter(mt, variable=="Total.dissolved.inorganic.nitrogen" |
                  variable=="Soluble.reactive.phosphorus" |
                  variable=="Bacteria" |
                  variable=="Ostracods" |
                  variable=="Harpacticoids")
## and to not detrend
ww.ntd <- filter(mt, variable!="Total.dissolved.inorganic.nitrogen" &
                  variable!="Soluble.reactive.phosphorus" &
                  variable!="Bacteria" &
                  variable!="Ostracods" &
                  variable!="Harpacticoids")
## detrend:
ww1 <- group_by(ww.td, variable) %>%
  mutate(trend=ksmooth(Day.number,fr.value,bandwidth=300,kernel="normal")$y)
ww1$dt.value <- ww1$fr.value-ww1$trend

## don't detrend
ww2 <- ww.ntd
ww2$trend <- 0
ww2$dt.value <- ww2$fr.value

## rejoin
detrended <- rbind(ww1, ww2)
```





## Rescale (zeros retained)

> Finally, the time series were linearly rescaled to have zero mean and a standard deviation of 1 (Fig. S1c).

(Note that this standardisation is not done in the code sent by Stephen. Probably shouldn't make a difference in the GAMs?)


```r
## standardise
#scaled_trimmed_final <- group_by(trimmed_detrended, variable) %>%
#  mutate(stand.y=as.numeric(scale(dt.value)))
#summarise(scaled_trimmed_final, mean=mean(stand.y), sd=sd(stand.y))

## or don't standardise
final <- detrended
final$y <- final$dt.value
#summarise(final, mean=mean(y), sd=sd(y))
```


```r
glimpse(final)
```

```
## Observations: 8,292
## Variables: 7
## $ Day.number (dbl) 343.35, 346.70, 350.05, 353.40, 356.75, 360.10, 363...
## $ variable   (fctr) Ostracods, Ostracods, Ostracods, Ostracods, Ostrac...
## $ value      (dbl) 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000...
## $ fr.value   (dbl) 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.00000...
## $ trend      (dbl) 0.03615422, 0.03640559, 0.03666087, 0.03690811, 0.0...
## $ dt.value   (dbl) -0.03615422, -0.03640559, -0.03666087, -0.03690811,...
## $ y          (dbl) -0.03615422, -0.03640559, -0.03666087, -0.03690811,...
```




## Figure S1 (visualising the transformation)

Choose a species to plot:

```r
soi <- "Bacteria"
```

Raw and interpolated data:


```r
g1 <- ggplot(filter(all.data, variable==soi), aes(x=Day.number, y=value)) +
  facet_wrap(~variable, ncol=2, scales="free_y") +
  geom_point(size=1, col="black") + geom_line(size=0.1) +
  scale_colour_manual(values = species.colour.mapping) + ggtitle("Raw and interpolated data")
g2 <- geom_line(data=filter(final, variable==soi), aes(x=Day.number, y=value), size=0.25, col="blue")
g1 + g2
```

![](report_files/figure-html/unnamed-chunk-60-1.png) 

Fourth root transformed with trend:


```r
g1 <- ggplot(filter(final, variable==soi), aes(x=Day.number, y=fr.value)) +
  facet_wrap(~variable, ncol=2, scales="free_y") +
  geom_point(size=0.5, col="black") + geom_line(size=0.1) +
  scale_colour_manual(values = species.colour.mapping) + ggtitle("Quarter power trans. and trend")
g2 <- geom_line(data=filter(final, variable==soi), aes(x=Day.number, y=trend), size=0.25, col="blue")
g1 + g2
```

![](report_files/figure-html/unnamed-chunk-61-1.png) 

Detrended and normalised:


```r
g1 <- ggplot(filter(final, variable==soi), aes(x=Day.number, y=y)) +
  facet_wrap(~variable, ncol=2, scales="free_y") +
  geom_point(size=0.5, col="black") + geom_line(size=0.1) +
  scale_colour_manual(values = species.colour.mapping) + ggtitle("Detrended and normalised")
g1
```

![](report_files/figure-html/unnamed-chunk-62-1.png) 


The functions used in the code below are based on code received from Stephen Ellner. The modifications have been tested, and produce the same results as Ellner's original code.

The function needs a matrix, X, with species abundances in wide format. Be careful to work on the unstandardised data (or standardised, if you wish) (comment out appropriate lines here).


```r
## use next line to work on unstandardised data
final.to.melt <- final[, c("variable", "dt.value", "Day.number")]
## use next line to work on standardised
#final.to.melt <- final[, c("variable", "y", "Day.number")]
names(final.to.melt)[1] <- "Species"
melted <- melt(final.to.melt, id=c("Species", "Day.number"))
X <- acast(melted, formula= Day.number ~ Species)
#str(X)
X <- as.data.frame(X)
```


Restrict range of data appropriately:

```r
## Select the time period to use 
start.longest=334; start.longer=808; start.shorter=1035; 
e=as.numeric(row.names(X)) > start.longer; X=X[e,];
e=as.numeric(row.names(X)) < 2654; X=X[e,]; 
```




Load and run the functions, or read in from data file. Chunk not run by default, as take a minute or so. At last run, gave GLE of 0.03730891


```r
# read script lines from website
script <- getURL(paste0(repo,"/code/indirect_method_functions.R"), ssl.verifypeer = FALSE)
eval(parse(text = script))

## Fit the GAMs
gams <- Fit_GAMs(X, gval=1.4, T_lag=1)

## calculate and print the GLE
Get_LE_from_fit2(gams, X1, epsval=0.01) 
```

This is quite far from the number using code and data from Steve (0.08415112). This is caused by small differences in the data (and not by the functions), as we've found in using the code in the folder "GAM_test".

# Predictability (Figure 2)

The procedure is to fit GAMs (neural nets in the original publication) with response variable a function of explanatory variables at lag T, then to record r-squared of observed and fit data. Chunk not run by default, due to time required.







And plot our version of figure 2:

```r
## Load the r-squared values that were previously saved (above)
## load the already saved data from github (this can take some time depending on the internet connection)
source_data(paste0(repo, "/data/reproduction/rsq_vals.Rdata?raw=True"))
```

```
## Downloading data from: https://raw.githubusercontent.com/opetchey/ReScience-submission/petchey-plebani-pennekamp-2016/data/reproduction/rsq_vals.Rdata?raw=True 
## 
## SHA-1 hash of the downloaded data file is:
## 8805fa71c9d401416058bf2b5ad490520d80d2a5
```

```
## [1] "rsq_vals"
```

```r
cors <- data.frame(dist=1:12, rsq_vals)
cors.long <- gather(cors, key=dist)
names(cors.long) <- c("Prediction_distance", "Variable", "Correlation")


cors.long$Variable <- factor(cors.long$Variable, c("Cyclopoids", "Rotifers", "Calanoid.copepods", "Picophytoplankton",
                             "Nanophytoplankton", "Filamentous.diatoms", "Soluble.reactive.phosphorus",
                             "Total.dissolved.inorganic.nitrogen", "Bacteria", "Ostracods", "Harpacticoids",
                             "Protozoa"))

## read in the data of the original article that came direct from Elisa
original_preds <- read.csv(text=getURL(paste0(repo,"/data/original/original_rsquared.csv")), skip=0, header=T)
opw <- gather(original_preds, Species, r_squared, 2:25)
opw <- separate(opw, Species, c("Species", "Model"))

name.mapping <- data.frame(Species=c("cyclo", "roti", "cala", "pico", "nano",
                                   "dia", "P", "N", "bact", "ostra", "harpa",
                                   "proto"),
                           Variable=c("Cyclopoids", "Rotifers", "Calanoid.copepods",
                                  "Picophytoplankton", "Nanophytoplankton",
                                  "Filamentous.diatoms", "Soluble.reactive.phosphorus",
                                  "Total.dissolved.inorganic.nitrogen", "Bacteria",
                                  "Ostracods", "Harpacticoids", "Protozoa"))
opw <- full_join(opw, name.mapping)
```

```
## Joining by: "Species"
```

```
## Warning in outer_join_impl(x, y, by$x, by$y): joining factor and character
## vector, coercing into character vector
```

```r
ggplot(filter(cors.long, Prediction_distance<13),
       aes(x=Prediction_distance*3.35, y=Correlation)) +
  xlab("Prediction time (days)") +
  geom_point(col="red") +
  facet_wrap(~Variable, ncol=3 ) +
  geom_point(data=filter(opw, Model=="nonlinear"),
             aes(x=prediction_time, y=r_squared)) +
  ylim(0,1)
```

![](report_files/figure-html/unnamed-chunk-68-1.png) 

Some species very similar to figure 2 of the nature paper, some more different. Difference must be due to data in figure 2 coming from neural net models, whereas here it comes from the GAMs. Qualitative result the same as in the original paper, however.
