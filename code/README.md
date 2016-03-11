### Code repository

Don't forget to choose a license. You're free to use one from a set of
well-understood licenses including BSD, GPL or Apache.

See [the Debian Free Software Guidelines](https://www.debian.org/social_contract#guidelines)
for a [list of licenses](https://www.debian.org/legal/licenses/).

**report.rmd**
An rmarkdown file that produced all elements of the ReScience article (including analyses, figures), and lots more. The easiest approach to running this yourself is to open the file `report.rmd` in Rstudio, and click the "Knit HTML" button (though see [here](http://rmarkdown.rstudio.com/authoring_quick_tour.html) for instructions about Rmarkdown). You will, however, need to make sure you have all the required packages installed.

R packages required to run `report.rmd` are:
```R
install.packages(c("tidyr","dplyr", "lubridate", "stringr",
                   "ggplot2", "RCurl", "pracma", "oce",
                   "tseriesChaos", "reshape2", "mgcv", "repmis",
                   "magrittr", "knitr"))
```

**indirect_method_functions.R**
Code for fitting the GAMs required for estimating the global Lyapunov exponent, and code for making the estimate.

**FoodWebGAM=AllSpeciesLongest.R**
Code provided by Stephen Ellner by email, for fitting the GAMs and estimating the global Lyapunov exponent.

**GAM_test folder**
Some scripts for testing the GAM fitting and estimation of global Lyapunov exponent.
