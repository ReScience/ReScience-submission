rm(list=ls())


##  Computing global LE with Owen's version of Steve's code
library(RCurl)
repo <- "https://raw.githubusercontent.com/opetchey/ReScience-submission/beninca_dev"
script <- getURL(paste0(repo,"/code/indirect_method_functions.R"), ssl.verifypeer = FALSE)
eval(parse(text = script))


# follow Elisa's detrending 
# green flag = nano 
# melosira = net 

### Before running this code, you must do at the R command line 
#	rm(list=ls(all=TRUE)); source("c:/lenns/R/lennsBIC.R"); 

require(mgcv); 

######################################################################
#  Load the data into a properly structured data frame 
#######################################################################

## load dataset that allows conversion of names
name_table <- read.csv("../../data/reproduction/repro_steve_name_table.csv")

## choose which dataset to use
dataset_to_use <- "reproduced" # 0.03748704
#dataset_to_use <- "from_steve" # 0.08415112

if(dataset_to_use=="from_steve") {
  X=read.csv("../../data/original/Steve_data_for_gam_test_transformed.csv", row.names=1)
  names(X) <- name_table$repro_names[match(names(X), name_table$steve_posttrans_names)]
  }

if(dataset_to_use =="reproduced") {
  X=read.csv("../../data/reproduction/Repro_data_for_gam_test_transformed.csv", row.names=1)
}

## Fit the GAMs
gams <- Fit_GAMs(X, gval=1.4, T_lag=1)

## calculate and print the GLE
Get_LE_from_fit2(gams, X, epsval=0.01) 


