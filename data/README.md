### Data repository

**species_abundances_original.csv**, **nutrients_original.csv** and **transformed_data_Nature2008.csv**

The data for the Benincà et al (2008) Nature article were subsequently released in a supplement to an [an Ecology Letters article](http://onlinelibrary.wiley.com/doi/10.1111/j.1461-0248.2009.01391.x/abstract). The released data were in the form of an Excel spreadsheet with multiple worksheets. We saved three of these worksheets as csv files with file name the same as the worksheet name. No changes were made in Excel prior to saving as csv files. 

**interp_short_allsystem_newnames.csv**
This dataset was received direct from Stephen Ellner by email, and contains interpolated variable values. Short refers to the period of the dataset being the short version used in the Nature paper analyses (rather than the longer version that was used otherwise).

**table1_original_article.csv**
We created this table of correlations among species abundances by manually typing in values from the original article.

**original_direct_LEs.csv
We created this list of Lyapunov exponents by manually typing in values from the original article.

**GLE_estimate.Rdata**
This is a Rdata file containing a list of three things:

1. The estimate of the global Lyapunov exponent.
2. A list of 12 GAMs, one for each variable.
3. A dataframe of the variables used in the GAMs.



