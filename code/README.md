### Code repository

Executing the main script `simulation.R` script will produce all the data, but beware, simulations will take some time - on a 16 core computer it takes about 8h. For this reason, the data is already provided in `../data` directory.


### Details on dependencies

The code has been tested on Linux machine:  
~~~
Distributor ID: Ubuntu
Description:    Ubuntu 17.04
Release:    17.04
Codename:   zesty
Kernel: Linux 4.12.7
~~~

R `sessionInfo()` output, with all necesary packages loaded:  
~~~
R version 3.4.0 (2017-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 17.04

Matrix products: default
BLAS: /usr/lib/openblas-base/libblas.so.3
LAPACK: /usr/lib/libopenblasp-r0.2.19.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] doRNG_1.6.6       rngtools_1.2.4    pkgmaker_0.22     registry_0.3     
 [5] doParallel_1.0.10 iterators_1.0.8   foreach_1.4.3     ggrepel_0.6.5    
 [9] reshape2_1.4.2    dplyr_0.5.0       ggplot2_2.2.1    

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.10     magrittr_1.5     munsell_0.4.3    colorspace_1.3-2
 [5] xtable_1.8-2     R6_2.2.0         stringr_1.2.0    plyr_1.8.4      
 [9] tools_3.4.0      grid_3.4.0       gtable_0.2.0     DBI_0.6-1       
[13] lazyeval_0.2.0   assertthat_0.1   digest_0.6.12    tibble_1.3.0    
[17] codetools_0.2-15 stringi_1.1.5    compiler_3.4.0   scales_0.4.1          
~~~


