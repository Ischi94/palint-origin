
<!-- add all libraries before knitting -->

# Paleoclimate Interaction Origination

<!-- badges: start -->

[![Generic
badge](https://img.shields.io/badge/Status-Published-green.svg)](https://shields.io/)
<!-- badges: end -->

This the R project to the second (potential) publication for PastKey
covering paleoclimate interactions and origination probability. This
repository contains all data files and code to reproduce results. All
figures are made with `ggplot` and can be therefore reproduced with this
repository as well. The raw fossil occurrence data set downloaded from
the Paleoiology Database is unfortunately too large for this repository
but is publically available on our [Figshare.com
archive](https://figshare.com/articles/dataset/Raw_fossil_data/14528925).

# File structure

Throughout our analysis, we use the `here` package to structure our data
files, which is based on the root directory of the *Rproj* file. We
recommend to either download the entire project or to fork the
repository to keep the directories tidy. If done so, it is not needed to
set working directories manually. All files will work correctly without
changing anything in the code. Every folder contains a README with
further instructions and explanations on individual files.

## R

This folder contains all R code as scripts.

## data

This folder contains both raw and processed files, as well as results
saved as *RData* or *csv* files.

## figures

This folder contains all figures from the main text and the
supplementary information file, as well as supplementary tables.

# Version info

Information about the operating system used to run the analysis and all
R packages including their dependencies and version info used throughout
the R-scripts:

``` r
sessionInfo()
```

    ## R version 4.0.4 (2021-02-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Linux Mint 20.1
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
    ##  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
    ##  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] effsize_0.8.1    infer_0.5.4      geiger_2.0.7     ape_5.4-1       
    ##  [5] lme4_1.1-26      Matrix_1.3-2     brms_2.15.0      Rcpp_1.0.6      
    ##  [9] divDyn_0.8.0     deeptime_0.0.5.2 patchwork_1.1.1  officer_0.3.16  
    ## [13] flextable_0.6.2  here_1.0.1       forcats_0.5.0    stringr_1.4.0   
    ## [17] dplyr_1.0.6      purrr_0.3.4      readr_1.4.0      tidyr_1.1.3     
    ## [21] tibble_3.1.2     ggplot2_3.3.5    tidyverse_1.3.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.3.1            uuid_0.1-4              backports_1.2.1        
    ##   [4] fastmatch_1.1-0         systemfonts_1.0.1       plyr_1.8.6             
    ##   [7] igraph_1.2.6            splines_4.0.4           crosstalk_1.1.1        
    ##  [10] rstantools_2.1.1        inline_0.3.17           digest_0.6.27          
    ##  [13] htmltools_0.5.1.1       rsconnect_0.8.16        fansi_0.5.0            
    ##  [16] phytools_0.7-70         magrittr_2.0.1          modelr_0.1.8           
    ##  [19] RcppParallel_5.0.2      matrixStats_0.57.0      xts_0.12.1             
    ##  [22] prettyunits_1.1.1       colorspace_2.0-1        rvest_0.3.6            
    ##  [25] haven_2.3.1             xfun_0.22               callr_3.5.1            
    ##  [28] crayon_1.4.1            jsonlite_1.7.2          phangorn_2.5.5         
    ##  [31] zoo_1.8-8               glue_1.4.2              gtable_0.3.0           
    ##  [34] V8_3.4.0                pkgbuild_1.2.0          rstan_2.21.2           
    ##  [37] maps_3.3.0              abind_1.4-5             scales_1.1.1           
    ##  [40] mvtnorm_1.1-1           DBI_1.1.1               miniUI_0.1.1.1         
    ##  [43] plotrix_3.8-1           xtable_1.8-4            tmvnsim_1.0-2          
    ##  [46] subplex_1.6             deSolve_1.28            stats4_4.0.4           
    ##  [49] StanHeaders_2.21.0-7    DT_0.17                 htmlwidgets_1.5.3      
    ##  [52] httr_1.4.2              threejs_0.3.3           ellipsis_0.3.2         
    ##  [55] pkgconfig_2.0.3         loo_2.4.1               dbplyr_2.0.0           
    ##  [58] utf8_1.2.1              tidyselect_1.1.1        rlang_0.4.11           
    ##  [61] reshape2_1.4.4          later_1.1.0.1           munsell_0.5.0          
    ##  [64] cellranger_1.1.0        tools_4.0.4             cli_2.5.0              
    ##  [67] generics_0.1.0          broom_0.7.7             ggridges_0.5.3         
    ##  [70] evaluate_0.14           fastmap_1.1.0           yaml_2.2.1             
    ##  [73] processx_3.4.5          knitr_1.33              fs_1.5.0               
    ##  [76] zip_2.1.1               nlme_3.1-152            mime_0.10              
    ##  [79] projpred_2.0.2          xml2_1.3.2              compiler_4.0.4         
    ##  [82] bayesplot_1.8.0         shinythemes_1.2.0       rstudioapi_0.13        
    ##  [85] gamm4_0.2-6             curl_4.3                clusterGeneration_1.3.7
    ##  [88] reprex_1.0.0            statmod_1.4.35          stringi_1.6.2          
    ##  [91] ps_1.5.0                Brobdingnag_1.2-6       gdtools_0.2.3          
    ##  [94] lattice_0.20-41         nloptr_1.2.2.2          markdown_1.1           
    ##  [97] shinyjs_2.0.0           vctrs_0.3.8             pillar_1.6.1           
    ## [100] lifecycle_1.0.0         combinat_0.0-8          bridgesampling_1.0-0   
    ## [103] data.table_1.13.6       httpuv_1.5.5            R6_2.5.0               
    ## [106] promises_1.1.1          gridExtra_2.3           codetools_0.2-18       
    ## [109] boot_1.3-27             colourpicker_1.1.0      MASS_7.3-53.1          
    ## [112] gtools_3.8.2            assertthat_0.2.1        rprojroot_2.0.2        
    ## [115] withr_2.4.2             mnormt_2.0.2            shinystan_2.5.0        
    ## [118] expm_0.999-6            mgcv_1.8-33             parallel_4.0.4         
    ## [121] hms_1.0.0               quadprog_1.5-8          grid_4.0.4             
    ## [124] coda_0.19-4             minqa_1.2.4             rmarkdown_2.6          
    ## [127] ggnewscale_0.4.5        scatterplot3d_0.3-41    numDeriv_2016.8-1.1    
    ## [130] shiny_1.6.0             lubridate_1.7.9.2       base64enc_0.1-3        
    ## [133] dygraphs_1.1.1.6
