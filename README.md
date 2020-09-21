# ocular_transcriptome_shiny

  <!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vinay-swamy/ocular_transcriptomes_shiny/master?urlpath=shiny)
  <!-- badges: end -->

R Shiny visualization for *de novo* transcriptomes.

To use run the following locally in R, or in the binder repo

```
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")# ignore R version mismatches
devtools::install_github('vinay-swamy/ocular_transcriptomes_shiny')
OcularTxome::downloadAppData()# get data for app
OcularTxome::runOcularTxomeVis()# run app
```
