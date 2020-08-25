# ocular_transcriptome_shiny

  <!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vinay-swamy/ocular_transcriptomes_shiny/master?urlpath=shiny)
  <!-- badges: end -->

R Shiny visualization for *de novo* transcriptomes.

To use locally(recommended) run the following

```
git clone https://github.com/vinay-swamy/ocular_transcriptomes_shiny.git
cd ocular_transcriptomes_shiny
wget http://hpc.nih.gov/~mcgaugheyd/ocular_transcriptomes_shiny/dl_data.tar.gz  && tar -xzf dl_data.tar.gz && rm dl_data.tar.gz 
wget http://hpc.nih.gov/~mcgaugheyd/ocular_transcriptomes_shiny/app_data.tar.gz   && tar -xzf app_data.tar.gz  && rm app_data.tar.gz  
```

Then open in Rstudio

The PanEye transcript annotation(gtf + fasta) is can be directly downloaded by
```
wget http://hpc.nih.gov/~mcgaugheyd/ocular_transcriptomes_shiny/paneye.zip
```


