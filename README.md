# ocular_transcriptome_shiny

  <!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/vinay-swamy/ocular_transcriptomes_shiny/master?urlpath=shiny)
  <!-- badges: end -->

R Shiny visualization for *de novo* transcriptomes.

To use locally(recommended) run the following

```
git clone https://github.com/vinay-swamy/ocular_transcriptomes_shiny.git
cd ocular_transcriptomes_shiny
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1mGmsWnnMo4Ll6x6AkSw5OWEMzlUhsIsc' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1mGmsWnnMo4Ll6x6AkSw5OWEMzlUhsIsc" -O  dl_data.tar.gz  && tar -xzf dl_data.tar.gz && rm /tmp/cookies.txt && rm dl_data.tar.gz 

wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17Z3qrd2bKlpOaJrGnWJ71VkQkv_QM7G4' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17Z3qrd2bKlpOaJrGnWJ71VkQkv_QM7G4" -O  app_data.tar.gz  && tar -xzf app_data.tar.gz && rm /tmp/cookies.txt && rm app_data.tar.gz 
```

Then open in Rstudio

The above command requires the GNU version of sed. If you are using macOS, install gnu-sed with [homebrew](brew.sh) (brew install gnu-sed) and swap out `sed` with `gsed`


