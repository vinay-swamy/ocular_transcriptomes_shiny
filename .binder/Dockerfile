FROM rocker/binder:3.5.3
LABEL maintainer='vinay_swamy'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}



RUN wget https://github.com/vinay-swamy/ocular_transcriptomes_shiny/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-04-24/')); devtools::install_deps()"

# get app_data and dl_data

RUN wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1mGmsWnnMo4Ll6x6AkSw5OWEMzlUhsIsc' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1mGmsWnnMo4Ll6x6AkSw5OWEMzlUhsIsc" -O  dl_data.tar.gz  && tar -xzf dl_data.tar.gz && rm /tmp/cookies.txt && rm dl_data.tar.gz 

RUN wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=17Z3qrd2bKlpOaJrGnWJ71VkQkv_QM7G4' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=17Z3qrd2bKlpOaJrGnWJ71VkQkv_QM7G4" -O  app_data.tar.gz  && tar -xzf app_data.tar.gz && rm /tmp/cookies.txt && rm app_data.tar.gz 

RUN rm DESCRIPTION.1; exit 0

