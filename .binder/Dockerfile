FROM rocker/binder:3.5.3
LABEL maintainer='vinay_swamy'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}



RUN wget https://github.com/vinay-swamy/ocular_transcriptomes_shiny/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-04-24/')); devtools::install_deps()"

RUN wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=19o_YpW3Egg42sD1WpOqb6gH-Mrt0wMHw' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=19o_YpW3Egg42sD1WpOqb6gH-Mrt0wMHw" -O  ocular_transcriptomes_shiny/app_data/DNTX_db.sql.gz  && gunzip ocular_transcriptomes_shiny/app_data/DNTX_db.sql.gz && rm /tmp/cookies.txt

RUN rm DESCRIPTION.1; exit 0

