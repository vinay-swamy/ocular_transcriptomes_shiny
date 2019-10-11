library(tidyverse)
gtf_raw <- rtracklayer::readGFF('dl_data/all_tissues.combined.gtf.gz')
psi_tab <- read_tsv('dl_data/all_tissues_psi.tsv.gz')
trans_tab <- read_tsv('dl_data/gfc_TCONS_to_st_MSTRG.tsv.gz')
colnames(trans_tab) <- c('transcript_id', paste0('exp_',colnames(trans_tab)[-1]))
colnames(psi_tab) <-  colnames(psi_tab) %>% str_split('_psi') %>% sapply(function(x) x[1])
psi_tab[is.na(psi_tab)] <- 0
sample_table <- read_tsv('dl_data/sampleTableV5.tsv', col_names = c('sample', 'run', 'paired', 'tissue', 'subtissue', 'origin')) %>% 
    filter(sample %in% colnames(psi_tab))

average_tissue_psi <- function(t_tissue){
    samples <- filter(sample_table, subtissue == t_tissue) %>% pull(sample)
    avg_psi <- psi_tab[,samples] %>% rowMeans
    return(avg_psi)
}
tissues <- unique(sample_table$subtissue)
#dont forget to adjust rmats coordinates
psi_by_tissue <-  lapply(tissues, average_tissue_psi) %>% bind_cols(psi_tab[,1:4],.) %>% mutate(start=start+1)
colnames(psi_by_tissue) <- c(colnames(psi_tab)[1:4], paste0(tissues, '_psi'))

#first, format transcript names and gene names so they are all the same length
#ENST is 17 chars, TCONS is 14 chars
max_char=max(nchar(gtf_raw$transcript_id), nchar(gtf_raw$oId %>% na.omit %>% .[grepl('ENST', .)]))


enst_tx <- filter(gtf_raw, grepl('ENST', oId)) %>% select(transcript_id, ensid=oId) %>% distinct
gtf <- gtf_raw %>% select(seqid, type, start, end, strand, transcript_id, gene_name, exon_number) %>% 
    left_join(psi_by_tissue) %>% left_join(trans_tab) %>% 
    left_join(enst_tx) %>% 
    mutate(transcript_id=replace(transcript_id,!is.na(ensid), ensid[!is.na(ensid)]))
gene_names <- unique(gtf$gene_name)
save(gtf, gene_names, tissues, file = 'data/shiny_data.Rdata')

set.seed(233234)
gtf <- sample(unique(gtf$gene_name), 1000) %>% {filter(gtf, gene_name %in% .)}
gene_names <- unique(gtf$gene_name)
save(gtf, gene_names, tissues, file = 'data/dev_shiny_data.Rdata')








