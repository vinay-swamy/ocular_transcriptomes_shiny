library(tidyverse)
getOverlappingGTF <- function(q_seq, q_start,  s_rang){
    q_start <- as.numeric(q_start)
    q_end=q_start
    res <- tibble(seqnames=q_seq, start=q_start, end=q_end) %>% 
        {plyranges::as_granges(.) } %>% 
        plyranges::join_overlap_inner(s_rang, .) %>% 
        as_tibble
    return(res$transcript_id[1])
}
getOverlappingGTF_no_range <- function(gtf,q_seq, q_start, q_end=q_start+1000, q_strand = '+'){
    stopifnot(q_start<q_end)
    df <- gtf %>% filter(seqid == q_seq, start >= q_start, start <= q_end, strand == q_strand)
    return(df$transcript_id %>% paste0(collapse = ':'))
}

