library(tidyverse)
getOverlappingGTF <- function(q_seq, q_start, q_end, s_rang){
    tibble(seqnames=q_seq, start=q_start, end=q_end) %>% 
        {plyranges::as_granges(.) } %>% 
        plyranges::join_overlap_inner(s_rang, .) %>% 
        as_tibble
    
}

