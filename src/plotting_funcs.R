library(tidyverse)
library(gridExtra)
library(ggiraph)



# piu <- dev_piu
# gene <- s_gene
# tissues <- s_tissues


make_piu_bargraph <- function(gene, piu, tissues,keep_tx){
      exp_piu <- filter(piu, transcript_id %in%keep_tx, gene_name == gene) %>% select(transcript_id, gene_name, tissues) %>% 
          gather(key='subtissue', value = 'piu', -transcript_id, -gene_name)
      plot <- ggplot(data = exp_piu) +
          geom_col(aes(x=transcript_id, y=piu, fill=subtissue), position = 'dodge')+
          #ggtitle('fraction of total gene expression attributed to each transcript in selected tissues')+
          ylim(c(0,1)) +
          theme_minimal() +
          theme(axis.text.x=element_text(angle=45, hjust = 1))
      return(plot)
    
}


make_num_det_bargraph <- function(gene, num_det, tissues, keep_tx){
    exp_numdet <- num_det %>% filter(transcript_id %in% keep_tx, gene_name == gene) %>% 
        select(transcript_id, gene_name, tissues) %>% 
        gather(key='subtissue', value = 'Fraction of Samples', -transcript_id, -gene_name)
    plot <- ggplot(data = exp_numdet) +
        geom_col(aes(x=transcript_id, y=`Fraction of Samples`, fill=subtissue), position = 'dodge')+
        #ggtitle('fraction of samples transcript was contructed in for selected tissues')+
        ylim(c(0,1)) +
        theme_minimal() +
        theme(axis.text.x=element_text(angle=45, hjust = 1))
    return(plot)
}




