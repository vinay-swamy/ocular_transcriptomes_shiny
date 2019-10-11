library(tidyverse)



draw_super_transcript <- function(gene, gtf){
  gtf_gene <- filter(gtf, gene_name == gene)
  maxchar <- max(nchar(gtf_gene$transcript_id), nchar(gtf_gene$gene_name), 
                 nchar(gtf_gene$oId %>% na.omit %>% .[grepl('ENST', .)] ))
  gtf_exons <- filter(gtf_gene, type == 'exon') %>% 
    select(seqid, strand, start, end, gene_name) %>% 
    distinct %>% 
    arrange(start) %>% 
    mutate(length=end-start, length=sqrt(length), Xmin=0, Xmax=0, Ymin=-1, Ymax=1,
           #transcript_id=str_pad(transcript_id, 'left', maxchar),
           gene_name=str_pad(string = gene_name,width = maxchar, side =  'left')) 
  gap=mean(gtf_exons$length)
  gtf_exons[1,'Xmax'] <- gtf_exons[1,'Xmin'] + gtf_exons[1,'length']
  for(i in 2:nrow(gtf_exons)){
      gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
      gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
    }
  plot <- ggplot(data = gtf_exons, 
                 aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax))+
   geom_rect(color='green')+
   facet_wrap(~gene_name, ncol=1, strip.position = 'left')+
   #ylim(c(-20,20))+
   #ggtitle('SuperTranscript')+
   theme_void() +
    theme(strip.text = element_text(angle = 180))
  #print(nchar(plot$data$gene_name))
  return(plot)
}

draw_all_transcripts <- function(gene, gtf){
  gtf_gene <- filter(gtf, gene_name == gene)
  maxchar <- max(nchar(gtf_gene$transcript_id), nchar(gtf_gene$gene_name), 
                 nchar(gtf_gene$oId %>% na.omit %>% .[grepl('ENST', .)] ))
  unique_tx <- length(unique(gtf_gene$transcript_id))
  gtf_exons <- filter(gtf_gene, type == 'exon') %>% 
    select(seqid, strand, start, end) %>% 
    distinct %>% 
    arrange(start) %>% 
    mutate(length=end-start, length=sqrt(length), Xmin=0, Xmax=0, Ymin=-1, Ymax=1) %>% 
    mutate(lab1=paste0('sdfs',1:nrow(.)), 
           lab2=paste0('wef',1:nrow(.)), 
           lab3=paste0('gfr',1:nrow(.)) )
  gap=mean(gtf_exons$length)
  gtf_exons[1,'Xmax'] <- gtf_exons[1,'Xmin'] + gtf_exons[1,'length']
  for(i in 2:nrow(gtf_exons)){
    gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
    gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
  }
  plot <- filter(gtf_gene, type == 'exon') %>% inner_join(gtf_exons) %>% 
          mutate(transcript_id=str_pad(transcript_id, maxchar, 'left'),
                 gene_name=str_pad(string = gene_name,width = maxchar, side =  'left')) %>% 
          ggplot(aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax, text=lab1, label=lab2, smooby=lab3))+
            geom_rect(color='blue') +
            facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
            theme_void() +
            theme(strip.text = element_text(angle = 180))
  #print(nchar(plot$data$transcript_id))
  return(plot)
  
}

get_num_uniq_tx <- function(gene, gtf){
  return(gtf %>% filter(type=='transcript', gene_name == gene) %>% pull(transcript_id) %>% unique %>% length)
}
# 
# 






