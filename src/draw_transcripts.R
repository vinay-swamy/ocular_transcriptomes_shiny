library(tidyverse)
library(gridExtra)
library(ggiraph)

draw_all_transcripts_static <- function(gene, gtf, tissues,keep_tx, wsvg=64, hsvg=1){
    main_cols <- c("seqid", "type",  "start", "end", "strand", "transcript_id", "gene_name", "exon_number" )
    gtf_gene <- filter(gtf, gene_name == gene, transcript_id %in% keep_tx)
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
    if(nrow(gtf_exons) > 1){
        for(i in 2:nrow(gtf_exons)){
            gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
            gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
        }
    }
    
    plot_data <- filter(gtf_gene, type == 'exon') %>%  select(main_cols,novel_exon_id) %>% inner_join(gtf_exons) %>% 
        mutate(transcript_id=str_pad(transcript_id, maxchar, 'left'),
               gene_name=str_pad(string = gene_name,width = maxchar, side =  'left'),
               `exon type`=ifelse(is.na(novel_exon_id), 'ref', 'novel')) 
    color_list <- c('black', 'red')
    names(color_list) <- c('ref', 'novel')
    plot <- ggplot(data = plot_data) +
        geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`))+
        scale_fill_manual(values = color_list) +
        facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
        ggtitle(gene) +
        theme_void() +
        theme(strip.text = element_text(angle = 180))
    #print(nchar(plot$data$transcript_id))
    return(plot)
    #return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))
    
}
get_num_uniq_tx <- function(gene, gtf, keep_tx){
    return(gtf %>% filter(type=='transcript', gene_name == gene, transcript_id %in% keep_tx) %>% pull(transcript_id) %>% unique %>% length)
}
