library(tidyverse)
library(gridExtra)
library(ggiraph)


TEXT_HEIGHT=50
draw_super_transcript_static <- function(gene, gtf, tissue, wsvg=64, hsvg=1){
    psi_var <- paste0(tissue, '_psi')
    main_cols <- c("seqid", "type",  "start", "end", "strand", "transcript_id", "gene_name")
    gtf_gene <- filter(gtf, gene_name == gene)
    maxchar <- max(nchar(gtf_gene$transcript_id), nchar(gtf_gene$gene_name), 
                 nchar(gtf_gene$oId %>% na.omit %>% .[grepl('ENST', .)] ))
    gtf_exons <- filter(gtf_gene, type == 'exon') %>% 
    select(seqid, strand, start, end, gene_name,  PSI= !!psi_var) %>% 
    distinct %>% 
    arrange(start) %>% 
    mutate(length=end-start, length=sqrt(length), Xmin=0, Xmax=0, Ymin=-1, Ymax=1, exon_number=NA) 
    gap=mean(gtf_exons$length)
    gtf_exons[1,'Xmax'] <- gtf_exons[1,'Xmin'] + gtf_exons[1,'length']
    if(nrow(gtf_exons) >1){
        for(i in 2:nrow(gtf_exons)){
        gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
        gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
        }
    }
    
    plot <- ggplot(data = gtf_exons)+
    geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax), color='green')+
    facet_wrap(~gene_name, ncol=1, strip.position = 'left')+
    #ylim(c(-20,20))+
    #ggtitle('SuperTranscript')+
    theme_void() +
    theme(strip.text = element_text(angle = 180))
    #print(nchar(plot$data$gene_name))
    return(plot)
    #return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))
}

draw_super_transcript <- function(gene, gtf, tissue, wsvg=64, hsvg=1){
  psi_var <- paste0(tissue, '_psi')
  main_cols <- c("seqid", "type",  "start", "end", "strand", "transcript_id", "gene_name")
  gtf_gene <- filter(gtf, gene_name == gene)
  maxchar <- max(nchar(gtf_gene$transcript_id), nchar(gtf_gene$gene_name), 
                 nchar(gtf_gene$oId %>% na.omit %>% .[grepl('ENST', .)] ))
  gtf_exons <- filter(gtf_gene, type == 'exon') %>% 
    select(seqid, strand, start, end, gene_name,  PSI= !!psi_var) %>% 
    distinct %>% 
    arrange(start) %>% 
    mutate(length=end-start, length=sqrt(length), Xmin=0, Xmax=0, Ymin=-1, Ymax=1, exon_number=NA) 
  gap=mean(gtf_exons$length)
  gtf_exons[1,'Xmax'] <- gtf_exons[1,'Xmin'] + gtf_exons[1,'length']
  if(nrow(gtf_exons) >1){
    for(i in 2:nrow(gtf_exons)){
      gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
      gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
    }
  }
  
  plot <- ggplot(data = gtf_exons)+
    geom_rect_interactive( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax, tooltip= PSI), color='green')+
    facet_wrap(~gene_name, ncol=1, strip.position = 'left')+
    #ylim(c(-20,20))+
    #ggtitle('SuperTranscript')+
    theme_void() +
    theme(strip.text = element_text(angle = 180, size = TEXT_HEIGHT))
  #print(nchar(plot$data$gene_name))

  return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))
}




draw_all_transcripts <- function(gene, gtf, tissue, wsvg=64, hsvg=1){
    psi_var <- paste0(tissue, '_psi')
    main_cols <- c("seqid", "type",  "start", "end", "strand", "transcript_id", "gene_name", "exon_number" )
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
    if(nrow(gtf_exons) > 1){
        for(i in 2:nrow(gtf_exons)){
            gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
            gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
            }
    }
    
    data <- filter(gtf_gene, type == 'exon') %>%  select(main_cols, PSI= !!psi_var) %>% inner_join(gtf_exons) %>% 
    mutate(transcript_id=str_pad(transcript_id, maxchar, 'left'),
           gene_name=str_pad(string = gene_name,width = maxchar, side =  'left')) 
    plot_data <- data %>% 
        select(seqid, strand, start, start, end, transcript_id=gene_name, Xmin, Xmax, Ymin, Ymax, exon_number) %>% 
        mutate(exon_number=NA) %>% 
        distinct %>% bind_rows(data)
    
    plot <- ggplot(data = plot_data) +
        geom_rect_interactive( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax, tooltip= PSI), color='blue')+
        facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
        theme_void() +
        theme(strip.text = element_text(angle = 180, size =TEXT_HEIGHT/2))
    #print(nchar(plot$data$transcript_id))
    #return(plot)
    return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))
  
}

draw_all_transcripts_static <- function(gene, gtf, tissue, wsvg=64, hsvg=1){
  main_cols <- c("seqid", "type",  "start", "end", "strand", "transcript_id", "gene_name", "exon_number" )
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
  if(nrow(gtf_exons) > 1){
    for(i in 2:nrow(gtf_exons)){
      gtf_exons[i,'Xmin'] <- gtf_exons[(i-1),'Xmax'] + gap 
      gtf_exons[i,'Xmax'] <- gtf_exons[i,'Xmin'] + gtf_exons[i,'length']
    }
  }
  
  data <- filter(gtf_gene, type == 'exon') %>%  select(main_cols) %>% inner_join(gtf_exons) %>% 
    mutate(transcript_id=str_pad(transcript_id, maxchar, 'left'),
           gene_name=str_pad(string = gene_name,width = maxchar, side =  'left')) 
  plot_data <- data %>% 
    select(seqid, strand, start, start, end, transcript_id=gene_name, Xmin, Xmax, Ymin, Ymax, exon_number) %>% 
    mutate(exon_number=NA) %>% 
    distinct %>% bind_rows(data)
  
  plot <- ggplot(data = plot_data) +
    geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax, ), color='blue')+
    facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
    theme_void() +
    theme(strip.text = element_text(angle = 180))
  #print(nchar(plot$data$transcript_id))
  return(plot)
  #return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))
  
}
get_num_uniq_tx <- function(gene, gtf){
  return(gtf %>% filter(type=='transcript', gene_name == gene) %>% pull(transcript_id) %>% unique %>% length)
}




#draw_super_transcript('BEX2', gtf, 'Retina_Adult.Tissue')

