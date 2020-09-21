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
draw_all_transcripts_interactive_v1 <- function(gene, gtf, tissues,keep_tx, wsvg=64, hsvg=1,k=1 ){
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
        geom_rect_interactive( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`, tooltip=start))+
        scale_fill_manual(values = color_list) +
        facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
        #ggtitle(gene) +
        theme_void() +
        theme(strip.text = element_text(angle = 180, size = 50))
    #print(nchar(plot$data$transcript_id))
    #return(plot)
    girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg)
    #return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))

}

draw_all_transcripts_interactive_v3 <- function(gene, gtf, tissues,keep_tx, wsvg=64, hsvg=1){
    gtf_gene <- filter(gtf, gene_name == gene, transcript_id %in% keep_tx)
    maxchar <- max(nchar(gtf_gene$transcript_id), nchar(gtf_gene$gene_name),
                   nchar(gtf_gene$oId %>% na.omit %>% .[grepl('ENST', .)] ))

    plot_data <- gtf_gene %>% filter( type == 'exon') %>%
        mutate(transcript_id=str_pad(transcript_id, maxchar, 'left'),
               `exon type`=ifelse(is.na(novel_exon_id), 'ref', 'novel'))
    color_list <- c('black', 'red')
    names(color_list) <- c('ref', 'novel')
    intron_gaps <- tibble(Xmin=plot_data$Xmax[1:(nrow(plot_data) -1)], Xmax=plot_data$Xmin[2:(nrow(plot_data))]) %>%
        distinct %>%
        mutate(midpoint=rowMeans(.)) %>%
        filter(midpoint >0)

    point_type = ifelse(plot_data$strand[1] == '+',62, 60)
    plot <- ggplot(data = plot_data) +
        geom_hline(yintercept = 0, colour='black', size=2)+
        geom_point(data=intron_gaps, aes(x=midpoint),y=0,size=20, fill='black', shape=point_type)+
        geom_rect_interactive(aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`, tooltip=ttip))+
        #geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`))+
        scale_fill_manual(values = color_list) +
        facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
        #ggtitle(gene) +
        theme_void() +
        theme(strip.text = element_text(angle = 180, size = 50),
              legend.title = element_text(size=80), legend.text =element_text(size = 70))
    #print(nchar(plot$data$transcript_id))
    #return(plot)
    #girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg)
    return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))

}



get_num_uniq_tx <- function(gene, gtf, keep_tx){
    return(gtf %>% filter(type=='transcript', gene_name == gene, transcript_id %in% keep_tx) %>% pull(transcript_id) %>% unique %>% length)
}

c_merge <- function(mat){
    #mat must bed sorted by start and interval length
    res <- list()
    cur <- mat[1, ]
    j <- 1
    for( i in 2:nrow(mat)){
        #1=start, 2=end#
        if(mat[i, 1]<= cur[1,  2]){
            #merge the interval
            if(mat[i, 2] > cur[1, 2]){
                cur[1, 2] <- mat[i, 2]
                cur[1, 3] <- paste(cur[1, 3],mat[i, 3], sep='-')
            }else{
                cur[1, 3] <- paste(cur[1, 3],mat[i, 3], sep='-')
            }

        }else{
            # interval is good
            if(any(cur %in% res)) {cur <- mat[i, ] }
            res[[j]] <- cur
            cur <- mat[i, ]
            j <- j+1
        }

    }
    if(!any(cur %in% res)) {
        #cur[1, 3] <- paste(cur[1, 3],mat[i, 3], sep='-')
        res[[j+1]] <- cur
    }
    return(do.call(rbind, res))

}



draw_all_transcripts_interactive_v4 <- function(gene, gtf, keep_tx, g, wsvg=64, hsvg=1){
    #dynamically scale introns, because it doesnt work well precopmuting it for genes with a lot of exons
    gtf_gene <- filter(gtf, gene_name == gene, transcript_id %in% keep_tx)
    distinct_exons <- gtf_gene %>% filter(type == 'exon') %>%
        select(seqid, strand, start, end, Xmin, Xmax) %>%
        distinct %>%
        mutate(length = sqrt(end-start),) %>%
        arrange(start, length) %>%
        mutate(id= as.character(1:nrow(.))) %>%
        select(-length)

    fin <- distinct_exons %>% select(Xmin, Xmax, id) %>% c_merge %>%
        mutate(length = mean(Xmax-Xmin))


    gap=mean(fin$length)/g
    min_igap <- {fin$Xmin[2:nrow(fin)] - fin$Xmax[1:(nrow(fin )-1)]}
    min_igap[min_igap<0] <- -Inf
    min_igap_fail <-  which(min_igap>gap)
    for(idx in min_igap_fail){
        fin_idx= idx+1
        #if(fin[fin_idx,'length'] >=gap){ gap <- fin[fin_idx,'length']+gap  }
        delta = gap- min_igap[idx]
        if(delta>0) print('MOOOOOO')
        nfin <- nrow(fin)
        fin[fin_idx:nfin,'Xmin'] <- fin[fin_idx:nfin,'Xmin']+delta
        fin[fin_idx:nfin,'Xmax'] <- fin[fin_idx:nfin,'Xmax']+delta

    }
    correct <-  fin %>%
        filter(!grepl('-',id))
    to_correct <-  fin %>%
        filter(grepl('-',id)) %>%
        mutate(id= str_split(id, '-')) %>%
        unnest(id) %>%
        rename(new_xmin = Xmin, new_xmax = Xmax) %>%
        inner_join(distinct_exons)

    res <- lapply(unique(to_correct$new_xmin), function(x) filter(to_correct, new_xmin == x) %>%
                      mutate(d=(min(Xmin)-x), Xmin =Xmin - d , Xmax =Xmax -d)) %>% bind_rows %>%
        select(all_of(colnames(correct)))
    all_correct <- bind_rows(correct, res, ) %>% arrange(Xmin) %>% select(-length)


    all_plot_data <- distinct_exons %>%
        select(-Xmax, -Xmin) %>%
        inner_join(all_correct) %>%
        inner_join(gtf_gene %>% select(-Xmin, -Xmax, -length)) %>%
        mutate(`exon type`=ifelse(is.na(novel_exon_id), 'ref', 'novel'))
    plot_data <- all_plot_data %>% filter(type == 'exon')
    color_list <- c('black', 'red')
    names(color_list) <- c('ref', 'novel')
    plot <- ggplot(data = plot_data) +
        geom_hline(yintercept = 0, colour='black', size=2)+
        geom_rect_interactive(aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`, tooltip=ttip))+
        #geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`))+
        scale_fill_manual(values = color_list) +
        facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
        #ggtitle(gene) +
        theme_void() +
        theme(strip.text.y.left = element_text(angle = 0, size = 50),
              legend.title = element_text(size=80), legend.text =element_text(size = 70))
    #print(nchar(plot$data$transcript_id))
    #return(plot)
    #girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg)
    return(girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg))

}

draw_all_transcripts_static_v4 <- function(gene, gtf, keep_tx, g){
    #dynamically scale introns, because it doesnt work well precopmuting it for genes with a lot of exons
    gtf <- new_plotting_gtf
    gene <- 'ABCA4'
    keep_tx <- keep_tx
    gtf_gene <- filter(gtf, gene_name == gene, transcript_id %in% keep_tx)
    distinct_exons <- gtf_gene %>% filter(type == 'exon') %>%
        select(seqid, strand, start, end, Xmin, Xmax) %>%
        distinct %>%
        mutate(length = sqrt(end-start),) %>%
        arrange(start, length) %>%
        mutate(id= as.character(1:nrow(.))) %>%
        select(-length)

    fin <- distinct_exons %>% select(Xmin, Xmax, id) %>% c_merge %>%
        mutate(length = mean(Xmax-Xmin))


    gap=mean(fin$length)/g
    min_igap <- {fin$Xmin[2:nrow(fin)] - fin$Xmax[1:(nrow(fin )-1)]}
    min_igap[min_igap<0] <- -Inf
    min_igap_fail <-  which(min_igap>gap)
    for(idx in min_igap_fail){
        fin_idx= idx+1
        #if(fin[fin_idx,'length'] >=gap){ gap <- fin[fin_idx,'length']+gap  }
        delta = gap- min_igap[idx]
        if(delta>0) print('MOOOOOO')
        nfin <- nrow(fin)
        fin[fin_idx:nfin,'Xmin'] <- fin[fin_idx:nfin,'Xmin']+delta
        fin[fin_idx:nfin,'Xmax'] <- fin[fin_idx:nfin,'Xmax']+delta

    }
    correct <-  fin %>%
        filter(!grepl('-',id))
    to_correct <-  fin %>%
        filter(grepl('-',id)) %>%
        mutate(id= str_split(id, '-')) %>%
        unnest(id) %>%
        rename(new_xmin = Xmin, new_xmax = Xmax) %>%
        inner_join(distinct_exons)

    res <- lapply(unique(to_correct$new_xmin), function(x) filter(to_correct, new_xmin == x) %>%
                      mutate(d=(min(Xmin)-x), Xmin =Xmin - d , Xmax =Xmax -d)) %>% bind_rows %>%
        select(all_of(colnames(correct)))
    all_correct <- bind_rows(correct, res, ) %>% arrange(Xmin) %>% select(-length)


    all_plot_data <- distinct_exons %>%
        select(-Xmax, -Xmin) %>%
        inner_join(all_correct) %>%
        inner_join(gtf_gene %>% select(-Xmin, -Xmax, -length)) %>%
        mutate(`exon type`=ifelse(is.na(novel_exon_id), 'ref', 'novel'))
    plot_data <- plot_data %>% filter(type == 'exon')
    color_list <- c('black', 'red')
    names(color_list) <- c('ref', 'novel')

    point_type = ifelse(plot_data$strand[1] == '+',62, 60)
    plot <- ggplot(data = plot_data) +
        geom_hline(yintercept = 0, colour='black', size=2)+
        geom_rect(aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`))+
        #geom_rect( aes(xmin=Xmin, xmax=Xmax, ymin=Ymin, ymax=Ymax,fill=`exon type`))+
        scale_fill_manual(values = color_list) +
        facet_wrap(~transcript_id, ncol=1, strip.position = 'left')+
        #ggtitle(gene) +
        theme_void() +
        theme(strip.text = element_text(angle = 180, size = 50),
              legend.title = element_text(size=80), legend.text =element_text(size = 70))
        #theme(strip.text = element_text(angle = 180))
    #print(nchar(plot$data$transcript_id))
    #return(plot)
    #girafe(ggobj = plot, width_svg = wsvg, height_svg = hsvg)
    return(plot)

}
