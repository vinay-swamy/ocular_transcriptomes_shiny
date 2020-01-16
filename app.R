library(shiny)
library(tidyverse, quietly = T)
library(DBI)

#load('~/NIH/ocular_transcriptomes_shiny/data/new_data_dev.Rdata')
#load('~/NIH/ocular_transcriptomes_shiny/data/V2_shinydata.Rdata')
load('~/NIH/ocular_transcriptomes_shiny/data/V3_shinydata.Rdata')
source('~/NIH/ocular_transcriptomes_shiny/src/draw_transcripts.R')
source('~/NIH/ocular_transcriptomes_shiny/src/plotting_funcs.R')
#subtissues <- colnames(dev_piu)[-(1:2)]
db_name <- '~/NIH/ocular_transcriptomes_pipeline/testing/testdb.sql'
con <- dbConnect(RSQLite::SQLite(),db_name)
ui <- fluidPage(
    titlePanel("De Novo Transcriptomes"),
    fluidRow(
        column(4,#----
               selectInput(inputId = 'u_gene',
                           label = 'Enter Gene Name',
                           choices = all_gene_names,
                           #selected = 'CENPS',
                           selected = "USP1",
                           multiple = F)
               
        ),
        column(4,#----
               selectizeInput(inputId = 'u_tissues',
                           label = 'Select up to 5 tissues',
                           options = list(maxItems = 5),
                           choices = all_tissues,
                           selected = c('Retina_Adult.Tissue', 'RPE_Adult.Tissue'),
                           multiple = T)
        )
    ),
    h2('Percentage of total gene expression attributed to its transcripts expressed in selected tissues'),
    plotOutput('piu_plot', height = 500),
    h3('Exon Diagram of Transcripts for selected gene'),
    uiOutput('tx_diag'),
    h4('Fraction of samples each transcript was constructed in'),
    plotOutput('num_samp_det')
    
    
)#----


server <- function(input, output, session) {
    #load('~/personal/shiny_practice/data/gtf.Rdata')
    
    output$tx_diag <- renderUI({
        s_gene <- input$u_gene
        s_tissues <- input$u_tissues
        #keep_tx <- dev_tissue_det %>% select(s_tissues) %>% {rowSums(.) >0} %>% {dev_tissue_det[.,]} %>% pull(transcript_id)
        s_keep_tx <- con %>% tbl('tissue_det') %>% 
             select(transcript_id, s_tissues) %>% collect %>% {.[rowSums(.[,-1]) >0 ,]}%>% pull(transcript_id)
        s_gtf <- con %>% tbl('gtf') %>% filter(gene_name == s_gene, transcript_id %in% s_keep_tx) %>% collect
        s_piu <- con %>% tbl('piu') %>% filter(gene_name == s_gene, transcript_id %in% s_keep_tx) %>% collect
        s_frac_det <- con %>% tbl('frac_samp_det') %>% filter(gene_name == s_gene, transcript_id %in% s_keep_tx) %>% collect
        if(length(s_keep_tx) == 0){
            output$warning <- renderText('Gene not expressed in selected tissues')
            textOutput('warning')
        }else {
            height <- get_num_uniq_tx(gene = s_gene, gtf = s_gtf,keep_tx = s_keep_tx)
            
            title <- 'All Transcripts'
            
            output$piu_plot <- renderPlot({
                make_piu_bargraph(gene = s_gene, tissues = s_tissues, keep_tx = s_keep_tx,piu = s_piu)
            }, height = 500)
            
            output$transcript_diag <- renderPlot({
                draw_all_transcripts_static(gene = s_gene, tissues=s_tissues, gtf=s_gtf,keep_tx = s_keep_tx) 
                
            }, height = 25*(height+1) )
            
            output$num_samp_det <- renderPlot({
                make_num_det_bargraph(gene = s_gene, tissues = s_tissues, num_det = s_frac_det, keep_tx = s_keep_tx)
            })
            plotOutput('transcript_diag')
        }
    })
    
    
}


onStop(function() dbDisconnect(con))
shinyApp(ui, server)


#example gene USP1
