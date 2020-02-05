library(shiny)
suppressMessages(library(tidyverse))
library(DBI)

load('app_data/shiny_data.Rdata')
source('src/draw_transcripts.R')
source('src/plotting_funcs.R')

db_name <- 'app_data/DNTX_db.sql'
con <- dbConnect(RSQLite::SQLite(),db_name)

ui <- fluidPage(
    titlePanel("De Novo Transcriptomes"),
    tabsetPanel(
        tabPanel('plots', #----
            fluidRow(
                column(4,
                       selectInput(inputId = 'u_gene',
                                   label = 'Enter Gene Name',
                                   choices = all_gene_names,
                                   selected = "USP1",
                                   multiple = F)
                       
                ),
                column(4,
                       selectizeInput(inputId = 'u_tissues',
                                   label = 'Select up to 5 tissues',
                                   options = list(maxItems = 5),
                                   choices = all_tissues,
                                   selected = c('Retina_Adult.Tissue', 'RPE_Adult.Tissue'),
                                   multiple = T)
                ),
                column(2, actionButton('draw', label = 'click to redraw plot' ))
                
            ),
            h2('Percentage of total gene expression attributed to each transcript expressed in selected tissues'),
            plotOutput('piu_plot', height = 500),
            h3('Exon Diagram of Transcripts for selected gene'),
            uiOutput('tx_diag'),
            h4('Fraction of samples each transcript was constructed in'),
            plotOutput('num_samp_det')
        ),#----
        tabPanel('download', 
                    h2('Download De Novo Transcriptomes'),
                    column(3, radioButtons(inputId = 'dl_base_choice',
                                           label = 'select an option', 
                                           c('pan-body(54 subtissues)' = 'panbody', 
                                             'pan-eye(6 subtissues)'  = 'paneye',
                                             'by subtissue(eye only)' = 'subtissue')
                                           ),
                           downloadButton(outputId = 'dl_download',
                                          label = 'download')
                            ),
                    column(3, conditionalPanel("input.dl_base_choice == 'subtissue'",
                                               selectInput(inputId = 'dl_tis_choice', 
                                                           label = 'select a tissue',
                                                           choices = c('Retina', 'RPE', 'Cornea'),
                                                           selected = NULL
                                                           ),
                                               checkboxGroupInput(inputId = 'dl_dev_choice', 
                                                            label = NULL, 
                                                            choices = c('Adult' = 'Adult', 
                                                                        'Fetal' = 'Fetal'),
                                                            selected = 'Adult')
                                                           
                                               )
                           
                           )

                 )
        
    )
    
)#----


server <- function(input, output, session) {
    output$tx_diag <- renderUI({
        get_gene <- eventReactive(input$draw, {
            return(input$u_gene)
        })
        get_tissues <-eventReactive(input$draw, {
            return(input$u_tissues)
        })
        s_gene <- get_gene()
        s_tissues <- get_tissues()
        
        s_keep_tx <- con %>% tbl('tissue_det') %>% filter(gene_name ==s_gene) %>% 
             select(transcript_id, s_tissues) %>% collect %>% {.[rowSums(.[,-1]) >0 ,]}%>% pull(transcript_id)
        if(length(s_keep_tx) == 0){
            output$warning <- renderText('Gene not expressed in any selected tissues')
            textOutput('warning')
        }else {
            s_gtf <- con %>% tbl('plotting_gtf') %>% filter(gene_name == s_gene) %>% collect %>% filter(transcript_id %in% s_keep_tx)
            s_piu <- con %>% tbl('piu') %>% filter(gene_name == s_gene) %>% collect %>% filter(transcript_id %in% s_keep_tx)
            s_frac_det <- con %>% tbl('frac_samp_det') %>% filter(gene_name == s_gene) %>% collect %>% filter(transcript_id %in% s_keep_tx)

            height <- get_num_uniq_tx(gene = s_gene, gtf = s_gtf,keep_tx = s_keep_tx)
            
            output$piu_plot <- renderPlot({
                make_piu_bargraph(gene = s_gene, tissues = s_tissues, keep_tx = s_keep_tx,piu = s_piu)
            }, height = 500)
            
            output$transcript_diag <- renderGirafe({
                draw_all_transcripts_interactive_v3(gene = s_gene, 
                                                    tissues=s_tissues, 
                                                    gtf=s_gtf,
                                                    #cds_df = s_cds_df,
                                                    keep_tx = s_keep_tx, 
                                                    hsvg = height) 
                
            })
            
            output$num_samp_det <- renderPlot({
                make_num_det_bargraph(gene = s_gene, tissues = s_tissues, num_det = s_frac_det, keep_tx = s_keep_tx)
            })
            girafeOutput('transcript_diag')
    
        }
    })
    output$dl_download <- downloadHandler(
        filename = function(){
            if(input$dl_base_choice == 'subtissue'){
                paste0(paste0(input$dl_tis_choice,'_', input$dl_dev_choice, '.Tissue', collapse = '-'), '.tar.gz', collapse = '-')
            } else {
                paste0(input$dl_base_choice, '.tar.gz')
            }
        },
        content = function(file){
            
            if(input$dl_base_choice == 'subtissue'){
                fn=paste0('dl_data/', paste0(input$dl_tis_choice,'_', input$dl_dev_choice, '.Tissue', collapse = '-'), '.tar.gz', collapse = '-')
            } else {
                fn=paste0('dl_data/', input$dl_base_choice, '.tar.gz')
            }
            stopifnot(file.exists(fn))
            file.copy(from = fn, to = file)
        }
        
    )
}

onStop(function() dbDisconnect(con))
shinyApp(ui, server)



