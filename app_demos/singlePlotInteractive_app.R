library(shiny)
library(ggiraph)

load('~/personal/shiny_practice/data/dev_shiny_data.Rdata')
source('~/personal/shiny_practice/src/draw_transcripts_singleplot.R')

ui <- fluidPage(
    titlePanel("TestApp"),
    fluidRow(
        column(4,
               selectInput(inputId = 'gene_entry',
                           label = 'Enter Gene Name',
                           choices = gene_names,
                           multiple = F),
               radioButtons('radio', 'Display All Transcripts?', choices=c('ALL', 'NONE'), selected='NONE')
               #textInput('p_height', label = 'Pick Plot Height',value = 5),
               #textInput('p_width', label='Pick Plot Width', value=6)
        ),
        column(4,
               selectInput(inputId = 'tissue_entry', 
                           label = 'Select a Tissue', 
                           choices=tissues, 
                           multiple = F,
                           selected='Retina_Adult.Tissue')
               )
    ),
    
    uiOutput('ui')
    
)#----


server <- function(input, output) {
    #load('~/personal/shiny_practice/data/gtf.Rdata')
    
    output$ui <- renderUI({
        t_gene <- input$gene_entry
        t_tissue <- input$tissue_entry
        height <- get_num_uniq_tx(gene = t_gene, gtf = gtf)
        # output$plot <- renderggiraph({
        #             t_gene <- input$gene_entry
        #             #height = 25
        #             draw_super_transcript(gene = t_gene, gtf=gtf, tissue = t_tissue, 
        #                                   wsvg=as.numeric(input$p_width), hsvg=as.numeric(input$p_height))
        #                 #layout(xaxis=ax, yaxis=ax)
        #         })
        # 
        if(input$radio == 'NONE'){
            #h4('Super Transcript')
            title <- 'SuperTranscript'
            height = 25
            output$plot <- renderGirafe({
                t_gene <- input$gene_entry
                draw_super_transcript(gene = t_gene, gtf=gtf, tissue = t_tissue)#, 
                                      #wsvg = 4, hsvg = 4)

            })
        } else {# if input$radio == ALL
            #h4('All Transcripts')
            title <- 'All Transcripts'
            height = 1*height
            output$plot <- renderGirafe({
                draw_all_transcripts(gene = t_gene, gtf=gtf, tissue = t_tissue, 
                                    hsvg =height)

            } )

        }
        
        girafeOutput('plot')
    })
}

shinyApp(ui, server)
