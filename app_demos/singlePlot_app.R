library(shiny)

load('~/personal/shiny_practice/data/gene_names.Rdata')
load('~/personal/shiny_practice/data/gtf.Rdata')
source('~/personal/shiny_practice/src/draw_transcripts_singleplot.R')

ui <- fluidPage(
  titlePanel("Smooby App"),
  fluidRow(
    column(4,#----
           selectInput(inputId = 'textbox',
                       label = 'Enter Gene Name',
                       choices = gene_names,
                       multiple = F),
           radioButtons('radio', 'Display All Transcripts?', choices=c('ALL', 'NONE'), selected='NONE')
      )
    ),
    
   uiOutput('ui')
    
)#----


server <- function(input, output) {
  #load('~/personal/shiny_practice/data/gtf.Rdata')

  output$ui <- renderUI({
    t_gene <- input$textbox
    height <- get_num_uniq_tx(gene = t_gene, gtf = gtf)
    if(input$radio == 'NONE'){
        #h4('Super Transcript')
        title <- 'SuperTranscript'
        output$plot <- renderPlot({
            t_gene <- input$textbox
            draw_super_transcript(gene = t_gene,tissue='Retina_Fetal.Tissue', gtf=gtf) 
            }, height = 25) 
        } else {
        #h4('All Transcripts')
            title <- 'All Transcripts'
            output$plot <- renderPlot({
            draw_all_transcripts(gene = t_gene,tissue='Retina_Fetal.Tissue', gtf=gtf) 
        
        }, height = 25*(height+1) )

     }
    
    plotOutput('plot')
  })
}

shinyApp(ui, server)
