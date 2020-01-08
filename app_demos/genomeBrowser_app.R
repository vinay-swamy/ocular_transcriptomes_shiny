library(shiny)
library(ggiraph)
setwd('~/NIH/ocular_transcriptomes_shiny/')
load('~/NIH/ocular_transcriptomes_shiny/data/chrom_names.Rdata')
load('~/NIH/ocular_transcriptomes_shiny/data/gtf_range.Rdata')
load('data/dev_gtf.Rdata')
source('~/NIH/ocular_transcriptomes_shiny/src/genomeBrowserTools.R')
chrom_names <- unique(gtf$seqid)
ui <- fluidPage(
   
    titlePanel("Genome Bowser"),
    fluidRow(
        column(4,
               selectInput(inputId = 'gb_chrom',
                           label = 'select a chromosome',
                           choices = chrom_names,
                           selected = 'chr1',
                           multiple = F)
               #textInput('p_height', label = 'Pick Plot Height',value = 5),
               #textInput('p_width', label='Pick Plot Width', value=6)
        ),
        column(4,
               textInput(inputId = 'gb_loc', 
                         label = 'select a location/coordinates.', 
                         value='15819122',
                         placeholder = "100000 or 50000-60000"
               )
        )
    ),
    
    textOutput(outputId = 'test')
    
)#----


server <- function(input, output) {
    output$test <- renderText({
        gb_chrom <- input$gb_chrom
        gb_loc <- input$gb_loc
        getOverlappingGTF_no_range(gtf, gb_chrom, as.numeric(gb_loc))}
    )
}

shinyApp(ui, server)