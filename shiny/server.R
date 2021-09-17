#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readxl)
library(DT)
library(ggplot2)
library(plotly)
library(data.table)
library(viridis)

load("data/mydata.rdata")

#####

#####



# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    #corr table
    gene.lipid.corr_data <- eventReactive(
        input$search_gene.lipid.corr,{
            #validate(
            #    need(nrow(gene_lipid.corr[Gene==input$gene,]) != 0, 'No correlation data fund.'),
            #    errorClass = "red"
            #)
            
            #converted to data.table, faster
            gene_lipid.corr[Gene==input$gene,]
        }
    )
    
    #output DT table with default buttons, server=T link to filtered data
    output$corr.table <- 
            DT::renderDT(server = FALSE, {
                validate(
                    need(nrow(gene.lipid.corr_data()) != 0, 'No correlation data fund.'),
                    errorClass = "red"
                )
                DT::datatable(gene.lipid.corr_data(), selection = 'single',
                              extensions=c("Buttons",'Scroller'),
                              options = list(dom = 'Bfrtip', 
                                    buttons = list(list(extend = 'copy'), 
                                                list(extend = "csv", 
                                                text = "Download csv", 
                                                filename = "correlation",
                                                exportOptions = list(
                                                                #modifier = list(page = "all")
                                                                    )
                                                )),
                               scrollY = 350,
                               scroller = TRUE
                               )
                             )
                
                
        })
    
    #output corr scatter plot
    scat.lipid <- eventReactive(input$corr.table_rows_selected, {
        gene.lipid.corr_data()$Lipid[input$corr.table_rows_selected]
    })
    scat.gene <- eventReactive(input$corr.table_rows_selected, {
        gene.lipid.corr_data()$Gene[input$corr.table_rows_selected]
    })
    
    genelip.scat <- eventReactive(input$corr.table_rows_selected, {
        gs.p <- rna_livlipid_merge[c(scat.gene(), scat.lipid(), "Group")] %>%
            ggplot(aes(x=!!as.symbol(scat.gene()), y=!!as.symbol(scat.lipid()),
                       color=Group) ) +
            geom_point(size=2.5, na.rm=TRUE, alpha=0.7)+
            scale_color_viridis(discrete = T)+
            theme_classic()
        gs.p <- ggplotly(gs.p, height = 500, width = 680)
    })
    
    output$genelip.scatplot <- renderPlotly({
        genelip.scat()
    })
    
    

    #gene boxplot
    gene.box <- eventReactive(input$search_gene.lipid.corr,{
        validate(
            need(input$gene %in% names(rnaseq), 'No gene name fund.'),
            errorClass = "red"
        )
        p <- rnaseq %>% select("Group", input$gene, "AB", "Index") %>%
            ggplot(aes(x = Group, y = !!as.symbol(input$gene), color=Group,
                       text=paste(AB, Index))) +
            geom_boxplot(color="black") +
            geom_jitter(size = 2, alpha = 0.7,
                        position=position_jitter(width=0.1, height=0),na.rm=TRUE)+
            scale_colour_viridis(discrete = T)+
            theme_classic()+
            theme(legend.position="bottom")
        p <- ggplotly(p, height = 650, width = 900) %>% 
            layout(legend = list(orientation = "h", x = 0.25, y = -0.1))
        p$x$data[1] <- lapply(p$x$data[1], 
                              FUN = function(x){
                                  x$marker = list(opacity = 0)
                                  return(x)
                              })
        p
    })
    output$gene.boxplot <- renderPlotly({
        #validate(
        #    need(input$gene %in% names(rnaseq), 'No gene name fund.'),
        #    errorClass = "red"
        #)
            gene.box()
        })
    
    #warning message if no gene match
    observeEvent(input$search_gene.lipid.corr, {
        if (!input$gene %in% names(rnaseq))
            showNotification(
                "No name matches!", 
                type = "warning", closeButton = T)
    })
    
    #gene list
    genelist <- data.frame(Gene=colnames(rnaseq)[-c(1:3)])
    output$gene.list <- 
        DT::renderDT(server = F, 
                     genelist,
                     extensions="Buttons",
                     options = list(dom = 'Blfrtip',
                                    buttons = list("copy", 
                                                list(extend="csv",
                                                     text = "Download csv", 
                                                     filename = "genelist")
                                    )
                        )
        )
    
    
 
    
    
    
    #read SpList
#    spList <- reactive({
#        spList.loc <- input$spList
#        if(is.null(spList.loc))
#            return(NULL)
#        spList <- read_excel(spList.loc$datapath, sheet = "SpList")
#        return(spList)
#    })
    

#    output$table.SpList <- reactive({
#        renderTable({
#        spList.loc <- input$spList
#        if(is.null(spList.loc))
#            return(NULL)
#        read_excel(spList.loc$datapath, sheet = "SpList")})
#    })
    
    #read species and compute labels
#    species <- reactive({
#        return(input$sp)
#    })
#    sp.full <- spList
    
}
)
