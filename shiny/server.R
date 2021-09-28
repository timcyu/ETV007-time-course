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
library(shinybusy)


load("data/mydata.rdata")

#####

#####



# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    
    #server-side loading input list
    updateSelectizeInput(
        session, 'gene', 
        choices = genelist, 
        server = TRUE, selected = "Apoh")
    
    
    #gene boxplot
    gene.box <- eventReactive(input$search_gene.lipid.corr, ignoreNULL = FALSE,
                              {
                                  validate(
                                      need(input$gene %in% genelist, 'No gene name fund.'),
                                      errorClass = "red"
                                  )
                                  p <- rna_livlipid_merge %>% select("Group", input$gene, "sampleID") %>%
                                      ggplot(aes(x = Group, y = !!as.symbol(input$gene), color=Group,
                                                 text=paste("ID ",sampleID))) +
                                      geom_boxplot(color="black") +
                                      geom_jitter(size = 2, alpha = 0.7,
                                                  position=position_jitter(width=0.1, height=0,seed=123),
                                                  na.rm=TRUE)+
                                      scale_colour_viridis(discrete = T)+
                                      theme_classic()+
                                      theme(legend.position="bottom")
                                  p <- ggplotly(p, height = 600, width = 900) %>% 
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
    
    
    #corr table
    gene.lipid.corr_data <- eventReactive(
        input$search_gene.lipid.corr, 
        ignoreNULL = FALSE, #ignoreNULL will initialze without event trigger input$search_gene.lipid.corr
        {
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
                              options = list(dom = 'Bfrti', 
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
    scat.lipid <- eventReactive(
        input$corr.table_rows_selected, {
        gene.lipid.corr_data()$Lipid[input$corr.table_rows_selected]
    })
    scat.gene <- eventReactive(input$corr.table_rows_selected, {
        gene.lipid.corr_data()$Gene[input$corr.table_rows_selected]
    })
    
    genelip.scat <- eventReactive(input$corr.table_rows_selected, {
        gs.p <- rna_livlipid_merge[c(scat.gene(), scat.lipid(), "Group","sampleID")] %>%
            ggplot(aes(x=!!as.symbol(scat.gene()), y=!!as.symbol(scat.lipid()),
                       color=Group, text = paste("ID ",sampleID)) ) +
            geom_point(size=2.5, na.rm=TRUE, alpha=0.7)+
            scale_color_viridis(discrete = T)+
            theme_classic()
        gs.p <- ggplotly(gs.p, height = 500, width = 700)
    })
    
    output$genelip.scatplot <- renderPlotly({
        genelip.scat()
    })


    
    #warning message if no gene match
    observeEvent(input$search_gene.lipid.corr, {
        if (!input$gene %in% genelist)
            showNotification(
                "No name matches!", 
                type = "warning", closeButton = T)
    })
    
    
    #tab2
    #rna-seq data
    rnadata <- eventReactive(input$search_gene.lipid.corr, ignoreNULL = FALSE,
        {
            validate(
                need(input$gene %in% genelist, 'No gene name fund.'),
                errorClass = "red"
            )
            
            rna_livlipid_merge[c("Group","sampleID",input$gene)]
    })
    
    output$rnaseq <- 
        DT::renderDT(server = FALSE, 
                     rnadata(),
                     selection=list(mode="single", target="cell"),
                     extensions="Buttons",
                     options = list(pageLength = 20, 
                                    dom = 'Blfrtip',
                                    buttons = list("copy", 
                                                list(extend="csv",
                                                     text = "Download csv", 
                                                     filename = "RNAseqData")
                                    )
                        )
        )
    
    #tab3
    #lipid data
    output$lipidtable <- DT::renderDT(server = FALSE, colnames = c(Species = 1),
                                 selection=list(mode="single", target="row"),
                                 setNames(data.frame(t(rna_livlipid_merge[c(22731:23306)])), rna_livlipid_merge[,2]),
                                 extensions= c("Buttons", "ColReorder"),
                                 options = list(pageLength = 5,
                                                dom = 'BRlfrtip',
                                                colReorder = list(realtime = FALSE),
                                                buttons = list("copy", 
                                                               list(extend="csv",
                                                                    text = "Download csv", 
                                                                    filename = "RNAseqData")
                                                              ),
                                                scrollX = TRUE
                                                )
    )
    
    #lipid plot
    lipidplotdata <- eventReactive(
        input$lipidtable_rows_selected, {
            ilipid <- names(rna_livlipid_merge)[input$lipidtable_rows_selected+22730]
            pLip <- rna_livlipid_merge %>% 
                select("Group", "sampleID", all_of(ilipid)) %>%
                ggplot(aes(x = Group, y = !!as.symbol(ilipid), color=Group,
                           text = paste("ID ",sampleID))) +
                geom_boxplot(outlier.shape = NA, color="black", outlier.color = NA,outlier.size = 0)+
                geom_jitter(size = 2, alpha = 0.7,
                            position=position_jitter(width=0.1, height=0,seed=123),
                            na.rm=TRUE)+
                scale_colour_viridis(discrete = T)+
                theme_classic()
            pLip <- ggplotly(pLip, height = 500, width = 900)
            pLip$x$data[1] <- lapply(pLip$x$data[1], 
                                  FUN = function(x){
                                      x$marker = list(opacity = 0)
                                      return(x)
                                  })
            pLip
        }
            
    )
    
    output$Lipid.boxplot <- renderPlotly({
        lipidplotdata()
    })
    
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
