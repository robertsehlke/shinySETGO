
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(shinythemes)
library(RMySQL)
library(DBI)
library(SETHRO)
library(Matrix)
library(GO.db)
library(cluster)
library(DT)
library(shinyBS)


# thanks to SBista (https://stackoverflow.com/questions/44319664/r-shiny-condition-a-tab-in-the-navbar-based-on-previous-tabs-condition)
jstabs = 
"shinyjs.disable_results =function(name){
$('ul li:has(a[data-value= \"Results\"])').addClass('disabled');
$('.nav li.disabled a').prop('disabled',true)
}
shinyjs.enable_results =function(name){
$('.nav li.disabled a').prop('disabled',false)
$('ul li:has(a[data-value= \"Results\"])').removeClass('disabled');
}"


############################## UI ############################## 

ui = navbarPage(title = "SETHRO", theme = shinytheme("flatly"), id="tabset",
                tabPanel("Analysis",
                         fluidPage( fluidRow(
                           column(1),
                           
                           column(4, wellPanel( style = "overflow-y:scroll; min-height: 900px; max-height: 900px",
                                                selectInput("species", label = "Select species", 
                                                            choices = list("Drosophila"="dro", 
                                                                           "C.elegans"="cel")),
                                                radioGroupButtons("evidence", label = "GO evidence codes", 
                                                                  choices = list("Curated and electronically inferred (IEA)"="iea",
                                                                                 "Curated only"="curated") ),
                                                checkboxGroupButtons("ontologies", label="GO ontologies", 
                                                                     choices = list("Cellular component"="CC",
                                                                                    "Biological process"="BP",
                                                                                    "Molecular Function"="MF"), 
                                                                     selected = c("CC","BP","MF")),
                                                HTML("<b>GO term size</b>"),
                                                HTML("</br>Between "),
                                                div(style="display: inline-block;vertical-align:middle", numericInput(inputId = "genes_min", value = 10, min = 1, max = 15000, label = "")),
                                                div(style="display: inline-block;vertical-align:middle", HTML(" and ")),
                                                div(style="display: inline-block;vertical-align:middle", numericInput(inputId = "genes_max", value = 300, min = 1, max = 15001, label = "")),
                                                div(style="display: inline-block;vertical-align:middle", HTML(" genes.")),
                                                
                                                HTML("</br><b>Genes of interest</b>"),
                                                HTML("</br>Require at least "),
                                                div(style="display: inline-block;vertical-align:middle", numericInput(inputId = "genes_vip", value = 5, min = 1, max = 15000, label = "")),
                                                div(style="display: inline-block;vertical-align:middle", HTML(" input genes in term."))
                           )),
                           
                           column(3, wellPanel( style = "overflow-y:scroll; min-height: 900px; max-height: 900px",
                                                fileInput("geneset_file", "Upload geneset (txt file, one gene per line):", multiple = FALSE, accept = c("text", "text/plain", ".txt")), 
                                                textAreaInput("geneset","OR - paste your geneset (one gene per line):", height = "260px"),
                                                fileInput("background_file", "Upload custom background set (txt file, one gene per line):", multiple = FALSE, accept = c("text", "text/plain", ".txt")),
                                                textAreaInput("background","OR - paste custom background (one gene per line):", height = "260px"),
                                                div( style="text-align: center;", actionLink("bt_example","[example data]", align="center") )
                             )
                           ),
                           
                           column(3, 
                                  wellPanel( style = "overflow-y:scroll; min-height: 900px; max-height: 900px",
                                             actionButton("load_data","Gather input data", style="width:100%"),
                                             htmlOutput("datastats"),
                                             tags$div(id="submit",
                                                      actionButton("submitdata", "Submit", icon("refresh"), style="width:100%") ) )
                           )
                           )
                         )
                ),
                
                tabPanel("Results", {
                  fluidPage( fluidRow(
                    column(1),
                    column(10, wellPanel( div(  dataTableOutput("results"), style = "font-size:80%" ),
                                          downloadButton("downloadResults", "Download results") ) )
                  ) )
                }),
                useShinyjs(),
                extendShinyjs(text = jstabs)
)




############################## SERVER ############################## 

server <- function(input, output,session) {
  output$datastats = renderText(HTML("</br>Upload the set of genes to analyse and click the button above to proceed!"))
  js$disable_results()
  shinyjs::hide("submit")
                                    
                                    
  
  dbcon = reactiveValues("c"=NULL)
  tmp = reactiveValues("queryresult"=NULL,
                       "idx"=NULL,
                       "id"=NULL,
                       "viable"=NULL,
                       "geneset"=NULL,
                       "res"=NULL,
                       "bg"=NULL,
                       "params"=NULL)
  
  
  
  #### > source functions ####
 source("functions.R")
  
  
  #### > Example data ####
  observeEvent(input$bt_example, {
    updateTextAreaInput(session, "geneset", value = paste0(read.delim("example_geneset.txt", header = F, stringsAsFactors = F)[,1], collapse = "\n") )
    updateTextAreaInput(session, "background", value = paste0(read.delim("example_background.txt", header = F, stringsAsFactors = F)[,1], collapse = "\n") )
  }, ignoreInit = T)
  
  
  #### > process file uploads (paste them into input fields) ####
  observeEvent(input$geneset_file, {
    this = read.delim(input$geneset_file$datapath, header = F, stringsAsFactors = F)[,1]
    updateTextAreaInput(session, "geneset", value = paste0(this,collapse = "\n") )
  })
  observeEvent(input$background_file, {
    this = read.delim(input$background_file$datapath, header = F, stringsAsFactors = F)[,1]
    updateTextAreaInput(session, "background", value = paste0(this,collapse = "\n") )
  })
  
  
  #### > gather data and reveal submit button ####
  observeEvent(input$load_data, {
    withProgress(message = 'Querying database...', value = 0, {
      incProgress(1/4)
      
      # construct database query
      dbquery = paste0("SELECT anno_",input$species,".goid, anno_",input$species,".entrez, genename, geneid FROM anno_", input$species, 
                       " JOIN gon_", input$species, " ON anno_", input$species, ".goid = gon_", input$species, ".goid",
                       " JOIN genes_", input$species, " ON anno_", input$species, ".entrez = genes_",input$species,".entrez",
                       " WHERE n BETWEEN ", input$genes_min, " AND ", input$genes_max,
                       " AND ontology IN ", paste0( "('", paste0( input$ontologies, collapse = "','"), "')" ) )
      
      # if curated only: exclude IEA evidence code
      if (input$evidence == "curated") { dbquery = paste0(dbquery, " AND evidence !='IEA'") }
      
      
      if (is.null(dbcon$c)) {
        output$datastats = renderText(NA)
        dbcon$c = dbConnect(RMySQL::MySQL(), default.file="~/.my.cnf", group="localSQL", user="shinylocal")
        incProgress(2/4)
      }
      
      if (!is.null(dbcon$c)) {
        # get query, then reset database connection
        tmp$queryresult = dbGetQuery(dbcon$c, dbquery)
        dbDisconnect(dbcon$c)
        dbcon$c = NULL
        incProgress(3/4)
        
        # select the geneset
        tmp$geneset = unlist( strsplit( input$geneset, split = "\n" ) )
        tmp$background = unlist( strsplit( input$background, split = "\n" ) )
        
        if( length(tmp$background) > 0 ) {
          tmp$background = union(tmp$geneset, tmp$background)    # by definition, the background needs to be a superset of the geneset of interest
        }
        
        # find out if gene set input is gene names, ids, entrez ...
        tmp$id = colnames(tmp$queryresult)[2:4][ which.max( apply(tmp$queryresult[,2:4], 2, function(x) sum(tmp$geneset %in% x) ) ) ]
        tmp$idx = which( tmp$queryresult[,tmp$id] %in% tmp$geneset )
        
        tmp$gomatrix = annotationMatrixFromLD( ldframe = tmp$queryresult[,c("goid",tmp$id)] )
        
        viable = sum( apply(tmp$gomatrix[,which(colnames(tmp$gomatrix) %in% tmp$geneset)],1,sum,na.rm=T) >= input$genes_vip )
        bg_go = sum( apply(tmp$gomatrix[,which(colnames(tmp$gomatrix) %in% tmp$background)],1,sum,na.rm=T) > 0 )
        bg_genes = length(intersect(tmp$background, colnames(tmp$gomatrix)))
        
        this = paste0("</br><b>",length(unique(tmp$queryresult$goid)), "</b> unique GO terms.</br>",
                      "<b>",nrow(tmp$queryresult), "</b> associations to genes.</br>",
                      "<b>",length(tmp$geneset), "</b> genes in input set.</br>",
                      if ( length(tmp$background) > 0 ) {
                        paste0( "<b>", bg_go,"</b> GO terms in custom background.</br>",
                                "<b>", bg_genes,"</b> annotated genes in custom background.</br>" )
                      } else { "" },
                      "<b>",length(unique(tmp$queryresult[tmp$idx,tmp$id])), "</b> genes annotated.</br>",
                      "<b>",viable, "</b> GO terms with more than ", input$genes_vip, " interesting genes.</br></br>")
        
        output$datastats =
          renderText( HTML( this ) )
      }
      
      if (length(tmp$geneset) > 0) { shinyjs::show(id="submit") }
      
    })
})
  
  #### SUBMIT: carry out analysis ####
  observeEvent(input$submitdata, {
    
    withProgress( message = "Calculating GO enrichments...", value=0, {
      incProgress(1/2)
      tmp$bg = if ( length(isolate(tmp$background)) > 0 ) { isolate(tmp$background) } else { colnames( isolate(tmp$gomatrix) ) }
      
      this_res = clusteredEnrichmentGO(
        tmp$gomatrix, tmp$geneset, 
        universe = tmp$bg, 
        min.genes = input$genes_vip)$results
      
      this_res = rs.signifDataframe(this_res)
      this_res$secondary_terms[this_res$Is.primary] = sapply( this_res$Primary[this_res$Is.primary], 
                                                            function(x) { paste0(this_res$Term[!this_res$Is.primary & this_res$Primary == x], collapse=", ") } )
      tmp$res = this_res
      js$enable_results()
      
      tmp$params = paste0("Species\t",isolate(input$species),"\n",
                          "Min_genes_annotated\t",isolate(input$genes_min),"\n",
                          "Max_genes_annotated\t",isolate(input$genes_max),"\n",
                          "Min_genes_significant\t",isolate(input$genes_vip),"\n",
                          "Ontologies\t",paste0( isolate(input$ontologies), collapse = ","),"\n",
                          "Evidence\t", switch(isolate(input$evidence), "iea"="curated and IEA", "curated"="curated only") )
    } )
    
    updateTabsetPanel(session, "tabset", selected="Results")
  })
  
  
  ##### OUPUT: Results table ####
  output$results = renderDataTable( cbind(' ' = '<i class=\"fa fa-plus\"></i>', tmp$res[tmp$res$Is.primary,]), 
                                    options = list( columnDefs=list( list("visible"=FALSE, targets=c(1,2,3,4,11,14,15,16)),
                                                                     list("searchable" = FALSE, className = 'details-control', targets = 0),
                                                                     list("width"="1%",targets=0),
                                                                     list("width"="8%",targets=c(6:10,12,13))),
                                                    pageLength = 20, lengthMenu = c(20, 50, 100) ),
                                    selection = 'single', style="bootstrap", rownames=FALSE, escape=0,
                                    callback = JS("
                                                  table.column(0).nodes().to$().css({cursor: 'pointer'});
                                                  var format = function(d) {
                                                  return '<div style=\"background-color:#eee; padding: .1em;\"><b>Sig. genes:</b> ' +
                                                  d[14] + '</br><b>Secondary terms:</b> ' + d[16] + '</div>';
                                                  };
                                                  table.on('click', 'td.details-control', function() {
                                                  var td = $(this), row = table.row(td.closest('tr'));
                                                  if (row.child.isShown()) {
                                                  row.child.hide();
                                                  td.html('<i class=\"fa fa-plus\"></i>');
                                                  } else {
                                                  row.child(format(row.data())).show();
                                                  td.html('<i class=\"fa fa-minus\"></i>');
                                                  }
                                                  });"))
  # ( adapted from DT example code (https://rstudio.github.io/DT/002-rowdetails.html) )
  
  
  #### Download handler ####
  output$downloadResults = downloadHandler(
    filename="results.zip",
    content=function(file) {
      write.table( tmp$res, file = "results.tsv", row.names = FALSE, sep = "\t")
      this_bg = data.frame( "Background_genes"=tmp$bg,
                            "In_annotated"=tmp$bg %in% colnames(tmp$gomatrix),
                            "In_geneset"=tmp$bg %in% tmp$geneset)
      write.table( this_bg, "background.tsv", sep="\t", row.names = FALSE )
      write(tmp$params, "parameters.tsv")
      
      zip(zipfile = file, files = c("results.tsv", "background.tsv","parameters.tsv"))
    }
  )
  
  
 
  
  
}


shinyApp(ui = ui, server = server)