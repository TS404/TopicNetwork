#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

shinyServer(function(input, output, session) {
#Setup --------
#> Packages  -------------

  library('WikidataR')
  library('WikidataQueryServiceR')
  library('tibble')
  library('devtools')
  library('stringr')
  library('dplyr')
  library('igraph')
  library('networkD3')
  library('dbscan')
  library('Hmisc')
  library('chorddiag')
  # library('lubridate')#
  # library('WikipediR')#
  # library('tidytext')#
  # library('htmltidy')#
  # library('readr')#
  # library('xml2')#
  
#> Functions ------
initials <- function(x,type="FLast"){
  if (type=="FLast"){
    gsub("^([A-Za-z]).* ([A-Za-z]*)", "\\1 \\2", x)
  }else{
    gsub("(.)\\S* *", "\\1", x)
  }
}
unspecial <- function(x){
  out <- x
  for(i in 1:ncol(x)){
    out[[i]] <- iconv(x[[i]],to = 'ASCII//TRANSLIT')
    if(Hmisc::all.is.numeric(x[[i]])){
      out[[i]] <- as.numeric(out[[i]])
    }else{
      out[[i]] <- as.factor(out[[i]])
    } 
  }
  return(as_tibble(out))
}
qid_from_name <- function(name   = 'Thomas Shafee',
                          limit  = 100,
                          unlist = TRUE){
  qid_from_name_nest1 <- function(x){lapply(x,"[[","id")}
  item.qs             <- lapply(name,find_item, limit=limit)
  item.qid            <- lapply(item.qs,qid_from_name_nest1)
  names(item.qid)     <- name
  if(unlist){item.qid <- unlist(item.qid)}
  return(item.qid)
}

colMax <- function(x) {if(is.null(dim(x))){x}else{apply(x, 2, max)}}

qdesc <- function(x){
  y   <- WikidataR::find_item(x)
  out <- paste0(y[[1]]$label,
                ' (<a target="_blank" href="https://wikidata.org/wiki/',
                y[[1]]$id,
                '">',
                y[[1]]$id,
                '</a>) - ',
                y[[1]]$description)
  out
}

qcite <- function(x){
  out <- paste0(1:nrow(x),
                ". ",
                x$authorLabel,', et al. "',
                '<a target="_blank" href="https://doi.org/',
                x$doi,
                '">',
                x$workLabel,'</a>" <i>',
                x$publisherLabel,
                "</i> (",x$year,") ",
                
                '<a href="https://doi.org/',
                x$doi,
                '">',
                x$doi,
                '</a> (<a target="_blank" href="',
                x$work,
                '">',
                gsub(".*/","",x$work),
                '</a>)',
                collapse = "</br>")
  out
}

# Calculations --------

input.button <- eventReactive(input$button.view,{
  nucleation.Qs <- qid_from_name(switch(input$rb,
                               example1 = c('Covid-19','SARS-CoV-2','2019–20 COVID-19 pandemic'),
                               example2 = c('antiviral drug','defensin', 'protease', 'cyclic peptide'),
                               custom   = unlist(strsplit(input$query,","))),
                        limit  = 1,
                        unlist = TRUE)
  
  sparql_query  <- paste0('SELECT ?node1 ?node2 ?node1Label ?node2Label ?pdate ?topicLabel ?authorLabel (COUNT(?work) AS ?count)
                         WHERE {
                           VALUES ?nucleation_topics { ',paste(paste0('wd:',nucleation.Qs),collapse = " "),' }
                           ?work wdt:P921 ?nucleation_topics, ?node1, ?node2. # Find co-nodes
                           FILTER   (?node1 != ?node2)                  # Exclude self-links
                         #  OPTIONAL { ?work wdt:P577 ?pdate. }
                         #  OPTIONAL { ?work wdt:P50  ?author. }
                           SERVICE  wikibase:label { bd:serviceParam wikibase:language "en,fr,de,ru,es,zh,jp". }
                         }
                         GROUP BY ?node1 ?node2 ?node1Label ?node2Label ?pdate ?topicLabel ?authorLabel
                         ORDER BY DESC(?count)
                         LIMIT 10000')
  links.qr   <- suppressMessages(unspecial(query_wikidata(sparql_query)))
  links.qr$node1Label <- ordered(links.qr$node1Label,         # define factor order by most linked
                                 unique(names(sort(table(links.qr$node1Label),TRUE))))
  links.qr   <- links.qr[with(links.qr, order(node1Label)),] # order all by that factor order
  links.qr.r <- t(apply(links.qr[,3:4], 1, sort))                 # find reciprocal links and flip them
  links.qr   <- links.qr[!duplicated(links.qr.r),]             # remove flipped reciprocal links (now show up as duplicates)

  #> igraph --------------
  links.qr.o <- links.qr[order(links.qr$count),]
  g    <- graph_from_edgelist(as.matrix(links.qr.o[,3:4]),
                              directed = 0) %>%
          set_edge_attr("weight", value = links.qr.o$count)
  clp  <- cluster_label_prop(g,weights = E(g)$weight)
  hcl  <- hclust(dist(g[]))
  hdb  <- dbscan::hdbscan(dist(g[]),minPts = 3)
  l    <- layout.graphopt(g,spring.length = E(g)$weight)
  
  groups           <- hdb$cluster
  groups           <- clp$membership
  max.groups       <- max(groups)
  colours          <- RColorBrewer::brewer.pal(max.groups,'Dark2')
  colours          <- viridis::plasma(max.groups)
  colours          <- colorRampPalette(c('#ee4400','#ff0099','#440055','#000055','#006688','#005533'))(max.groups)
  if(min(groups)==0){colours[1]<-"white"}
  barplot(rep(1,length(colours)),col=colours,space=0,border=if(length(colours)>20){NA},main="group colour palette used")
  group.order      <- order(groups,hdb$cluster)
  group.colours    <- colours[groups]
  group.colours[1:length(nucleation.Qs)] <- "black"
  group.colours    <- group.colours[group.order]
  data.grouped     <- log10(1+100*as.matrix(get.adjacency(g,attr = "weight",sparse = T))[group.order,group.order]/max(g[]))
  
  #> D3 net ------------------
  
  links.D3 <- data.frame(node1Q     = links.qr$node1,
                         node2Q     = links.qr$node2,
                         node1Label = links.qr$node1Label,
                         node2Label = links.qr$node2Label,
                         count      = links.qr$count)
  
  node_names <- iconv(clp$names,to = 'ASCII//TRANSLIT')
  
  links.D3$node1Label <- ordered(iconv(links.D3$node1Label,to='ASCII//TRANSLIT'),node_names)
  links.D3$node2Label <- ordered(iconv(links.D3$node2Label,to='ASCII//TRANSLIT'),node_names)
  
  links.D3$node1Label <- as.numeric(links.D3$node1Label)-1 # zero-indexed
  links.D3$node2Label <- as.numeric(links.D3$node2Label)-1 # zero-indexed
  links.D3$count      <- as.numeric(links.D3$count)
  links.D3$colour     <- colorRampPalette(c(rgb(0,0,0,0.3),rgb(0,0,0,0.9)), alpha=TRUE)(max(E(g)$weight))[E(g)$weight]
  
  # node_sizes <- c(max(links.D3$count),links.D3$count[links.D3$node1Label==0|1|2][order(links.D3$node2Label[links.D3$node1Label==0])]) #old version
  if(length(nucleation.Qs)==1){
    node_sizes <- g[unique(links.qr$node1Label[grep(paste(nucleation.Qs[1]),links.qr$node1)]),][node_names] 
  }else{
    node_sizes <- colMax(as.matrix(g[unique(links.qr$node1Label[grep(paste(nucleation.Qs,collapse="|"),links.qr$node1)]),]))[node_names]
  }
  node_sizes[1:length(nucleation.Qs)] <- max(node_sizes)
  
  normalise <- 1
  if(normalise){
    links.D3$count     <- 1+9*sqrt((links.qr$count-1)/max(links.qr$count-1))
    node_sizes         <- 1+3*(log10(1+1000*(node_sizes-1)/max(node_sizes-1)))
  }
  
  nodes.D3 <- data.frame(name  = node_names,
                         group = groups,
                         size  = node_sizes)
  nodes.D3$group[1:length(nucleation.Qs)] <- 0

  list(nucleation.Qs = nucleation.Qs,
       sparql_query  = sparql_query,
       links.qr      = links.qr,
       g             = g,
       links.D3      = data.frame(links.D3),
       nodes.D3      = data.frame(nodes.D3),
       data.grouped  = data.grouped,
       group.colours = group.colours
       )
})



# Visualisations ---------
#> renderForceNetwork ---------
output$force <- renderForceNetwork({
  
  forceNetwork(Links        = input.button()$links.D3,
               Nodes        = input.button()$nodes.D3,
               linkColour   = input.button()$links.D3$colour,
               charge       = input$charge*-3,
               Source       = "node1Label",
               Target       = "node2Label",
               Value        = "count",
               NodeID       = "name",
               Group        = "group", 
               Nodesize     = "size",
               zoom         = 1,
               opacity      = 0.8,
               fontFamily   = "sans-serif",
               fontSize     = 20,
               height       = 5000,
               clickAction  = 'Shiny.onInputChange("nodeid", d.name)',
               linkWidth    = JS('function(d){return 1+(d.value)}'),
               linkDistance = JS('function(d){return 2/(1+d.value)}'),
               radiusCalculation = JS("(d.nodesize)+3")
  )
})

#> renderChordDiag ---------
output$chord <- renderChorddiag({
  chorddiag(data              = as.matrix(input.button()$data.grouped),
            groupColors       = as.matrix(input.button()$group.colours),
            groupPadding      = if(ncol(input.button()$g[])<100){1}else{0},
            showTicks         = F,
            margin            = 80,
            groupnameFontsize = 12,
            groupnamePadding  = 10,
            chordedgeColor    = '#00000022',
            groupedgeColor    = '#00000099')
})

# Text -------------


observe({
  input.button()$nodes.D3$name <- input.button()$nodes.D3$name
  choices.df <- data.frame(names  = sort(input.button()$nodes.D3$name),
                         toggle = rep(FALSE,length(input.button()$nodes.D3$name)))
  
  choices.df$toggle[which(choices.df$names==input$nodeid)] <- !choices.df$toggle[which(choices.df$names==input$nodeid)]
  
  updateCheckboxGroupInput(session, "inCheckboxGroup2",
                           label    = NULL,
                           choices  = choices.df$names,
                           selected = choices.df$toggle
  )
  output$nodes.sel <- renderUI({ HTML(paste(selected)) })
})


output$qids  <- renderUI({ HTML(paste(lapply(input.button()$nucleation.Qs,qdesc),collapse = "<br/>")) })

output$cites <- renderUI({
  subset.Qs <- qid_from_name(input$inCheckboxGroup2,limit=1)
  subset.sparql <- paste0('SELECT ?work ?workLabel ?authorLabel ?publisherLabel ?year ?pdate ?doi
                         WHERE { ?work ',
                          paste0("wdt:P921? wd:",subset.Qs,";",collapse=""),
                          '(p:P50|p:P2093) [ (ps:P50|ps:P2093) ?author; pq:P1545 ?ordinal ].  
                           OPTIONAL { ?work wdt:P1433 ?publisher. }
                           OPTIONAL { ?work wdt:P577 ?pdate. }
                           OPTIONAL { ?work wdt:P356 ?doi. }
                           bind(str(year(?pdate)) as ?year)
                           FILTER (xsd:integer(?ordinal) = 1)
                           SERVICE  wikibase:label { bd:serviceParam wikibase:language "en,fr,de,ru,es,zh,jp". }
                         }
                         LIMIT 10000')
  subset.qr   <- suppressMessages(unspecial(query_wikidata(subset.sparql)))
  if(length(subset.Qs)>0){
    if(nrow(subset.qr)>0){
      subset.cite <- paste0("<b>",
                            nrow(subset.qr),
                            " publication",
                            if(nrow(subset.qr)>1){"s"},
                            " found on that ",
                            if(length(subset.Qs)>1){"combination of topics"}else{"topic"},
                            "</b></br>",
                            qcite(subset.qr))
      HTML(subset.cite)
    }else{
      HTML("No publications found on that combination of topics")
    }    
  }else{
    HTML("Select topics from list on the left")
  }


})




}) # end of shinyServer






