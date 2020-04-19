#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(networkD3)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Co-topic network from Wikidata"),
    wellPanel(
        div(class="header",
            p("This tool calculates a network of topics.
               Each 'node' is a topic, and 'edges' are shown whenever a scientific publication is about both topics."),
            p("The data is drawn from Wikidata by searching outwards from the starting topic(s).
               If clear clusters can be identified, they'll be assigned group colours."),
            p("The more publications that are about a pair of topics, the thicker the link between them.
               Node size indicates topics that commly co-occur in publications  with the starting topic(s).")
            )
    ),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            textAreaInput("query",resize="horizontal",rows=1,
                          label       = "Starting topic",
                          value       = "COVID-19",
                          placeholder = "Search term to start the graph"),
            radioButtons("rb", "Choose one:",
                         choiceNames = list(
                             "custom",
                             "example 1 (Covid-19, SARS-CoV-2, 2019–20 COVID-19 pandemic)",
                             "example 2 (Protein, defensin, protease, cyclic peptide)"
                         ),
                         choiceValues = list(
                             "custom", "example1", "example2"
                         )),
            actionButton("button.view", "Calculate network"),
            sliderInput("charge",
                        "node separation",
                        min = 1,
                        max = 50,
                        value = 30)
        ),
        # Show the network
        mainPanel(
            forceNetworkOutput("force"),
            renderText("query")
        )
    )
))
