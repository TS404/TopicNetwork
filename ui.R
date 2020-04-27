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
library(chorddiag)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

    tags$head(
        tags$style(HTML("hr {border-top: 1px solid #999;}"))
    ),
    # Application title
    titlePanel("Co-topic network from Wikidata"),
    wellPanel(
        div(class="header",
            p("This tool calculates a network of topics.
               Each 'node' is a topic, and 'edges' are shown whenever a scientific publication is about both topics.
               The data is drawn from Wikidata by searching outwards from the starting topic(s).
               Searches can take >30 seconds.
               If clear clusters can be identified, they'll be assigned group colours.
               The more publications that are about a pair of topics, the thicker the link between them.
               Node size indicates topics that commly co-occur in publications  with the starting topic(s)."
              )
            )
    ),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            h4("Starting topic for network"),
            conditionalPanel(
                condition = "input.rb == 'custom'",
                textInput("query",
                          label       = NULL,
                          value       = "COVID-19",
                          placeholder = "Search term to start the graph"),
            ),
            htmlOutput("qids"),
            hr(),
            radioButtons("rb",
                         label = NULL,
                         choiceNames = list(
                             "custom",
                             "example 1 (Covid-19, SARS-CoV-2, 2019–20 COVID-19 pandemic)",
                             "example 2 (antiviral drug, defensin, protease, cyclic peptide)"
                         ),
                         choiceValues = list(
                             "custom",
                             "example1",
                             "example2"
                         )),
            actionButton("button.view", "Calculate network"),
            hr(),
            sliderInput("charge",
                        "node separation",
                        min = 1,
                        max = 50,
                        value = 30),
            hr(),
            helpText(p("Code for this page can be found at",
                       a(href="https://github.com/TS404/TopicNetwork","github.com/TS404/TopicNetwork"),
            )),
        ),
        # Show the network
        mainPanel(
            forceNetworkOutput("force", height = "800px"),
            # chorddiagOutput("chord", height = "800px")
        )
    )
))
