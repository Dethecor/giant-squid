# ui.R

shinyUI(fluidPage(
  titlePanel("Metagenomic Data Browser - powered by HDF5"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("GenusSelect", label = "Genus", min = 1, max = floor(nEntities(cohortFile, studyID, entity = "Genera")), value = 1),
      sliderInput("GeneSelect", label = "Gene", min = 1, max = 1000, value = 1, step = 10),
      tableOutput("CohortInfo")
      ),
    mainPanel(
      "Correlation Overview Plot",
      plotOutput("CorBars"),
      plotOutput("CorViz")
      )
  )
))
