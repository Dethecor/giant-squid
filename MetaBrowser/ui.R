# ui.R

shinyUI(fluidPage(
  titlePanel("Metagenomic Data Browser - powered by HDF5"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("SampleSelect", label = "Samples", min = 1, max = floor(nSamples(cohortFile, studyID)), value = c(1,nSamples(cohortFile, studyID))),
      sliderInput("GeneSelect", label = "Genes", min = 1, max = 1000, value = c(1,100), step = 10),
      tableOutput("CohortInfo")
      ),
    mainPanel(
      "Abundance Overview Plot",
      plotOutput("AbundanceHeat")
      )
  )
))
