shinyServer(function(input, output) {
  output$CohortInfo <- renderTable({
    ret <- data.frame(
    "Value" = c(
      cohortFile,
      studyID,
      nSamples(cohortFile, studyID),
      nGenes(cohortFile, studyID),
      nEntities(file = cohortFile, studyID = studyID, entity = "Genera")
      )
    )
    rownames(ret) <- c("Active File", "Study ID", "#Samples", "#Genes", "#Genera")
    ret
    })
  output$AbundanceHeat <- renderPlot({
    idx <- list(
      seq(input$GeneSelect[1], input$GeneSelect[2]),
      seq(input$SampleSelect[1], input$SampleSelect[2])
      )
    heatmap(fetchAbundance(cohortFile, studyID, type = "Relative", entity = "Genes", index = idx))
  })
})
