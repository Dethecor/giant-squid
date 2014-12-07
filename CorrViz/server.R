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
  output$CorViz <- renderPlot({
    genus <- fetchAbundance(cohortFile, studyID, type = "Relative", entity = "Genera", index = list(input$GenusSelect,NULL))
    gene  <- fetchAbundance(cohortFile, studyID, type = "Relative", entity = "Genes", index = list(input$GenusSelect,NULL))
    stopifnot(colnames(genus) == colnames(gene))
    plotData <- data.frame("Genus" = genus[1,], "Gene" = gene[1,])
    p <- ggplot(plotData, aes(x=Genus,y=Gene)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle(paste("Pearson Correlation", cor(plotData$Genus, plotData$Gene), sep = " - "))
    print(p)
  })
  output$CorBars <- renderPlot({
    stopifnot(colnames(genus) == colnames(gene))
    plotData <- data.frame("Genus" = genus[1,], "Gene" = gene[1,])
    p <- ggplot(plotData, aes(x=Genus,y=Gene)) + geom_point() + scale_x_log10() + scale_y_log10() + ggtitle(paste("Pearson Correlation", cor(plotData$Genus, plotData$Gene), sep = " - "))
    print(p)
  })
})