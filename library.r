#metaHDF5
require(rhdf5)
addTaxon <- function( file, studyID, taxonAnnotation, taxon = "Species", taxonChunk = 10e3){
  nTaxa     <- nrow(taxonAnnotation)
  for(aType in c("Raw", "Normalised", "Relative")){
    group <- paste( studyID, "Abundance", aType, sep = "/" )
   h5createDataset(file, dataset = paste(group, taxon, sep = "/"), dims = c(nTaxa, nSamples), storage.mode = "double", chunk = c(min(taxonChunk, nTaxa),min(sampleChunk, nSamples)), level = compressionLevel)
  }
  h5createDataset(file, dataset = paste( studyID, "Annotation", taxon, sep = "/"), dims = c(nTaxa, ncol(taxonAnnotation)), storage.mode = "character", size = 100, level = compressionLevel)
  h5write(as.matrix(taxonAnnotation), file, name = paste( studyID, "Annotation", taxon, sep = "/"))
  f <- H5Fopen(file)
  g <- H5Gopen(f, paste( studyID, "Annotation", sep = "/"))
  d <- H5Dopen(g, name = taxon)
  h5writeAttribute.character(colnames(taxonAnnotation), d, "ColNames")
  h5writeAttribute.character(sapply(taxonAnnotation[1,], class), d, "Types")
  H5Dclose(d)
  H5Gclose(g)
  H5Fclose(f)
}

createHDF5CohortFile <- function(
  file = tempfile(), studyID = "/MetaGenomics",
  sampleAnnotation, geneAnnotation, taxonAnnotation, taxon = "Genera",
  compressionLevel = 7, geneChunk = 10e3,
  sampleChunk = 100, taxonChunk = 10e3
){
  nSamples  <- nrow(sampleAnnotation)
  nGenes    <- nrow(geneAnnotation)
  nTaxa     <- nrow(taxonAnnotation)
  h5createFile(file)
  h5createGroup(file, studyID)
  h5createGroup(file, paste( studyID, "Abundance", sep = "/" ) )
  H5close()
  for(aType in c("Raw", "Normalised", "Relative")){
    group <- paste( studyID, "Abundance", aType, sep = "/" )
    h5createGroup(file, group )
    h5createDataset(file, dataset = paste(group, "Genes", sep = "/"), dims = c(nGenes, nSamples), storage.mode = "double", chunk = c(min(geneChunk, nGenes),min(sampleChunk,nSamples)), level = compressionLevel)
    h5createDataset(file, dataset = paste(group, taxon, sep = "/"), dims = c(nTaxa, nSamples), storage.mode = "double", chunk = c(min(taxonChunk, nTaxa),min(sampleChunk, nSamples)), level = compressionLevel)
  }
  h5createGroup(file, paste(studyID, "Annotation", sep = "/"))
  h5createDataset(file, dataset = paste( studyID, "Annotation", "Samples", sep = "/"), dims = c(nSamples, ncol(sampleAnnotation)), storage.mode = "character", size = 100, level = compressionLevel)
  h5write(as.matrix(sampleAnnotation), file, name = paste( studyID, "Annotation", "Samples", sep = "/"))
  h5createDataset(file, dataset = paste( studyID, "Annotation", "Genes", sep = "/"), dims = c(nGenes, ncol(geneAnnotation)), storage.mode = "character", size = 100, level = compressionLevel, chunk = c(geneChunk, ncol(geneAnnotation)))
  h5write(as.matrix(geneAnnotation), file, name = paste( studyID, "Annotation", "Genes", sep = "/"))
  h5createDataset(file, dataset = paste( studyID, "Annotation", taxon, sep = "/"), dims = c(nTaxa, ncol(taxonAnnotation)), storage.mode = "character", size = 100, level = compressionLevel)
  h5write(as.matrix(taxonAnnotation), file, name = paste( studyID, "Annotation", taxon, sep = "/"))
  H5close()
  f <- H5Fopen(file)
  g <- H5Gopen(f, paste( studyID, "Annotation", sep = "/"))
  d <- H5Dopen(g, name = "Samples")
  h5writeAttribute.character(colnames(sampleAnnotation), d, "ColNames")
  h5writeAttribute.character(sapply(sampleAnnotation[1,], class), d, "Types")
  H5Dclose(d)
  d <- H5Dopen(g, name = "Genes")
  h5writeAttribute.character(colnames(geneAnnotation), d, "ColNames")
  h5writeAttribute.character(sapply(geneAnnotation[1,], class), d, "Types")
  H5Dclose(d)
  d <- H5Dopen(g, name = taxon)
  h5writeAttribute.character(colnames(taxonAnnotation), d, "ColNames")
  h5writeAttribute.character(sapply(taxonAnnotation[1,], class), d, "Types")
  H5Dclose(d)
  H5Gclose(g)
  H5Fclose(f)
}

prepareForCorrelationMatrix <- function(
  file, studyID, rows = "Samples", columns = "Taxa",
  rowChunk = 100, colChunk = 10e3, compressionLevel = 7
){
  nRows    <- nEntities(file, studyID, rows)
  nColumns <- nEntities(file, studyID, columns)
  if(!(paste("/", studyID, "Correlation", sep = "/") %in% h5ls("Data/cohort.h5")$group)){
    h5createGroup(file, paste(studyID, "Correlation", sep = "/"))
  }
  h5createDataset(file, dataset = paste(studyID, "Correlation", paste(rows, columns, sep = "_"), sep = "/"), dims = c(nRows, nColumns), storage.mode = "double", chunk = c(min(rowChunk, nRows),min(colChunk,nColumns)), level = compressionLevel)
}

correlateEntities <- function(file, studyID, type = "Relative", rows, columns, nRows = 1000, nCols = 1000, method = "pearson", samples = NULL, verbose = FALSE){
  rowLoc <- paste(studyID, type, rows)
  nRow   <- nEntities(file, studyID, entity = rows)
  colLoc <- paste(studyID, type, columns)
  nCol   <- nEntities(file, studyID, entity = columns)
  corDest <- paste(studyID, "Correlation", paste(rows, columns, sep = "_"), sep = "/")
  stopifnot(fetchDimensions(file,location = corDest) == c(nRow, nCol)) #will the output fit?
  f <- H5Fopen(file)
  g <- H5Gopen(f, paste( studyID, "Correlation", sep = "/"))
  d <- H5Dopen(g, name = paste(rows, columns, sep = "_"))
  h5writeAttribute.character(method, d, "Method")
  H5Dclose(d)
  H5Gclose(g)
  H5Fclose(f)
  for(startCol in seq(1,nCol,nCols)){
    colAbundance <- fetchAbundance(file, studyID, type = type, entity = columns, index = list(startCol:min(startCol + nCols - 1, nCol), samples))
    for(startRow in seq(1,nRow,nRows)){
      rowAbundance <- fetchAbundance(file, studyID, type = type, entity = rows, index = list(startRow:min(startRow + nRows - 1, nRow), samples))
      h5write(obj = cor(t(rowAbundance), t(colAbundance), method = method), file = file, name = corDest, index = list(startRow:min(startRow+nRows-1,nRow), startCol:min(startCol+nCols-1, nCol)))
      if(verbose){
        print(paste("Correlation of block ", startRow, ":", min(startRow+nRows-1,nRow), " x ", startCol, ":", min(startCol+nCols-1,nCol), " written!", sep=""))
      }
    }
  }
}

fetchCorrelations <- function( file, studyID, rows, columns, index = list(NULL, NULL)){
  corLoc <- paste(studyID, "Correlation", paste(rows, columns, sep = "_"), sep = "/")
  ret <- h5read(file, corLoc, index)
  if(nrow(ret) <= 1000){
    rownames(ret) <- fetchNames(file = file, studyID = studyID, entity = rows, index = index[[1]])
  }
  if(ncol(ret) <= 1000){
    colnames(ret) <- fetchNames(file = file, studyID = studyID, entity = columns, index = index[[2]])
  }
  return(ret)
}

fetchAnnotation <- function( file, studyID, entity = "Samples", index = list(NULL, NULL)){
  ret <- h5read(file, paste(studyID, "Annotation", entity, sep = "/"), index)
  attr <- h5readAttributes(file, paste(studyID, "Annotation", entity, sep = "/"))
  colnames(ret) <- attr$ColNames[ifelse(is.null(index[[2]]), TRUE, index[[2]])]
  ret <- as.data.frame(ret)
  for(i in seq_along(attr$Types)){
    ret[[i]] <- as(object = ret[[i]], Class = attr$Types[i])
  }
  return(ret)
}

fetchNames <- function( file, studyID, entity = "Samples", index = NULL ){
  attr <- h5readAttributes(file, paste(studyID, "Annotation", entity, sep = "/"))
  if("Name" %in% attr$ColNames){
    nameCol = which(attr$ColNames == "Name")
    ret <- h5read(file, paste(studyID, "Annotation", entity, sep = "/"), list(index, nameCol))[,1]
    return(ret)
  }else{
    stop(paste("No 'Name' annotation found at location:", paste(studyID, "Annotation", entity, sep = "/"), "in file:", file))
  }
}

writeAbundance <- function( abundanceFile, destination, studyID, type = "Relative", entity = "Genes", chunkSize = 20e3, verbose = TRUE, skipBlocks = 0 ){
  f <- file(abundanceFile, open = "r")
  sampleIDs <- scan(f, nlines = 1, what = "")
  sampleAnnotation <- fetchAnnotation(destination, studyID, "Samples")
  nSamples <- length(sampleIDs)
  stopifnot(nSamples == nrow(sampleAnnotation))
  nEntities <- nEntities(destination, studyID, entity)
  idx <- 1 + skipBlocks
  if(skipBlocks > 0){
    if(verbose){
      print(paste("Skipping", skipBlocks, "blocks!"))
    }
    tmp <- scan(f, skip = skipBlocks * chunkSize - 1, nlines = 1) #reposition the connection ... wish I could use seek here
  }
  theWhat <- c(list(""), as.list(as.double(1:nSamples)))
  for(i in (1+skipBlocks):max(ceiling(nEntities/chunkSize), 1)){
    data <- scan(f, nlines = chunkSize, what = theWhat)
    data <- do.call(cbind, data[-1])
    h5write(data, destination, paste(studyID, "Abundance", type, entity, sep = "/"), index = list(idx:(idx+nrow(data)-1), NULL))
    if(verbose){
      print(paste("[", Sys.time(), "] Block #", i, " written at index ", idx, ":", idx+nrow(data)-1, sep = ""))
    }
    idx <- idx + nrow(data)
  }
  close(con = f)
}

fetchAbundance <- function(
  file, studyID, type = "Relative",
  entity = "Taxa", index = list(NULL, NULL),
  namecols = 10e3, namerows = 10e3
){
  ret <- h5read(file, paste(studyID, "Abundance", type, entity, sep = "/"), index = index)
  if(nrow(ret) <= namerows){
    rownames(ret) <- fetchNames(file = file, studyID = studyID, entity = entity, index = index[[1]])
  }
  if(ncol(ret) <= namecols){
    colnames(ret) <- fetchNames(file = file, studyID = studyID, entity = "Samples",  index = index[[2]])
  }
  return(ret)
}

fetchDimensions <- function(file, location){
  f = H5Fopen(file)
  d = H5Dopen(f, location)
  s = H5Dget_space(d)
  ret = H5Sget_simple_extent_dims(s)$size
  H5Sclose(s)
  H5Dclose(d)
  H5Fclose(f)
  return(ret)
}

nGenes <- function(file, studyID){
  nEntities(file, studyID, "Genes")
}

nSamples <- function(file, studyID){
  nEntities(file, studyID, "Samples")
}

nTaxa <- function(file, studyID){
  nEntities(file, studyID, "Taxa")
}

nEntities <- function(file, studyID, entity = "Samples"){
  return( fetchDimensions(file, location = paste(studyID, "Annotation", entity, sep = "/"))[1] )
}

happly <- function(file, location, blockSize, dim, FUN, ...){
  targetDims <- fetchDimensions(file = file, location = location)
  if(dim > length(targetDims)){
    stop(paste("'dim' argument must map to one of the", length(targetDims), "dimensions of the target dataset."))
  }
  blockStarts <- seq(1, targetDims[dim], blockSize)
  ret <- lapply(blockStarts, function(blockStart){
    idx <- vector(mode="list", length=length(targetDims))
    blockEnd <- min(blockStart + blockSize - 1, targetDims[dim])
    idx[[dim]] <- blockStart:blockEnd
    block <- h5read(file = file, name = location, index = idx)
    FUN(block, index = idx, ...)
  })
  return(ret)
}

#apply a function for all combinations of subject and query
blockApply <- function(file, studyID, type = "Relative", subjectEntity = "Genes", queryEntity = "Taxa"){
  query = fetchAbundance(file, studyID, type, queryEntity)
  nSubjects = nEntities(file, studyID, subjectEntity)
}