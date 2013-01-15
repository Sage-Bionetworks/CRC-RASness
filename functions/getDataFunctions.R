## FUNCTIONS TO EXTRACT RAW DATA OBJECTS FROM SYNAPSE
#####
## ANALYST: BRIAN M. BOT
#####


#####
## KFSYSCC DATA
#####
getKFSYSCCdata <- function(){
  
  require(synapseClient)
  require(affy)
  
  ## GRAB THE ARCHIVE FILE FROM KOO - UNTAR IN TEMP DIRECTORY
  kfEnt <- downloadEntity("syn1528362")
  kfDir <- tempfile(pattern="kfDir")
  dir.create(kfDir)
  untar(file.path(kfEnt$cacheDir, kfEnt$files), exdir = kfDir)
  theseFiles <- list.celfiles(kfDir, recursive=T, full.names=T)
  
  ## READ IN CEL FILES
  kooExpr <- ReadAffy(filenames=theseFiles)
  
  ## READ IN THE CLINICAL DATA
  clinEnt <- downloadEntity("syn1588162")
  kooClin <- read.csv(file.path(clinEnt$cacheDir, clinEnt$files), as.is=T)
  
  stopifnot(all(sapply(strsplit(sampleNames(kooExpr), ".", fixed=T), "[[", 1) == kooClin$SN))
  rownames(kooClin) <- sampleNames(kooExpr)
  pData(kooExpr) <- kooClin
  
  return(kooExpr)
}


#####
## TCGA DATA -- COAD AND READ
#####
loadTCGAFileFromEntity <- function(synId){
  require(synapseClient)
  
  ent <- downloadEntity(synId)
  df <- read.delim(file.path(ent$cacheDir, ent$files), header=F, as.is=T)
  colnames(df) <- as.character(df[1, ])
  df <- df[-1, ]
  rownames(df) <- as.character(df[, 1])
  df <- df[, -1]
  return(df)
}

#####
## AGILENT
getTCGAcrcAgilent <- function(){
  coadAgilent <- loadTCGAFileFromEntity("syn417828")
  readAgilent <- loadTCGAFileFromEntity("syn418082")
  
  if( all(rownames(coadAgilent) == rownames(readAgilent)) ){
    theseFeatures <- rownames(coadAgilent)
    crcAgilent <- cbind(coadAgilent, readAgilent)
  } else{
    stop("rownames do not match")
  }
  
  thesePatients <- sapply(strsplit(colnames(crcAgilent), "-", fixed=T), function(x){
    paste(x[1:3], collapse="-")
  })
  
  if( all(duplicated(thesePatients) == FALSE) ){
    colnames(crcAgilent) <- thesePatients
  } else{
    stop("duplicated patients")
  }
  
  ## CONVERT TO NUMERIC MATRIX
  crcAgilent <- apply(crcAgilent, 2, as.numeric)
  rownames(crcAgilent) <- theseFeatures
  
  return(crcAgilent)
}

#####
## RNAseq
getTCGAcrcRNAseq <- function(){
  coadRNAseq <- loadTCGAFileFromEntity("syn417839")
  readRNAseq <- loadTCGAFileFromEntity("syn418090")
  
  if( all(rownames(coadRNAseq) == rownames(readRNAseq)) ){
    theseFeatures <- rownames(coadRNAseq)
    crcRNAseq <- cbind(coadRNAseq, readRNAseq)
  } else{
    stop("rownames do not match")
  }
  
  thesePatients <- sapply(strsplit(colnames(crcRNAseq), "-", fixed=T), function(x){
    paste(x[1:3], collapse="-")
  })
  
  if( all(duplicated(thesePatients) == FALSE) ){
    colnames(crcRNAseq) <- thesePatients
  } else{
    stop("duplicated patients")
  }
  
  ## CONVERT TO NUMERIC MATRIX
  crcRNAseq <- apply(crcRNAseq, 2, as.numeric)
  rownames(crcRNAseq) <- theseFeatures
  
  return(crcRNAseq)
}



#####
## CLINICAL DATA AT THE PATIENT LEVEL
getTCGAcrcClinical <- function(){
  coadClin <- loadTCGAFileFromEntity("syn1446080")
  readClin <- loadTCGAFileFromEntity("syn1446153")
  coadClin$rns <- rownames(coadClin)
  readClin$rns <- rownames(readClin)
  
  crcClin <- merge(x=coadClin, y=readClin, all=T)
  rownames(crcClin) <- crcClin$rns
  crcClin$rns <- NULL
  
  return(crcClin)
}




#####
## GAEDCKE DATASET FROM GEO
## Agilent-014850 Whole Human Genome Microarray 4x44K G4112F 
#####
getGaedckeFromGEO <- function(){
  require(GEOquery)
  
  geoFiles <- getGEO("GSE20842", GSEMatrix=T, AnnotGPL=T)
  gaedckeEset <- geoFiles$GSE20842_series_matrix.txt.gz
  # NOTE: KRAS STATUS - pData(gaedckeEset)$characteristics_ch1.5 IN pData
  
  ## SUBSET TO ONLY TUMORS
  clin <- pData(gaedckeEset)
  gaedckeEset <- gaedckeEset[, clin$characteristics_ch1.4 == "tissue: tumor" ]
  
  return(gaedckeEset)
}


#####
## KHAMBATA-FORD DATASET FROM GEO
## 
#####
getKhambataFromGEO <- function(){
  require(GEOquery)
  
  geoFiles <- getGEO("GSE5851", GSEMatrix=T, AnnotGPL=T)
  khambataEset <- geoFiles$GSE5851_series_matrix.txt.gz
  
  return(khambataEset)
}


#####
## EMEXP3557
#####
## NEED INFORMATION ON THIS CHIP -- ADXCRCG2a520319
## NEED CLINICAL INFORMATION
#####
getEMEXP3557 <- function(){
  require(synapseClient)
  require(affy)
  
  ent <- downloadEntity("syn372544")
  exprSet <- ReadAffy(celfile.path=ent$cacheDir)
  
}


#####
## EMEXP991
#####
## NEED CLINICAL INFORMATION
#####
getEMEXP991 <- function(){
  require(synapseClient)
  
  ent <- downloadEntity("syn202357")
  exprSet <- ReadAffy(celfile.path=ent$cacheDir)
  ## exprSet@protocolData@data$ScanDate
  ## ALL SCANED ON SAME DAY
  
}



