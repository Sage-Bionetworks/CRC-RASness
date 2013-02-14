## A COLLECTION OF FUNCTIONS USED IN MODELING
#####
## DEPENDS ON THE utilityFunctions.R PROGRAM BEING SOURCED
#####
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####


## PREDICTION USING ELASTIC NET -- TRAINING AND TEST ExpressionSets PASSED
#####
## FIRST DETERMINES OVERLAP IN FEATURE SETS SO MODELS
## FIT ON THE SAME TRAINING SET BUT WITH DIFFERENT
## TESTING SETS MAY PRODUCE DIFFERENT MODELS
#####
binomialPredictEN <- function(trainEset, trainResponse, testEsets, alpha=.1, seed=2012, quantileNormalize=F){
  require(glmnet)
  require(caret)
  require(Biobase)
  
  # find common set of features across train and validation data
  commonFeatures <- featureNames(trainEset)
  for(eset in testEsets){
    commonFeatures <- intersect(commonFeatures, featureNames(eset))
  }
  commonFeatures <- sort(commonFeatures)
  
  idxs <- match(commonFeatures, featureNames(trainEset))
  trainEset <- trainEset[na.omit(idxs),]
  for(i in seq_along(testEsets)){
    eset <- testEsets[[i]]
    idxs <- match(commonFeatures, featureNames(eset))
    testEsets[[i]] <- eset[na.omit(idxs),]
  }
  
  if(quantileNormalize){
    cat("quantile normalizing...\n")
    ref <- exprs(trainEset)
    for(i in 1:length(testEsets)){
      eset <- testEsets[[i]]
      qndata <- normalize2Reference(exprs(eset), rowMeans(ref))
      exprs(eset) <- qndata
      testEsets[[i]] <- eset
    }
  }
  
  # train EN model
  set.seed(seed)
  cvFit <- cv.glmnet(t(exprs(trainEset)),
                     factor(trainResponse),
                     nfolds=5,
                     alpha=alpha,
                     family="binomial")
  fitM <- glmnet(t(exprs(trainEset)),
                 factor(trainResponse),
                 alpha=alpha, 
                 lambda=cvFit$lambda.1se, 
                 family="binomial")
  
  # predict on validation esets
  yhats <- lapply(testEsets, function(eset){
    testX <- normalizeToX(rowMeans(exprs(trainEset)), 
                          apply(exprs(trainEset), 1, sd), 
                          exprs(eset))
    
    yHat <- predict(fitM, t(testX),type="response")
    yHat
  })
  
  list(yhats=yhats,
       model=fitM,
       featureVec=featureNames(trainEset))
}


## KHAMBATA-FORD SIGNATURE
khambataEnrichmentMethod <- function(eset){
  require(Biobase)
  
  posIdxs <- featureNames(eset) %in% c("AREG","EREG","DUSP6")
  negIdxs <- featureNames(eset) %in% c("SLC26A3")
  E <- apply(exprs(eset), 2, function(x){
    v <- mean(x[posIdxs]) - mean(x[negIdxs])
    v
  })
  E
}

## LOBODA SIGNATURE
lobodaEnrichmentMethod <- function(eset){
  require(synapseClient)
  require(Biobase)
  
  sigsEnt <- downloadEntity("syn1680777")
  sigs <- load.gmt.data(file.path(sigsEnt$cacheDir, sigsEnt$file))
  positiveGeneset <- sigs[["loboda_ras_up"]]
  negativeGeneset <- sigs[["loboda_ras_down"]]
  
  posIdxs <- which(featureNames(eset) %in% positiveGeneset)
  negIdxs <- which(featureNames(eset) %in% negativeGeneset)	
  E <- apply(exprs(eset), 2, function(x){
    v <- mean(x[pos_idxs]) - mean(x[negIdxs])
    v
  })
  E
}

## GO ENRICHMENT SCORES
GOenrichment <- function(geneSet, backgroundGenes, minSz=10, maxSz=300){
  require(org.Hs.eg.db)
  require(GO.db)
  
  geneSetEg <- na.omit(unique(mget(geneSet, org.Hs.egSYMBOL2EG, ifnotfound=NA)))
  backgroundEg <- na.omit(unique(mget(backgroundGenes, org.Hs.egSYMBOL2EG, ifnotfound=NA)))
  
  gokeys <- keys(org.Hs.egGO2EG)
  GOMap <- as.list(org.Hs.egGO2EG)
  sizes <- lapply(gokeys, function(x){ length(GOMap[[x]])})
  
  gokeys <- gokeys[ (sizes >= minSz) & (sizes <= maxSz) ]
  ni <- length(gokeys)
  terms <- Term(GOTERM)
  out <- data.frame(Overlap=rep(-1, ni), TermSize=rep(-1, ni), Expected=rep(-1.0, ni), PvalueUpper=rep(-1.0, ni), PvalueLower=rep(-1.0, ni), Term=NA, stringsAsFactors=FALSE)
  rownames(out) <- gokeys[1:ni]
  for( i in 1:ni ){
    if( i%%1000 == 1 ) print(paste(i, "/", ni))
    gm <- GOMap[[gokeys[i]]]
    
    gm <- intersect(gm, backgroundEg)
    overlap <- intersect(gm, geneSetEg)
    p1 <- phyper(as.numeric(length(overlap)-1), as.numeric(length(gm)), as.numeric(length(backgroundEg)-length(gm)), as.numeric(length(geneSetEg)), lower.tail=FALSE)
    p2 <- phyper(as.numeric(length(overlap)), as.numeric(length(gm)), as.numeric(length(backgroundEg)-length(gm)), as.numeric(length(geneSetEg)), lower.tail=TRUE)
    
    out[i, 1] <- length(overlap)
    out[i, 2] <- length(gm)
    out[i, 3] <- as.numeric(length(geneSetEg)) / length(backgroundEg) * length(gm)
    out[i, 4] <- p1
    out[i, 5] <- p2
    out[i, 6] <- terms[[gokeys[i]]]
  }
  
  out <- out[order(as.numeric(out[,4])),]
  
  return(out)
}

