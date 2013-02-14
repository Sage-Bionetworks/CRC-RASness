## A COLLECTION OF FUNCTIONS USED IN MODELING
#####
## DEPENDS ON THE utilityFunctions.R PROGRAM BEING SOURCED
#####
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####


## PREDICTION USING ELASTIC NET -- TRAINING AND TEST
binomial_predict_EN <- function(trainEset, trainResponse, testEsets, alpha=.1, seed=2012, quantile.normalize=F){
  require(glmnet)
  require(Biobase)
  
  # find common set of features across train and validation data
  common_features <- featureNames(trainEset)
  for(eset in testEsets){
    common_features <- intersect(common_features, featureNames(eset))
  }
  common_features <- sort(common_features)
  
  idxs <- match(common_features, featureNames(trainEset))
  trainEset <- trainEset[na.omit(idxs),]
  for(i in seq_along(testEsets)){
    eset <- testEsets[[i]]
    idxs <- match(common_features, featureNames(eset))
    testEsets[[i]] <- eset[na.omit(idxs),]
    
  }
  
  if(quantile.normalize){
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
  cv.fit <- cv.glmnet(t(exprs(trainEset)),
                      factor(trainResponse),
                      nfolds=5,
                      alpha=alpha,
                      family="binomial")
  fit.m <- glmnet(t(exprs(trainEset)),
                  factor(trainResponse),
                  alpha=alpha, 
                  lambda=cv.fit$lambda.1se, 
                  family="binomial")
  
  # predict on validation esets
  yhats <- lapply(testEsets, function(eset){
    testX <- normalize_to_X(rowMeans(exprs(trainEset)), 
                            apply(exprs(trainEset),1,sd), 
                            exprs(eset))
    
    y_hat <- predict(fit.m, t(testX),type="response")
    y_hat
  })
  
  list(yhats=yhats,model=fit.m,featureVec=featureNames(trainEset))
}


## KHAMBATA-FORD SIGNATURE
khambata_enrichment_method <- function(eset){
  require(Biobase)
  
  pos_idxs <- featureNames(eset) %in% c("AREG","EREG","DUSP6")
  neg_idxs <- featureNames(eset) %in% c("SLC26A3")
  E <- apply(exprs(eset), 2, function(x){
    v <- mean(x[pos_idxs]) - mean(x[neg_idxs])
    v
  })
  E
}

## LOBODA SIGNATURE
loboda_enrichment_method <- function(eset){
  require(synapseClient)
  require(Biobase)
  
  sigsEnt <- downloadEntity("syn1680777")
  sigs <- load.gmt.data(file.path(sigsEnt$cacheDir, sigsEnt$file))
  positive_geneset <- sigs[["loboda_ras_up"]]
  negative_geneset <- sigs[["loboda_ras_down"]]
  
  pos_idxs <- which(featureNames(eset) %in% positive_geneset)
  neg_idxs <- which(featureNames(eset) %in% negative_geneset)	
  E <- apply(exprs(eset), 2, function(x){
    v <- mean(x[pos_idxs]) - mean(x[neg_idxs])
    v
  })
  E
}

## GO ENRICHMENT SCORES
GOenrichment <- function(geneSet, backgroundGenes,minSz=10,maxSz=300){
  
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

