## PROGRAM PULLING TOGETHER ALL ANALYSIS STEPS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

require(glmnet)
require(Biobase)
require(ROCR)

findActiveMethylationSites <- function(methEset, geneEset, outlierPercent=.05, rhoThreshold=-.25){
	
	require(IlluminaHumanMethylation27k.db)
	cpgSymbols <- as.list(IlluminaHumanMethylation27kSYMBOL)
	
	cpgSymbols <- cpgSymbols[!is.na(unlist(cpgSymbols))]
	cpgSymbols <- cpgSymbols[unlist(cpgSymbols) %in% featureNames(geneEset)]
	
	idxs <- match(extractTcgaPatientIds(sampleNames(methEset)), 
			extractTcgaPatientIds(sampleNames(geneEset)))
	
	methEsetM <- methEset[, !is.na(idxs)]
	geneEsetM <- geneEset[, na.omit(idxs)]
	
	idxs <- match(featureNames(methEsetM), names(cpgSymbols))
	methEsetM <- methEsetM[!is.na(idxs),]
	cpgSymbolsM <- cpgSymbols[na.omit(idxs)]
	
	# test 1: find outliers
	ranks <- apply(exprs(methEsetM), 1, function(x){
			x <- x - mean(x)
			sd <- sd(x)
			max(sum(x > 2 * sd), sum(x < -2 * sd))
	})
	outlierTest <- ranks > ncol(methEsetM) * outlierPercent
	
	# test 2: find anti-correlation with gene expression
	M <- exprs(methEsetM)
	G <- exprs(geneEsetM)
	geneIdxs <- match(unlist(cpgSymbolsM), featureNames(geneEsetM))
	rhos <- sapply(1:nrow(methEsetM), function(i) {
		cor(M[i,], G[geneIdxs[i],],method="spearman")
	})

	return(list(outlierCpgs=featureNames(methEsetM)[outlierTest],
				 rhoCpgs=featureNames(methEsetM)[ outlierTest & (rhos < rhoThreshold) ]))
}

khambataEnrichmentMethod <- function(eset){
  posIdxs <- featureNames(eset) %in% c("AREG","EREG","DUSP6")
  negIdxs <- featureNames(eset) %in% c("SLC26A3")
  E <- apply(exprs(eset), 2, function(x){
    v <- mean(x[posIdxs]) - mean(x[negIdxs])
    v
  })
  return(E)
}


lobodaEnrichmentMethod <- function(eset){
	
	sigs <- getRASSigs()
	positiveGeneset <- sigs[["loboda_ras_up"]]
	negativeGeneset <- sigs[["loboda_ras_down"]]
	
	posIdxs <- which(featureNames(eset) %in% positiveGeneset)
	negIdxs <- which(featureNames(eset) %in% negativeGeneset)
	E <- apply(exprs(eset), 2, function(x){
		v <- mean(x[posIdxs]) - mean(x[negIdxs])
		v
	})
	return(E)
}

binomialPredictEN <- function(trainEset, trainResponse, testEsets, alpha=.1, seed=2012, quantileNormalize=F){
	
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

	list(yhats=yhats, model=fitM, featureVec=featureNames(trainEset))
}


