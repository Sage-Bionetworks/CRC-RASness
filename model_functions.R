require(glmnet)
require(Biobase)
require(ROCR)

find.active.methylation.sites <- function(meth.eset, gene.eset, outlier.percent=.05, rho.threshold=-.25){
	
	library(IlluminaHumanMethylation27k.db)
	cpgSymbols <- as.list(IlluminaHumanMethylation27kSYMBOL)
	
	cpgSymbols <- cpgSymbols[!is.na(unlist(cpgSymbols))]
	cpgSymbols <- cpgSymbols[unlist(cpgSymbols) %in% featureNames(gene.eset)]
	
	idxs <- match(extract.tcga.patientIds(sampleNames(meth.eset)), 
			extract.tcga.patientIds(sampleNames(gene.eset)))
	
	meth.eset.m <- meth.eset[, !is.na(idxs)]
	gene.eset.m <- gene.eset[, na.omit(idxs)]
	
	idxs <- match(featureNames(meth.eset.m), names(cpgSymbols))
	meth.eset.m <- meth.eset.m[!is.na(idxs),]
	cpgSymbols.m <- cpgSymbols[na.omit(idxs)]
	
	# test 1: find outliers
	ranks <- apply(exprs(meth.eset.m), 1, function(x){
			x <- x - mean(x)
			sd <- sd(x)
			max(sum(x > 2 * sd), sum(x < -2 * sd))
	})
	outlier.test <- ranks > ncol(meth.eset.m) * outlier.percent
	
	# test 2: find anti-correlation with gene expression
	M <- exprs(meth.eset.m)
	G <- exprs(gene.eset.m)
	gene.idxs <- match(unlist(cpgSymbols.m), featureNames(gene.eset.m))
	rhos <- sapply(1:nrow(meth.eset.m), function(i) {
		cor(M[i,], G[gene.idxs[i],],method="spearman")
	})

	return (list(outlier.cpgs=featureNames(meth.eset.m)[outlier.test],
				 rho.cpgs=featureNames(meth.eset.m)[outlier.test & rhos < rho.threshold]))
}

khambata_enrichment_method <- function(eset){
  pos_idxs <- featureNames(eset) %in% c("AREG","EREG","DUSP6")
  neg_idxs <- featureNames(eset) %in% c("SLC26A3")
  E <- apply(exprs(eset), 2, function(x){
    v <- mean(x[pos_idxs]) - mean(x[neg_idxs])
    v
  })
  E
}


loboda_enrichment_method <- function(eset){
	
	sigs <- getRASSigs()
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

binomial_predict_EN <- function(trainEset, trainResponse, testEsets, alpha=.1, seed=2012,quantile.normalize=F){
	
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


