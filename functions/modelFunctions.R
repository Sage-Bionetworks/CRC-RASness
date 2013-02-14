## A COLLECTION OF FUNCTIONS USED IN MODELING
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####


binomial_predict_EN <- function(trainEset, trainResponse, testEsets, alpha=.1, seed=2012, quantile.normalize=F){
  
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
