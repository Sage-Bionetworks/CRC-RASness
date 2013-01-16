## SCRIPT TO QC AND PROCESS THE KOO FOUNDATION DATASET
#####
## ANALYST: BRIAN M. BOT
#####

kooFndnQC <- function(){
  
  require(rGithubClient)
  require(affyPLM)
  require(limma)
  require(chron)
  require(corpcor)
  require(lattice)
  require(hgu133plus2.db)
  
  ## SOURCE IN THE DATA EXTRACTION FUNCTIONS
  myRepo <- getRepo(repository="Sage-Bionetworks/CRC-RASness",
                    ref="branch", refName="dev")
  sourceRepoFile(myRepo, "functions/getDataFunctions.R")
  sourceRepoFile(myRepo, "functions/utilityFunctions.R")
  
  ## PULL IN THE DATA AND SUBSET TO COHORT OF INTEREST (TUMOR SAMPLES)
  kooRaw <- getKFSYSCCdata()
  
  scanDate <- protocolData(kooRaw)$ScanDate
  scanDate <- chron(gsub(" .*", "", scanDate))
  
  ## ORDER SAMPLES IN CHRONOLOGICAL ORDER OF SCAN DATE
  kooRaw <- kooRaw[, order(scanDate)]
  
  ## PERFORM SOME PRELIMINARY QC
  pset <- fitPLM(kooRaw)
  
  kooDir <- file.path(tempdir(), "kooQC")
  dir.create(kooDir)
  
  png(file.path(kooDir, "QA_CRC322.png"))
  par(mfrow=c(2,1))
  RLE(pset, main="RLE for CRC322")
  NUSE(pset, main="NUSE for CRC322")
  dev.off()
  rm(pset)
  
  # REMOVE EARLY SAMPLES; LARGE SE FROM NUSE PLOT
  kooRaw <- kooRaw[ , scanDate > "01/01/07" ]
  exprSet <- rma(kooRaw)
  expr <- exprs(exprSet)
  
  ## CONVENIENCE FUNCTION FOR SVD EVALUATIONS
  fs <- function(x){
    u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
    u$d <- u$d^2/sum(u$d^2)
    return(u)
  }
  
  ## SVD ON EXPRESSION MATRIX -- ASSESS OVERALL STRUCTURE AND POSSIBLE LATENT STRUCTURE
  s <- fs(expr)
  png(filename=file.path(kooDir, "koo_probesetRma_pctVar.png"))
  xyplot(s$d ~ 1:length(s$d),
         xaxt="n",
         xlab="singular value",
         ylab="% variation explained")
  dev.off()
  png(filename=file.path(kooDir, "koo_probesetRma_svd1v2.png"))
  xyplot(s$v[,2] ~ s$v[,1],
         xlab="1st svd",
         ylab="2nd svd")
  dev.off()
  xyplot(s$v[,2] ~ s$v[,1], groups=factor(pData(exprSet)$Grade),
         xlab="1st svd",
         ylab="2nd svd")
  
  fit <- lmFit(exprSet, model.matrix(~pData(exprSet)$AGE))
  bfit <- eBayes(fit)
  png(file.path(kooDir, "ageEffectHist.png"))
  hist(bfit$p.value[,2], main="Limma: Effect of age on gene expression", xlab="p-value")
  dev.off()
  exprSetCorrected <- exprSet
  exprs(exprSetCorrected) <- exprs(exprSet) - fitted.values(fit)
  
  ## COMBINE TO ONE PROBE PER GENE
  mapDat <- unlist(as.list(hgu133plus2ENSEMBL))[ rownames(expr) ]
  exprGene <- combineProbesToGene(exprs(exprSetCorrected), mapDat)
  exprs(exprSetCorrected) <- exprGene
  
  ## MAKE SURE NO LARGE AMOUNTS OF VARIATION HAVE BEEN ADDED BY PROCESS
  sGene <- fs(exprGene)
  png(filename=file.path(kooDir, "koo_probesetRmaGeneSvd_pctVar.png"))
  xyplot(sGene$d ~ 1:length(sGene$d),
         xaxt="n",
         xlab="eigen gene",
         ylab="% variance explained")
  dev.off()
  png(filename=file.path(kooDir, "koo_probesetRmaGeneSvd_svd1v2.png"))
  xyplot(sGene$v[,2] ~ sGene$v[,1],
         xlab="1st svd",
         ylab="2nd svd")
  dev.off()
  
  
  ## FINALIZE THE SAMPLES THAT WILL BE KEPT FOR FURTHER ANALYSIS
  ## REMOVE UNANNOTATED KRAS
  exprSetCorrected <- exprSetCorrected[ , pData(exprSetCorrected)$KRAS_Mutation != "NaN"]
  kras_status <- pData(exprSetCorrected)$KRAS_Mutation
  kras_status[kras_status=="None"] <- "WT"
  kras_status[kras_status=="Yes"] <- "MUT"
  pData(exprSetCorrected)$kras_status <- kras_status
  
  kfsysccQcedExpressionSet <- exprSetCorrected
  
  kooEnt <- Data(name="kfsysccQcedExpressionSet.rda", parentId="syn1528361")
  kooEnt <- createEntity(kooEnt)
  kooEnt <- addObject(kooEnt, kfsysccQcedExpressionSet)
  attachThese <- list.files(kooDir, full.names=T)
  attachThese <- attachThese[ grep(".png", attachThese, fixed=T) ]
  for(i in attachThese){
    kooEnt <- synapseClient:::addAttachment(kooEnt, i)
  }
  
  kooEnt$properties$description <- 
    paste("KFSYSCC generously provided 322 CRC samples with expression profiling on ",
          "the Affymetrix U133 Plus 2 platform. These samples were QCed initially ",
          "using the `affyPLM` package. The following plot shows the Relative Log ",
          "Expression (RLE) and Normalized Unscaled Standard Errors (NUSE) plots ",
          "for all 322 samples.\n\n",
          "${image$fileName=QA_CRC322.png}\n\n",
          "Given the high NUSE for samples collected before 01/01/2007, those 15 ",
          "samples were removed from all further analyses. We then ran `rma` on ",
          "the resulting 307 samples which consolidates the feature space down to 54,675 ",
          "probesets. Singular Value Decomposition (SVD) was performed on the resulting ",
          "expression matrix to assess its structure.\n\n",
          "${image?fileName=koo_probesetRma_pctVar.png} ",
          "${image?fileName=koo_probesetRma_svd1v2.png}\n\n",
          "While the first SVD appears to explain about 20% of the variation, however it was ",
          "not strongly associated with any clinical covariates and there appears to be no ",
          "single outlier nor distinct group of outliers that was the cause.\n\n",
          "The effect of age of the patient was explored further for association with ",
          "`rma` summarized expression levels. The `eBayes`-corrected p-value histogram ",
          "of the age effect is presented below.\n\n",
          "${image?fileName=ageEffectHist.png}\n\n",
          "The gene-specific fitted values of age were removed from the expression matrix ",
          "and the first SVD was used to collapse multiple probesets to a single ",
          "expression value per gene. SVD plots were repeated to ensure no artifacts ",
          "were introduced in this process.\n\n",
          "${image?fileName=koo_probesetRmaGeneSvd_pctVar.png}  ",
          "${image?fileName=koo_probesetRmaGeneSvd_svd1v2.png}\n\n",
          "Once this final processing was complete, all patients who did not have KRAS ",
          "status evaluated were removed from the dataset and will not be used for further ",
          "analysis. This final R object is an `ExpressionSet` containing gene-level ",
          "expression values as well as clinical information on the samples within ",
          "the `phenoData` slot of the object.", sep="")
  kooEnt <- storeEntity(kooEnt)
  kooEnt <- synapseClient:::storeAttachment(kooEnt)
  act <- Activity(list(name="QC of KFSYSCC expression data", used=list(list(entity="syn1528362", wasExecuted=F))))
  act <- createEntity(act)
  generatedBy(kooEnt) <- act
  kooEnt <- storeEntity(kooEnt)
  
}
