## SCRIPT TO QC AND PROCESS THE GAEDCKE DATASET
#####
## ANALYST: BRIAN M. BOT
#####

require(rGithubClient)
require(limma)
require(Biobase)
require(corpcor)
require(lattice)

## SOURCE IN THE DATA EXTRACTION FUNCTIONS
myRepo <- getRepo(repository="Sage-Bionetworks/CRC-RASness",
                  ref="branch", refName="dev")
sourceRepoFile(myRepo, "functions/getDataFunctions.R")
sourceRepoFile(myRepo, "functions/utilityFunctions.R")

## CREATE A DIRECTORY FOR PLOTS TO BE UPLOADED TO SYNAPSE
gaedckeDir <- file.path(tempdir(), "gaedckeQC")
dir.create(gaedckeDir)

## PULL IN THE DATA AND SUBSET TO COHORT OF INTEREST (TUMOR SAMPLES)
exprSet <- getGaedckeFromGEO()
expr <- exprs(exprSet)
clin <- pData(exprSet)
feat <- fData(exprSet)

## FIX UP CLINICAL DATA A BIT
charNames <- "characteristics_ch1"
charNames <- c(charNames, paste(charNames, 1:5, sep="."))
for(header in charNames){
  tmp <- t(sapply(as.character(clin[[header]]), function(x){ strsplit(x, ": ")[[1]] } ))
  stopifnot(length(unique(tmp[, 1])) == 1)
  name <- unique(tmp[, 1])
  clin[ name ] <- tmp[, 2]
}
clin$kras_status <- ifelse(clin$"genome/variation" == "mutated KRAS", "MUT", "WT")


## GET RID OF GENES FOR WHICH WE HAVE NO ANNOTATION
mapDat <- as.character(feat$"Gene symbol")
mapDat[mapDat == ""] = NA
expr <- expr[ !is.na(mapDat), ]
feat <- feat[!is.na(mapDat), ]
mapDat <- mapDat[ !is.na(mapDat) ]

## CONVENIENCE FUNCTION FOR SVD EVALUATIONS
fs <- function(x){
  u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
  u$d <- u$d^2/sum(u$d^2)
  return(u)
}

## MAKE SURE NO LARGE AMOUNTS OF VARIATION HAVE BEEN ADDED BY PROCESS
s <- fs(expr)
png(filename=file.path(gaedckeDir, "gaedcke_probesetRma_pctVar.png"))
xyplot(s$d ~ 1:length(s$d),
       xaxt="n",
       xlab="eigen gene",
       ylab="% variance explained")
dev.off()
png(filename=file.path(gaedckeDir, "gaedcke_probesetRma_svd1v2.png"))
xyplot(s$v[,2] ~ s$v[,1],
       xlab="1st svd",
       ylab="2nd svd")
dev.off()

## ADJUST DATA FOR AGE AND GENDER
fit <- lmFit(expr, model.matrix(~as.numeric(clin$age) + factor(clin$gender)))
exprAdj <- residuals(fit, expr)

## COMBINE TO ONE MEASUREMENT PER GENE
exprGene <- combineProbesToGene(exprAdj, mapDat)
exprSetGene <- ExpressionSet(exprGene)
pData(exprSetGene) <- clin

## PLOTS FOR ADJUSTED EXPRESSION
sAdj <- fs(exprs(exprSetGene))
png(filename=file.path(gaedckeDir, "gaedcke_probesetRmaGeneSvdAdjusted_pctVar.png"))
xyplot(sAdj$d ~ 1:length(sAdj$d),
       xaxt="n",
       xlab="eigen gene",
       ylab="% variance explained")
dev.off()
png(filename=file.path(gaedckeDir, "gaedcke_probesetRmaGeneSvdAdjusted_svd1v2.png"))
xyplot(sAdj$v[,2] ~ sAdj$v[,1],
       xlab="1st svd",
       ylab="2nd svd")
dev.off()

gaedckeQcedExpressionSet <- exprSetGene



#http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20842
## UPLOAD INFORMATION TO SYNAPSE
gaedckeEnt <- Data(name="gaedckeQcedExpressionSet.rda", parentId="syn1586647")
gaedckeEnt <- createEntity(gaedckeEnt)
gaedckeEnt$properties$description <- 
  paste("The Gaedcke gene expression data was accessed from the Gene Expression Omnibus (GEO)",
        "using the R package `GEOquery`. These data were generated on the Agilent-014850 ",
        "Whole Human Genome Microarray 4x44K G4112F and were reported as being previously ",
        "quantile normalized using the `limma` package in R. We took the subset of samples ",
        "which were designated as tumor tissue and further consolidated the feature space",
        "by collapsing ",
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
        "${image?fileName=gaedcke_probesetRmaGeneSvd_pctVar.png}  ",
        "${image?fileName=gaedcke_probesetRmaGeneSvd_svd1v2.png}\n\n",
        "Once this final processing was complete, all patients who did not have KRAS ",
        "status evaluated were removed from the dataset and will not be used for further ",
        "analysis. This final R object is an `ExpressionSet` containing gene-level ",
        "expression values as well as clinical information on the samples within ",
        "the `phenoData` slot of the object.", sep="")
gaedckeEnt <- addObject(gaedckeEnt, gaedckeQcedExpressionSet)
act <- Activity(list(name="KFSYSCC expression QC", used=list(list(entity="syn1528362", wasExecuted=F))))
act <- createEntity(act)
generatedBy(gaedckeEnt) <- act
gaedckeEnt <- storeEntity(gaedckeEnt)

attachThese <- list.files(gaedckeDir, full.names=T)
attachThese <- attachThese[ grep(".png", attachThese, fixed=T) ]
for(i in attachThese){
  gaedckeEnt <- synapseClient:::addAttachment(gaedckeEnt, i)
}
gaedckeEnt <- synapseClient:::storeAttachment(gaedckeEnt)



