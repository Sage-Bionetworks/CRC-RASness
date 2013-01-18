## SCRIPT TO QC AND PROCESS THE KOO FOUNDATION DATASET
#####
## ANALYST: BRIAN M. BOT
#####

require(rGithubClient)
require(affyPLM)
require(limma)
require(chron)
require(corpcor)
require(lattice)
require(hgu133plus2.db)

## CREATE A DIRECTORY FOR PLOTS TO BE UPLOADED TO SYNAPSE
kooDir <- file.path(tempdir(), "kooQC")
dir.create(kooDir)

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

png(file.path(kooDir, "QA_CRC322.png"), width=600, height=600)
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

fit <- lmFit(exprSet, model.matrix(~as.numeric(pData(exprSet)$AGE)))
bfit <- eBayes(fit)
png(file.path(kooDir, "ageEffectHist.png"))
hist(bfit$p.value[,2], main="Limma: Effect of age on gene expression", xlab="p-value")
dev.off()
exprSetCorrected <- exprSet
exprs(exprSetCorrected) <- exprs(exprSet) - fitted.values(fit)

## COMBINE TO ONE PROBE PER GENE
mapDat <- unlist(as.list(hgu133plus2SYMBOL))[ rownames(expr) ]
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



## UPLOAD INFORMATION TO SYNAPSE
require(synapseClient)
kooEnt <- Data(name="kfsysccQcedExpressionSet.rda", parentId="syn1528361")
kooEnt <- createEntity(kooEnt)
kooEnt <- addObject(kooEnt, kfsysccQcedExpressionSet)
act <- Activity(list(name="KFSYSCC expression QC", used=list(list(entity="syn1528362", wasExecuted=F))))
act <- createEntity(act)
generatedBy(kooEnt) <- act
kooEnt <- storeEntity(kooEnt)
## ADD DESCRIPTION ON WEB

attachThese <- list.files(kooDir, full.names=T)
attachThese <- attachThese[ grep(".png", attachThese, fixed=T) ]
for(i in attachThese){
  kooEnt <- synapseClient:::addAttachment(kooEnt, i)
}
kooEnt <- synapseClient:::storeAttachment(kooEnt)
