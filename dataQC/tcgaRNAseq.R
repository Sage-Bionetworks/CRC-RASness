## SVD PLOTS LOOK WORRISOME -- WILL NEED TO CONSIDER REMOVAL OF LATENT STRUCTURE


## FUNCTIONS TO EXTRACT DATA OBJECTS FROM SYNAPSE AND QC
#####
## ANALYST: BRIAN M. BOT
#####

require(rGithubClient)
require(Biobase)
require(corpcor)
require(lattice)

## SOURCE IN THE DATA EXTRACTION FUNCTIONS
myRepo <- getRepo(repository="Sage-Bionetworks/CRC-RASness",
                  ref="branch", refName="dev")
sourceRepoFile(myRepo, "functions/getDataFunctions.R")
sourceRepoFile(myRepo, "functions/utilityFunctions.R")

## CREATE A DIRECTORY FOR PLOTS TO BE UPLOADED TO SYNAPSE
tcgaRNAseqDir <- file.path(tempdir(), "tcgaRNAseqQC")
dir.create(tcgaRNAseqDir)

## PULL IN THE RNAseq TCGA DATA AND SUBSET TO COHORT OF INTEREST
expr <- getTCGAcrcRNAseq()
clin <- getTCGAcrcClinical()
clin <- clin[colnames(expr), ]

## REMOVE GENES WHERE THERE ARE NO COUNTS AVAILABLE FOR ANY PATIENTS
expr <- expr[ rowSums(expr) != 0, ]

## GET RID OF GENES WITH NO GENE SYMBOL
rns <- sapply(strsplit(rownames(expr), "|", fixed=T), "[[", 1)
idx <- rns != "?"
expr <- expr[idx, ]

## ONE DUPLICATED GENE - COLLAPSE BY USING SVD
expr <- combineProbesToGene(expr, rns[idx])

## CONVENIENCE FUNCTION FOR SVD EVALUATIONS
fs <- function(x){
  u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
  u$d <- u$d^2/sum(u$d^2)
  return(u)
}

## SVD ON EXPRESSION MATRIX -- ASSESS OVERALL STRUCTURE AND POSSIBLE LATENT STRUCTURE
s <- fs(expr)
xyplot(s$d ~ 1:length(s$d),
       xaxt="n",
       xlab="eigen gene",
       ylab="% variance explained")
xyplot(s$v[,2] ~ s$v[,1], groups=factor(grepl("RECTUM", clin$tumor_tissue_site)),
       xlab="1st svd",
       ylab="2nd svd")

# factor(grepl("RECTUM", clin$tumor_tissue_site))

tcgaRNAseqQcedExpressionSet <- ExpressionSet(expr)
pData(tcgaRNAseqQcedExpressionSet) <- clin

## UPLOAD INFORMATION TO SYNAPSE
require(synapseClient)
tcgaRNAseqEnt <- Data(name="tcgaRNAseqQcedExpressionSet.rda", parentId="syn1585011")
tcgaRNAseqEnt <- createEntity(tcgaRNAseqEnt)
tcgaRNAseqEnt <- addObject(tcgaRNAseqEnt, tcgaRNAseqQcedExpressionSet)
act <- Activity(list(name="tcgaRNAseq expression QC",
                     used=list(list(entity="syn417828", wasExecuted=F),
                               list(entity="syn418082", wasExecuted=F),
                               list(entity="syn1446080", wasExecuted=F),
                               list(entity="syn1446153", wasExecuted=F))))
act <- createEntity(act)
generatedBy(tcgaRNAseqEnt) <- act
tcgaRNAseqEnt <- storeEntity(tcgaRNAseqEnt)
## ADD DESCRIPTION ON WEB

attachThese <- list.files(tcgaRNAseqDir, full.names=T)
attachThese <- attachThese[ grep(".png", attachThese, fixed=T) ]
for(i in attachThese){
  tcgaRNAseqEnt <- synapseClient:::addAttachment(tcgaRNAseqEnt, i)
}
tcgaRNAseqEnt <- synapseClient:::storeAttachment(tcgaRNAseqEnt)

