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

## CREATE A DIRECTORY FOR PLOTS TO BE UPLOADED TO SYNAPSE
tcgaAgilentDir <- file.path(tempdir(), "tcgaAgilentQC")
dir.create(tcgaAgilentDir)

## PULL IN THE AGILENT TCGA DATA AND SUBSET TO COHORT OF INTEREST
expr <- getTCGAcrcAgilent()
clin <- getTCGAcrcClinical()
clin <- clin[colnames(expr), ]

## CONVENIENCE FUNCTION FOR SVD EVALUATIONS
fs <- function(x){
  u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
  u$d <- u$d^2/sum(u$d^2)
  return(u)
}

## ONE PATIENT HAS 112 MISSING VALUES
# > colnames(expr)[colSums(is.na(expr))==112]
# [1] "TCGA-AY-4071"

## REMOVE GENES WHERE AT LEAST ONE PATIENT HAS A MISSING VALUE
expr <- expr[ rowSums(is.na(expr)) == 0, ]

## SVD ON EXPRESSION MATRIX -- ASSESS OVERALL STRUCTURE AND POSSIBLE LATENT STRUCTURE
s <- fs(expr)
png(file.path(tcgaAgilentDir, "tcgaAgilent_tcgaFreeze_pctVar.png"))
xyplot(s$d ~ 1:length(s$d),
       xaxt="n",
       xlab="eigen gene",
       ylab="% variance explained")
dev.off()
png(file.path(tcgaAgilentDir, "tcgaAgilent_tcgaFreeze_svd1v2.png"))
xyplot(s$v[,2] ~ s$v[,1], groups=factor(ifelse(grepl("RECTUM", clin$tumor_tissue_site), "rectum", "colon")),
       xlab="1st svd",
       ylab="2nd svd", auto.key=list(space="right"))
dev.off()

## PATIENT NOTED ABOVE DOES NOT SEEM LIKE MAJOR OUTLIER WHEN GENES WITH MISSING VALUES ARE REMOVED
## LIKELY NOT STRUCTURAL SO WILL LEAVE IN FOR ANALYSIS
# xyplot(s$v[,2] ~ s$v[,1], groups=factor(colnames(expr)=="TCGA-AY-4071"))


tcgaAgilentQcedExpressionSet <- ExpressionSet(expr)
pData(tcgaAgilentQcedExpressionSet) <- clin

## UPLOAD INFORMATION TO SYNAPSE
require(synapseClient)
tcgaAgilentEnt <- Data(name="tcgaAgilentQcedExpressionSet.rda", parentId="syn1585011")
tcgaAgilentEnt <- createEntity(tcgaAgilentEnt)
tcgaAgilentEnt <- addObject(tcgaAgilentEnt, tcgaAgilentQcedExpressionSet)
act <- Activity(list(name="tcgaAgilent expression QC",
                     used=list(list(entity="syn417828", wasExecuted=F),
                               list(entity="syn418082", wasExecuted=F),
                               list(entity="syn1446080", wasExecuted=F),
                               list(entity="syn1446153", wasExecuted=F))))
act <- createEntity(act)
generatedBy(tcgaAgilentEnt) <- act
tcgaAgilentEnt <- storeEntity(tcgaAgilentEnt)
## ADD DESCRIPTION ON WEB

attachThese <- list.files(tcgaAgilentDir, full.names=T)
attachThese <- attachThese[ grep(".png", attachThese, fixed=T) ]
for(i in attachThese){
  tcgaAgilentEnt <- synapseClient:::addAttachment(tcgaAgilentEnt, i)
}
tcgaAgilentEnt <- synapseClient:::storeAttachment(tcgaAgilentEnt)



