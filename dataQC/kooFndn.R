## SCRIPT TO QC AND PROCESS THE KOO FOUNDATION DATASET
#####
## ANALYST: BRIAN M. BOT
#####

require(rGithubClient)
require(corpcor)
require(lattice)
require(hgu133plus2.db)

## SOURCE IN THE DATA EXTRACTION FUNCTIONS
myRepo <- getRepo(repository="Sage-Bionetworks/CRC-RASness",
                  ref="branch", refName="dev")
sourceRepoFile(myRepo, "functions/getDataFunctions.R")
sourceRepoFile(myRepo, "functions/utilityFunctions.R")

## PULL IN THE DATA AND SUBSET TO COHORT OF INTEREST (TUMOR SAMPLES)
exprSet <- getKFSYSCCdata()
scanDate <- as.Date(substr(exprSet@protocolData@data$ScanDate, 1, 8), format="%m/%d/%y")
exprSet <- exprSet[ , order(scanDate) ]
expr <- exprs(exprSet)


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
xyplot(s$v[,2] ~ s$v[,1],
       xlab="1st svd",
       ylab="2nd svd")
# xyplot(s$v[,4] ~ 1:nrow(s$v),
#        xaxt="n",
#        ylab="eigen gene")

## COMBINE TO ONE PROBE PER GENE
mapDat <- unlist(as.list(hgu133plus2ENSEMBL))[rownames(expr)]
exprGene <- combineProbesToGene(expr, mapDat)

## MAKE SURE NO LARGE AMOUNTS OF VARIATION HAVE BEEN ADDED BY PROCESS
sGene <- fs(exprGene)
xyplot(sGene$d ~ 1:length(sGene$d),
       xaxt="n",
       xlab="eigen gene",
       ylab="% variance explained")
xyplot(sGene$v[,2] ~ sGene$v[,1],
       xlab="1st svd",
       ylab="2nd svd")
xyplot(s$v[,1] ~ 1:nrow(s$v),
       xaxt="n",
       ylab="eigen gene")


