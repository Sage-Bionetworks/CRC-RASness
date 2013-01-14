## SCRIPT TO QC AND PROCESS THE GAEDCKE DATASET
#####
## ANALYST: BRIAN M. BOT
#####

require(rGithubClient)
require(corpcor)
require(lattice)

## SOURCE IN THE DATA EXTRACTION FUNCTIONS
myRepo <- getRepo(repository="Sage-Bionetworks/CRC-RASness",
                  ref="branch", refName="dev")
sourceRepoFile(myRepo, "functions/getDataFunctions.R")
sourceRepoFile(myRepo, "functions/utilityFunctions.R")

## PULL IN THE DATA AND SUBSET TO COHORT OF INTEREST (TUMOR SAMPLES)
exprSet <- getGaedckeFromGEO()
expr <- exprs(exprSet)
clin <- pData(exprSet)
feat <- exprSet@featureData@data

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

## COMBINE TO ONE PROBE PER GENE
idx <- (is.na(feat$"Gene ID")) | (feat$"Gene ID" == "")
expr <- expr[ !idx, ]
feat <- feat[ !idx, ]

exprGene <- combineProbesToGene(expr, feat$"Gene ID")

## MAKE SURE NO LARGE AMOUNTS OF VARIATION HAVE BEEN ADDED BY PROCESS
sGene <- fs(exprGene)
xyplot(sGene$d ~ 1:length(sGene$d),
       xaxt="n",
       xlab="eigen gene",
       ylab="% variance explained")
xyplot(sGene$v[,2] ~ sGene$v[,1],
       xlab="1st svd",
       ylab="2nd svd")
