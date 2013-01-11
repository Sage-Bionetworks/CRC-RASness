## NEED TO CONSOLIDATE DOWN TO ONE SINGLE GENE MEASUREMENT


## SCRIPT TO QC AND PROCESS THE KHAMBATA-FORD DATASET
#####
## ANALYST: BRIAN M. BOT
#####

require(rGithubClient)
require(corpcor)
require(lattice)

## SOURCE IN THE DATA EXTRACTION FUNCTIONS
sourceRepoFile(repository="Sage-Bionetworks/CRC-RASness",
               repositoryPath="functions/getDataFunctions.R",
               ref="branch", refName="dev")

## PULL IN THE DATA AND SUBSET TO COHORT OF INTEREST (TUMOR SAMPLES)
exprSet <- getKhambataFromGEO()
expr <- exprs(exprSet)
clin <- pData(exprSet)

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



## ONE MAJOR OUTLIER
# > colnames(expr)[ which(s$v[,2] < -.7) ]
# [1] "GSM136626"
# Tissue: Adrenal gland; Kidney

