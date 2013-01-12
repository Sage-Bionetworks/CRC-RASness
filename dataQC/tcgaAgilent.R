## FUNCTIONS TO EXTRACT DATA OBJECTS FROM SYNAPSE AND QC
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
xyplot(s$d ~ 1:length(s$d),
       xaxt="n",
       xlab="eigen gene",
       ylab="% variance explained")
xyplot(s$v[,2] ~ s$v[,1], groups=factor(grepl("COLON", clin$tumor_tissue_site)),
       xlab="1st svd",
       ylab="2nd svd")

## PATIENT NOTED ABOVE DOES NOT SEEM LIKE MAJOR OUTLIER WHEN GENES WITH MISSING VALUES ARE REMOVED
## LIKELY NOT STRUCTURAL SO WILL LEAVE IN FOR ANALYSIS
# xyplot(s$v[,2] ~ s$v[,1], groups=factor(colnames(expr)=="TCGA-AY-4071"))

