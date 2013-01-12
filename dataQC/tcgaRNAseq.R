## SVD PLOTS LOOK WORRISOME -- WILL NEED TO CONSIDER REMOVAL OF LATENT STRUCTURE


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
expr <- getTCGAcrcRNAseq()
clin <- getTCGAcrcClinical()
clin <- clin[colnames(expr), ]

## REMOVE GENES WHERE THERE ARE NO COUNTS AVAILABLE FOR ANY PATIENTS
expr <- expr[ rowSums(expr) != 0, ]


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
