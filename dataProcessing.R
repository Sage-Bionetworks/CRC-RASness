## DATA PREPROCESSING
## WHICH STORE INTERMEDIATE DATA IN SYNAPSE
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

require(synapseClient)
require(Biobase)
require(GEOquery)

#####
## PULL TOGETHER TCGA COLON (COAD) AND RECTUM (READ) DATA
#####

## AGILENT LEVEL 3 DATA
tcgaCRCAgilent <- do.call(lapply(as.list(c("syn1460914", "syn1465735")), mergeAgilentLevel3), cbind)

## RNA-SEQ LEVEL 3 DATA
tcgaCRCAgilent <- do.call(lapply(as.list(c("syn1460925", "syn1465752")), mergeRNAseqLevel3), cbind)
#tcgaCRCAgilent <- do.call(lapply(as.list(c("syn1460925", "syn1465752", "syn1460944", "syn1465772")), mergeRNAseqLevel3), cbind)

## MUTATION INFORMATION -- BOTH ILLUMINAGA AND SOLID PLATFORMS FOR EACH TISSUE TYPE
tcgaCRCMut <- do.call(lapply(as.list(c("syn1460948", "syn1460953", "syn1465778", "syn1465784")), mergeMafLevel2), rbind)


#####
## PULL GAEDCKE DATA FROM GEO
#####
geoFiles <- getGEO("GSE20842", GSEMatrix=T, AnnotGPL=T)
gaedckeEset <- geoFiles$GSE20842_series_matrix.txt.gz
# NOTE: KRAS STATUS - pData(gaedckeEset)$characteristics_ch1.5 IN pData


