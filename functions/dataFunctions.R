## COLLECTION OF DATA RETRIEVAL FUNCTIONS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

require(synapseClient)
require(Biobase)
require(affxparser)
require(org.Hs.eg.db)

#synapseLogin()

#####
## FUNCTION WHICH LOOKS IN GLOBAL ENVIRONMENT TO SEE
## IF VARIABLE ALREADY EXISTS BEFORE FETCHING IT
getOrFetch <- function(var, fetchExpr){
	tmp <- .GlobalEnv[[var]]
	if(is.null(tmp)){
		tmp <- eval(fetchExpr)
		.GlobalEnv[[var]] <- tmp
	}
	tmp
}

## not fully functional until support for external locations is available
#####
## READ IN TCGA MAF FILES FOR COAD AND READ
# getCleanTcgaMaf <- function(){
#   
#   ## GET TCGA COAD
#   coadIllFolder <- synapseQuery("SELECT id, name FROM entity WHERE entity.parentId=='syn1460948'")
#   coadIllFiles <- sapply(as.list(coadIllFolder$entity.id), function(synId){
#     tmp <- downloadEntity(synId)
#     return(file.path(tmp$cacheDir, tmp$files))
#   })
#   
#   coadSolidFolder <- synapseQuery("SELECT id, name FROM entity WHERE entity.parentId=='syn1460953'")
#   coadSolidFiles <- sapply(as.list(coadSolidFolder$entity.id), function(synId){
#     tmp <- downloadEntity(synId)
#     return(file.path(tmp$cacheDir, tmp$files))
#   })
#   
#   ## GET TCGA READ
#   readIllFolder <- synapseQuery("SELECT id, name FROM entity WHERE entity.parentId==''")
#   readIllFiles <- sapply(as.list(readIllFolder$entity.id), function(synId){
#     tmp <- downloadEntity(synId)
#     return(file.path(tmp$cacheDir, tmp$files))
#   })
#   
#   readSolidFolder <- synapseQuery("SELECT id, name FROM entity WHERE entity.parentId==''")
#   readSolidFiles <- sapply(as.list(readSolidFolder$entity.id), function(synId){
#     tmp <- downloadEntity(synId)
#     return(file.path(tmp$cacheDir, tmp$files))
#   })
#   
#   ## PICK THE FILES OF INTEREST
#   allFiles <- c(coadIllFiles, coadSolidFiles, readIllFiles, readSolidFiles)
#   mafFiles <- grep(".maf", allFiles, fixed=T)
#   
#   allMafs <- lapply(as.list(mafFiles), function(maf){
#     read.delim(maf, quote="", as.is=T, header=T, comment.char="", skip=1, fill=T)
#   })
#   masterMaf <- do.call(rbind, allMafs)
#   
#   return(masterMaf)
# }

## rename to getTCGAAgilent()
# getTCGACRC <- function(){
# 	getOrFetch("TCGACRC",loadEntity(getEntity(SYN_TCGA_CRCEXPR_ADJUSTED_ID))$objects$eset)
# }

# getTCGARasMuts <- function(){
# 	getOrFetch("TCGAMUTS",loadEntity(getEntity(SYN_TCGA_RASMUT_ID))$objects$eset)
# }


## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE
# getRASSigs <- function(){
# 	loadGmtData("resources/ras_signatures.gmt")
# }


# getGaedcke <- function(){
# 	getOrFetch("GAEDCKE",loadEntity(getEntity(SYN_GAEDCKE_RECTAL_ID))$objects$adjusted_eset)
# }

# getKhambata <- function(){
# 	getOrFetch("KHAMBATA",loadEntity(getEntity(SYN_KHAMBATA_CRC_ID))$objects$eset)
# }

# getKFSYSCC <- function(){
# 	getOrFetch("KFSYSCC", loadEntity(getEntity(SYN_KFSYSCC_ADJUSTED_ID))$objects$eset)
# }

getCCLE <- function(){
	getOrFetch("CCLE", makeCcleEset())
}


## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE - syn1490292 - 3 zip files
# getEMEXP3549 <- function(){
#   env <- new.env()
#   load("data~/MEXP3549/MEXP3549_eset.rda", env)
#   return(env$eset)
# }

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE - syn1490290
# getEMEXP3557 <- function(){
#   env <- new.env()
#   load("data~/MEXP3557/MEXP3557_eset.rda", env)
#   return(env$eset)
# }

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE - syn1503591
# getEMEXP991 <- function(){
#   annot <- read.table("data~/MEXP991/MEXP991.annotation.txt",sep="\t",as.is=TRUE,header=TRUE)
#   env <- new.env()
#   load(file="data~/MEXP991/MEXP991_eset.rda", env)
#   return(list(eset=env$gene.eset, annot=annot))
# }

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE - no study found
# getEMTAB333 <- function(){
# 	f <- list.files("data~/EMTAB333",pattern="LEM")
# 	tmp <- read.table(paste("data~/EMTAB333/",f[1],sep=""),sep="\t",header=T,as.is=T)
# 	M <- matrix(NA,nrow=nrow(tmp),ncol=length(f),dimnames=list(tmp$Probe_ID,f))
# 	for(fname in f){
# 		tmp <- read.table(paste("data~/EMTAB333/",fname,sep=""),sep="\t",header=T,as.is=T)
# 		M[tmp$Probe_ID,fname] <- tmp[,8]
# 	}
# 	M <- log(M - min(M)+1)
# 	sdrf <- read.table("data~/EMTAB333/E-MTAB-333.sdrf.txt",header=T,sep="\t",as.is=T,comment="")
# 	rownames(sdrf) <- sdrf$Array.Data.File
# 	idxs <- match(colnames(M),rownames(sdrf))
# 	Mm <- M[,!is.na(idxs)]
# 	sdrfM <- sdrf[na.omit(idxs),]
# 	eset <- new("ExpressionSet",exprs=Mm)
# 	pData(eset) <- sdrfM
# 	
# 	adf <- read.table("data~/EMTAB333/A-MEXP-503.adf.txt",sep="\t",header=T,as.is=T,skip=21,quote="")
# 	
# 	idxs <- match(featureNames(eset), adf$Reporter.Name)
# 	esetM <- eset[!is.na(idxs),]
# 	adfM <- adf[na.omit(idxs),]
# 	
# 	genes <- adfM$Comment.AEReporterName.
# 	tmp <- combineProbesToGene(exprs(esetM), genes)
# 	exprs(esetM) <- tmp
# 	esetM
# }

makeCcleEset <- function(){
	ccleEntityExpr <- loadEntity(getEntity(SYN_CCLE_EXPR_ID))$objects$exprSet
	ccleEntityResponse <- loadEntity(getEntity(SYN_CCLE_DRUGRESPONSE_ID))$objects$responseADF
	ccleEntityMut <- loadEntity(getEntity(SYN_CCLE_MUT_ID))$objects$oncomapSet
	
	# merge expr with mut
	idxs <- match(sampleNames(ccleEntityExpr), sampleNames(ccleEntityMut))
	
	tmp <- ccleEntityExpr[,!is.na(idxs)]
	pData(tmp) <- data.frame(t(exprs(ccleEntityMut[,na.omit(idxs)])))
	
	ccleEntityExpr <- tmp
	
	idxs <- match(sampleNames(ccleEntityExpr), sampleNames(ccleEntityResponse))
	
	ccleEset <- ccleEntityExpr[,!is.na(idxs)]
	ccleResponse <- ccleEntityResponse[na.omit(idxs)]
	(list(ccleEset, ccleResponse))
}

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE
# getSangerCellLineData <- function(){
# 	
# 	expr <- getOrFetch("sanger.celline.expr",loadEntity("syn210931")$objects$eSet_expr)
# 	drug <- getOrFetch("sanger.celline.drug",pData(loadEntity(getEntity("syn220680"))$objects$adf_drug))
# 	idxs <- match(rownames(drug), sampleNames(expr))
# 	drugM <- drug[!is.na(idxs),]
# 	exprM <- expr[, na.omit(idxs)]
# 	colnames(drugM) <- paste("DRUG.",colnames(drugM),sep="")
# 	
# 	meta <- read.table("resources/Sanger_affy_n798_sample_info_published.txt",sep="\t",header=T,fill=TRUE,as.is=TRUE,quote="")
# 	idxs <-  match(rownames(drugM), gsub("-","",meta$SampleName))
# 	clinical.m <- cbind(drugM, meta[idxs,])
# 	pData(exprM) <- clinical.m
# 	return (exprM)
# }

getCCLEmetaGenomics <- function(){
  pharmaDat <- getCCLEPharmaMetaGenomics()
  exprDat <- getCCLEExprMetaGenomics()
  
  idxs <- match(rownames(pharmaDat), sampleNames(exprDat))
  pharmaDatM <- as.data.frame(pharmaDat[!is.na(idxs),])
  exprDatM <- exprDat[, na.omit(idxs)]
  eset <- new("ExpressionSet", expr=exprs(exprDatM), phenoData=new("AnnotatedDataFrame",data=pharmaDatM))
  eset
}

getCCLEPharmaMetaGenomics <- function(){
  cellLines <- loadEntity('syn1417611')$objects$cellLines
  
  # Load the information on the pharmacologic profiling data
  cclePharma <- loadEntity('syn1412540')$objects$pharma
  
  compounds <- read.delim("~bmecham/compounds.txt",sep="\t",header=FALSE, stringsAsFactors=FALSE)
  compounds[which(is.na(compounds[,2])),2] <- compounds[which(is.na(compounds[,2])),1]
  compounds[ compounds[,2] == "",2] <- compounds[ compounds[,2] == "",1]
  compounds[,1] <- toupper(compounds[,1])
  compounds[,2] <- toupper(compounds[,2])
  sym2name <- c(split(compounds[,2], compounds[,1]), split(compounds[,2], compounds[,2]))
    
  # CCLE
  inCommon <- intersect(cclePharma$DRUG, names(sym2name))
  for(i in 1:length(inCommon)){
    cclePharma$DRUG[which(cclePharma$DRUG == inCommon[i])] <- sym2name[[inCommon[i]]]
  }

  ccleMap <- matchColumns(cclePharma, cellLines)
  cclePharma$NAME <- ccleMap
  
  # Build IC50 matrices
  ccleMat <- matrix(NA, nr=length(unique(cclePharma$NAME)), nc=length(unique(cclePharma$DRUG)))
  rownames(ccleMat) <- unique(cclePharma$NAME); colnames(ccleMat) <- unique(cclePharma$DRUG)
  
  for(i in 1:nrow(cclePharma)){
    ccleMat[cclePharma$NAME[i], cclePharma$DRUG[i]] <- cclePharma$activityArea[i]
  }
  
  return(ccleMat)
}

getCCLEExprMetaGenomics <- function(){
  
  averageReplicates <- function(dat,map){
    ids <- which(map != "notAvailable")
    dat <- dat[,ids]
    map <- map[ids]
    map2row <- split(1:ncol(dat), map)
    tmp <- matrix(NA, nr=nrow(dat),nc=length(map2row))
    colnames(tmp) <- names(map2row)
    rownames(tmp) <- rownames(dat)
    for(i in 1:length(map2row)){
      cat("\r",i)
      if(length(map2row[[i]]) == 1){
        tmp[,i] <- dat[,map2row[[i]]]
      }else{
        tmp[,i] <- rowMeans(dat[,map2row[[i]]])
      }
    }
    return(tmp)
  }
  
  getExpr <- function(mapId, dataId){
    cellLines <- loadEntity('syn1417611')$objects$cellLines
    
    # Load the various data matrices
    mapDat <- loadEntity(mapId)$objects$mapDat
    data <- loadEntity(dataId)
    dat <- exprs(data$objects$eset)
    
    # Find the columns in the mapDat objects corresponding to the column names
    mapDat <- mapDat[match(colnames(dat), mapDat$technologyID),]
    
    map <- matchColumns(mapDat, cellLines)
    
    averageExpr <- averageReplicates(dat, map)
    return(averageExpr)
  }
  
  ccleExpr <- getExpr(mapId='syn1417166',dataId='syn425002')
  entrezIds <- gsub("(\\d*)_mt","\\1", rownames(ccleExpr))
  genes <- unlist(mget(entrezIds, org.Hs.egSYMBOL, ifnotfound=NA))
  
  tmp <- combineProbesToGene(ccleExpr, genes)
  eset <- new("ExpressionSet", expr=tmp)
  eset
}

matchColumns <- function(dat, cellLines){
  ids1 <- match(toupper(dat$cellLineName), toupper(cellLines$Cell.line.name))
  ids2 <- match(toupper(dat$cellLineName), toupper(cellLines$Simplified.cell.name))
  ids3 <- match(toupper(dat$NAME), toupper(cellLines$Cell.line.name))
  ids4 <- match(toupper(dat$NAME), toupper(cellLines$Simplified.cell.name))
  mat <- cbind(ids1, ids2, ids3, ids4)
  newNames <- apply(mat, 1, function(x){ 
    if(sum(!is.na(x)) == 0){ 
      return("notAvailable")
    }else{ 
      cellLines$TIP_name[x[which(!is.na(x))[1]]]
    }
  })
 
  return(newNames)
}

