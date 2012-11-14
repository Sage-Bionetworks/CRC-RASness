## PROGRAM PULLING TOGETHER ALL ANALYSIS STEPS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

# global synapse variables
SYN_KFSYSCC_ADJUSTED_ID <- "163118"
SYN_TCGA_CRCEXPR_ADJUSTED_ID <- "140743"
SYN_TCGA_RASMUT_ID <- "162410"
SYN_TCGA_27KMETHYLATION <- "167612"
SYN_CCLE_EXPR_ID <- "48344"
SYN_CCLE_DRUGRESPONSE_ID <- "48359"
SYN_CCLE_MUT_ID <- "48341"
SYN_GAEDCKE_RECTAL_ID <- "140741"
SYN_KHAMBATA_CRC_ID <- "140742"

require(synapseClient)
require(Biobase)
require(affxparser)
require(org.Hs.eg.db)

#synapseLogin()

getOrFetch <- function(var, fetchExpr){
	tmp <- .GlobalEnv[[var]]
	if(is.null(tmp)){
		tmp <- eval(fetchExpr)
		.GlobalEnv[[var]] <- tmp
	}
	tmp
}

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE
# getCleanTcgaMaf <- function(){
#   coad.illumina <- read.table("data~/tcga_coad_read/hgsc.bcm.edu_COAD.IlluminaGA_DNASeq.1.maf",
#                               sep="\t",quote="",as.is=TRUE,header=TRUE,comment.char="",skip=1,fill=TRUE)
#   read.illumina <- read.table("data~/tcga_coad_read/hgsc.bcm.edu_READ.IlluminaGA_DNASeq.1.maf",
#                               sep="\t",quote="",as.is=TRUE,header=TRUE,comment.char="",skip=1,fill=TRUE)
#   all.illumina <- rbind(coad.illumina, read.illumina)
#   
#   coad.solid <- read.table("data~/tcga_coad_read/hgsc.bcm.edu_COAD.SOLiD_DNASeq.1.maf",
#                            sep="\t",quote="",as.is=TRUE,header=TRUE,comment.char="",skip=1,fill=TRUE)
#   read.solid <- read.table("data~/tcga_coad_read/hgsc.bcm.edu_READ.SOLiD_DNASeq.1.maf",
#                            sep="\t",quote="",as.is=TRUE,header=TRUE,comment.char="",skip=1,fill=TRUE)
#   all.solid <- rbind(coad.solid, read.solid)
#   
#   stopifnot(all(colnames(all.illumina) == colnames(all.solid)))
#   
#   maf <- rbind(all.illumina, all.solid)
#   
#   return (maf)
# }

getTCGAEpiReg <- function(){
	getOrFetch("TCGAEPIREG", loadEntity("SYN317477")$objects$epigRes)
}

getTCGAMeth <- function(){
	getOrFetch("TCGAMETH",loadEntity(getEntity(SYN_TCGA_27KMETHYLATION))$objects$eset)
}

getTCGACRC <- function(){
	getOrFetch("TCGACRC",loadEntity(getEntity(SYN_TCGA_CRCEXPR_ADJUSTED_ID))$objects$eset)
}

getTCGARasMuts <- function(){
	getOrFetch("TCGAMUTS",loadEntity(getEntity(SYN_TCGA_RASMUT_ID))$objects$eset)
}

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE
# getRASSigs <- function(){
# 	loadGmtData("resources/ras_signatures.gmt")
# }

getGaedcke <- function(){
	getOrFetch("GAEDCKE",loadEntity(getEntity(SYN_GAEDCKE_RECTAL_ID))$objects$adjusted_eset)
}

getKhambata <- function(){
	getOrFetch("KHAMBATA",loadEntity(getEntity(SYN_KHAMBATA_CRC_ID))$objects$eset)
}

getKFSYSCC <- function(){
	getOrFetch("KFSYSCC", loadEntity(getEntity(SYN_KFSYSCC_ADJUSTED_ID))$objects$eset)
}

getCCLE <- function(){
	getOrFetch("CCLE", make_ccle_eset())
}


## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE - study found -- no expression data
# getEMEXP3549 <- function(){
#   env <- new.env()
#   load("data~/MEXP3549/MEXP3549_eset.rda", env)
#   return (env$eset)
# }

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE - syn372544
# getEMEXP3557 <- function(){
#   env <- new.env()
#   load("data~/MEXP3557/MEXP3557_eset.rda", env)
#   return (env$eset)
# }

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE - syn202357
# getEMEXP991 <- function(){
#   annot <- read.table("data~/MEXP991/MEXP991.annotation.txt",sep="\t",as.is=TRUE,header=TRUE)
#   env <- new.env()
#   load(file="data~/MEXP991/MEXP991_eset.rda", env)
#   return (list(eset=env$gene.eset, annot=annot))
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
# 	M.m <- M[,!is.na(idxs)]
# 	sdrf.m <- sdrf[na.omit(idxs),]
# 	eset <- new("ExpressionSet",exprs=M.m)
# 	pData(eset) <- sdrf.m
# 	
# 	adf <- read.table("data~/EMTAB333/A-MEXP-503.adf.txt",sep="\t",header=T,as.is=T,skip=21,quote="")
# 	
# 	idxs <- match(featureNames(eset), adf$Reporter.Name)
# 	eset.m <- eset[!is.na(idxs),]
# 	adf.m <- adf[na.omit(idxs),]
# 	
# 	genes <- adf.m$Comment.AEReporterName.
# 	tmp <- combine_probes_2_gene(exprs(eset.m), genes)
# 	exprs(eset.m) <- tmp
# 	eset.m
# }

make_ccle_eset <- function(){
	ccle_entity_expr <- loadEntity(getEntity(SYN_CCLE_EXPR_ID))$objects$exprSet
	ccle_entity_response <- loadEntity(getEntity(SYN_CCLE_DRUGRESPONSE_ID))$objects$responseADF
	ccle_entity_mut <- loadEntity(getEntity(SYN_CCLE_MUT_ID))$objects$oncomapSet
	
	# merge expr with mut
	idxs <- match(sampleNames(ccle_entity_expr), sampleNames(ccle_entity_mut))
	
	tmp <- ccle_entity_expr[,!is.na(idxs)]
	pData(tmp) <- data.frame(t(exprs(ccle_entity_mut[,na.omit(idxs)])))
	
	ccle_entity_expr <- tmp
	
	idxs <- match(sampleNames(ccle_entity_expr), sampleNames(ccle_entity_response))
	
	ccle_eset <- ccle_entity_expr[,!is.na(idxs)]
	ccle_response <- ccle_entity_response[na.omit(idxs)]
	(list(ccle_eset, ccle_response))
}

## brian-bot TO DO: CHANGE TO PULL FROM SYNAPSE
# getSangerCellLineData <- function(){
# 	
# 	expr <- getOrFetch("sanger.celline.expr",loadEntity("syn210931")$objects$eSet_expr)
# 	drug <- getOrFetch("sanger.celline.drug",pData(loadEntity(getEntity("syn220680"))$objects$adf_drug))
# 	idxs <- match(rownames(drug), sampleNames(expr))
# 	drug.m <- drug[!is.na(idxs),]
# 	expr.m <- expr[, na.omit(idxs)]
# 	colnames(drug.m) <- paste("DRUG.",colnames(drug.m),sep="")
# 	
# 	meta <- read.table("resources/Sanger_affy_n798_sample_info_published.txt",sep="\t",header=T,fill=TRUE,as.is=TRUE,quote="")
# 	idxs <-  match(rownames(drug.m), gsub("-","",meta$SampleName))
# 	clinical.m <- cbind(drug.m, meta[idxs,])
# 	pData(expr.m) <- clinical.m
# 	return (expr.m)
# }

getCCLE_MetaGenomics <- function(){
  pharmaDat <- getCCLEPharma_MetaGenomics()
  exprDat <- getCCLEExpr_MetaGenomics()
  
  idxs <- match(rownames(pharmaDat), sampleNames(exprDat))
  pharmaDat.m <- as.data.frame(pharmaDat[!is.na(idxs),])
  exprDat.m <- exprDat[, na.omit(idxs)]
  eset <- new("ExpressionSet", expr=exprs(exprDat.m), phenoData=new("AnnotatedDataFrame",data=pharmaDat.m))
  eset
}

getCCLEPharma_MetaGenomics <- function(){
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
  
  ccleMat
}

getCCLEExpr_MetaGenomics <- function(){
  
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
    tmp
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
    averageExpr
  }
  
  ccleExpr <- getExpr(mapId='syn1417166',dataId='syn425002')
  entrez.ids <- gsub("(\\d*)_mt","\\1", rownames(ccleExpr))
  genes <- unlist(mget(entrez.ids, org.Hs.egSYMBOL, ifnotfound=NA))
  
  tmp <- combine_probes_2_gene(ccleExpr, genes)
  eset <- new("ExpressionSet",expr=tmp)
  eset
}

matchColumns <- function(dat, cellLines){
  ids1 <- match(toupper(dat$cellLineName), toupper(cellLines$Cell.line.name))
  ids2 <- match(toupper(dat$cellLineName), toupper(cellLines$Simplified.cell.name))
  ids3 <- match(toupper(dat$NAME), toupper(cellLines$Cell.line.name))
  ids4 <- match(toupper(dat$NAME), toupper(cellLines$Simplified.cell.name))
  mat <- cbind(ids1,ids2,ids3,ids4)
  newNames <- apply(mat, 1, function(x){ 
    if(sum(!is.na(x)) == 0){ 
      return("notAvailable")
    }else{ 
      cellLines$TIP_name[x[which(!is.na(x))[1]]]
    }
  })
 
  newNames
}


############################################################
### utility functions to initially upload the data

# synapseUpload <- function(){
# 	data.dir <- "~/projects/az_ras/resources/colon/"
# 	e <- new.env()
# 	load(paste(data.dir,"kfsyscc/eset_kfsyscc_colon_processed.rda",sep=""),e)
# 	upload.eset(SYN_KFSYSCC_ADJUSTED_ID, e$corrected_eset_genes)
# 	
# 	e <- new.env()
# 	load(paste(data.dir,"khambata/eset_khambata_colon_processed.rda",sep=""),e)
# 	upload.eset(SYN_KHAMBATA_CRC_ID, e$eset_genes)
# 	
# 	e <- new.env()
# 	load(paste(data.dir,"tcga_crc_adjusted_eset_new.rda",sep=""),e)
# 	upload.eset(SYN_TCGA_CRCEXPR_ADJUSTED_ID, e$tcga_crc_adjusted_eset)
# 	
# 	e <- new.env()
# 	load("/gluster/work/DAT_088__TCGA_COAD_READ/Results/ColonCancerTCGA_DNAmethylation_237patients.RData", e)
# 	upload.eset(SYN_TCGA_27KMETHYLATION, e$methSet)
# 	
# 	e <- new.env()
# 	load("~/projects/az_ras/analysis/jguinney/tcga_crc/kras_braf_pik3ca_nras_hras.rda",e)
# 	upload.eset(SYN_TCGA_RASMUT_ID, e$tcga_aa_mat)
# }

# upload.eset <- function(entity_id, eset){
# 	layer <- getEntity(entity_id)
# 	layer <- addObject(entity=layer, object=eset)
# 	storeEntity(layer)
# }
