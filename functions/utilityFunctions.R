## A COLLECTION OF UTILILTY FUNCTIONS USED IN ANALYSIS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

## COMBINE PROBES TO GENES BY FIRST SV
combineProbesToGene <- function(expr, genes, method="svd"){
  
  if(is.list(genes)) genes <- unlist(genes)
  
  stopifnot(dim(expr)[1] ==  length(genes))
  ugenes <- unique(genes)
  ugenes <- sort(ugenes[!is.na(ugenes)])
  M <- matrix(NaN, ncol=dim(expr)[2], nrow=length(ugenes),
              dimnames=list(ugenes, colnames(expr)))
  
  for(gene in ugenes){
    subExpr <- as.matrix(expr[which(genes == gene),])
    if(dim(subExpr)[2] == 1){
      M[gene, ] <- subExpr
    }else{
      tmp <- svd(subExpr - rowMeans(subExpr))$v[,1]
      tmpC <- mean(cor(tmp, t(subExpr)))
      multiplier <- ifelse(tmpC < 0, -1, 1)
      M[gene,] <- tmp * multiplier
    }
  }
  return(M)
}

## CONVENIENCE FUNCTION FOR SVD EVALUATIONS
fs <- function(x){
  require(corpcor)
  u <- fast.svd(t(scale(t(x), scale = FALSE)), tol = 0)
  u$d <- u$d^2/sum(u$d^2)
  return(u)
}

## READ IN FILES AS FORMATED BY TCGA PAN CANCER GROUP
loadTCGAFileFromEntity <- function(synId){
  require(synapseClient)
  
  ent <- downloadEntity(synId)
  df <- read.delim(file.path(ent$cacheDir, ent$files), header=F, as.is=T)
  colnames(df) <- as.character(df[1, ])
  df <- df[-1, ]
  rownames(df) <- as.character(df[, 1])
  df <- df[, -1]
  return(df)
}


buildTopTableFromCorrelationMatrix <- function(C, idxs, top=20){
	lblsCol <- colnames(C)
	lblsRow <- rownames(C)
	nRow <- nrow(C)
	nCol <- ncol(C)
	
	c <- floor(idxs[1:top] / nRow)
	r <- idxs[1:top] %% nRow
	tbl <- cbind(lblsCol[c],lblsRow[r],C[idxs[1:top]])
	return(tbl)
}

## LOAD GMT-TYPE DATA
loadGmtData <- function(gmtFilePath){
	tmp <- readLines(gmtFilePath)
	gsets <- list()
	for(i in 1:length(tmp)){
		t <- strsplit(tmp[i],'\t')[[1]]
		t2 <- unlist(sapply(t[3:length(t)], function(x){ 
							strsplit(x,"///")
						}),use.names=FALSE)
		gsets[[t[1]]] <- gsub(" ","",t2)
	}
	return(gsets)
}

## EXTRACT TCGA PATIENT IDS FROM LONGER TCGA IDS
extractTcgaPatientIds <- function(tcgaIds){
	
	fixIds <- function(tcgaIds){
		return (gsub("\\.","-", as.matrix(tcgaIds)))
	}
	
	parts <- strsplit(fixIds(tcgaIds),"-",fixed=TRUE)
	patientIds <- sapply(parts,
			function(x){ paste(x[1],x[2],x[3],sep="-") },
			simplify=TRUE)
	return(patientIds)
}

## CREATE INDICES FOR MERGING ACROSS MULTIPLE DATA SETS
mergeAcross <- function(...){
  tmp <- list(...)
  xIds <- tmp[[1]]
  for(i in 2:length(tmp)){
    xIds <- intersect(xIds, tmp[[i]])
  }
  these <- sapply(1:length(tmp), function(i){ match(xIds, tmp[[i]]) })
  return(these)
}

## NORMALIZE A MATRIX Y TO ANOTHER X
normalizeToX <- function(meanX, sdX, Y){
	mY <- rowMeans(Y)
	sdY <- apply(Y, 1, sd)
	adjY <- (Y - mY) * sdX / sdY  + meanX 
	adjY[sdY == 0] <- meanX[sdY==0]
	return(adjY)
}

## PLOTTING FOR GENOMIC FEATURES -- CHANGE IN FUTURE?
displayGenomicFeatures <- function(featureList, colorSchemes=NULL, maxSampleWidth=100){
	
	if(is.null(colorSchemes)){
		colorSchemes <- list(c("red","blue"),c("green","yellow"),c("purple","orange"))
	}
	
	makePlot <- function(FL, sampleWidthCount){
		Nsample <- length(FL[[1]])
		Nfeatures <- length(FL)	
		sampleSpace <- 5
		featureSpace <- 10
		iconWidth <- 15
		
		maxWidth <- sampleWidthCount * (sampleSpace + iconWidth)
		
		startX <- seq(0,(Nsample-1) * (iconWidth+sampleSpace), by=(iconWidth+sampleSpace))
		endX <- startX + iconWidth
		startY <- seq(0,(Nfeatures-1) * (iconWidth+featureSpace), by=(iconWidth+featureSpace))
		endY <- startY + iconWidth
		op <- par(mar = c(1,5,1,1))
		plot(c(-2, maxWidth + 2), c(0, Nfeatures * (featureSpace + iconWidth)), bty="n",
				type = "n", xlab="", ylab="",main="",yaxt="n",xaxt="n")
		for(i in 1:Nfeatures){
			colorScheme <- colorSchemes[[i]]
			cols <- rep(colorScheme[1], Nsample)
			cols[FL[[i]] > 0] <- colorScheme[2]
			rect(startX, rep(startY[i], Nsample), endX, rep(endY[i], Nsample),col=cols,lwd=2)
		}
		axis(side=2, at=(startY + iconWidth/2), labels=names(FL),las=2)
	}		
	
	N <- length(featureList[[1]])
	
	nGroups <- ceiling(N / maxSampleWidth)
	sampleWidthCount <- min(N, maxSampleWidth)
	par(mfrow=c(nGroups,1))
	
	for(i in 1:nGroups){
		tmpFeatureList <- list()
		for(j in 1:length(featureList)){
			start <- maxSampleWidth * (i-1) + 1
			end <- start + (maxSampleWidth-1)
			if(end > N){ end <- N }
			tmpFeatureList[[j]] <- featureList[[j]][start:end]
		}
		names(tmpFeatureList) <- names(featureList)
		makePlot(tmpFeatureList, sampleWidthCount)
	}
}

## GO ENRICHMENT SCORES
GOenrichment <- function(geneSet, backgroundGenes,minSz=10,maxSz=300){
	
	require(org.Hs.eg.db)
	require(GO.db)
	
	geneSetEg <- na.omit(unique(mget(geneSet, org.Hs.egSYMBOL2EG, ifnotfound=NA)))
	backgroundEg <- na.omit(unique(mget(backgroundGenes, org.Hs.egSYMBOL2EG, ifnotfound=NA)))
	
	gokeys <- keys(org.Hs.egGO2EG)
	GOMap <- as.list(org.Hs.egGO2EG)
	sizes <- lapply(gokeys, function(x){ length(GOMap[[x]])})
	
	gokeys <- gokeys[ (sizes >= minSz) & (sizes <= maxSz) ]
	ni <- length(gokeys)
	terms <- Term(GOTERM)
	out <- data.frame(Overlap=rep(-1, ni), TermSize=rep(-1, ni), Expected=rep(-1.0, ni), PvalueUpper=rep(-1.0, ni), PvalueLower=rep(-1.0, ni), Term=NA, stringsAsFactors=FALSE)
	rownames(out) <- gokeys[1:ni]
	for( i in 1:ni ){
		if( i%%1000 == 1 ) print(paste(i, "/", ni))
		gm <- GOMap[[gokeys[i]]]
		
		gm <- intersect(gm, backgroundEg)
		overlap <- intersect(gm, geneSetEg)
		p1 <- phyper(as.numeric(length(overlap)-1), as.numeric(length(gm)), as.numeric(length(backgroundEg)-length(gm)), as.numeric(length(geneSetEg)), lower.tail=FALSE)
		p2 <- phyper(as.numeric(length(overlap)), as.numeric(length(gm)), as.numeric(length(backgroundEg)-length(gm)), as.numeric(length(geneSetEg)), lower.tail=TRUE)
		
		out[i, 1] <- length(overlap)
		out[i, 2] <- length(gm)
		out[i, 3] <- as.numeric(length(geneSetEg)) / length(backgroundEg) * length(gm)
		out[i, 4] <- p1
		out[i, 5] <- p2
		out[i, 6] <- terms[[gokeys[i]]]
	}
	
	out <- out[order(as.numeric(out[,4])),]
	
	return(out)
}

