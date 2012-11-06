
# create indices for merging across multiple data sets
merge.across <- function(...){
  tmp <- list(...)
  x.ids <- tmp[[1]]
  for(i in 2:length(tmp)){
    x.ids <- intersect(x.ids, tmp[[i]])
  }
  sapply(1:length(tmp), function(i){ match(x.ids, tmp[[i]]) })
}



buildTopTableFromCorrelationMatrix <- function(C, idxs, top=20){
	lbls.col <- colnames(C)
	lbls.row <- rownames(C)
	N.row <- nrow(C)
	N.col <- ncol(C)
	
	
	c <- floor(idxs[1:top] / N.row)
	r <- idxs[1:top] %% N.row
	tbl <- cbind(lbls.col[c],lbls.row[r],C[idxs[1:top]])
	tbl

}

load.gmt.data <- function(gmt.file.path){
	tmp <- readLines(gmt.file.path)
	gsets <- list()
	for(i in 1:length(tmp)){
		t <- strsplit(tmp[i],'\t')[[1]]
		t2 <- unlist(sapply(t[3:length(t)], function(x){ 
							strsplit(x,"///")
						}),use.names=FALSE)
		gsets[[t[1]]] <- gsub(" ","",t2)
	}
	return (gsets)
}

extract.tcga.patientIds <- function(tcga_ids){
	
	fixIds <- function(tcga_ids){
		return (gsub("\\.","-", as.matrix(tcga_ids)))
	}
	
	parts = strsplit(fixIds(tcga_ids),"-",fixed=TRUE)
	patient_ids = sapply(parts,
			function(x){ paste(x[1],x[2],x[3],sep="-") },
			simplify=TRUE)
}


normalize_to_X <- function(mean.x, sd.x, Y){
	m.y <- rowMeans(Y)
	sd.y <- apply(Y, 1, sd)
	Y.adj <- (Y - m.y) * sd.x / sd.y  + mean.x 
	Y.adj[sd.y == 0] <- mean.x[sd.y==0]
	Y.adj
}

custom.xtable <- function(df, custom.function, digits=NULL){
		
	doDigit <- function(x,col){
		if(is.numeric(x) & !is.null(digits)){
			dig <- ifelse(length(digits) > 1, digits[col], digits)
			return (format(x, digits=abs(dig),scientific=dig<0 ))
		}else{
			return (x)
		}	
	}
	tmp <- df
	for(r in 1:nrow(tmp)){
		for(c in 1:ncol(tmp)){
			tmp[r,c] <- custom.function(doDigit(df[r,c],c),r,c)
		}
	}
	tmp
}

displayGenomicFeatures <- function(featureList,colorSchemes=NULL,max.sample.width=100){
	
	if(is.null(colorSchemes)){
		colorSchemes <- list(c("red","blue"),c("green","yellow"),c("purple","orange"))
	}
	
	makePlot <- function(FL, sampleWidthCount){
		N.sample <- length(FL[[1]])
		N.features <- length(FL)	
		sample.space=5
		feature.space=10
		iconWidth=15
		
		maxWidth <- sampleWidthCount * (sample.space + iconWidth)
		
		startX <- seq(0,(N.sample-1) * (iconWidth+sample.space), by=(iconWidth+sample.space))
		endX <- startX + iconWidth
		startY <- seq(0,(N.features-1) * (iconWidth+feature.space), by=(iconWidth+feature.space))
		endY <- startY + iconWidth
		op <- par(mar = c(1,5,1,1))
		plot(c(-2, maxWidth + 2), c(0, N.features * (feature.space + iconWidth)), bty="n",
				type = "n", xlab="", ylab="",main="",yaxt="n",xaxt="n")
		for(i in 1:N.features){
			colorScheme <- colorSchemes[[i]]
			cols <- rep(colorScheme[1], N.sample)
			cols[FL[[i]] > 0] <- colorScheme[2]
			rect(startX, rep(startY[i], N.sample), endX, rep(endY[i], N.sample),col=cols,lwd=2)
		}
		axis(side=2, at=(startY + iconWidth/2), labels=names(FL),las=2)
	}		
	
	
	
	N <- length(featureList[[1]])
	
	nGroups <- ceiling(N / max.sample.width)
	sampleWidthCount <- min(N, max.sample.width)
	par(mfrow=c(nGroups,1))
	
	for(i in 1:nGroups){
		tmpFeatureList <- list()
		for(j in 1:length(featureList)){
			start <- max.sample.width * (i-1) + 1
			end <- start + (max.sample.width-1)
			if(end > N){ end <- N }
			tmpFeatureList[[j]] <- featureList[[j]][start:end]
		}
		names(tmpFeatureList) <- names(featureList)
		makePlot(tmpFeatureList, sampleWidthCount)
	}
}

GOenrichment <- function(geneSet, backgroundGenes,minSz=10,maxSz=300){
	
	library(org.Hs.eg.db)
	library(GO.db)
	
	geneSetEg <- na.omit(unique(mget(geneSet,org.Hs.egSYMBOL2EG,ifnotfound=NA)))
	backgroundEg <- na.omit(unique(mget(backgroundGenes, org.Hs.egSYMBOL2EG,ifnotfound=NA)))
	
	gokeys <- keys(org.Hs.egGO2EG)
	GOMap <- as.list(org.Hs.egGO2EG)
	sizes <- lapply(gokeys, function(x){ length(GOMap[[x]])})
	
	gokeys <- gokeys[sizes >= minSz & sizes <= maxSz]
	ni <- length(gokeys)
	terms <- Term(GOTERM)
	out <- data.frame(Overlap=rep(-1,ni), TermSize=rep(-1,ni), Expected=rep(-1.0,ni), PvalueUpper=rep(-1.0,ni),PvalueLower=rep(-1.0,ni), Term=NA, stringsAsFactors=FALSE)
	rownames(out) <- gokeys[1:ni]
	for (i in 1:ni)
	{
		if (i%%1000==1) print(paste(i,"/",ni))
		#gm <-  unique(get(gokeys[i],IlluminaHumanMethylation450kGO2ALLPROBES))
		gm <- GOMap[[gokeys[i]]]
		
		gm <- intersect(gm, backgroundEg)
		overlap <- intersect(gm,geneSetEg)
		p1 <- phyper(as.numeric(length(overlap)-1), as.numeric(length(gm)), as.numeric(length(backgroundEg)-length(gm)), as.numeric(length(geneSetEg)), lower.tail=FALSE)
		p2 <- phyper(as.numeric(length(overlap)), as.numeric(length(gm)), as.numeric(length(backgroundEg)-length(gm)), as.numeric(length(geneSetEg)), lower.tail=TRUE)
		
		out[i,1] <- length(overlap)
		out[i,2] <- length(gm)
		out[i,3] <- as.numeric(length(geneSetEg)) / length(backgroundEg) * length(gm)
		out[i,4] <- p1
		out[i,5] <- p2
		out[i,6] <- terms[[gokeys[i]]]
	}
	
	out <- out[order(as.numeric(out[,4])),]
	
	return(out)
}


combine_probes_2_gene <- function(expr, genes, method="svd"){
	
	if(is.list(genes)) genes <- unlist(genes)
	
	stopifnot(dim(expr)[1] ==  length(genes))
	ugenes <- unique(genes)
	ugenes <- sort(ugenes[!is.na(ugenes)])
	M <- matrix(NaN, ncol=dim(expr)[2],nrow=length(ugenes),
			dimnames=list(ugenes, colnames(expr)))
	
	for(gene in ugenes){
		sub.expr <- as.matrix(expr[which(genes == gene),])
		if(dim(sub.expr)[2] == 1){
			M[gene,] <- sub.expr
		}else{
			tmp <- svd(sub.expr - rowMeans(sub.expr))$v[,1]
			tmp.c <- mean(cor(tmp, t(sub.expr)))
			#cat(gene," ", tmp.c, "\n")
			multiplier <- ifelse(tmp.c < 0, -1, 1)
			M[gene,] <- tmp * multiplier
		}
	}
	M
}

