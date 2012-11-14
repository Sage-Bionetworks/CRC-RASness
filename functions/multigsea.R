## PROGRAM PULLING TOGETHER ALL ANALYSIS STEPS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

ES <- function(scoreList, gsets){
	allGenes <- sort(unique(c(unlist(sapply(scoreList, names)))))
	
	M <- matrix(NA, nrow=length(allGenes), ncol=length(scoreList),
			dimnames=list(allGenes,1:4))
	for(i in 1:length(scoreList)){
		idxs <- match(names(scoreList[[i]]), allGenes)
		M[idxs, i] <- scoreList[[i]]	
	}
	
	scoresMean <- rowMeans(M,na.rm=TRUE)
	ranksSort <- sort(scoresMean, decreasing=TRUE)
	
	genesSort <- names(ranksSort)
	E <- sapply(gsets, function(gset){
		mask <- genesSort %in% gset
		tmp <- length(ranksSort):1
		dec <- sum(abs(tmp[mask])) / (length(genesSort) - sum(mask))
		tmp[!mask] <- -dec
		rw <- cumsum(tmp)
		max(rw)
	})

	return(E)
}
