## PROGRAM PULLING TOGETHER ALL ANALYSIS STEPS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

ES <- function(scoreList, gsets){
	all.genes <- sort(unique(c(unlist(sapply(scoreList, names)))))
	
	M <- matrix(NA, nrow=length(all.genes),ncol=length(scoreList),
			dimnames=list(all.genes,1:4))
	for(i in 1:length(scoreList)){
		idxs <- match(names(scoreList[[i]]), all.genes)
		M[idxs, i] <- scoreList[[i]]	
	}
	
	scores.mean <- rowMeans(M,na.rm=TRUE)
	ranks.sort <- sort(scores.mean, decreasing=TRUE)
	
	genes.sort <- names(ranks.sort)
	E <- sapply(gsets, function(gset){
		mask <- genes.sort %in% gset
		tmp <- length(ranks.sort):1
		dec <- sum(abs(tmp[mask])) / (length(genes.sort) - sum(mask))
		tmp[!mask] <- -dec
		rw <- cumsum(tmp)
		max(rw)
	})

	
}
