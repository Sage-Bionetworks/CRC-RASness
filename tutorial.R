## PROGRAM PULLING TOGETHER ALL ANALYSIS STEPS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

library(multicore)
library(GSVA)
library(IlluminaHumanMethylation27k.db)
library(limma)
library(edgeR)
library(cqn)
source("code/data_functions.R")
source("code/util_functions.R")

eset <- getKFSYSCC()
fit <- eBayes(lmFit(eset, model.matrix(~factor(eset$kras))))

sigs <- load.gmt.data("resources/c2.all.v3.0.symbols.gmt")
hist(fit$p.value[,2])
#genes_001 <- featureNames(eset)[fit$p.value[,2] < .001]
genes <- featureNames(eset)[p.adjust(fit$p.value[,2],method="BH") < .05]
genes <- featureNames(eset)[fit$p.value[,2] < .001]

sz <- sapply(sigs, function(x){
	length(intersect(x,genes))
})


sig <- sigs[["SWEET_LUNG_CANCER_KRAS_UP"]]
sig.x <- intersect(sig, featureNames(eset))

percent.1 <- length(genes) / dim(eset)
E.1 <- length(sig.x) * .001
O.1 <- length(intersect(sig.x, genes))
pbinom(O.1, length(sig.x), .001,lower.tail=FALSE)


percent <- length(genes) / dim(eset)[1]
E <- percent * length(sig.x)
O <- length(intersect(genes, sig.x))
tbl <- table(featureNames(eset) %in% genes, featureNames(eset) %in% sig.x)
fisher.test(tbl,alternative="greater")

pvals <- sapply(sigs, function(x,y){
			fisher.test(y, featureNames(eset) %in% x,alternative="greater")$p.value
		},featureNames(eset) %in% genes)
