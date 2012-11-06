library(ROCR)
library(glmnet)
library(caret)
library(GSVA)
library(parallel)
library(limma)
library(survival)
library(GenomeGraphs)
source("code/lung_data_functions.R")
source("code/model_functions.R")

chemores <- function(){
	eset <- getCHEMORES()
	surv <- Surv(as.numeric(eset$Time), as.numeric(eset$Status_dcd))
	survdiff(surv ~ eset$Chemo)
	survdiff(surv ~ eset$P_Stage)
	
	sapply(eset$Dat)
}

build_model_gse26939 <- function(){
	load("data~/lung_gse26939.rda")
	N <- dim(gene.eset)[2]
	aucs <- sapply(1:20, function(i){
				mask <- sample(N, round(N*.7))
				cv.fit <- cv.glmnet(t(exprs(gene.eset[, mask])), factor(gene.eset$kras[mask]), 
						alpha=.1,nfolds=5,family="binomial")
				
				fit <- glmnet(t(exprs(gene.eset[, mask])), factor(gene.eset$kras[mask]), 
						alpha=.1, 
						lambda=cv.fit$lambda.min,
						family="binomial")
				y_hat <- predict(fit, t(exprs(gene.eset[, -mask])),type="response")
				pred <- prediction(y_hat, factor(gene.eset$kras[-mask]))
				auc <- performance(pred, 'auc')@y.values
			})
}

test_model_gse26939 <- function(){
	env <- new.env()
	load("data~/GSE26939/lung_gse26939.rda", env)
	lungA <- env$gene.eset

	#load("data~/tcga.luad_rnaseq.rda",env)
	#tcga.luad <- env$eset
	tcga.luad <- get.luad.RNAseq(tumor.only=TRUE)
	tmp <- get.luad.exome(tcga.luad)
	tcga.luad.eset <- tmp$eset
	exome <- tmp$exome
	
	kras.snps <- list(gly12a="rs121913530",gly12="rs121913529",ala146thr="rs121913527",gly13="rs112445441",gln61="rs121913240")
	kras.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "KRAS" & exome$Variant_Classification=="Missense_Mutation"]))
	braf.v600e.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "BRAF" & exome$Start_Position==140453136]))
	braf.missense.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "BRAF" & exome$Variant_Classification=="Missense_Mutation"]))
	nras.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "NRAS" & exome$Variant_Classification=="Missense_Mutation"]))
	hras.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "HRAS"]))
	egfr.l858r.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol=="EGFR" & exome$dbSNP_RS=="rs121434568"]))
	egfr.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol=="EGFR" & grepl("Missense|In_Frame_Del", exome$Variant_Classification)]))
	
	patids <- extract.tcga.patientIds(sampleNames(tcga.luad.eset))
	
	# CCLE
	ccle <- getCCLE()
	ccle_eset <- ccle[[1]]
	ccle_drug <- ccle[[2]]
	idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
	ccle_eset <- ccle_eset[,!is.na(idxs)]
	ccle_response <- ccle_drug[na.omit(idxs),]
	lung.mask <- grepl("LUNG", sampleNames(ccle_eset))
	ccle_lung <- ccle_eset[, lung.mask]
			
	fit <- binomial_predict_EN(lungA,lungA$egfr, list(lungA, tcga.luad.eset, ccle_lung), alpha=.1, seed=20)
	pred <- prediction(fit$yhats[[1]], factor(lungA$egfr))
	auc1 <- performance(pred, 'auc')@y.values
	pred <- prediction(fit$yhats[[2]], factor(patids %in% egfr.l858r.mut.patients))
	auc2 <- performance(pred, 'auc')@y.values		
	
	summary(lm(ccle_response$Erlotinib[lung.mask] ~ ccle_lung$EGFR  + fit$yhats[[3]]))
	summary(lm(ccle_response$Erlotinib[lung.mask][ccle_lung$EGFR==0] ~ fit$yhats[[3]][ccle_lung$EGFR==0]))

	# colon ras adds nothing to prediction of erlotinib in lung
	kfsyscc_eset <- getKFSYSCC()
	y_hats <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status), 
			list(ccle_lung))$yhats
	summary(lm(ccle_response$Erlotinib[lung.mask][ccle_lung$EGFR==0] ~ y_hats[ccle_lung$EGFR==0]))
	
}

p53_in_tcga <- function(){
  tcga.luad <- get.luad.RNAseq(tumor.only=TRUE)
  tmp <- get.luad.exome(tcga.luad)
  tcga.luad.eset <- tmp$eset
  exome <- tmp$exome
  
  kras.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "KRAS" & exome$Variant_Classification=="Missense_Mutation"]))
  atm.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "ATM" & exome$Variant_Classification!="INTRON" & exome$Variant_Classification!="Silent"] ))
  mdm2.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "MDM2" & exome$Variant_Classification!="INTRON" & exome$Variant_Classification!="Silent"] ))
  tp53.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "TP53" & exome$Variant_Classification!="INTRON" & exome$Variant_Classification!="Silent"] ))
  tp53.wt.patients <- extract.tcga.patientIds(sampleNames(tcga.luad.eset))[!(extract.tcga.patientIds(sampleNames(tcga.luad.eset)) %in% tp53.mut.patients)]
  p53_factor <- factor(extract.tcga.patientIds(sampleNames(tcga.luad.eset)) %in% tp53.mut.patients)
  N <- length(p53_factor)
  fits <- mclapply(1:50, function(i){
    mask <- sample(N, N/2)
    yhats <- binomial_predict_EN(tcga.luad.eset[,mask], p53_factor[mask], list(tcga.luad.eset[,-mask]))$yhats[[1]]
    auc <- performance(prediction(yhats, p53_factor[-mask]), 'auc')@y.values
    return (list(auc=auc,yhats=yhats))
  },mc.cores=10,mc.preschedule = F)
  
  scores.mut <- matrix(NA, nrow=50,ncol=length(tp53.mut.patients),dimnames=list(1:50, tp53.mut.patients))
  for(i in 1:50){
     idxs <- match(tp53.mut.patients, extract.tcga.patientIds(rownames(fits[[i]]$yhats)))
     scores.mut[i, ] <- fits[[i]]$yhats[idxs]
  }
  scores.wt <- matrix(NA, nrow=50,ncol=length(tp53.wt.patients),dimnames=list(1:50, tp53.wt.patients))
  for(i in 1:50){
    idxs <- match(tp53.wt.patients, extract.tcga.patientIds(rownames(fits[[i]]$yhats)))
    scores.wt[i, ] <- fits[[i]]$yhats[idxs]
  }
  
  med.scores <- c(apply(scores.mut, 2, function(x) median(na.omit(x))),
                  apply(scores.wt, 2, function(x) median(na.omit(x))))
  mut.factor <- rep("WT",length(med.scores))
  mut.factor[names(med.scores) %in% tp53.mut.patients & !(names(med.scores) %in% kras.mut.patients) ] <- "p53" 
  mut.factor[!(names(med.scores) %in% tp53.mut.patients) & names(med.scores) %in% kras.mut.patients ] <- "kras" 
  mut.factor[names(med.scores) %in% tp53.mut.patients & names(med.scores) %in% kras.mut.patients ] <- "p53+kras" 
  
  pdf("plots_lung/p53_kras_boxplot_p53model.pdf",width=9,height=4)
  par(mfrow=c(1,2))
  boxplot(med.scores ~ factor(names(med.scores) %in% tp53.mut.patients),names=c("WT","p53"))
  mtext(side=3, paste("p=",format(kruskal.test(med.scores ~ factor(names(med.scores) %in% tp53.mut.patients))$p.value,digits=2)))
  boxplot(med.scores ~ factor(mut.factor))
  mtext(side=3, paste("p=",format(kruskal.test(med.scores ~ factor(mut.factor))$p.value,digits=2)))
  dev.off()
  
  
  mut.idxs <- order(apply(scores.mut, 2, function(x) median(na.omit(x))))
  wt.idxs <- order(apply(scores.wt, 2, function(x) median(na.omit(x))))
  pdf("plots_lung/p53_ris_distcompare.pdf",width=8,height=8)
  par(mfrow=c(2,1))
  boxplot(scores.mut[,mut.idxs],col=(colnames(scores.mut)[mut.idxs] %in% kras.mut.patients)+1, ylim=c(0,1),main="Dist p53 mut")
  abline(h=.5,col="red",lty=2)
  boxplot(scores.wt[,wt.idxs],col=(colnames(scores.wt)[wt.idxs] %in% kras.mut.patients)+1, ylim=c(0,1),main="Dist p53 WT")
  abline(h=.5,col="red",lty=2)
  dev.off()
   
  pdf("plots_lung/p53_mut_vs_wt_auc_boxplot.pdf")
  boxplot(unlist(sapply(fits, function(x) x$auc)))
  dev.off()

  # make genomic plot of RIS vs mutation location
  tp53.exome <- exome[exome$Hugo_Symbol == "TP53",]
  idxs <- match(extract.tcga.patientIds(tp53.exome$Tumor_Sample_Barcode), names(med.scores))
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  gr <- makeGeneRegion(start=min(tp53.exome$Start_Position)-500, end=max(tp53.exome$Start_Position)+500,chromosome=17,strand = "-",biomart=mart)
  tp53ness <- makeBaseTrack(value = med.scores[idxs],
                             base = tp53.exome$Start_Position, dp = DisplayPars(lwd=.4,color = "darkred"))
  #gene <- makeGene(id = "ENSG00000141510", type = "ensembl_gene_id", biomart = mart)
  pdf("plots_lung/tp53_genome_ris_plot.pdf",width=10,height=4)
  gdPlot(list(makeGenomeAxis(add53 = TRUE), RIS=tp53ness, gr))
  dev.off()
  
  
  ## cell line sensitivity
  ccle <- getCCLE()
  ccle_eset <- ccle[[1]]
  ccle_drug <- ccle[[2]]
  idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
  ccle_eset <- ccle_eset[,!is.na(idxs)]
  ccle_response <- ccle_drug[na.omit(idxs),]
  lung.mask <- grepl("LUNG", sampleNames(ccle_eset))
  ccle_lung <- ccle_eset[, lung.mask]
  ccle_response_lung <- ccle_response[lung.mask, ]
  
  yhats <- binomial_predict_EN(tcga.luad.eset, p53_factor, list(ccle_lung))$yhats[[1]]
  drug.pvals <- apply(pData(ccle_response_lung), 2, function(x){
     cor.test(x, yhats,method="spearman")$p.value
  })
}

ccle_lung_featureselection <- function(){
  ccle <- getCCLE()
  ccle_eset <- ccle[[1]]
  ccle_drug <- ccle[[2]]
  idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
  ccle_eset <- ccle_eset[,!is.na(idxs)]
  ccle_response <- ccle_drug[na.omit(idxs),]
  lung.mask <- grepl("LUNG", sampleNames(ccle_eset))
  ccle_lung <- ccle_eset[, lung.mask]
  
  N <- dim(ccle_lung)[2]
  
  features <- mclapply(1:1000, function(i){
    mask <- sample(N,replace=TRUE)
    fit <- binomial_predict_EN(ccle_lung[, mask],ccle_lung$KRAS[mask], 
                               list(ccle_lung[, -mask]), alpha=.1)
    fvec <- fit$featureVec[as.numeric(fit$model$beta) != 0]
    fvec
  },mc.cores=10,mc.preschedule=FALSE)
  
  count.vec <- rep(0, dim(ccle_lung)[1])
  names(count.vec) <- featureNames(ccle_lung)
  for(i in 1:length(features)){
     mask <- featureNames(ccle_lung) %in% features[[i]]
     count.vec[mask] <- count.vec[mask] + 1
  }
  
  s.count.vec <- sort(count.vec, decreasing=TRUE) / length(features)
  
  sigs <- getRASSigs()
  loboda_geneset <- union(sigs[["loboda_ras_up"]],sigs[["loboda_ras_down"]])
  r <- sapply(seq(.95,.5, by=-.05), function(x)){
     genes <- names(s.count.vec[s.count.vec > x])
     c(length(genes), length(intersect(genes, loboda_geneset)))
  })

 write.table(s.count.vec, file="lung_kras_ccle_model_bootstrapped.txt",sep="\t",quote=FALSE,col.names=FALSE)
}

ccle_model <- function(){
	ccle <- getCCLE()
	ccle_eset <- ccle[[1]]
	ccle_drug <- ccle[[2]]
	idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
	ccle_eset <- ccle_eset[,!is.na(idxs)]
	ccle_response <- ccle_drug[na.omit(idxs),]
	lung.mask <- grepl("LUNG", sampleNames(ccle_eset))
	ccle_lung <- ccle_eset[, lung.mask]
	
	tcga.luad <- get.luad.RNAseq(tumor.only=TRUE)
	tmp <- get.luad.exome(tcga.luad)
	tcga.luad.eset <- tmp$eset
	exome <- tmp$exome

	kras.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "KRAS" & exome$Variant_Classification=="Missense_Mutation"]))
	egfr.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol=="EGFR" & grepl("Missense|In_Frame_Del", exome$Variant_Classification)]))
	nf1.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "NF1" & exome$Variant_Classification=="Missense_Mutation"]))
	tp53.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "TP53" & exome$Variant_Classification!="INTRON" & exome$Variant_Classification!="Silent"] ))
	stk11.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "STK11" & exome$Variant_Classification!="INTRON" & exome$Variant_Classification!="Silent"] ))
	patids <- extract.tcga.patientIds(sampleNames(tcga.luad.eset))

	battle <- getBattle()
	battle.m <- rowMeans(exprs(battle))
	pc <- svd(exprs(battle) - battle.m)
	pc$d[1] <- 0
	battle.pc1 <- battle
	exprs(battle.pc1) <- pc$u %*% diag(pc$d) %*% t(pc$v) + battle.m 
	featureNames(battle.pc1) <- featureNames(battle)
	
	env <- new.env()
	load("data~/GSE26939/lung_gse26939.rda", env)
	lungA <- env$gene.eset
	
	chemores <- getCHEMORES()
	
	
	fit <- binomial_predict_EN(ccle_lung,ccle_lung$KRAS, list(lungA,tcga.luad.eset,battle,chemores), alpha=.1, seed=20,quantile.normalize=TRUE)
	pred1 <- prediction(fit$yhats[[1]], factor(lungA$kras))
	auc1 <- performance(pred1, 'auc')@y.values
	pred2 <- prediction(fit$yhats[[2]], factor(patids %in% kras.mut.patients))
	auc2 <- performance(pred2, 'auc')@y.values	
	pred3 <- prediction(fit$yhats[[3]], factor(battle$KRAS=="Mutant"))
	auc3 <- performance(pred3, 'auc')@y.values	
	pred4 <- prediction(fit$yhats[[4]][chemores$Histology=="AC"], factor(chemores$KRAS==1)[chemores$Histology=="AC"])
	auc4 <- performance(pred4, 'auc')@y.values	
	model.genes <- fit$featureVec[as.numeric(fit$model$beta) != 0]
	
	pdf("plots_lung/ccle_auc_model.pdf")
	plot(performance(pred1, 'tpr', 'fpr'), lty=2, col="blue",main="Predictive Model: CCLE Lung Trained")
	plot(performance(pred2, 'tpr', 'fpr'), lty=2, col="red",add=TRUE)
	plot(performance(pred3, 'tpr', 'fpr'), lty=2, col="darkgreen",add=TRUE)
	plot(performance(pred4, 'tpr', 'fpr'), lty=2, col="orange",add=TRUE)
	legend(0.5,0.6,c(paste("AUC=",format(auc1,digits=2),', GSE26939'),
					paste("AUC=",format(auc2,digits=2),', TCGA LUAD'),
					paste("AUC=",format(auc3,digits=2),', BATTLE'),
					paste("AUC=",format(auc4,digits=2),', CHEMORES')),
			col=c('blue','red','darkgreen','orange'),lwd=3)
	dev.off()
  
  # loboda compare
	loboda.1 <- loboda_enrichment_method(lungA)
	loboda.2 <- loboda_enrichment_method(tcga.luad.eset)
	loboda.3 <- loboda_enrichment_method(battle)
	loboda.4 <- loboda_enrichment_method(chemores)
	pred1 <- prediction(loboda.1, factor(lungA$kras))
	auc1 <- performance(pred1, 'auc')@y.values
	pred2 <- prediction(loboda.2, factor(patids %in% kras.mut.patients))
	auc2 <- performance(pred2, 'auc')@y.values	
	pred3 <- prediction(loboda.3, factor(battle$KRAS=="Mutant"))
	auc3 <- performance(pred3, 'auc')@y.values	
	pred4 <- prediction(loboda.4[chemores$Histology=="AC"], factor(chemores$KRAS==1)[chemores$Histology=="AC"])
	auc4 <- performance(pred4, 'auc')@y.values	
	pdf("plots_lung/loboda_auc_model.pdf")
  plot(performance(pred1, 'tpr', 'fpr'), lty=2, col="blue",main="Loboda on Lung AD")
	plot(performance(pred2, 'tpr', 'fpr'), lty=2, col="red",add=TRUE)
	plot(performance(pred3, 'tpr', 'fpr'), lty=2, col="darkgreen",add=TRUE)
	plot(performance(pred4, 'tpr', 'fpr'), lty=2, col="orange",add=TRUE)
	legend(0.5,0.6,c(paste("AUC=",format(auc1,digits=2),', GSE26939'),
	                 paste("AUC=",format(auc2,digits=2),', TCGA LUAD'),
	                 paste("AUC=",format(auc3,digits=2),', BATTLE'),
	                 paste("AUC=",format(auc4,digits=2),', CHEMORES')),
	       col=c('blue','red','darkgreen','orange'),lwd=3)
	dev.off()
  #######################################
	# see where EGFR falls on the spectrum
	tcga.factor <- rep("WT", length(patids))
	tcga.factor[patids %in% kras.mut.patients & !(patids %in% tp53.mut.patients)] <- "KRAS"
	tcga.factor[patids %in% kras.mut.patients & patids %in% tp53.mut.patients] <- "KRAS+TP53"
	#tcga.factor[patids %in% egfr.mut.patients] <- "EGFR"
	tcga.factor[patids %in% tp53.mut.patients & !(patids %in% kras.mut.patients)] <- "TP53"
	#tcga.factor[patids %in% stk11.mut.patients & !(patids %in% kras.mut.patients)] <- "STK1"

	#tcga.factor[patids %in% nf1.mut.patients] <- "NF1"
  pdf("plots_lung/boxplot_kras_p53_wt.pdf")
  boxplot(fit$yhats[[2]] ~ factor(tcga.factor),ylab="RIS",main="LUAD: RIS vs mutation status")
  dev.off()
	
  run.kruskal.compare <- function(factors){
    filter <- tcga.factor %in% c(factors)
    kruskal.test(fit$yhats[[2]][filter], factor(tcga.factor)[filter])$p.value
  }
  
  kras_vs_wt <- run.kruskal.compare(c("KRAS","WT"))
  p53_vs_wt <- run.kruskal.compare(c("TP53","WT"))
	p53_vs_kras <- run.kruskal.compare(c("TP53","KRAS"))
	mnames <- c("KRAS","TP53","WT")
  tmp <- matrix("",nrow=3,ncol=3, dimnames=list(mnames, mnames))
	tmp[1,2] <- format(p53_vs_kras,digits=2)
	tmp[1,3] <- format(kras_vs_wt,digits=2)
	tmp[2,3] <- format(p53_vs_wt,digits=2)
  tmp[1,1] <- paste("n=",sum(tcga.factor=="KRAS"))
	tmp[2,2] <- paste("n=",sum(tcga.factor=="TP53"))
	tmp[3,3] <- paste("n=",sum(tcga.factor=="WT"))
  pdf("plots_lung/kras_vs_tp53_vs_wt_kruskal.pdf")
  textplot(tmp)
  title("LUAD: RIS comparison (Mann-Whitney)")
  dev.off()
  
  tmp <- exprs(tcga.luad.eset[featureNames(tcga.luad.eset) %in% model.genes, ])
	pc <- svd(tmp - rowMeans(tmp))
	plot(pc$v[,1], pc$v[,2],col=as.numeric(factor(tcga.factor)))
	
	
	chemores.factor <- rep("WT", dim(chemores)[2])
	chemores.factor[chemores$KRAS==1] <- "KRAS"
	chemores.factor[chemores$EGFR!="WT"] <- "EGFR"
	boxplot(fit$yhats[[4]]  ~ factor(chemores.factor))
	
	#######################################
	# gsva
	gsets <- load.gmt.data("resources/c2.all.v3.0.symbols.gmt")
	tcga.es <- gsva(tcga.luad.eset,gsets,min.sz=10,max.sz=500,parallel.sz=10)$es
	tcga.kras.factor <- factor(patids %in% kras.mut.patients)
	tcga.pvals <- apply(exprs(tcga.es), 1, function(x){
				coef(summary(lm(fit$yhats[[2]] ~ x + tcga.kras.factor)))[2,4]
			})
	
	chemores.es <- gsva(chemores,gsets,min.sz=10,max.sz=500,parallel.sz=10)$es
	chemores.kras.factor <- factor(chemores$KRAS==1)
	chemores.pvals <- apply(exprs(chemores.es), 1, function(x){
				coef(summary(lm(fit$yhats[[4]] ~ x + chemores.kras.factor)))[2,4]
			})
	
}

combined.gsea.analysis <- function(){
	
	tcga.luad <- get.luad.RNAseq(tumor.only=TRUE)
	tmp <- get.luad.exome(tcga.luad)
	tcga.luad.eset <- tmp$eset
	exome <- tmp$exome
	kras.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "KRAS" & exome$Variant_Classification=="Missense_Mutation"]))
	patids <- extract.tcga.patientIds(sampleNames(tcga.luad.eset))
	
	battle <- getBattle()
	battle.m <- rowMeans(exprs(battle))
	pc <- svd(exprs(battle) - battle.m)
	pc$d[1] <- 0
	battle.pc1 <- battle
	exprs(battle.pc1) <- pc$u %*% diag(pc$d) %*% t(pc$v) + battle.m 
	featureNames(battle.pc1) <- featureNames(battle)
	
	env <- new.env()
	load("data~/GSE26939/lung_gse26939.rda", env)
	lungA <- env$gene.eset
	
	chemores <- getCHEMORES()
	source("code/multigsea.R")	
	fit.1 <- eBayes(lmFit(lungA,model.matrix(~factor(lungA$kras))))
	fit.2 <- eBayes(lmFit(tcga.luad.eset, model.matrix(~factor(patids %in% kras.mut.patients))))
	fit.3 <- eBayes(lmFit(battle.pc1, model.matrix(~factor(battle$KRAS=="Mutant"))))
	fit.4 <- eBayes(lmFit(chemores, model.matrix(~factor(chemores$KRAS==1))))
	gsets <- load.gmt.data("resources/c2.all.v3.0.symbols.gmt")
	
	scores <- list(-log10(fit.1$p.value[,2]), -log10(fit.2$p.value[,2]), -log10(fit.3$p.value[,2]), -log10(fit.4$p.value[,2]))
	es <- ES(scores, gsets)
	
	null.es <- mclapply(1:1000, function(i){
		fit.1 <- eBayes(lmFit(lungA,model.matrix(~factor(lungA$kras[sample(dim(lungA)[2])]))))
		fit.2 <- eBayes(lmFit(tcga.luad.eset, model.matrix(~factor(patids %in% kras.mut.patients)[sample(length(patids))])))
		fit.3 <- eBayes(lmFit(battle.pc1, model.matrix(~factor(battle$KRAS=="Mutant")[sample(dim(battle)[2])])))
		fit.4 <- eBayes(lmFit(chemores, model.matrix(~factor(chemores$KRAS==1)[sample(dim(chemores)[2])])))
		scores <- list(-log10(fit.1$p.value[,2]), -log10(fit.2$p.value[,2]), -log10(fit.3$p.value[,2]), -log10(fit.4$p.value[,2]))
		es <- ES(scores, gsets)
		es
	},mc.set.seed=TRUE)
	null.m <- matrix(NA, nrow=length(gsets), ncol=length(null.es),dimnames=list(names(gsets),1:length(null.es)))
	for(i in 1:length(null.es)){null.m[,i] <- null.es[[i]] }
	pvals <- sapply(1:length(gsets), function(i){
		sum(null.m[i,] > es[i]) / length(null.es)
	})
	names(pvals) <- names(gsets)
	

	
	
	
}
