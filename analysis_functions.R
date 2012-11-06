
library(parallel)
library(GSVA)
library(IlluminaHumanMethylation27k.db)
library(limma)
library(edgeR)
library(cqn)
library(survival)
library(hgu133plus2.db)
library(org.Hs.eg.db)

source("code/util_functions.R")
source("code/data_functions.R")
source("code/model_functions.R")

# global variables
canonical.kras <- c("kras.12","kras.13","kras.61","kras.146")
canonical.braf <- c("braf.600")
canonical.hras <- c("hras.G12","hras.G13")
canonical.nras <- c("nras.G12","nras.G13", "nras.Q61")
canonical.pik3ca <- c("pik3ca.542","pik3ca.545","pik3ca.546","pik3ca.1047")

raf_resistance <- function(){
  eset <- getGEO("GSE34299")
  eset <- eset$GSE34299_series_matrix.txt.gz
  pheno <- eset$characteristics_ch1
  genes <- unlist(mget(featureNames(eset), hgu133plus2SYMBOL,ifnotfound=NA))
  tmp <- combine_probes_2_gene(exprs(eset), genes)
  has.na <- apply(tmp, 1, function(x){any(is.na(x))})
  raf_eset <- new("ExpressionSet", exprs=tmp[!has.na,])
  
  kfsyscc_eset <- getKFSYSCC()
  
  yhat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status == "MUT"), 
                              list(raf_eset))$yhats[[1]]
  
  plot(1:4, yhat)
}

kfsyscc_survival <- function(){
 
  eset <- getKFSYSCC()
  N <- ncol(eset)
  
  splitAndTest <- function(i){
    idxs <- sample(N, round(N/1.4))
    kras <- factor(eset$kras_status)
    r <- binomial_predict_EN(eset[,idxs], kras[idxs], list(eset[,-idxs]))
    y_hat <- r$yhats[[1]]
    y_hat
  }
  rlist <- mclapply(1:20, splitAndTest,mc.cores=5,mc.set.seed=TRUE,mc.preschedule =FALSE)
  M <- matrix(NA, nrow=N,ncol=20,dimnames = list(sampleNames(eset), 1:20))
  for(i in 1:length(rlist)){
     idxs <- match(rownames(rlist[[i]]), rownames(M))
     M[idxs, i] <- rlist[[i]]
  }
  yhat <- rowMedians(M,na.rm=TRUE)
  surv <- Surv(kfsyscc_eset$FT, kfsyscc_eset$Death)
  rasMask <- kfsyscc_eset$kras_status =="MUT"
  metastatic <- kfsyscc_eset$StageM == 1
  
  summary(cox.fit <- coxph(surv ~ factor(kras_status),data=pData(kfsyscc_eset)))
  cox.fit <- coxph(surv ~ yhat + factor(kfsyscc_eset$kras_status))
  cox.fit <- coxph(surv[!metastatic ] ~ yhat[!metastatic ] )
  summary(cox.fit)
  plot(survfit(surv[metastatic] ~ rasMask[metastatic]))
}


uterineAnalysis <- function(){
  uterine.eset.eg <- loadEntity('syn585089')$objects$eset
  entrez.ids <- gsub("(\\d*)_eg","\\1", featureNames(uterine.eset.eg))
  genes <- unlist(mget(entrez.ids, org.Hs.egSYMBOL, ifnotfound=NA))
  
  tmp <- combine_probes_2_gene(exprs(uterine.eset.eg), genes)
  uterine.eset <- new("ExpressionSet",exprs=tmp)
  
  maf.entity <- loadEntity('syn350455')
  maf <- read.table(paste(maf.entity$cacheDir ,"/",maf.entity$files,sep=""),header=TRUE,sep="\t",quote="")
  all.patients <- extract.tcga.patientIds(unique(maf$Tumor_Sample_Barcode))
  kras.patients <- extract.tcga.patientIds(unique(maf[maf$Hugo_Symbol=="KRAS",]$Tumor_Sample_Barcode))
  
  uterine.eset.m <- uterine.eset[, extract.tcga.patientIds(sampleNames(uterine.eset)) %in% all.patients]
  
  yhat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status == "MUT"), 
                              list(uterine.eset.m))$yhats[[1]]
  pred <- prediction(yhat, extract.tcga.patientIds(sampleNames(uterine.eset.m)) %in% kras.patients)
  auc <- performance(pred, 'auc')@y.values[[1]]
}

crcCellLineToTumor <- function(){
  
  ccle <- getCCLE()
  ccle_eset <- ccle[[1]]
  colon_mask <- grepl("LARGE_INTESTINE", sampleNames(ccle_eset))
  
  kfsyscc_eset <- getKFSYSCC()
  khambata_eset <- getKhambata()
  gaedcke_eset <- getGaedcke()
  tcga_eset <- getTCGACRC()
  tcga_eset <- tcga_eset[,!is.na(tcga_eset$kras)]
  
  
  y_hats_all <- binomial_predict_EN(ccle_eset, factor(ccle_eset$KRAS==1), 
                                list(kfsyscc_eset, tcga_eset, khambata_eset, gaedcke_eset))$yhats
  y_hats_crc <- binomial_predict_EN(ccle_eset[,colon_mask], factor(ccle_eset$KRAS[colon_mask]==1), 
                                    list(kfsyscc_eset, tcga_eset, khambata_eset, gaedcke_eset))$yhats
  
 
  auc_and_perf <- function(yhat, actual){
    pred <- prediction(yhat, actual)
    perf <- performance(pred, 'tpr', 'fpr')
    auc <- performance(pred, 'auc')@y.values[[1]]
    
    list(pred=pred, perf=perf, auc=auc)
  }
  
  kfsyscc_perf <- auc_and_perf(y_hats_all[[1]], factor(kfsyscc_eset$kras_status=="MUT"))
  tcga_perf <- auc_and_perf(y_hats_all[[2]], factor(tcga_eset$kras==1))
  khambata_perf <- auc_and_perf(y_hats_all[[3]], factor(khambata_eset$kras_status=="MUT"))
  gaedcke_perf <- auc_and_perf(y_hats_all[[4]], factor(gaedcke_eset$kras_status=="MUT"))
  
  kfsyscc_perf <- auc_and_perf(y_hats_crc[[1]], factor(kfsyscc_eset$kras_status=="MUT"))
  tcga_perf <- auc_and_perf(y_hats_crc[[2]], factor(tcga_eset$kras==1))
  khambata_perf <- auc_and_perf(y_hats_crc[[3]], factor(khambata_eset$kras_status=="MUT"))
  gaedcke_perf <- auc_and_perf(y_hats_crc[[4]], factor(gaedcke_eset$kras_status=="MUT"))
  
  
  
}

crcXenoGraft <- function(){
  xenoData <- getEMEXP991()
  kfsyscc_eset <- getKFSYSCC()
  
  eset <- xenoData$eset
  annot <- xenoData$annot
  fisher.test(factor(annot$Cetuximab==3), factor(annot$Kras == "MUT" | annot$Braf=="MUT"))
  
  # QC eset
  model.names <- gsub("([A-Z]{2}?-[A-Z]{2,3}?-[0-9A-Z]{4}).*","\\1",eset$Source.Name)
  xeno.passage <- as.numeric(gsub("[A-Z]{2}?-[A-Z]{2,3}?-[0-9A-Z]{4}-P([0-9]).*","\\1",eset$Source.Name))
  is.xeno <- !is.na(xeno.passage)
  
  
  pc <- svd(exprs(eset) - rowMeans(exprs(eset)))
  
  # remove top PC
  adj.eset <- eset
  exprs(adj.eset) <- pc$u %*% diag(c(0,pc$d[-c(1)])) %*% t(pc$v)
  featureNames(adj.eset) <- featureNames(eset)
  pc.adj <- svd(exprs(adj.eset) - rowMeans(exprs(adj.eset)))
  
  pdf("plots/xeno_correction.pdf")
  par(mfrow=c(2,2))
  plot(pc$d^2 / sum(pc$d^2),ylab="%var",xlab="eigenrank",main="Pre adjust")
  plot(pc$v[,1], pc$v[,2], col=(as.numeric(is.xeno)+1),xlab="PC1",ylab='PC2')
  
  plot(pc.adj$d^2 / sum(pc.adj$d^2),ylab="%var",xlab="eigenrank",main="Post adjust")
  plot(pc.adj$v[,1], pc.adj$v[,2], col=(as.numeric(is.xeno)+1),xlab="PC1",ylab='PC2')
  dev.off()
  
  yhat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status == "MUT"), 
                              list(adj.eset))$yhats[[1]]
  
  #yhat <- loboda_enrichment_method(adj.eset)

   
  ## all data
  idxs <- match(model.names, annot$Code)
  annot.model <- annot[na.omit(idxs),]
  xeno.model <- xeno.passage[!is.na(idxs)]
  early.model.mask <- xeno.model < 4 & !is.na(xeno.model)
  late.model.mask <- xeno.model >= 4 & !is.na(xeno.model)
  yhat.model <- yhat[!is.na(idxs)]
  
  colors = rep("azure4",nrow(annot.model))
  colors[annot.model$Kras == "MUT" & annot.model$Braf == "WT"] <- "red"
  colors[annot.model$Kras == "WT" & annot.model$Braf == "MUT"] <- "darkgoldenrod1"
  colors[annot.model$Kras == "MUT" & annot.model$Braf == "MUT"] <- "green"
  
  kras_braf_factor <- factor(annot.model$Kras == "MUT" | annot.model$Braf=="MUT")
  response_factor <- factor(annot.model$Cetuximab==3)
  fisher.test(factor(annot$Cetuximab==3), kras_braf_factor)
  
  pdf("plots/xenomodel_primaryandmouse_cetux.pdf",width=8,height=6)
  par(mfrow=c(2,2), oma=c(0,0,0,8))
  model <- yhat.model[is.na(xeno.model)] ~ factor(annot.model$Cetuximab==3)[is.na(xeno.model)]
  mut_pval <- fisher.test(response_factor[is.na(xeno.model)], kras_braf_factor[is.na(xeno.model)])$p.value
  boxplot(model,ylab="RIS",names=c("Non-response","Response"),
          pars=list(outpch=NA,medlwd=.5),main=paste("tumor, n=",sum(is.na(xeno.model)),sep=""))
  points(model,pch=19,cex=.9,col=colors[!is.na(xeno.model)])
  ptext <- ifelse(mut_pval > .5,  " (p>.5)",paste(" (p=",format(mut_pval,digits=2),")",sep=""))
  mtext(side=3,paste("p=",format(kruskal.test(model)$p.value,digits=2),ptext,sep=""),cex=.8)
  
  
  model <- yhat.model[early.model.mask] ~ factor(annot.model$Cetuximab==3)[early.model.mask]
  mut_pval <- fisher.test(response_factor[early.model.mask], kras_braf_factor[early.model.mask])$p.value
  boxplot(model,ylab="RIS",names=c("Non-response","Response"),
          pars=list(outpch=NA,medlwd=.5),main=paste("xeno early, n=",sum(early.model.mask),sep=""))
  points(model,pch=19,cex=.9,col=colors[early.model.mask])
  ptext <- ifelse(mut_pval > .5,  " (p>.5)",paste(" (p=",format(mut_pval,digits=2),")",sep=""))
  mtext(side=3,paste("p=",format(kruskal.test(model)$p.value,digits=2),ptext,sep=""),cex=.8)
  
  
  model <- yhat.model[late.model.mask] ~ factor(annot.model$Cetuximab==3)[late.model.mask]
  mut_pval <- fisher.test(response_factor[late.model.mask], kras_braf_factor[late.model.mask])$p.value
  boxplot(model,ylab="RIS",names=c("Non-response","Response"),
          pars=list(outpch=NA,medlwd=.5),main=paste("xeno late, n=",sum(late.model.mask),sep=""))
  points(model,pch=19,cex=.9,col=colors[late.model.mask])
  ptext <- ifelse(mut_pval > .5,  " (p>.5)",paste(" (p=",format(mut_pval,digits=2),")",sep=""))
  mtext(side=3,paste("p=",format(kruskal.test(model)$p.value,digits=2),ptext,sep=""),cex=.8)
  
  
  model <- yhat.model ~ factor(annot.model$Cetuximab==3)
  mut_pval <- fisher.test(response_factor, kras_braf_factor)$p.value
  boxplot(model,ylab="RIS",names=c("Non-response","Response"),
          pars=list(outpch=NA,medlwd=.5),ylim=c(0,1),main=paste("tumor + xeno, n=",length(yhat.model),sep=""))
  points(model,pch=19,cex=.9,col=colors)
  ptext <- ifelse(mut_pval > .5, " (p>.5)",paste(" (p=",format(mut_pval,digits=2),")",sep=""))
  mtext(side=3,paste("p=",format(kruskal.test(model)$p.value,digits=2),ptext,sep=""),cex=.8)
  
  
  
  par(xpd=NA)
  legend(par("usr")[2] + .2, mean(par("usr")[3:4]), legend=c("kras","braf","kras+braf","wt"),
         fill=c("red","darkgoldenrod1","green","azure4"))
  par(xpd=FALSE)
  dev.off()
  
  
 
  boxplot(yhat.model ~ factor(annot.model$Cetuximab),ylab="RIS",main="Cetuximab response",names=c(">42%","10-42%","-10-10%", "<-10%"))
  mtext(side=3,paste("p=",format(kruskal.test(yhat.model ~ factor(annot.model$Cetuximab))$p.value,digits=2),sep=""))
  
  summary(glm(factor(annot.model$Cetuximab==3) ~ factor(annot.model$Kras == "MUT" | annot.model$Braf),family="binomial"))
  summary(glm(factor(annot.model$Cetuximab==3) ~ factor(annot.model$Kras == "MUT" | annot.model$Braf=="MUT" | annot.model$Pik3ca=="MUT") + (yhat.model > .35),family="binomial"))
  
}

irontecanResponse2 <- function(){
	
  eset <- getEMEXP3549()
  kfsyscc_eset <- getKFSYSCC()
  
  yhat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status), 
                              list(eset))$yhats[[1]]
  
  mut.mask <- grepl("wild-type",eset$Factor.Value.ClinicalHistory.)
  boxplot(yhat ~ factor(eset$Factor.Value..ResponseTo5.FU.irinotecanTreatment. ))
  boxplot(yhat ~ mut.mask)
}

mekInhibition <- function(){
	
	tmp <- getEMEXP3557()
	kfsyscc_eset <- getKFSYSCC()
  
  # 1 array is major outlier
	pc <- svd(exprs(tmp) - rowMeans(exprs(tmp)))
	bad.idxs <- which(pc$v[,1] > .8)
  mekEset <- tmp[,-bad.idxs]
	
	yhat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status == "MUT"), 
	                            list(mekEset))$yhats[[1]]
  
	mut.mask <- mekEset$Factor.Value.GENOTYPE == "Kras mutant"
	cntrl.mask <- mekEset$Factor.Value.COMPOUND. == "control"
  plot(yhat, type="p",pch=as.numeric(factor(mekEset$Factor.Value.GENOTYPE.)),
       col=as.numeric(factor(mekEset$Factor.Value.COMPOUND.)))
  

  
  #leveneTest(yhat[!mut.mask] ~ cntrl.mask[!mut.mask])
	boxplot(yhat ~ factor(mekEset$Factor.Value.COMPOUND.),main="Kras mut",ylab="RIS")
	points(factor(mekEset$Factor.Value.COMPOUND.)[mut.mask], yhat,col=c(as.numeric(factor(mekEset$Factor.Value.COMPOUND.)) + 1))
	mtext(side=3, paste("p=",format(kruskal.test(yhat ~ factor(mekEset$Factor.Value.COMPOUND.))$p.value,digits=2),sep=""))
	
 
  
  
  pdf("plots/mek_mouse_xenograft.pdf",width=6,height=4)
  par(mfrow=c(1,2))
	boxplot(yhat[mut.mask] ~ factor(mekEset$Factor.Value.COMPOUND.)[mut.mask],main="Kras mut",ylab="RIS",ylim=c(.25, .55))
  points(factor(mekEset$Factor.Value.COMPOUND.)[mut.mask], 
         yhat[mut.mask],col=c(as.numeric(factor(mekEset$Factor.Value.COMPOUND.)[mut.mask]) + 1),pch=19)
  mtext(side=3, paste("p=",format(kruskal.test(yhat[mut.mask] ~ factor(mekEset$Factor.Value.COMPOUND.)[mut.mask])$p.value,digits=2),sep=""))
	boxplot(yhat[!mut.mask] ~ factor(mekEset$Factor.Value.COMPOUND.)[!mut.mask],main="Kras wt",ylab="RIS",ylim=c(.25, .55))
	points(factor(mekEset$Factor.Value.COMPOUND.)[!mut.mask], 
         yhat[!mut.mask],col=c(as.numeric(factor(mekEset$Factor.Value.COMPOUND.)[!mut.mask]) + 1),pch=19)
  mtext(side=3, paste("p=",format(kruskal.test(yhat[!mut.mask] ~ factor(mekEset$Factor.Value.COMPOUND.)[!mut.mask])$p.value,digits=2),sep=""))
	dev.off()
  
	#boxplot(yhat[cntrl.mask] ~ factor(mekEset$Factor.Value.GENOTYPE)[cntrl.mask])
}

irontecanResponse <- function(){
	eset <- getEMTAB333()
	
	kfsyscc_eset <- getKFSYSCC()
	
	yhat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status), 
			list(eset))$yhats[[1]]
	flofiri.mask <- eset$Factor.Value..treatment=="folfiri"
	table(yhat[flofiri.mask] < .4, factor(eset.m$Factor.Value..response)[flofiri.mask])
	table(yhat > .4, factor(eset$Factor.Value..response))
}

notchSignaling <- function(){

	eset <- getGEO("GSE37645",GSEMatrix=TRUE)
	eset <- eset$GSE37645_series_matrix.txt.gz
	genes <- mget(featureNames(eset), hgu133plus2SYMBOL,ifnotfound=NA)
	M <- combine_probes_2_gene(exprs(eset), genes)
	g.eset <- eset
	exprs(g.eset) <- M
	notch_sensitive <- !grepl("Non-sensitive", g.eset$title)
	
	kfsyscc_eset <- getKFSYSCC()
	
	yhat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status), 
			list(g.eset))$yhats[[1]]
	kruskal.test(yhat ~ notch_sensitive)
	
}

getRISandMut <- function(){
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
	
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rasness[,1]), rownames(tcga_aa_mat))
	rasScore.m <- rasness[!is.na(idxs),2]
	names(rasScore.m) <- rasness[!is.na(idxs),1]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	idxs <- pmatch(canonical, colnames(tcga_aa_mat.m))	
	has.mut <- rowSums(tcga_aa_mat.m[,idxs]) > 0
	return(list(ris=rasScore.m, mut=has.mut))
}

foo <- function(){
	genes <- read.table("resources/chr_symbols.txt",sep="\t",header=F,as.is=T)
	pvals <- read.table("resources/CNA_conditionalpval_rasness_coad_read.txt",sep="\t",header=F,as.is=T)
	genes.chr12 <- genes[genes[,1] == "hs12",]
	pvals.chr12 <- pvals[pvals[,1] == "hs12",]
	mask <- pvals.chr12[,4] > 7
	ir.1 <- IRanges(as.numeric(pvals.chr12[mask,2]), 
			as.numeric(pvals.chr12[mask,3]))
	ir.2 <- IRanges(as.numeric(genes.chr12[,2]), as.numeric(genes.chr12[,3]))
	
	idxs <- findOverlaps(ir.2, ir.1,select="first")
	sel.genes <- genes.chr12[!is.na(idxs),]
	write.table(sel.genes, file="resources/top_chr12_cna_ras_genes.txt",sep="\t",quote=F,col.names=F,row.names=F)
}

tcga.ras.mirna.diff.coexpr <- function(){
  r <- getRISandMut()
  load("data~/tcga_crc_mirnaseq.rda")
  eset <- getTCGACRC()
  idxs <- merge.across(extract.tcga.patientIds(names(r$ris)),
                extract.tcga.patientIds(colnames(tcga.crc.mirnaseq)),
                       extract.tcga.patientIds(sampleNames(eset)))
  rasScore.m <- r$ris[idxs[,1]]
  has.mut.m <- r$mut[idxs[,1]]
  miRna.m <- tcga.crc.mirnaseq[,idxs[,2]]
  eset.m <- eset[,idxs[,3]]
  
  
  # differential miRNA expression based on ras mut status
  pvals <- apply(log(miRna.m+1), 1, function(x){
    mask <- !is.na(x)
    kruskal.test(x[mask] ~ has.mut.m[mask])$p.value
  })	
  
  #top differentially expressed
  boxplot(log(miRna.m["hsa-mir-582", ]) ~ has.mut.m,main="hsa-mir-582",names=c("RAS WT","RAS MUT"))
  mtext(side=3, format(kruskal.test(log(miRna.m["hsa-mir-582", ]) ~ has.mut.m)$p.value,digits=2))
  
  C.mut <- cor(t(exprs(eset.m[, has.mut.m])), t(miRna.m[, has.mut.m]),method="spearman")
  C.wt <- cor(t(exprs(eset.m[, !has.mut.m])), t(miRna.m[, !has.mut.m]),method="spearman")
  C.diff <- C.mut^2 - C.wt^2
  foo <- C.mut - C.wt
  sort.idxs <- order(abs(C.diff),decreasing=TRUE)
  tbl <- buildTopTableFromCorrelationMatrix( C.mut - C.wt, sort.idxs,top=1000)
  gene.counts.by.mir <- sort(table(tbl[,1]),decreasing=T)
  gene.counts.by.mir <- gene.counts.by.mir[gene.counts.by.mir > 10]
  #hsa-mir-146a has most
  barplot(gene.counts.by.mir,las=3,ylab="Number of genes highly differentially correlated",
          main="Diff correlation analsyis: kras mut vs wt")

  mir146a.regulated.genes <- sort(C.mut[,"hsa-mir-146a"],decreasing=T)
  mir146b.regulated.genes <- sort(C.mut[,"hsa-mir-146b"],decreasing=F)
  
  # MINDY style test: Effector -> Modulator -> Target
}

tcga.ras.mirna <- function(){
	r <- getRISandMut()
	load("data~/tcga_crc_mirnaseq.rda")
	idxs <- match(extract.tcga.patientIds(names(r$ris)),
			extract.tcga.patientIds(colnames(tcga.crc.mirnaseq)))
	rasScore.m <- r$ris[!is.na(idxs)]
	has.mut.m <- r$mut[!is.na(idxs)]
	miRna <- tcga.crc.mirnaseq[, na.omit(idxs)]
	
	# differential miRNA expression based on mut status
	pvals <- apply(log(miRna+1), 1, function(x){
				mask <- !is.na(x)
				kruskal.test(x[mask] ~ has.mut.m[mask])$p.value
			})	
  
  
	pvals <- t(apply(miRna, 1, function(x){
		if(all(x == 0)){ return (c(NA,NA)) }
		p1 <- coef(summary(lm(rasScore.m ~ log(x+1))))[2,4]	
		p2 <- coef(summary(lm(rasScore.m ~ log(x+1) + has.mut.m)))[2,4]
		c(p1,p2)
	}))

	pvals.mut.only <- apply(miRna[, has.mut.m], 1, function(x){
					if(all(x == 0)){ return (NA) }
					p1 <- coef(summary(lm(rasScore.m[has.mut.m] ~ log(x+1))))[2,4]	
					p1
				})
		
	pvals.wt.only <- apply(miRna[, !has.mut.m], 1, function(x){
				if(all(x == 0)){ return (NA) }
				p1 <- coef(summary(lm(rasScore.m[!has.mut.m] ~ log(x+1))))[2,4]	
				p1
			})

  pdf("plots/mirna_146b_vs_ris.pdf",width=8,height=4)
	par(mfrow=c(1,2))
	plot(miRna["hsa-mir-146b",has.mut.m], rasScore.m[has.mut.m],main="hsa-mir-146b: kras mut",xlab="miRna",ylab="RIS")
	lines(lowess(miRna["hsa-mir-146b",has.mut.m], rasScore.m[has.mut.m]), col = 2,lty=2)
  mtext(paste("p=",format(
    cor.test(miRna["hsa-mir-146b",has.mut.m], rasScore.m[has.mut.m], method="spearman")$p.value,digits=2),"(spearman)"),3)
	plot(miRna["hsa-mir-146b",!has.mut.m], rasScore.m[!has.mut.m],main="hsa-mir-146b: kras wt",xlab="miRna",ylab="RIS")
	lines(lowess(miRna["hsa-mir-146b",!has.mut.m], rasScore.m[!has.mut.m]), col = 2,lty=2)
  mtext(paste("p=",format(
	  cor.test(miRna["hsa-mir-146b",!has.mut.m], rasScore.m[!has.mut.m], method="spearman")$p.value,digits=2),"(spearman)"),3)
	
  dev.off()
  
	eset <- getTCGACRC()
	idxs <- match(extract.tcga.patientIds(sampleNames(eset)),
					extract.tcga.patientIds(colnames(miRna)))
	
	eset.m <- eset[,!is.na(idxs)]
	tmp <- miRna[, na.omit(idxs)]
	tmp.mut.mask <- has.mut.m[na.omit(idxs)]
	
	C.mut <- cor(t(exprs(eset.m[, tmp.mut.mask])), t(tmp[, tmp.mut.mask]),method="spearman")
	C.wt <- cor(t(exprs(eset.m[, !tmp.mut.mask])), t(tmp[, !tmp.mut.mask]),method="spearman")
	
	
	mir146b.regulated.genes <- sort(C.mut[,"hsa-mir-146b"],decreasing=T)
	
	C <- cor(t(exprs(eset.m)), t(tmp),method="spearman")
	mir181.regulated.genes <- sort(C[,"hsa-mir-181a-2"],decreasing=F)
	mir181.regulated.genes <- sort(C[,"hsa-mir-181a-2"],decreasing=F)
	
	pdf("plots/hist_mir181a2_rnaseq.pdf",width=5,height=5)
	hist(tmp["hsa-mir-181a-2",],breaks=20,main="MIR-181A-2",xlab="RNAseq counts")
	dev.off()
	
	
	pdf("plots/density_mir181a2_geneexpr_correlation.pdf",width=5,height=5)
	plot(density(mir181.regulated.genes),main="MIR-181A-2",xlab="Correlation (spearman) with gene expression")
	abline(v=mir181.regulated.genes[1],col="red")
	text(-0.42, 3, "FOXO1",adj = c(0, 0.5))
	dev.off()
	
	colnames(pvals) <- c("marginal","conditional")
	pvals
}

tcga.ras.gsva <- function(){
	r <- getRISandMut()
	
	gsets <- load.gmt.data("resources/c2.all.v3.0.symbols.gmt")
	mir.gsets <- load.gmt.data("resources/c3.mir.v3.0.symbols.gmt")
	tft.gsets <- load.gmt.data("resources/c3.tft.v3.0.symbols.gmt")
	eset <- getTCGACRC()
	es <- gsva(eset,gsets,min.sz=10,max.sz=500,parallel.sz=10)$es
	

	
	idxs <- match(sampleNames(es), names(r$ris))
	es.m <- exprs(es[, !is.na(idxs)])
	rasScore.m <- r$ris[na.omit(idxs)]
	has.mut.m <- r$mut[na.omit(idxs)]
	
	diff.expr <- eBayes(lmFit(es.m, model.matrix(~has.mut.m)))
	#differential pathway activation, mut vs WT
	sort(diff.expr$p.value[,2])[1:20]
	
	pvals <- apply(es.m, 1, function(x){
		coef(summary(lm(rasScore.m ~ x + has.mut.m)))[2,4]
	})
	#
	pval.corrected <- p.adjust(pvals,method="bonferroni")
	# most significant geneset is SABATES_COLORECTAL_ADENOMA_UP
	# http://mcr.aacrjournals.org/content/5/12/1263.long
	# mostly corresponds to WNT signaling
	x <- es.m["SABATES_COLORECTAL_ADENOMA_UP",]
	pdf("plots/sabates_adenoma_pathway_vs_rasness.pdf")
	plot(x, rasScore.m, col=(has.mut.m+1),pch=19,cex=.8,xlab="Sabates CRC Adenoma Enrichment",ylab="RIS")
	lines(lowess(rasScore.m ~ x), col="green", lwd=3,lty=2)
	legend("topleft",c("RAS WT","RAS MUT"), col=c("black","red"), pch=19)
	dev.off()
	
	# mut only
	idxs <- match(sampleNames(eset), names(r$ris))
	eset.m <- exprs(eset[, !is.na(idxs)])
	rasScore.m <- r$ris[na.omit(idxs)]
	has.mut.m <- r$mut[na.omit(idxs)]
	
	
	
	es.mut.only <- gsva(eset.m[ ,has.mut.m],gsets,min.sz=10,max.sz=500,parallel.sz=10)$es
	pvals.mut.only <- apply(es.mut.only, 1, function(x){
				cor.test(rasScore.m[has.mut.m], x)$p.value
			})
	mir.mut.only <- gsva(eset.m[ ,has.mut.m],mir.gsets,min.sz=10,max.sz=500,parallel.sz=10)$es
	pvals.mir.mut.only <- apply(mir.mut.only, 1, function(x){
				cor.test(rasScore.m[has.mut.m], x)$p.value
			})
	tft.mut.only <- gsva(eset.m[ ,has.mut.m],tft.gsets,min.sz=10,max.sz=500,parallel.sz=10)$es
	pvals.tft.mut.only <- apply(tft.mut.only, 1, function(x){
				cor.test(rasScore.m[has.mut.m], x)$p.value
			})
	# wt only
	es.wt.only <- gsva(eset.m[ ,!has.mut.m],gsets,min.sz=10,max.sz=500,parallel.sz=10)$es
	pvals.wt.only <- apply(es.wt.only, 1, function(x){
				cor.test(rasScore.m[!has.mut.m], x)$p.value
			})
	
	#mask <- grepl("REACTOME|KEGG|BIOCARTA",rownames(es.mut.only))
	## differential correlation
	C.mut <- cor(t(es.mut.only), method="spearman")
	diag(C.mut) <- 0
	C.wt <- cor(t(es.wt.only), method="spearman")
	diag(C.wt) <- 0
	C.diff <- C.mut^2 - C.wt^2
	
	geneset.names <- rownames(es.mut.only)
	printTop <- function(C, idxs, top=20){
		N <- dim(C)[2]
		#idxs <- order(abs(C.mut),decreasing=TRUE)
		for(i in 1:top){
			c <- floor(idxs[i] / N)
			r <- idxs[i] %% N
			cat(geneset.names[c],"-",geneset.names[r],": ", C[idxs[i]],"\n")
		}
	}
	printTop(C.diff, order(abs(C.diff),decreasing=TRUE), 500)
	printTop(C.wt, order(abs(C.wt),decreasing=TRUE), 500)
}

tcga.ras.methylation <- function(){
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
	
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rasness[,1]), rownames(tcga_aa_mat))
	rasScore.m <- rasness[!is.na(idxs),2]
	names(rasScore.m) <- rasness[!is.na(idxs),1]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	idxs <- pmatch(canonical, colnames(tcga_aa_mat.m))	
	has.mut <- rowSums(tcga_aa_mat.m[,idxs]) > 0
	
	
	meth_eset <- getTCGAMeth()
	
	idxs <- match(extract.tcga.patientIds(sampleNames(meth_eset)),extract.tcga.patientIds(names(rasScore.m)))
	
	meth_eset.m <- meth_eset[,!is.na(idxs)]
	meth_rasScore.m <-rasScore.m[na.omit(idxs)]
	has.mut.m <- has.mut[na.omit(idxs)]
	
	
	epigRes <- getTCGAEpiReg()
	cutDiff <- epigRes$hyperMethDiff >= 0.3
	cutFC <- epigRes$hyperExprFC <= 0.5
	cutPval <- p.adjust(epigRes$hyperExprPval, method="bonferroni") <= 0.05
	
	cutAll <- cutDiff & cutFC & cutPval
	table(cutAll)
	
	epigResCut <- epigRes[which(cutAll), ]
	cpgs <- rownames(epigResCut)
	
	pvals <- sapply(cpgs, function(cpg){
		meth <- as.vector(exprs(meth_eset.m[which(featureNames(meth_eset.m) == cpg),]))
		coef(summary(lm(meth_rasScore.m~ as.vector(meth) + has.mut.m)))[2,4]		
	})
	idxs <- order(pvals)
	R <- cbind(pvals[idxs], epigResCut[idxs,])
	R[1:10,]
	#LRRC2, potential tumor suppressor; similarity to RAS suppressor protein
	# http://www.nature.com/ejhg/journal/v10/n1/full/5200758a.html
	# NAT2; enzyme which metabolizes carcinogens
	
	meth <- as.vector(exprs(meth_eset.m[which(featureNames(meth_eset.m) == "cg02831294"),]))
	plot(meth, meth_rasScore.m, col=(has.mut.m+1),pch=19,cex=.8, xlab="Methylation", ylab="RIS")
	lines(lowess(meth_rasScore.m ~ meth),col="gray",lwd=4,lty=2)
	
	#abline(fit.1, col="red")
	#abline(fit.2, col="black")
}

ccle.predictions.lung <- function(){
	ccle <- getCCLE()
	ccle_eset <- ccle[[1]]
	ccle_drug <- ccle[[2]]
	idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
	
	ccle_eset <- ccle_eset[,!is.na(idxs)]
	ccle_response <- ccle_drug[na.omit(idxs),]
	
	lung.mask <- grepl("LUNG", sampleNames(ccle_eset))
	summary(lm(ccle_response$Erlotinib[lung.mask] ~ ccle_eset$EGFR[lung.mask] + ccle_eset$PIK3CA[lung.mask]))
}

build.circos.data <- function(){
	load("data~/tcga.coad.read.CNA_hg19_nocnv.rda")
	rm <- apply(coad_and_read.CN.byregion[,-1:-3], 1, function(x){ mean(as.numeric(t(x))) } )
	rv <- apply(coad_and_read.CN.byregion[,-1:-3], 1, function(x){ var(as.numeric(t(x))) } )
	RM <- cbind(paste("hs",coad_and_read.CN.byregion[,1],sep=""),
			coad_and_read.CN.byregion[,2:3], rm)
	RV <- cbind(paste("hs",coad_and_read.CN.byregion[,1],sep=""),
			coad_and_read.CN.byregion[,2:3], rv)
	write.table(RM, file="resources/CNA_mean_coad_read.txt", sep="\t",quote=F,row.names=F,col.names=F)
	write.table(RV, file="resources/CNA_var_coad_read.txt", sep="\t",quote=F,row.names=F,col.names=F)
	
	
#	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
#	tcga_aa_mat <- getTCGARasMuts()
	
#	idxs <- match(extract.tcga.patientIds(rasness[,1]), rownames(tcga_aa_mat))
#	rasScore.m <- rasness[!is.na(idxs),2]
#	names(rasScore.m) <- rasness[!is.na(idxs),1]
#	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
#	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	
#	has.mut <- rowSums(tcga_aa_mat.m[,pmatch(canonical, colnames(tcga_aa_mat.m))]) > 0
	
	r <- getRISandMut()
	
	# load CN data
	load("data~/tcga.coad.read.CNA_hg19_nocnv.rda")
	
	cna.annot <- coad_and_read.CN.byregion[, 1:3]
	cna <- coad_and_read.CN.byregion[, -1:-3]
	
	idxs <- match(extract.tcga.patientIds(names(r$ris)),
			extract.tcga.patientIds(colnames(cna)))
	rasScore.m <- r$ris[!is.na(idxs)]
	has.mut.m <- r$mut[!is.na(idxs)]
	cna.m <- cna[, na.omit(idxs)]
	cna.m <- apply(cna.m, 2, function(x){ as.numeric(as.character(x))})
	
	cna.ras.cor <- cor(t(cna.m), rasScore.m,method="spearman")
	cna.ras.cor.wt <- cor(t(cna.m[,!has.mut.m]), rasScore.m[!has.mut.m],method="spearman")
	cna.ras.cor.mut <- cor(t(cna.m[,has.mut.m]), rasScore.m[has.mut.m],method="spearman")
	
	pvals <- apply(cna.m, 1, function(x) {
		coef(summary(lm(rasScore.m ~ x + has.mut.m)))[2,4]
	})
	top.cna <- cna.m[-log10(pvals) > 6,]
	top.cna.ann <- cna.annot[-log10(pvals) > 6,]
	idxs <- order(apply(top.cna, 1, var), decreasing=T)
	top.cna.ann[idxs[1],]
	
	tmp <- cbind(cbind(paste("hs",cna.annot[,1],sep=""),
					cna.annot[,2:3], cna.ras.cor))
	write.table(tmp, file="resources/CNA_rasness_coad_read.txt", sep="\t",quote=F,row.names=F,col.names=F)
	tmp <- cbind(cbind(paste("hs",cna.annot[,1],sep=""),
					cna.annot[,2:3], cna.ras.cor.wt))
	write.table(tmp, file="resources/CNAwt_rasness_coad_read.txt", sep="\t",quote=F,row.names=F,col.names=F)
	tmp <- cbind(cbind(paste("hs",cna.annot[,1],sep=""),
					cna.annot[,2:3], cna.ras.cor.mut))
	write.table(tmp, file="resources/CNAmut_rasness_coad_read.txt", sep="\t",quote=F,row.names=F,col.names=F)
	tmp <- cbind(cbind(paste("hs",cna.annot[,1],sep=""),
					cna.annot[,2:3], -log10(pvals)))
	write.table(tmp, file="resources/CNA_conditionalpval_rasness_coad_read.txt", sep="\t",quote=F,row.names=F,col.names=F)
	
}

cna <- function(){
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rasness[,1]), rownames(tcga_aa_mat))
	rasScore.m <- rasness[!is.na(idxs),2]
	names(rasScore.m) <- rasness[!is.na(idxs),1]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	
	has.mut <- rowSums(tcga_aa_mat.m[,pmatch(canonical, colnames(tcga_aa_mat.m))]) > 0
	
	# load CN data
	load("data~/tcga.coad.read.CNA_hg19_nocnv.rda")
	
	cna.annot <- coad_and_read.CN.byregion[, 1:3]
	cna <- coad_and_read.CN.byregion[, -1:-3]
	
	idxs <- match(extract.tcga.patientIds(names(rasScore.m)),
			extract.tcga.patientIds(colnames(cna)))
	rasScore.m <- rasScore.m[!is.na(idxs)]
	has.mut.m <- has.mut[!is.na(idxs)]
	cna.m <- cna[, na.omit(idxs)]
	cna.m <- apply(cna.m, 2, function(x){ as.numeric(as.character(x))})

	# compute genomic instability using CN
	GI.m <- apply(cna.m, 2, function(x) { sum(x < -2) } ) 
	
	summary(lm(rasScore.m ~ log(GI.m+1) + has.mut.m))
	
	getGeneSummary <- function(chr, start, end){
		chr.mask <- cna.annot$chrom == chr
		ir <- IRanges(as.numeric(as.character(cna.annot[chr.mask,]$start)), 
				as.numeric(as.character(cna.annot[chr.mask,]$end)))
		idxs <- as.matrix(findOverlaps(ir, IRanges(start, end)))[,"query"]
		gene.cn <- cna.m[chr.mask,][idxs,]
		if(length(idxs) > 1){
			gene.cn <- colMeans(apply(gene.cn, 2, function(x){ as.numeric(as.character(x)) } ))
		}else{
			gene.cn <- as.numeric(t(gene.cn))
		}
		gene.cn
	}
	kras.cn <- getGeneSummary(12, start=25358180, end=25403854)
	braf.cn <- getGeneSummary(7, start=140433813, end=140624564)
	nras.cn <- getGeneSummary(1, start=115247085, end=115259515)
	pik3ca.cn <- getGeneSummary(3, 178866311, 178952497)

	rasScore.m[kras.cn > .5 & !has.mut.m]
	rasScore.m[braf.cn > 1 & !has.mut.m]
	rasScore.m[nras.cn > 1 & !has.mut.m]
}


EMT <- function(){
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
	
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rasness[,1]), rownames(tcga_aa_mat))
	rasScore.m <- rasness[!is.na(idxs),2]
	names(rasScore.m) <- rasness[!is.na(idxs),1]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	
	has.mut <- rowSums(tcga_aa_mat.m[,pmatch(canonical, colnames(tcga_aa_mat.m))]) > 0
	has.kras_nras <- rowSums(tcga_aa_mat.m[,pmatch(c(canonical.kras, canonical.nras), colnames(tcga_aa_mat.m))]) > 0
	has.braf <- tcga_aa_mat.m[,pmatch(canonical.braf, colnames(tcga_aa_mat.m))] > 0
	
	
	eset <- getTCGACRC()
	idxs <- match(extract.tcga.patientIds(names(rasScore.m)), extract.tcga.patientIds(sampleNames(eset)))
	rasScore.m <- rasScore.m[!is.na(idxs)]
	has.mut.m <- has.mut[!is.na(idxs)]
	eset.m <- eset[, na.omit(idxs)]
	
	vim_rasmut_pval <- kruskal.test(exprs(eset.m)["VIM",] ~ has.mut.m)$p.value
	
	#summary(lm(rasScore.m ~ exprs(eset.m)["VIM",] + has.mut.m))
	pdf("plots/vimentin_ras_association.pdf", height=5, width=10)
	par(mfrow=c(1,2))
	fit.1 <- lm(rasScore.m[has.mut.m] ~ exprs(eset.m)["VIM",][has.mut.m])
	vim_pvalue_mut <- summary(fit.1)$coef[2,4]
	plot( exprs(eset.m)["VIM",][has.mut.m], rasScore.m[has.mut.m],
			xlab="VIM expression", 
			ylab="RASness",
			main="RAS mutant",
			pch=16, col="black")
	abline(fit.1,col="blue")
	mtext(paste("p=",format(vim_pvalue_mut,digits=2),sep=""))
	
	fit.2 <- lm(rasScore.m[!has.mut.m] ~ exprs(eset.m)["VIM",][!has.mut.m])
	vim_pvalue_wt <- summary(fit.2)$coef[2,4]
	plot( exprs(eset.m)["VIM",][!has.mut.m], rasScore.m[!has.mut.m],
			xlab="VIM expression", 
			ylab="RASness",
			pch=16, col="black",
			main="RAS wild-type")
	abline(fit.2,col="blue")
	mtext(paste("p=",format(vim_pvalue_wt,digits=2),sep=""))
	dev.off()
	# vimentin has a significant negative association with RASness in mut samples, but not WT
	
	# association is significant even when conditioning on braf
	vim_pvalue_mut_condbraf <- summary(lm(rasScore.m[has.mut.m] ~ exprs(eset.m)["VIM",][has.mut.m] + has.braf[has.mut.m]))$coef[2,4]
	
	list(vim_rasmut_pval, vim_pvalue_mut, vim_pvalue_wt, vim_pvalue_mut_condbraf)
}

cbioValidation <- function(){
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
	
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rasness[,1]), rownames(tcga_aa_mat))
	rasScore.m <- rasness[!is.na(idxs),2]
	names(rasScore.m) <- rasness[!is.na(idxs),1]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	idxs <- pmatch(canonical, colnames(tcga_aa_mat.m))	
	has.mut <- rowSums(tcga_aa_mat.m[,idxs]) > 0
	
	foo <- extract.tcga.patientIds(names(rasScore.m))
	# alternative codons (K117N and Q22K) appear to be ras-like; R68S and amplification doesn't
	rasScore.m[foo %in% c("TCGA-AA-3975","TCGA-AA-A01G","TCGA-AG-A025","TCGA-AG-A026")]
	
	# double nras+kras are not additive
	rasScore.m[foo %in% c("TCGA-AA-3558","TCGA-AF-2691")]
}

genomicInstability <- function(){
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
	load("data~/all_mut.RData")
	
	tcga_aa_mat <- getTCGARasMuts()
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	idxs <- pmatch(canonical, colnames(tcga_aa_mat))
	has.mut.mask <- rowSums(tcga_aa_mat[,idxs]) > 0
	idxs <- match(names(has.mut.mask), extract.tcga.patientIds(rasness[,1]))
	ras.mut.mask.m <- has.mut.mask[!is.na(idxs)]
	rasness.m <- rasness[na.omit(idxs),2]
	
	idxs <- match(names(ras.mut.mask.m), rownames(Mtx_mut))
	ras.mut.mask.m <- ras.mut.mask.m[!is.na(idxs)]
	rasness.m <- rasness.m[!is.na(idxs)]
	Mtx.mut.m <- Mtx_mut[na.omit(idxs),]
	
	# genomic instability in wild-type
	GI <- apply(Mtx.mut.m, 1, function(x) { sum(grepl("frame_shift", x)) } )
	
	GI.wt <- apply(Mtx.mut.m[!ras.mut.mask.m,], 1, function(x) { sum(grepl("frame_shift", x)) } )
	GI.by.gene <- apply(Mtx.mut.m[!ras.mut.mask.m,], 2, function(x) { sum(grepl("frame_shift", x)) } )
	pdf("plots/genomicInstability_hist_in_WTRAS.pdf")
	hist(GI.wt, xlab="# aberations per sample, all genes",main="Genomic Instability, RAS wild-type population")
	dev.off()
	GI.wt.factor <- GI.wt > 6000
	
	pdf("plots/genomicInstability_boxplot_in_WTRAS.pdf")
	boxplot(rasness.m[!ras.mut.mask.m] ~ GI.wt.factor, names=c("Stable","Unstable"), ylab="RIS",
			main="RIS by genomic instability")
	dev.off()
	
	pval.GI.wtras <- kruskal.test(rasness.m[!ras.mut.mask.m] ~ GI.wt.factor)$p.value
	summary(lm(rasness.m[!ras.mut.mask.m] ~ GI.wt.factor + pik3ca.factor[!ras.mut.mask.m]))
	
	# genomic instability in MUT
	GI.mut <- apply(Mtx.mut.m[ras.mut.mask.m,], 1, function(x) { sum(grepl("frame_shift", x)) } )
	hist(GI.mut)
	GI.mut.factor <- GI.mut > 6000
	pval.GI.mutras <- kruskal.test(rasness.m[ras.mut.mask.m] ~ GI.mut.factor)$p.value
	summary(lm(rasness.m[ras.mut.mask.m] ~ GI.mut.factor + pik3ca.factor[ras.mut.mask.m]))
	
	# evaluate all other combinations with wt
	pvals.cond.GI <- apply(Mtx.mut.m[!ras.mut.mask.m,], 2, function(x) { 
			fs <- grepl("frame_shift", x)
			if(sum(fs)==0) return (NA)
			summary(lm((rasness.m[!ras.mut.mask.m] ~ fs + GI.wt.factor)))$coef[2,4]
		})

	sort(pvals.cond.GI)
	braf.fs <- grepl("frame_shift", Mtx.mut.m[,"BRAF"])
	fit.1 <- lm(rasness.m[!ras.mut.mask.m] ~ braf.fs[!ras.mut.mask.m] + 
											pik3ca.factor[!ras.mut.mask.m] + GI.wt.factor)
	fit.2 <- lm(rasness.m[!ras.mut.mask.m] ~  pik3ca.factor[!ras.mut.mask.m] + GI.wt.factor)
	anova(fit.1, fit.2)
							
}

evaluateAlternativeMutations <- function(){
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
	load("data~/all_mut.RData")
	
	tcga_aa_mat <- getTCGARasMuts()
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	idxs <- pmatch(canonical, colnames(tcga_aa_mat))
	has.mut.mask <- rowSums(tcga_aa_mat[,idxs]) > 0
	idxs <- match(names(has.mut.mask), extract.tcga.patientIds(rasness[,1]))
	ras.mut.mask.m <- has.mut.mask[!is.na(idxs)]
	rasness.m <- rasness[na.omit(idxs),2]
	
	idxs <- match(names(ras.mut.mask.m), rownames(Mtx_mut))
	ras.mut.mask.m <- ras.mut.mask.m[!is.na(idxs)]
	rasness.m <- rasness.m[!is.na(idxs)]
	Mtx.mut.m <- Mtx_mut[na.omit(idxs),]
	
	#GI <- apply(Mtx.mut.m, 1, function(x) { sum(grepl("frame_shift", x)) } )
	find.add.mut <- function(searchTerm, ras.conditional=FALSE){
		pvals <- sapply(1:dim(Mtx.mut.m)[2], function(idx){
				mut <- factor(grepl(searchTerm, Mtx.mut.m[!ras.mut.mask.m,idx]))
				pval <- NA
				if(length(levels(mut)) > 1){
					if(ras.conditional){
						fit <- lm(rasness.m ~ mut + factor(ras.mut.mask.m))
					}else{
						fit <- lm(rasness.m[!ras.mut.mask.m] ~ mut)
					}
					pval <- summary(fit)$coef[2,4]
				}
				pval
		})
		s.idxs <- order(pvals)
		R <- cbind(colnames(Mtx.mut.m)[s.idxs], pvals[s.idxs])
	}
	
	LPAR1.idx <- which(colnames(Mtx.mut.m) == "LPAR1")
	FGFR4.idx <- which(colnames(Mtx.mut.m) == "FGFR4")
	mut.fgfr4 <- grepl("frame_shift", Mtx.mut.m[!ras.mut.mask.m,FGFR4.idx])
	mut.lpar1 <- grepl("frame_shift", Mtx.mut.m[!ras.mut.mask.m,LPAR1.idx])
	mut.braf <- grepl("frame_shift", Mtx.mut.m[!ras.mut.mask.m, which(colnames(Mtx.mut.m) == "BRAF")])
	
	R.cn <- find.add.mut("Cn",FALSE)
	R.fs <- find.add.mut("Cn|frame_shift",FALSE)
	R.cn.cond <- find.add.mut("Cn",TRUE)
	R.fs.cond <- find.add.mut("Cn|frame_shift",TRUE)
	
	## evaluation of ZGPAT
	zgpat.fs <- grepl("frame_shift", Mtx.mut.m[,"ZGPAT"])
	names(zgpat.fs) <- rownames(Mtx.mut.m)
	
	hoxd13.fs <- grepl("Cn|frame_shift", Mtx.mut.m[,"HOXD13"])
	
	# fishter test with other traits
	pvals <- apply(Mtx.mut.m, 2, function(x){
		mask <- grepl("Cn", x)
		pval = NA
		if(length(levels(factor(mask))) > 1) { pval <- chisq.test(zgpat.fs, mask)$p.value }
		pval
	})
	
	tbl.counts <- table(ras.mut.mask.m, zgpat.fs)
	mutex.pval <- fisher.test(mutex.pval)
	pval.wt <- kruskal.test(rasness.m[!ras.mut.mask.m] ~ zgpat.fs[!ras.mut.mask.m])
	pval.mut <- kruskal.test(rasness.m[ras.mut.mask.m] ~ zgpat.fs[ras.mut.mask.m])
	
	## combination of pik3ca and ZGPAT
	idxs <- pmatch(canonical.pik3ca, colnames(tcga_aa_mat))
	pik3ca.mut.pats <- rownames(tcga_aa_mat)[rowSums(tcga_aa_mat[,idxs]) > 0]
	pik3ca.factor <- names(ras.mut.mask.m) %in% pik3ca.mut.pats
	
	pval.wt <- kruskal.test(rasness.m[!ras.mut.mask.m] ~ pik3ca.factor[!ras.mut.mask.m])
	pten.factor <- grepl("frame_shift", Mtx.mut.m[,"PTEN"])
	egfr.factor.cn <- grepl("Cn", Mtx.mut.m[,"EGFR"])
	egfr.factor.fs <- grepl("frame_shift", Mtx.mut.m[,"EGFR"])
	rasness.m[egfr.factor]
	
	fit.1 <- lm(rasness.m ~ factor(ras.mut.mask.m) + zgpat.fs + pik3ca.factor + egfr.factor)
	fit.2 <- lm(rasness.m[!ras.mut.mask.m] ~ zgpat.fs[!ras.mut.mask.m] + pik3ca.factor[!ras.mut.mask.m] + egfr.factor[!ras.mut.mask.m])
	
	fit.3 <- lm(rasness.m[ras.mut.mask.m] ~ zgpat.fs[ras.mut.mask.m] + pik3ca.factor[ras.mut.mask.m])
	
	#fit.2 <- lm(rasness.m[!ras.mut.mask.m] ~ zgpat.fs[!ras.mut.mask.m] + pten.factor[!ras.mut.mask.m])
	tcga_eset <- getTCGACRC()
	idxs <- match(names(zgpat.fs), extract.tcga.patientIds(sampleNames(tcga_eset)))
	tcga_eset.m <- tcga_eset[,idxs]
	egfr_gene <- exprs(tcga_eset.m)[featureNames(tcga_eset.m) == "EGFR"]
	zgpat_gene <- exprs(tcga_eset.m)[featureNames(tcga_eset.m) == "ZGPAT"]
	hoxd13_gene <- exprs(tcga_eset.m)[featureNames(tcga_eset.m) == "HOXD13"]
	
	zapat.expr.pval <- apply(exprs(tcga_eset.m), 1, function(x){
				summary(lm(x ~ zgpat.fs))$coef[2,4]
			})
	#positive correlation .24
	cor.test(egfr_gene[!zgpat.fs], zgpat_gene[!zgpat.fs])
	# neg correlation; -.18
	cor.test(egfr_gene[zgpat.fs], zgpat_gene[zgpat.fs],method="spearman")
	
	cor.test(egfr_gene[!egfr.factor], zgpat_gene[!egfr.factor],method="spearman")
	cor.test(egfr_gene[egfr.factor], zgpat_gene[egfr.factor])
	
}

test_rasness_in_luad <- function(){
	kfsyscc_eset <- getKFSYSCC()
	env <- new.env()
	load("../AZLung/tcga_luad_eset.rda", env)
	luad_eset <- env$eset
	
	hgnc_symbols <- sapply(featureNames(luad_eset), function(x){ strsplit(x,"|",fixed=TRUE)[[1]][1] } )
	uniq_symbols <- setdiff(unique(hgnc_symbols), "?")
	idxs <- match(uniq_symbols, hgnc_symbols)
	luad_eset <- luad_eset[idxs,]
	featureNames(luad_eset) <- hgnc_symbols[idxs]
	y_hats <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(luad_eset))$yhats
	
	pred <- prediction(y_hats[[1]], list(as.logical(luad_eset$kras)))
	auc.kras <- performance(pred, 'auc')@y.values
	
	pred <- prediction(y_hats[[1]], list(as.logical(luad_eset$ras)))
	auc.ras <- performance(pred, 'auc')@y.values
	
	pdf("plots/raspredict_luad.pdf",width=5,height=5)
	plot(1:4, c(NULL, auc.kras, auc.ras, NULL),xaxt="n",col="red",pch=19,main="Ras prediction in LUAD", ylab="AUC", xlab="Platform/Normalization")
	axis(1, 1:4, labels=c(NULL, "KRAS", "KRAS/BRAF/NRAS",NULL))
	dev.off()
}

compute_rnaseq_rasness <- function(){
	kfsyscc_eset <- getKFSYSCC()
	tcga_eset <- getTCGACRC()
	tcga_aa_mat <- getTCGARasMuts()
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	
	idxs <- pmatch(canonical, colnames(tcga_aa_mat))
	has.mut.mask <- rowSums(tcga_aa_mat[,idxs]) > 0
	
	coad.env <- new.env()
	load("data~/tcga_coad_rnaseq.rda",env=coad.env)
	read.env <- new.env()
	load("data~/tcga_read_rnaseq.rda",env=read.env)
	stopifnot(rownames(coad.env$RPKM) == rownames(read.env$RPKM))
	
	RPKM <- cbind(coad.env$RPKM, read.env$RPKM)
	Counts <- cbind(coad.env$Counts, read.env$Counts)
	
	#idxs <- match(extract.tcga.patientIds(rownames(tcga_aa_mat)), extract.tcga.patientIds(colnames(RPKM)))
	#tcga_aa_mat.m <- tcga_aa_mat[!is.na(idxs),]
	#RPKM <- RPKM[, na.omit(idxs)]
	#Counts <- Counts[, na.omit(idxs)]
	hgnc_symbols <- sapply(rownames(RPKM), function(x){ strsplit(x,"|",fixed=TRUE)[[1]][1] } )
	uniq_symbols <- setdiff(unique(hgnc_symbols), "?")
	idxs <- match(uniq_symbols, hgnc_symbols)
	RPKM <- RPKM[idxs,]
	Counts <- Counts[idxs,]
	rownames(RPKM) <- uniq_symbols
	rownames(Counts) <- uniq_symbols
	
	RPKM.quant <- normalizeBetweenArrays(RPKM,method="quantile")
	
	# cqn normalization
	gene.annot <- read.table("data~/gene.annot.txt",sep="\t",header=TRUE)
	idxs <- match(rownames(Counts), gene.annot$hgnc_symbol)
	Counts.m <- Counts[!is.na(idxs),]
	gene.annot.m <- gene.annot[na.omit(idxs),]
	
	cqn <- cqn(Counts.m, lengths=gene.annot.m$length, x=gene.annot.m$gcpercent,sizeFactors=colSums(Counts),verbose=TRUE)
	RPKM.cqn <- cqn$y + cqn$offset	
	
	make.ras.factor <- function(matrix){
		if(class(matrix)=="ExpressionSet"){
			eset <- matrix
		}else{
			# remove duplicate rows
			idxs <- match(unique(extract.tcga.patientIds(colnames(matrix))),
							extract.tcga.patientIds(colnames(matrix)))
			eset <- new("ExpressionSet",exprs=matrix[, idxs])
		}
		
		
		idxs <- match(extract.tcga.patientIds(rownames(tcga_aa_mat)), extract.tcga.patientIds(sampleNames(eset)))
		tcga_aa_mat.m <- tcga_aa_mat[!is.na(idxs),]
		eset.m <- eset[, na.omit(idxs)]
		idxs <- pmatch(canonical, colnames(tcga_aa_mat.m))
		y <- rowSums(tcga_aa_mat.m[,idxs]) > 0
		return (list(eset=eset.m, y=y))
	}
	
	rpkm <- make.ras.factor(RPKM)
	rpkm.quant <- make.ras.factor(RPKM.quant)
	rpkm.cqn <- make.ras.factor(RPKM.cqn)
	agilent <- make.ras.factor(tcga_eset)
	
	factors <- list(rpkm=factor(!rpkm$y),rpkm.quant=factor(!rpkm.quant$y),rpkm.cqn=factor(!rpkm.cqn$y),agilent=factor(!agilent$y))
	aucs <- mclapply(1:10, function(idx){
		idx <- idx
		N <- dim(kfsyscc_eset)[2]
		mask <- sample(N, round(.66 * N), replace=TRUE)
		y_hats <- binomial_predict_EN(kfsyscc_eset[,mask], factor(kfsyscc_eset$kras_status)[mask], 
			testEsets=list(rpkm$eset, rpkm.quant$eset, rpkm.cqn$eset, agilent$eset),seed=(idx+1),alpha=.1)$yhats
		aucs <- NULL
		for(i in 1:length(y_hats)){
			pred <- prediction(y_hats[[i]], factors[[i]])	
			auc <- performance(pred, 'auc')@y.values[[1]]
			aucs <- c(aucs, auc)
		}
		aucs
	},mc.cores=10,mc.set.seed=TRUE,mc.preschedule =FALSE)
	aucs <- t(sapply(aucs, function(x) x))	
	auc.list <- lapply(1:4, function(x) { aucs[,x] } )		
	names(auc.list) <- c("rpkm", "rpkm.quant", "rpkm.cqn", "agilent")

	pdf("plots/raspredict_rnaseq_vs_agilent.pdf",width=6,height=6)
	boxplot(auc.list, main="Platform AUCs", ylab="AUC", xlab="Platform/Normalization",
          col=c(rep("red",3),"cornflowerblue"),
          ylim=c(.81, .90))
  legend("topright",col=c("red","cornflowerblue"),legend=c("RNAseq","microarray"),lty=1,lwd=3)
	dev.off()
	
	(list(pval.diff=kruskal.test(auc.list)$p.value))
}


# compute enrichment of RAS signatures in colon cancer data sets
# results: colon_gsva_loboda_ras_mut_status.pdf
compute_ras_signature_enrichment <- function(){
	require(GSVA)
	
	ras_gsets <- getRASSigs()
	kfsyscc_eset <- getKFSYSCC()
	khambata_eset <- getKhambata()
	gaedcke_eset <- getGaedcke()
	tcga_eset <- getTCGACRC()
	
	kfsyscc_es_ras <- gsva(kfsyscc_eset, ras_gsets,min.sz=10,max.sz=500)$es
	kfsyscc_es <- rbind(exprs(kfsyscc_es_ras), khambata_enrichment_method(kfsyscc_eset))
	khambata_es_ras <- gsva(khambata_eset,ras_gsets,min.sz=10,max.sz=500)$es
	khambata_es <- rbind(exprs(khambata_es_ras), khambata_enrichment_method(khambata_eset))
	gaedcke_es_ras <- gsva(gaedcke_eset,ras_gsets,min.sz=10,max.sz=500)$es
	gaedcke_es <- rbind(exprs(gaedcke_es_ras), khambata_enrichment_method(gaedcke_eset))
	tcga_es_ras <- gsva(tcga_eset, ras_gsets, min.sz=10,max.sz=500)$es
	tcga_es <- rbind(exprs(tcga_es_ras), khambata_enrichment_method(tcga_eset))
	
	
	
	makeSigPlot <- function(sigIdx){
		names <- c("MUT","WT")
		f1 <- formula(as.vector(kfsyscc_es[sigIdx,]) ~ factor(pData(kfsyscc_eset)$kras_status))
		f2 <- formula(as.vector(tcga_es[sigIdx,]) ~ factor(pData(tcga_eset)$kras==0))
		f3 <- formula(as.vector(khambata_es[sigIdx,]) ~ factor(pData(khambata_eset)$kras_status))
		f4 <- formula(as.vector(gaedcke_es[sigIdx,]) ~ factor(pData(gaedcke_eset)$kras_status))
		par(mfrow=c(1,4))
		boxplot(f1,names=names,main="KFSYSCC")
		mtext(paste("p=",format(kruskal.test(f1)$p.value,digits=2),sep=""),side=1,line=3)
		boxplot(f2,names=names, main="TCGA-CRC")
		mtext(paste("p=",format(kruskal.test(f2)$p.value,digits=2),sep=""),side=1,line=3)
		boxplot(f3,names=names, main="Khambata-Ford")
		mtext(paste("p=",format(kruskal.test(f3)$p.value,digits=2),sep=""),side=1,line=3)
		boxplot(f4,names=names, main="Gaedcke")
		mtext(paste("p=",format(kruskal.test(f4)$p.value,digits=2),sep=""),side=1,line=3)
	}
	
	
	pdf("plots/colon_gsva_loboda_ras_mut_status.pdf",width=7,height=2.5)
	makeSigPlot(5)
	dev.off()
	
	pdf("plots/colon_gsva_settleman_ras_mut_status.pdf",width=7,height=2.5)
	makeSigPlot(3)
	dev.off()
  
	pdf("plots/colon_gsva_bild_ras_mut_status.pdf",width=7,height=2.5)
	makeSigPlot(7)
	dev.off()
  
	pdf("plots/colon_gsva_bild_ras_mut_status.pdf",width=7,height=2.5)
	makeSigPlot(7)
	dev.off()
  
	pdf("plots/colon_gsva_baker_ras_mut_status.pdf",width=7,height=2.5)
	makeSigPlot(8)
	dev.off()
  
  
}

compute_kfsyscc_bootstrapped_aucs <- function(nbootstraps=100,num.processors=5){
	require(ROCR)
	
	eset <- getKFSYSCC()
	N <- ncol(eset)
	
	splitAndTest <- function(i){
		idxs <- sample(N, N/2)
		kras <- factor(eset$kras_status)
		r <- binomial_predict_EN(eset[,idxs], kras[idxs], list(eset[,-idxs]))
		y_hat <- r$yhats[[1]]
		selectedFeatures <- r$featureVec[as.logical(coefficients(r$model) != 0)]
		auc <- performance(prediction(y_hat, factor(kras[-idxs])), 'auc')@y.values[[1]]
		list(auc=auc, sf=selectedFeatures)
	}
	rlist <- mclapply(1:nbootstraps, splitAndTest,mc.cores=num.processors,mc.set.seed=TRUE,mc.preschedule =FALSE)
	aucs <- unlist(lapply(rlist, function(x){ x$auc}))
	modelSizes <- unlist(lapply(rlist, function(x){ length(x$sf)}))
	
	df <- data.frame(gene=featureNames(eset),freq=rep(0,nrow(eset)))
	
	for(r in rlist){
		mask <- df$gene %in% r$sf
		df$freq[mask] <- df$freq[mask] + 1
	}
	df$freq <- df$freq / nbootstraps
	df <- df[order(df$freq,decreasing=T),]
	rownames(df) <- NULL
	
	pdf("plots/kfsyscc_bootstrapped_aucs.pdf",width=8,height=4)
	par(mfrow=c(1,2))
	boxplot(aucs, ylab="AUC",notch=TRUE,col="gold",main="Bootstrapped AUC")
	boxplot(modelSizes, ylab="# Genes Selected",col="red", main="Bootstrapped Model Size")
	dev.off()
	list(mean.auc=mean(aucs), quantile.auc=quantile(aucs,probs=c(.05, .95)),featureFreqTbl=df)
}

compute_kfsyscc_validation <- function(){
	
	kfsyscc_eset <- getKFSYSCC()
	khambata_eset <- getKhambata()
	gaedcke_eset <- getGaedcke()
	tcga_eset <- getTCGACRC()
	tcga_eset <- tcga_eset[,!is.na(tcga_eset$kras)]
  
  all.lung <- getAllLung()
	
	y_hats <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(tcga_eset, khambata_eset, gaedcke_eset))$yhats
	
	y_hats_lung <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
	                                   list(all.lung$tcga.luad, 
                                          all.lung$battle,
                                          all.lung$gse26939,
                                          all.lung$chemores))$yhats
	
	auc_and_perf <- function(yhat, actual){
		pred <- prediction(yhat, actual)
		perf <- performance(pred, 'tpr', 'fpr')
		auc <- performance(pred, 'auc')@y.values[[1]]
    
		idx <- order(performance(pred,"acc")@y.values[[1]],decreasing=T)[1]
		sens <- pred@tp[[1]][idx] / (pred@tp[[1]][idx] + pred@fn[[1]][idx])
		spec <- pred@tn[[1]][idx] / (pred@tn[[1]][idx] + pred@fp[[1]][idx])
    
		list(pred=pred, perf=perf, auc=auc, sens=sens, spec=spec)
	}
	
	tcga_perf <- auc_and_perf(y_hats[[1]], factor(tcga_eset$kras==1))
	khambata_perf <- auc_and_perf(y_hats[[2]], factor(khambata_eset$kras_status=="MUT"))
	gaedcke_perf <- auc_and_perf(y_hats[[3]], factor(gaedcke_eset$kras_status=="MUT"))
  
  tcgaluad_perf <- auc_and_perf(y_hats_lung[[1]], factor(all.lung$tcga.luad$kras))
	battle_perf <- auc_and_perf(y_hats_lung[[2]], factor(all.lung$battle$KRAS=="Mutant"))
	gse26939_perf <- auc_and_perf(y_hats_lung[[3]], factor(all.lung$gse26939$kras))
	chemores_perf <- auc_and_perf(y_hats_lung[[4]], factor(all.lung$chemores$KRAS==1))
 
	pdf("plots/kfsyscc_kras_prediction_validtion.pdf",width=5,height=5)
  #par(bg="gray")
	#plot(tcga_perf$perf@x.values[[1]], tcga_perf$perf@y.values[[1]],type="l",lty=2,col="red")
  
	plot(tcga_perf$perf, lty=1, col="red",bg="gray")
	plot(khambata_perf$perf, lty=1, col="blue",add=TRUE)
	plot(gaedcke_perf$perf, lty=1, col="darkgoldenrod1",add=TRUE)
  
  tmp <- list(tcga_perf, khambata_perf,gaedcke_perf)
  tmp.lung <- list(tcgaluad_perf, battle_perf, gse26939_perf, chemores_perf)
  
  df <- data.frame(auc=unlist(sapply(tmp, function(x) format(x$auc,digits=3))), 
             sens=unlist(sapply(tmp, function(x) format(x$sens,digits=3))),
             spec=unlist(sapply(tmp, function(x) format(x$spec,digits=3))),
              row.names=c("TCGA CRC","Khambata-Ford","Gaedcke"))
  
	df.lung <- data.frame(auc=unlist(sapply(tmp.lung, function(x) format(x$auc,digits=3))), 
	                 sens=unlist(sapply(tmp.lung, function(x) format(x$sens,digits=3))),
	                 spec=unlist(sapply(tmp.lung, function(x) format(x$spec,digits=3))),
	                 row.names=c("TCGA LUAD","Battle","GSE26939","CHEMORES"))
	
	#addtable2plot(.4, .4, df, display.rownames=TRUE,hlines=TRUE)
	legend(.4, .6,
			legend=c("TCGA CRC",
        "Khambata-Ford", 
        "Gaedcke"),
			col=c("red","blue","darkgoldenrod1"),
			lty=1)
	
	dev.off()
	
	return(list(tcga=tcga_perf$auc,KF=khambata_perf$auc,gaedcke=gaedcke_perf$auc))
}

## function for generating mutual exclusivity plot for RAS defects
compute_tcga_ras_feature_summary <- function(){
	
	tcga_aa_mat <- getTCGARasMuts()
	bgColor <- "azure3"
	
	colorSchemes= list(c(bgColor,"green"),c(bgColor,"green"),c(bgColor,"green"),c(bgColor,"green"),
			c(bgColor,"red"),
			c(bgColor,"orange"),c(bgColor,"orange"),c(bgColor,"orange"),
			c(bgColor,"purple"),c(bgColor,"purple"),c(bgColor,"purple"),
			c(bgColor,"white"),c(bgColor,"white"),c(bgColor,"white"),c(bgColor,"white"))
	
	canonical <- c(canonical.kras, canonical.braf, canonical.nras, canonical.hras)
	
	idxs <- pmatch(canonical, colnames(tcga_aa_mat))
	has.mut.mask <- rowSums(tcga_aa_mat[,idxs]) > 0
	# 
	
	plotFeatures <- list()
	for(i in 1:length(canonical)){
		plotFeatures[[canonical[i]]] <- tcga_aa_mat[has.mut.mask,idxs[i]]	
	}
	for(i in rev(seq_along(plotFeatures))){
		idxs <- order(plotFeatures[[i]],decreasing=TRUE)
		for(j in seq_along(plotFeatures)){
			plotFeatures[[j]] <- plotFeatures[[j]][idxs]
		}
	}
	pdf("plots/kras_braf_hars_nras_freq_plot.pdf",width=12)
	displayGenomicFeatures(plotFeatures,colorSchemes=colorSchemes,max.sample.width=50)
	dev.off()
	
	df <- data.frame(mut.count=sapply(plotFeatures, sum))
	df$mut.freq <- df$mut.count / nrow(tcga_aa_mat)
	return (df)
}

compte_tcga_rasact_differences_for_ras_aminoacids <- function(){
  kfsyscc_eset <- getKFSYSCC()
  tcga_eset <- getTCGACRC()
  maf <- getCleanTcgaMaf()
  rasScore <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
                                  list(tcga_eset))$yhats[[1]]
   
  
  overlap.pats <- intersect(extract.tcga.patientIds(sampleNames(tcga_eset)), 
                    unique(extract.tcga.patientIds(unique(maf$Tumor_Sample_Barcode))))
  
  tcga_eset_m <- tcga_eset[,extract.tcga.patientIds(sampleNames(tcga_eset)) %in% overlap.pats]
  rasScore_m <- rasScore[extract.tcga.patientIds(sampleNames(tcga_eset)) %in% overlap.pats,1]
  maf_m <- maf[extract.tcga.patientIds(maf$Tumor_Sample_Barcode) %in% overlap.pats,]
  
  kras_pats <- extract.tcga.patientIds(maf_m[maf_m$Hugo_Symbol == "KRAS",]$Tumor_Sample_Barcode)
  nras_pats <- extract.tcga.patientIds(maf_m[maf_m$Hugo_Symbol == "NRAS",]$Tumor_Sample_Barcode)
  braf_pats <- extract.tcga.patientIds(maf_m[maf_m$Hugo_Symbol == "BRAF",]$Tumor_Sample_Barcode)
  ras_pats <- c(kras_pats, nras_pats, braf_pats)
  braf__kras_wt <- setdiff(braf_pats, kras_pats)
  pik3ca_pats <- extract.tcga.patientIds(maf_m[maf_m$Hugo_Symbol == "PIK3CA",]$Tumor_Sample_Barcode)
  pik3ca__ras_wt <- setdiff(pik3ca_pats, ras_pats)
  pik3ca__ras_mut <- intersect(pik3ca_pats, ras_pats)
  
  ##### KRAS,NRAS,BRAF
  
  kras.maf <- maf_m[maf_m$Hugo_Symbol == "KRAS" ,]
  kras.maf$AAChange <- paste("kras.",kras.maf$AAChange,sep="")
  braf.maf <- maf_m[maf_m$Hugo_Symbol == "BRAF", ]
  braf.maf$AAChange <- paste("braf.",braf.maf$AAChange,sep="")
  nras.maf <- maf_m[maf_m$Hugo_Symbol == "NRAS", ]
  nras.maf$AAChange <- paste("nras.",nras.maf$AAChange,sep="")
  
  kras.aa.changes <- sort(table(kras.maf$AAChange),decreasing=TRUE)
  kras.aa.changes <- kras.aa.changes[names(kras.aa.changes) != "kras.."]
  braf.aa.changes <- sort(table(braf.maf$AAChange),decreasing=TRUE)
  nras.aa.changes <- sort(table(nras.maf$AAChange),decreasing=TRUE)
  
  both.maf <- rbind(kras.maf,braf.maf,nras.maf)
  aa.changes <- c(kras.aa.changes, braf.aa.changes,nras.aa.changes)
  
  aa.factor <- rep("WT", ncol(tcga_eset_m))
  names(aa.factor) <- extract.tcga.patientIds(sampleNames(tcga_eset_m))
  for(aa in names(aa.changes)){
    pats <- extract.tcga.patientIds(both.maf[both.maf$AAChange == aa,]$Tumor_Sample_Barcode)
    aa.factor[pats] <- aa
  }

  aa.factor <- factor(aa.factor,levels=c(names(aa.changes),"WT"))
  
  
  colors <- rep("azure4",length(f))
  colors[as.character(aa.factor) %in% names(kras.aa.changes)] = "red"
  colors[as.character(aa.factor) %in% names(braf.aa.changes)] = "cornflowerblue"
  colors[as.character(aa.factor) %in% names(nras.aa.changes)] = "darkgoldenrod1"
  
  pdf("plots/tcga_aminoacid_ris_kras_braf.pdf",width=10,height=5)
  par(mar=c(8, 4, 2, 8))
  boxplot(rasScore_m ~ aa.factor,las=2, pars=list(outpch=NA,medlwd=.5),ylim=c(0,1), ylab="RIS")
  points(rasScore_m ~ aa.factor,pch=19,cex=.7,col=colors)
  abline(h= quantile(rasScore_m[aa.factor=="WT"])[4], lty=2)
  par(xpd=TRUE)
  legend(length(levels(aa.factor))+3, 1, legend=c("kras","braf","nras","wt"),
         fill=c("red","cornflowerblue","darkgoldenrod1","azure4"),xjust=0)
  par(xpd=FALSE)
  dev.off()
  
  ## PIK3CA
  pik3ca.maf <- maf_m[maf_m$Hugo_Symbol == "PIK3CA" ,]
  pik3ca.exon.factor <- rep("WT", ncol(tcga_eset_m))
  pik3ca.exon.changes <- sort(table(pik3ca.maf$Exon),decreasing=TRUE)
  names(pik3ca.exon.factor) <- extract.tcga.patientIds(sampleNames(tcga_eset_m))
  for(exon in names(pik3ca.exon.changes)){
    pats <- extract.tcga.patientIds(pik3ca.maf[pik3ca.maf$Exon == exon,]$Tumor_Sample_Barcode)
    pik3ca.exon.factor[pats] <- exon
  }
  
  pik3ca.exon.factor <- factor(pik3ca.exon.factor,levels=c(names(pik3ca.exon.changes),"WT"))
  pik3ca.exon.factor[pik3ca.exon.factor == "WT" & aa.factor != "WT"] <- NA
  ras_wt_mask <- pik3ca.exon.factor != "WT" & aa.factor == "WT"
  ras_mut_mask <- pik3ca.exon.factor != "WT" & aa.factor != "WT"
  names(pik3ca.exon.factor[ras_mut_mask]) <- paste("ras.mut_",names(pik3ca.exon.factor[ras_mut_mask]),sep="")
  names(pik3ca.exon.factor[ras_wt_mask]) <- paste("ras.wt_",names(pik3ca.exon.factor[ras_wt_mask]),sep="")
 
  model <- rasScore_m[pik3ca.exon.factor != "WT"] ~ (aa.factor != "WT")[pik3ca.exon.factor != "WT"]
  boxplot(model)
  kruskal.test(model)
  
  colors = rep("azure4",length(pik3ca.exon.factor))
  colors[pik3ca.exon.factor != "WT"] <- "cornflowerblue"
  colors[pik3ca.exon.factor != "WT" & aa.factor != "WT"] <- "red"
  pdf("plots/tcga_pik3ca_exon_ris.pdf",width=8,height=5)
  boxplot(rasScore_m ~ factor(pik3ca.exon.factor),las=2, pars=list(outpch=NA,medlwd=.5),ylim=c(0,1), ylab="RIS")
  points(rasScore_m ~ factor(pik3ca.exon.factor),pch=19,cex=.6,col=colors)
  abline(h= quantile(na.omit(rasScore_m[pik3ca.exon.factor=="WT"]))[4], lty=2)
  legend(6, .95, legend=c("kras/braf wt","kras/braf mut","all wt"),
         fill=c("cornflowerblue","red","azure4"))
  dev.off()
  
  tmp.factor <- rep("WT", ncol(tcga_eset_m))
  names(tmp.factor) <- extract.tcga.patientIds(sampleNames(tcga_eset_m))
  tmp.factor[names(tmp.factor) %in% kras_pats] <- "kras"
  tmp.factor[names(tmp.factor) %in% braf_pats] <- "braf"
  tmp.factor[names(tmp.factor) %in% nras_pats] <- "nras"
  tmp.factor[names(tmp.factor) %in% pik3ca__ras_mut] <- "pik3ca/ras+"
  tmp.factor[names(tmp.factor) %in% pik3ca__ras_wt] <- "pik3ca/ras-"
  colors = rep("azure4",length(tmp.factor))
  colors[tmp.factor != "WT"] <- "azure4"
  colors[tmp.factor == "kras"] <- "red"
  colors[tmp.factor == "braf"] <- "cornflowerblue"
  colors[tmp.factor == "nras"] <- "darkgoldenrod1"
  colors[tmp.factor == "pik3ca/ras+"] <- "aquamarine3"
  colors[tmp.factor == "pik3ca/ras-"] <- "aquamarine4"
  
  factor_color <- c("cornflowerblue","red","darkgoldenrod1","aquamarine3","aquamarine4","azure4")
  pdf("plots/tcga_kras_braf_hras_pik3ca.pdf",width=10,height=5)
  par(mar=c(8, 4, 2, 10))
  boxplot(rasScore_m ~ factor(tmp.factor),las=2, pars=list(outpch=NA,medlwd=.5),ylim=c(0,1), ylab="RIS")
  #points(rasScore_m ~ factor(tmp.factor),pch=19,cex=.8,col=colors)
  stripchart(rasScore_m ~ factor(tmp.factor),pch=19,cex=.8,col=factor_color,method="jitter",vertical=TRUE,add=TRUE)
  abline(h= quantile(na.omit(rasScore_m[tmp.factor=="WT"]))[4], lty=2)
  par(xpd=TRUE)
  legend(6.8, 1, legend=c("braf","kras","nras","pik3ca/ras-","pik3ca/ras+","wt"),
         fill=c("cornflowerblue","red","darkgoldenrod1","aquamarine3","aquamarine4","azure4"))
  par(xpd=FALSE)
  dev.off()  
  
  lvls <- levels(factor(tmp.factor))
  M <- matrix(NA, nrow=5,ncol=5)
  for(i in 1:4){
    for(j in (i+1):5){
      mask <- tmp.factor %in% lvls[c(i,j)]
      M[i,j] <- kruskal.test(rasScore_m[mask] ~ factor(tmp.factor[mask]))$p.value
    }
  }
  
  #################################
  # look for novel mutations associated in Ras WT samples
  wt_pats <- names(aa.factor)[aa.factor == "WT"]
  maf_wt <- maf_m[extract.tcga.patientIds(maf_m$Tumor_Sample_Barcode) %in% wt_pats 
                  & maf_m$Variant_Classification %in% c("Missense_Mutation","Nonsense_Mutation"),]
  maf_patids <- extract.tcga.patientIds(maf_wt$Tumor_Sample_Barcode)
  rasScore_wt <- rasScore_m[aa.factor == "WT"]
  rasScore_wt_patids <- extract.tcga.patientIds(names(rasScore_wt))
  all.genes <- unique(maf_wt$Hugo_Symbol)
  pvals <- rep(NA, length(all.genes))
  names(pvals) <- all.genes
  for(i in 1:length(all.genes)){
    gene <- all.genes[i]
    tmp <- unique(maf_patids[maf_wt$Hugo_Symbol == gene])
    f <- factor(rasScore_wt_patids %in% tmp)
    if(length(levels(f)) == 1) { next }
    pvals[i] <- kruskal.test(rasScore_wt ~ f)$p.value
    if(i %% 100 == 0){ cat(".")}
  }
  
  gene <- "CTNNB1"
  tmp <- unique(maf_patids[maf_wt$Hugo_Symbol == gene])
  f <- factor(rasScore_wt_patids %in% tmp)
  pdf("plots/ctnnb1_ris_assoc_in_WT.pdf",width=6,height=6)
  boxplot(rasScore_wt ~ f,las=2, pars=list(outpch=NA,medlwd=.5),ylim=c(0,1), ylab="RIS",names=c("WT","MUT"),main='CTNNB1')
  stripchart(rasScore_wt ~ f,pch=19,cex=.8,col=c("darkgoldenrod1","cornflowerblue"),method="jitter",vertical=TRUE,add=TRUE)
  dev.off()
  
  
  # look for novel mutations associated in Ras Mutant samples
  mut_pats <- names(aa.factor)[grepl("nras|kras",as.character(aa.factor),perl=T)]
  maf_mut <- maf_m[extract.tcga.patientIds(maf_m$Tumor_Sample_Barcode) %in% mut_pats 
                  & !(maf_m$Variant_Classification %in% c("Silent")),]
  maf_patids <- extract.tcga.patientIds(maf_mut$Tumor_Sample_Barcode)
  rasScore_mut <- rasScore_m[grepl("nras|kras",as.character(aa.factor),perl=T)]
  rasScore_mut_patids <- extract.tcga.patientIds(names(rasScore_mut))
  all.genes <- unique(maf_mut$Hugo_Symbol)
  pvals <- rep(NA, length(all.genes))
  names(pvals) <- all.genes
  for(i in 1:length(all.genes)){
    gene <- all.genes[i]
    tmp <- unique(maf_patids[maf_mut$Hugo_Symbol == gene])
    f1 <- factor(rasScore_mut_patids %in% tmp)
    f2 <- factor(rasScore_mut < .35)
    if(length(levels(f)) == 1) { next }
    pvals[i] <- fisher.test(f1,f2,alternative="less")$p.value
    if(i %% 100 == 0){ cat(".")}
  }
}

# Function for demonstrating RASness on a per-codon basis in TCGA
compute_tcga_rasact_differences_for_ras_codons <- function(){
	
	kfsyscc_eset <- getKFSYSCC()
	tcga_eset <- getTCGACRC()
	rasScore <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(tcga_eset))$yhats[[1]]
	
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rownames(rasScore)), rownames(tcga_aa_mat))
	rasScore.m <- rasScore[!is.na(idxs),]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	
	canonical <- c(canonical.kras, canonical.braf, canonical.nras, canonical.hras)
	# boxplot for each codon, vs WT
	codon_factor <- rep("WT", dim(tcga_aa_mat.m)[1])
	for(codon in canonical){
		mask <- rowSums(cbind(tcga_aa_mat.m[,pmatch(codon, colnames(tcga_aa_mat.m))])) > 0
		codon_factor[mask & codon_factor=="WT"] <- codon
	}
	lvls <- levels(factor(codon_factor))
	col <- factor(gsub("([A-Za-z]+).*", "\\1", lvls))
	wt_75_quantile <- quantile(rasScore.m[codon_factor=="WT"])[4]
	pdf("plots/boxplot_kras_braf_nras_hras_wt_rasScores.pdf",width=8,height=5)
	boxplot(rasScore.m ~ codon_factor,col=as.numeric(col)+1, 
			ylab="RASness Score",main="TCGA RASness by Codon",xaxt="n")
	axis(1, at=1:length(lvls), labels=lvls,las=2)
	abline(h=wt_75_quantile, lty=2)
	dev.off()
	
	p.vals <- sort(sapply(setdiff(codon_factor,"WT"), function(x){
			mask <- codon_factor %in% c(x, "WT")
			kruskal.test(rasScore.m[mask] ~ factor(codon_factor[mask]))$p.val
	}))
	p.vals.kras.12 <- sapply(setdiff(codon_factor,"kras.12"), function(x){
			mask <- codon_factor %in% c(x, "kras.12")
			kruskal.test(rasScore.m[mask] ~ factor(codon_factor[mask]))$p.val
		})
	
	######################################################
	## how do PI3K mutations contribute to Ras score?
	pik3ca.mut.mask <- rowSums(cbind(tcga_aa_mat.m[,pmatch(canonical.pik3ca, colnames(tcga_aa_mat.m))])) > 0
	ras.mut.mask <- codon_factor != "WT"
	f0 <- rasScore.m ~ factor(pik3ca.mut.mask)
	f1.ras_mut <- rasScore.m[ras.mut.mask] ~ factor(pik3ca.mut.mask[ras.mut.mask])
	f2.ras_wt <- rasScore.m[!ras.mut.mask] ~ factor(pik3ca.mut.mask[!ras.mut.mask])
	
	pik3ca.pval <- kruskal.test(f0)$p.value
	pik3ca_ras_mut.pval <- kruskal.test(f1.ras_mut)$p.value
	pik3ca_ras_wt.pval <- kruskal.test(f2.ras_wt)$p.value
	
	pik3ca_tbl <- cbind(pik3ca.pval,pik3ca_ras_mut.pval,pik3ca_ras_wt.pval)
	
	pdf("plots/pik3ca_ras_score_vs_ras.pdf",width=10,height=4)
	par(mfrow=c(1,3))
	boxplot(f0,notch=TRUE,col=c("gold","darkgreen"),names=c("PIK3CA WT","PIK3CA Mut"),
			ylab="Ras activiation score",
			main="All samples")
	mtext(paste("p=",format(pik3ca.pval,digits=2),sep=""),side=1,line=3)
	
	boxplot(f1.ras_mut,notch=TRUE,col=c("gold","darkgreen"),names=c("PIK3CA WT","PIK3CA Mut"),
			ylab="Ras activiation score",
			main="Ras mutant samples")
	mtext(paste("p=",format(pik3ca_ras_mut.pval,digits=2),sep=""),side=1,line=3)
	boxplot(f2.ras_wt,notch=TRUE,col=c("gold","darkgreen"),
			names=c("PIK3CA WT","PIK3CA Mut"),
			ylab="Ras activiation score",
			main="Ras wild type samples")
	mtext(paste("p=",format(pik3ca_ras_wt.pval,digits=2),sep=""),side=1,line=3)
	dev.off()
	pik3ca.mut.count.in.ras.wt <- table(pik3ca.mut.mask, ras.mut.mask)[2,1]
	pik3ca.mut.count.in.ras.mut <- table(pik3ca.mut.mask, ras.mut.mask)[2,2]
	############ done pik3ca ##
	(list(pvals.from.WT=data.frame(pvalue=p.vals),
		  pvals.from.kras12=p.vals.kras.12,
		  pik3ca.pvals=pik3ca_tbl,
		  pik3ca.mut.count=c(ras.wt=pik3ca.mut.count.in.ras.wt,
				  			 ras.mut=pik3ca.mut.count.in.ras.mut)))
}

#compute.alternative.mutations <- function(rasActivityScore, rasSampleMutants){
#	env <- new.env()
#	load("/gluster/home/qmeng/My_Experiments/Projects/ColonCa/TCGA_RE.CO_DNASeq_mutatoins/Mutations_ALL/3rd_round-Five_centers_TCGA_somatic_mutations/Step_f_TCGA_mut_summary##/All_genes/Results_Rdata/all_mut.RData",env)
#	mutTbl <- env$Mtx_mut
#	idxs <- match(extract.tcga.patientIds(names(rasActivityScore)), 
#			rownames(mutTbl))
#	rasActivityScore.m <- rasActivityScore[!is.na(idxs)]
#	mutTbl.m <- mutTbl[na.omit(idxs),]
#	cnOnly <- apply(mutTbl.m, 2, function(x){ grepl("Cn",x)})
#	cnOrFs <- apply(mutTbl.m, 2, function(x){ grepl("Cn|frame_shift",x)})
#	
#	ras_mut_mask <- factor(names(rasActivityScore.m) %in% rasSampleMutants)
#	#fit.1 <- lm(rasActivityScore.m ~ ras_mut_mask)
#	p.vals.cn.only <- apply(cnOnly, 2, function(test.mut) {
#				if(sum(test.mut)==0){ return (NA)}
#				fit.2 <- lm(rasActivityScore.m ~ ras_mut_mask + factor(test.mut))
#				anova(fit.1, fit.2)[2,6]
#			})
#	mut.cts <- apply(cnOrFs[!ras_mut_mask,], 2, sum)
#	p.vals.cn.or.fs <- apply(cnOrFs[!ras_mut_mask,mut.cts > 8], 2, function(test.mut) {
#				if(sum(test.mut)==0){ return (NA)}
#				fit.2 <- lm(rasActivityScore.m[!ras_mut_mask] ~ factor(test.mut))
#				anova(fit.2)[1,5]
#			})
#}

compute_sanger_drug_response <- function(){
	kfsyscc_eset <- getKFSYSCC()
	sanger_expr <- getSangerCellLineData()
	
	intestine.mask <- sanger_expr$PrimarySite == "large intestine"
		
	colon_eset <- sanger_expr[, intestine.mask]
	
	pred <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(sanger_expr, colon_eset))

	all.yhat <- pred$yhats[[1]]
	crc.yhat <- pred$yhats[[2]]
	
	drug.mask <- grepl("DRUG",colnames(pData(sanger_expr)))
	
	test.drug <- function(drug.tbl, yhat){
		tmp <- t(apply(drug.tbl, 2, function(x){
							t <- cor.test(x, yhat, method="spearman")
							c(stat=t$estimate, pval=t$p.value, n=sum(!is.na(x)))
						}))
		tmp <- cbind(tmp, pval.adj=p.adjust(tmp[,"pval"], method="bonferroni"))
		tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
		rownames(tmp) <- gsub("DRUG\\.","", rownames(tmp))
		tmp
	}
	
	idxs <- order(all.yhat)
	sanger_expr$PrimarySite[idxs]
	drug.all <- test.drug(pData(sanger_expr)[, drug.mask], all.yhat)
	drug.crc <- test.drug(pData(sanger_expr)[intestine.mask,drug.mask],crc.yhat)
	
	idxs <- order(stats)
	pdf(paste("plots/ccle_ic50_rasact_association_",sample.origin,".pdf",sep=""),width=8,height=8)
	par(mar=c(6, 8, 2, 2))
	col=rep("black",length(pvals))
	col[pvals[idxs] < .05] <- "red"
	bp <- barplot(stats[idxs],horiz=T,las=2,xlim=c(-1, 1),col=col,xlab="Spearman rho")
	text(.7, bp,labels=format(pvals[idxs],digits=2))
	dev.off()
	
	return (df)
}


compute_ccle_drug_response <- function(sample.origin="LARGE_INTESTINE"){
	kfsyscc_eset <- getKFSYSCC()
	ccle <- getCCLE()
	
	ccle_eset <- ccle[[1]]
	ccle_drug <- ccle[[2]]
	idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
	
	ccle_eset <- ccle_eset[,!is.na(idxs)]
	ccle_response <- ccle_drug[na.omit(idxs)]
	
	if(!is.null(sample.origin)){
		sub_eset <- ccle_eset[, grep(sample.origin,sampleNames(ccle_eset))]
		sub_response <- pData(ccle_response[grep(sample.origin,sampleNames(ccle_response))])
	}else{
		sub_eset <- ccle_eset
		sub_response <- pData(ccle_response)
	}
	RASact <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(sub_eset))$yhats[[1]]
	
	stats <- apply(sub_response[mask,], 2, function(x){ cor.test(RASact[mask], x,method="spearman")$estimate })
	pvals <- apply(sub_response[mask,], 2, function(x){ cor.test(RASact[mask], x,method="spearman")$p.value })
	has.ras <- length(levels(factor(sub_eset$KRAS))) > 1
	has.raf <- length(levels(factor(sub_eset$BRAF))) > 1
	
	pvals <- apply(sub_response[mask,], 2, function(x){ cor.test(RASact[mask], x,method="spearman")$p.value })
	pvals <- apply(sub_response[mask,], 2, function(x){ 
				coef(summary(lm(RASact[mask] ~ x + has.mut[mask])))[2,4] })
	
	if(has.ras){
		mut.kras.pval=apply(sub_response, 2, function(x) { kruskal.test(x ~ factor(sub_eset$KRAS))$p.value } )
	}else{
		mut.kras.pval=rep(NA, length(pvals))
	}
	
	if(has.raf){
		mut.braf.pval=apply(sub_response, 2, function(x) { kruskal.test(x ~ factor(sub_eset$BRAF))$p.value } )
	}else{
		mut.braf.pval <- rep(NA, length(pvals))
	}
	if(has.ras | has.raf){
		mut.kras_braf.pval=apply(sub_response, 2, function(x) { kruskal.test(x ~ factor(sub_eset$KRAS | sub_eset$BRAF))$p.value } )
	}else{
		mut.kras_braf.pval <- rep(NA, length(pvals))
	}
	df <- data.frame(model.rho=stats,
			model.pval=pvals,
			mut.kras.pval=mut.kras.pval,
			mut.braf.pval=mut.braf.pval,
			mut.kras_braf.pval=mut.kras_braf.pval)
	
	idxs <- order(stats)
	pdf(paste("plots/ccle_ic50_rasact_association_",sample.origin,".pdf",sep=""),width=8,height=8)
	par(mar=c(6, 8, 2, 2))
	col=rep("black",length(pvals))
	col[pvals[idxs] < .05] <- "red"
	bp <- barplot(stats[idxs],horiz=T,las=2,xlim=c(-1, 1),col=col,xlab="Spearman rho")
	text(.7, bp,labels=format(pvals[idxs],digits=2))
	dev.off()
	
	sig.drugs <- pvals < .01
	glyphs <- rep(19, length(pvals))
	glyphs[sig.drugs] <- 1:sum(sig.drugs)
	pdf(paste("plots/ccle_ic50_rasact_association_",sample.origin,"_volcano.pdf",sep=""),width=8,height=8)
	plot(stats, -log10(pvals),xlab="Spearman rho",pch=glyphs)
	abline(h=-log10(.01),lty=2,col="red")
	dev.off()
	
	return (df)
}

compute_ccle_drug_response_revised <- function(sample.origin="LARGE_INTESTINE"){
  kfsyscc_eset <- getKFSYSCC()
  ccle <- getCCLE_MetaGenomics()
  
    
  if(!is.null(sample.origin)){
    sub_eset <- ccle[, grep(sample.origin,sampleNames(ccle))]
  }else{
    sub_eset <- ccle
  }
  RASact <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
                                list(sub_eset))$yhats[[1]]
  
  stats <- apply(pData(sub_eset), 2, function(x){ cor.test(RASact, x,method="spearman")$estimate })
  pvals <- apply(pData(sub_eset), 2, function(x){ cor.test(RASact, x,method="spearman")$p.value })
  
  ##TODO .....
  
}

find_methylation_rasact_pairs <- function(wt.ras_act.threshold=.4, gene.plots=c("RASSF1","STK11")){
	
	library(IlluminaHumanMethylation27k.db)
	cpgSymbols <- as.list(IlluminaHumanMethylation27kSYMBOL)
	
	
	kfsyscc_eset <- getKFSYSCC()
	tcga_eset <- getTCGACRC()
	rasScore <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(tcga_eset))$yhats[[1]]
	
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rownames(rasScore)), rownames(tcga_aa_mat))
	rasScore.m <- rasScore[!is.na(idxs),]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	wt.mask <- rowSums(tcga_aa_mat.m[,pmatch(canonical,colnames(tcga_aa_mat.m))]) == 0
	
	meth_eset <- getTCGAMeth()
	idxs <- match(extract.tcga.patientIds(sampleNames(meth_eset)),extract.tcga.patientIds(names(rasScore.m)))
	
	meth_eset.m <- meth_eset[,!is.na(idxs)]
	meth_rasScore.m <- rasScore.m[na.omit(idxs)]
	meth_wt.mask <- wt.mask[na.omit(idxs)]
	
	grps_wt <- factor(meth_rasScore.m[meth_wt.mask] > wt.ras_act.threshold)
	grps_mut <- factor(meth_rasScore.m[!meth_wt.mask] > wt.ras_act.threshold)
	pvals.meth.fp <- apply(exprs(meth_eset.m)[,meth_wt.mask], 1, function(cpg) {		
				kruskal.test(cpg ~ grps_wt)$p.value	
			})
	pvals.meth.fp.sort <- sort(pvals.meth.fp,decreasing=FALSE)
	
	cpg.annot.diff <- setdiff(names(pvals.meth.fp.sort), names(cpgSymbols))
	pvals.meth.fp.sort <- pvals.meth.fp.sort[!(names(pvals.meth.fp.sort) %in% cpg.annot.diff)]
	cpgSymbols <- cpgSymbols[!(names(cpgSymbols) %in% cpg.annot.diff)]
	
	df <- data.frame(pvals.meth.fp.sort,gene=unlist(cpgSymbols[names(pvals.meth.fp.sort)]))
	df <- df[!is.na(df$gene),]
	
	makeplot_fnfp <- function(gene){
		idx <- which(df$gene == gene)[1]
		cpg <- rownames(df)[idx]
		
		par(mfrow=c(1,3))
		plot(density(exprs(meth_eset.m)[cpg,]),
				main=paste(cpg,"/",gene),
				xlab="Methylation")
		f1 <- as.vector(exprs(meth_eset.m[cpg,meth_wt.mask])) ~ grps_wt
		boxplot(f1,
				names=c(c("low Ras activity","high Ras activity")),
				main=paste("WT Ras samples",sep=""),
				ylab="Methylation",
				col=c("gold","darkgreen"),
				notch=TRUE)
		mtext(paste("p=",format(kruskal.test(f1)$p.value,digits=2),sep=""),side=1,line=3)
		
		f2 <- as.vector(exprs(meth_eset.m[cpg,!meth_wt.mask])) ~ grps_mut
		boxplot(f2,
				names=c("low Ras activity","high Ras activity"),
				main=paste("Mut Ras samples",sep=""),
				ylab="Methylation",
				col=c("gold","darkgreen"),
				notch=TRUE)
		mtext(paste("p=",format(kruskal.test(f2)$p.value,digits=2),sep=""),side=1,line=3)
	}
	for(i in 1:length(gene.plots)){
		pdf(paste("plots/methyl_rasact_",gene.plots[i],".pdf",sep=""),width=10,height=4)
		makeplot_fnfp(gene.plots[i])
		dev.off()
	}
	
	
	return (df)
}

find_methylation_rasact_pairs_3 <- function(){
	
	library(IlluminaHumanMethylation27k.db)
	cpgSymbols <- as.list(IlluminaHumanMethylation27kSYMBOL)
	
	rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
		
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rasness[,1]), rownames(tcga_aa_mat))
	rasScore.m <- rasness[!is.na(idxs),2]
	names(rasScore.m) <- rasness[!is.na(idxs),1]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	idxs <- pmatch(canonical, colnames(tcga_aa_mat.m))	
	has.mut <- rowSums(tcga_aa_mat.m[,idxs]) > 0
	
	
	meth_eset <- getTCGAMeth()
	
	idxs <- match(extract.tcga.patientIds(sampleNames(meth_eset)),extract.tcga.patientIds(names(rasScore.m)))
	
	meth_eset.m <- meth_eset[,!is.na(idxs)]
	meth_rasScore.m <-rasScore.m[na.omit(idxs)]
	has.mut.m <- has.mut[na.omit(idxs)]
	
	pvals <- apply(meth_eset.m, 1, function(x){
			mask <- (x > 2)
			pval = NA
			if(sum(mask) > 5){
				pval <- phyper(sum(mask & has.mut.m), sum(mask), length(mask)-sum(mask), sum(has.mut.m))
				#pval <- fisher.test(table(mask, has.mut.m))$p.value
			}
			pval
		})

	cpgSymbols[names(sort(pvals))][1:20]
	fisher.test(table(exprs(meth_eset.m["cg21460081", ]) > 2), cbind(has.mut.m))
	
	fit <- eBayes(lmFit(meth_eset.m, model.matrix(~meth_rasScore.m + factor(has.mut.m))))
	cpgs <- featureNames(meth_eset.m)[order(fit$p.value[,2])]
	cpgSymbols[cpgs[1:20]]
}

find_methylation_rasact_pairs_2 <- function(){
	
	library(IlluminaHumanMethylation27k.db)
	cpgSymbols <- as.list(IlluminaHumanMethylation27kSYMBOL)
	
	
	kfsyscc_eset <- getKFSYSCC()
	tcga_eset <- getTCGACRC()
	rasScore <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(tcga_eset))$yhats[[1]]
	
	tcga_aa_mat <- getTCGARasMuts()
	
	idxs <- match(extract.tcga.patientIds(rownames(rasScore)), rownames(tcga_aa_mat))
	rasScore.m <- rasScore[!is.na(idxs),]
	tcga_aa_mat.m <- tcga_aa_mat[na.omit(idxs),]
	
	canonical <- c(canonical.kras, canonical.braf, canonical.nras)
	wt.mask <- rowSums(tcga_aa_mat.m[,pmatch(canonical,colnames(tcga_aa_mat.m))]) == 0
	
	meth_eset <- getTCGAMeth()
	
	active.methyl.cpgs <- find.active.methylation.sites(meth_eset, tcga_eset)
	idxs <- match(extract.tcga.patientIds(sampleNames(meth_eset)),extract.tcga.patientIds(names(rasScore.m)))
	
	meth_eset.m <- meth_eset[,!is.na(idxs)]
	meth_rasScore.m <- rasScore.m[na.omit(idxs)]
	meth_wt.mask <- wt.mask[na.omit(idxs)]
	
	M <- exprs(meth_eset.m)
	R <- data.frame(t(sapply(active.methyl.cpgs$outlier.cpgs, function(cpg){
		#summary(lm(meth_rasScore.m ~ factor(meth_wt.mask) + M[cpg,]))$coefficients[3,4]
		p1 <- kruskal.test(M[cpg, meth_wt.mask] ~ meth_rasScore.m[meth_wt.mask] > .35)$p.value
		p2 <- kruskal.test(M[cpg, !meth_wt.mask] ~ meth_rasScore.m[!meth_wt.mask] > .35)$p.value
		c(logp.diff=-log10(p1) + log10(p2), wt=p1, mut=p2)
	})))
	R <- cbind(R, gene=unlist(cpgSymbols[rownames(R)]))
	R[order(R[,3])[1:20],]
	
	s.idxs <- order(abs(R[,1]),decreasing=TRUE)
	R[s.idxs[1:10],]
	
	cpg="cg04577715"
	cpg="cg08047457"
	par(mfrow=c(1,2))
	table(low_ras=M[cpg, meth_wt.mask] > 2 & meth_rasScore.m[meth_wt.mask] < .4,
			high_ras=M[cpg, meth_wt.mask] > 2 & meth_rasScore.m[meth_wt.mask] > .4)

	plot(M[cpg, !meth_wt.mask] ~ meth_rasScore.m[!meth_wt.mask] > .4,main="RAS MUT")
	plot(M[cpg, meth_wt.mask] ~ meth_rasScore.m[meth_wt.mask] > .4,main="RAS WT")
	
	R
}

predictCetuximabResponse <- function(){
	kfsyscc_eset <- getKFSYSCC()
	khambata_eset <- getKhambata()
	
	y_hat <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
			list(khambata_eset))$yhats[[1]]
	not_utd_mask <- khambata_eset$Best != "UTD"
	wt.mask <- khambata_eset$kras_status == "WT"

	grps <- factor(khambata_eset$Best[not_utd_mask]=="PD")
	by_ras_score <- kruskal.test(y_hat[not_utd_mask] ~ grps)$p.value
	by_kras_status <- fisher.test(factor(khambata_eset$kras_status[not_utd_mask]), grps)$p.value
	
	by_ras_score_WT <- kruskal.test(y_hat[not_utd_mask][wt.mask] ~ grps[wt.mask])$p.value
	by_ras_score_MUT <- kruskal.test(y_hat[not_utd_mask][!wt.mask] ~ grps[!wt.mask])$p.value
	
	pdf("plots/khambata_ford_cetuximab_response.pdf",width=10, height=4)
	par(mfrow=c(1,3))
  model <- y_hat[not_utd_mask] ~ factor(khambata_eset$Best[not_utd_mask])
	boxplot(model,
	    pars=list(outpch=NA,medlwd=.5),
			ylab="RIS",
			main="All Samples")
  stripchart(model, 
             vertical=TRUE,
             col=c("gold","cornflowerblue","cornflowerblue"),
             pch=19,cex=1,
             method = "jitter",add=TRUE)
	mtext(paste("p=",format(by_ras_score,digits=3),sep=""),side=1,line=3)
	
  model <- y_hat[not_utd_mask][wt.mask] ~ factor(khambata_eset$Best[not_utd_mask][wt.mask])
  boxplot(model,
    pars=list(outpch=NA,medlwd=.5),
		ylab="RIS",
		main="Kras wild-type")
	stripchart(model, 
	           vertical=TRUE,
	           col=c("gold","cornflowerblue","cornflowerblue"),
	           pch=19,cex=1,
	           method = "jitter",add=TRUE)
	mtext(paste("p=",format(by_ras_score_WT,digits=2),sep=""),side=1,line=3)
	
  model <- y_hat[not_utd_mask][!wt.mask] ~ factor(khambata_eset$Best[not_utd_mask][!wt.mask])
	boxplot(model,
	    pars=list(outpch=NA,medlwd=.5),
			ylab="RIS",
			main="Kras mutant")
	stripchart(model, 
	           vertical=TRUE,
	           col=c("gold","cornflowerblue","cornflowerblue"),
	           pch=19,cex=1,
	           method = "jitter",add=TRUE)
	mtext(paste("p=",format(by_ras_score_MUT,digits=2),sep=""),side=1,line=3)
	dev.off()
	
	return (list(pval_by_ras_score=by_ras_score, 
				pval_by_kras_status=by_kras_status, 
				pval_by_ras_score_WT=by_ras_score_WT))
}

checkLobodaCCLEPredictions <- function(sample.origin=NULL){
	ccle <- getCCLE()
	
	ccle_eset <- ccle[[1]]
	ccle_drug <- ccle[[2]]
	idxs <- match(sampleNames(ccle_eset), sampleNames(ccle_drug))
	
	ccle_eset <- ccle_eset[,!is.na(idxs)]
	ccle_response <- ccle_drug[na.omit(idxs)]
	
	if(!is.null(sample.origin)){
		sub_eset <- ccle_eset[, grep(sample.origin,sampleNames(ccle_eset))]
		sub_response <- pData(ccle_response[grep(sample.origin,sampleNames(ccle_response))])
	}else{
		sub_eset <- ccle_eset
		sub_response <- pData(ccle_response)
	}
	
	RASact <- loboda_enrichment_method(sub_eset)
	
	stats <- apply(sub_response, 2, function(x){ cor.test(RASact, x,method="spearman")$estimate })
	pvals <- apply(sub_response, 2, function(x){ cor.test(RASact, x,method="spearman")$p.value })
	
	df <- data.frame(model.rho=stats,
			model.pval=pvals)
	df
}
