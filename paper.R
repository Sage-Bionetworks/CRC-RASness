## PROGRAM PULLING TOGETHER ALL ANALYSIS STEPS
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

require(synapseClient)
require(rGithubClient)
require(parallel)
require(limma)
require(edgeR)
require(cqn)
require(IlluminaHumanMethylation27k.db)

## PULL IN THE PROJECT REPOSITORY IN ORDER TO SOURCE FUNCTIONS VIA GITHUB
rasRepo <- getRepo("Sage-Bionetworks/CRC-RASness")
sourceRepoFile(rasRepo, "util_functions.R")
sourceRepoFile(rasRepo, "data_functions.R")
sourceRepoFile(rasRepo, "model_functions.R")

# GLOBAL VARIABLES WHICH DEFINE CANONICAL MUTATIONS PER GENE
canonical.kras <- c("kras.12","kras.13","kras.61","kras.146")
canonical.braf <- c("braf.600")
canonical.hras <- c("hras.G12","hras.G13")
canonical.nras <- c("nras.G12","nras.G13", "nras.Q61")
canonical.pik3ca <- c("pik3ca.542","pik3ca.545","pik3ca.546","pik3ca.1047")



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
  
  pdf("figs/xeno_correction.pdf")
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
  
  pdf("figs/xenomodel_primaryandmouse_cetux.pdf",width=8,height=6)
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
  
  
  
  
  pdf("figs/mek_mouse_xenograft.pdf",width=6,height=4)
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
  
  pdf("figs/raspredict_luad.pdf",width=5,height=5)
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
  
  pdf("figs/raspredict_rnaseq_vs_agilent.pdf",width=6,height=6)
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
  tcga_eset <- tcga_eset[,!is.na(tcga_eset$kras)]
  
  kfsyscc_es_ras <- gsva(kfsyscc_eset, ras_gsets,min.sz=10,max.sz=500)$es
  kfsyscc_es <- t(scale(t(rbind(exprs(kfsyscc_es_ras), khambata_enrichment_method(kfsyscc_eset)))))
  khambata_es_ras <- gsva(khambata_eset,ras_gsets,min.sz=10,max.sz=500)$es
  khambata_es <- t(scale(t(rbind(exprs(khambata_es_ras), khambata_enrichment_method(khambata_eset)))))
  gaedcke_es_ras <- gsva(gaedcke_eset,ras_gsets,min.sz=10,max.sz=500)$es
  gaedcke_es <- t(scale(t(rbind(exprs(gaedcke_es_ras), khambata_enrichment_method(gaedcke_eset)))))
  tcga_es_ras <- gsva(tcga_eset, ras_gsets, min.sz=10,max.sz=500)$es
  tcga_es <- t(scale(t(rbind(exprs(tcga_es_ras), khambata_enrichment_method(tcga_eset)))))
  
  
  makeSigPlot <- function(sigIdx, title=FALSE,sig=""){
    names <- c("MU","WT")
    f1 <- formula(as.vector(kfsyscc_es[sigIdx,]) ~ factor(pData(kfsyscc_eset)$kras_status))
    f2 <- formula(as.vector(tcga_es[sigIdx,]) ~ factor(pData(tcga_eset)$kras==0))
    f3 <- formula(as.vector(khambata_es[sigIdx,]) ~ factor(pData(khambata_eset)$kras_status))
    f4 <- formula(as.vector(gaedcke_es[sigIdx,]) ~ factor(pData(gaedcke_eset)$kras_status))
    #par(mfrow=c(1,4))
    par(mar=c(0,0,0,0))
    plot(0:10,0:10, type="n",axes=FALSE)
    text(5,5,sig,cex=1.5)
    par(mar=c(3,2,2,2))
    par(cex.axis=.7)
    boxplot(f1,names=names,yaxt='n')
    if(title){ mtext("KFSYSCC",side=3,cex=.8)}
    mtext(paste("p=",format(kruskal.test(f1)$p.value,digits=2),sep=""),side=1,line=3,cex=.7)
    boxplot(f2,names=names,yaxt='n')
    if(title){ mtext("TCGA-CRC",side=3,cex=.8)}
    mtext(paste("p=",format(kruskal.test(f2)$p.value,digits=2),sep=""),side=1,line=3,cex=.7)
    boxplot(f3,names=names,yaxt='n')
    if(title){ mtext("Khambata-Ford",side=3,cex=.8)}
    mtext(paste("p=",format(kruskal.test(f3)$p.value,digits=2),sep=""),side=1,line=3,cex=.7)
    boxplot(f4,names=names,yaxt='n')
    if(title){ mtext("GAEDCKE",side=3,cex=.8)}
    mtext(paste("p=",format(kruskal.test(f4)$p.value,digits=2),sep=""),side=1,line=3,cex=.7)
  }
  
  
  
  pdf("figs/ras_sig_eval.pdf",width=5,5,useDingbats=FALSE)
  par(omi=c(.5,0,0,0))
  layout(matrix(1:20,4,5,byrow=TRUE))
  makeSigPlot(5,TRUE,sig="Loboda")
  makeSigPlot(3,sig="Meyerson")
  makeSigPlot(7,sig="Bild")
  makeSigPlot(8,sig="Baker")
  dev.off()
  
  
  make.roc.plot <- function(sigIdx,sig){
    
    kfsyscc.perf <- performance(prediction(as.vector(kfsyscc_es[sigIdx,]), 
                                           factor(kfsyscc_eset$kras_status == "MUT")),'tpr','fpr')
    tcga.perf <- performance(prediction(tcga_es[sigIdx,],
                                        factor(tcga_eset$kras==1)),'tpr','fpr')
    khambata.perf <- performance(prediction(khambata_es[sigIdx,],
                                            factor(khambata_eset$kras_status == "MUT")),'tpr','fpr')
    gaedcke.perf <- performance(prediction(gaedcke_es[sigIdx,],
                                           factor(gaedcke_eset$kras_status=="MUT")),'tpr','fpr')
    
    plot(tcga.perf, lty=1, lwd=1.5,col="red",main=sig)
    plot(khambata.perf, lty=1,lwd=1.5, col="cornflowerblue",add=TRUE)
    plot(gaedcke.perf, lty=1, lwd=1.5,col="darkgoldenrod1",add=TRUE)
    plot(kfsyscc.perf, lty=1, lwd=1.5,col="darkgreen",add=TRUE)
   
  }
  
  pdf("figs/ras_sig_eval_ROC.pdf",width=12,height=3,useDingbats=FALSE)
  par(mfrow=c(1,5))
  make.roc.plot(5,sig="Loboda")
  make.roc.plot(3,sig="Meyerson")
  make.roc.plot(7,sig="Bild")
  make.roc.plot(8,sig="Baker")
  plot(1:5,1:, type="n",axes=F,xlab="",ylab="")
  legend("topleft",
         legend=c("TCGA CRC",
                  "Khambata-Ford", 
                 "Gaedcke",
                  "KFSYSCC"),
         col=c("red","cornflowerblue","darkgoldenrod1","darkgreen"),
         lty=1,lwd=1.5)
  dev.off()
  
}

compute_kfsyscc_bootstrapped_aucs <- function(nbootstraps=100,num.processors=5){
  require(ROCR)
  
  set.seed(2012)
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
  
  pdf("figs/kfsyscc_bootstrapped_aucs.pdf",width=8,height=4,useDingbats=FALSE)
  par(mfrow=c(1,2))
  boxplot(aucs, ylab="AUC",notch=TRUE,col="gold",main="Bootstrapped AUC")
  boxplot(modelSizes, ylab="# Genes Selected",col="red", main="Bootstrapped Model Size")
  dev.off()
  #list(mean.auc=mean(aucs), quantile.auc=quantile(aucs,probs=c(.05, .95)),featureFreqTbl=df)
}

compute_kfsyscc_validation <- function(){
  
  kfsyscc_eset <- getKFSYSCC()
  khambata_eset <- getKhambata()
  gaedcke_eset <- getGaedcke()
  tcga_eset <- getTCGACRC()
  tcga_eset <- tcga_eset[,!is.na(tcga_eset$kras)]
  
  y_hats <- binomial_predict_EN(kfsyscc_eset, factor(kfsyscc_eset$kras_status=="MUT"), 
                                list(tcga_eset, khambata_eset, gaedcke_eset))$yhats
  
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
  
  pdf("figs/kfsyscc_kras_prediction_validtion.pdf",width=5,height=5,useDingbats=FALSE) 
  plot(tcga_perf$perf, lty=1, lwd=1.5,col="red")
  plot(khambata_perf$perf, lty=1,lwd=1.5, col="cornflowerblue",add=TRUE)
  plot(gaedcke_perf$perf, lty=1, lwd=1.5,col="darkgoldenrod1",add=TRUE)
  legend(.4, .6,
         legend=c("TCGA CRC",
                  "Khambata-Ford", 
                  "Gaedcke"),
         col=c("red","cornflowerblue","darkgoldenrod1"),
         lty=1,lwd=1.5)
  
  dev.off()
 
  tmp <- list(tcga_perf, khambata_perf,gaedcke_perf)
  
  df <- data.frame(auc=unlist(sapply(tmp, function(x) format(x$auc,digits=3))), 
                   sens=unlist(sapply(tmp, function(x) format(x$sens,digits=3))),
                   spec=unlist(sapply(tmp, function(x) format(x$spec,digits=3))),
                   row.names=c("TCGA CRC","Khambata-Ford","Gaedcke"))
  
  return(list(tcga=tcga_perf$auc,KF=khambata_perf$auc,gaedcke=gaedcke_perf$auc))
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
  
  
  colors <- rep("azure4",length(aa.factor))
  colors[as.character(aa.factor) %in% names(kras.aa.changes)] = "red"
  colors[as.character(aa.factor) %in% names(braf.aa.changes)] = "cornflowerblue"
  colors[as.character(aa.factor) %in% names(nras.aa.changes)] = "darkgoldenrod1"
  
  pdf("figs/tcga_aminoacid_ris_kras_braf.pdf",width=10,height=5,useDingbats=FALSE)
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
  pdf("figs/tcga_pik3ca_exon_ris.pdf",width=8,height=5,useDingbats=FALSE)
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
  pdf("figs/tcga_kras_braf_hras_pik3ca.pdf",width=10,height=5,useDingbats=FALSE)
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
  #all.genes <- unique(maf_wt$Hugo_Symbol)
  #pvals <- rep(NA, length(all.genes))
  #names(pvals) <- all.genes
  #for(i in 1:length(all.genes)){
  #  gene <- all.genes[i]
  #  tmp <- unique(maf_patids[maf_wt$Hugo_Symbol == gene])
  #  f <- factor(rasScore_wt_patids %in% tmp)
  #  if(length(levels(f)) == 1) { next }
  #  pvals[i] <- kruskal.test(rasScore_wt ~ f)$p.value
  #  if(i %% 100 == 0){ cat(".")}
  #}
  
  gene <- "CTNNB1"
  tmp <- unique(maf_patids[maf_wt$Hugo_Symbol == gene])
  f <- factor(rasScore_wt_patids %in% tmp)
  pdf("figs/ctnnb1_ris_assoc_in_WT.pdf",width=6,height=6,useDingbats=FALSE)
  boxplot(rasScore_wt ~ f,las=2, pars=list(outpch=NA,medlwd=.5),ylim=c(0,1), ylab="RIS",names=c("WT","MUT"),main='CTNNB1')
  stripchart(rasScore_wt ~ f,pch=19,cex=.8,col=c("darkgoldenrod1","cornflowerblue"),method="jitter",vertical=TRUE,add=TRUE)
  dev.off()
  
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
  pdf(paste("figs/ccle_ic50_rasact_association_",sample.origin,".pdf",sep=""),width=8,height=8)
  par(mar=c(6, 8, 2, 2))
  col=rep("black",length(pvals))
  col[pvals[idxs] < .05] <- "red"
  bp <- barplot(stats[idxs],horiz=T,las=2,xlim=c(-1, 1),col=col,xlab="Spearman rho")
  text(.7, bp,labels=format(pvals[idxs],digits=2))
  dev.off()
  
  sig.drugs <- pvals < .01
  glyphs <- rep(19, length(pvals))
  glyphs[sig.drugs] <- 1:sum(sig.drugs)
  pdf(paste("figs/ccle_ic50_rasact_association_",sample.origin,"_volcano.pdf",sep=""),width=8,height=8)
  plot(stats, -log10(pvals),xlab="Spearman rho",pch=glyphs)
  abline(h=-log10(.01),lty=2,col="red")
  dev.off()
  
  return (df)
}

compute_ccle_drug_response_metagenomics <- function(sample.origin="LARGE_INTESTINE"){
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
  
  pdf("figs/khambata_ford_cetuximab_response.pdf",width=10, height=4)
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
