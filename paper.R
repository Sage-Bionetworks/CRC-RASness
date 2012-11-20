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
require(ROCR)
require(IlluminaHumanMethylation27k.db)

## PULL IN THE PROJECT REPOSITORY IN ORDER TO SOURCE FUNCTIONS VIA GITHUB
rasRepo <- getRepo("Sage-Bionetworks/CRC-RASness")
sourceRepoFile(rasRepo, "functions/utilityFunctions.R")
sourceRepoFile(rasRepo, "functions/dataFunctions.R")
sourceRepoFile(rasRepo, "functions/modelFunctions.R")

# GLOBAL VARIABLES WHICH DEFINE CANONICAL MUTATIONS PER GENE
canonicalKras <- c("kras.12","kras.13","kras.61","kras.146")
canonicalBraf <- c("braf.600")
canonicalHras <- c("hras.G12","hras.G13")
canonicalNras <- c("nras.G12","nras.G13", "nras.Q61")
canonicalPik3ca <- c("pik3ca.542","pik3ca.545","pik3ca.546","pik3ca.1047")



uterineAnalysis <- function(){
  uterineEsetEg <- loadEntity('syn585089')$objects$eset
  entrezIds <- gsub("(\\d*)_eg","\\1", featureNames(uterineEsetEg))
  genes <- unlist(mget(entrezIds, org.Hs.egSYMBOL, ifnotfound=NA))
  
  tmp <- combineProbesToGene(exprs(uterineEsetEg), genes)
  uterineEset <- new("ExpressionSet",exprs=tmp)
  
  mafEntity <- loadEntity('syn350455')
  maf <- read.table(file.path(mafEntity$cacheDir , mafEntity$files), header=TRUE, sep="\t", quote="")
  allPatients <- extractTcgaPatientIds(unique(maf$Tumor_Sample_Barcode))
  krasPatients <- extractTcgaPatientIds(unique(maf[maf$Hugo_Symbol=="KRAS",]$Tumor_Sample_Barcode))
  
  uterineEsetM <- uterineEset[, extractTcgaPatientIds(sampleNames(uterineEset)) %in% allPatients]
  
  ## NEED TO LOAD THE KFSYSCC DATA??
  yhat <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status == "MUT"), 
                              list(uterineEsetM))$yhats[[1]]
  pred <- prediction(yhat, extractTcgaPatientIds(sampleNames(uterineEsetM)) %in% krasPatients)
  auc <- performance(pred, 'auc')@y.values[[1]]
}

crcXenoGraft <- function(){
  xenoData <- getEMEXP991()
  kfsysccEset <- getKFSYSCC()
  
  eset <- xenoData$eset
  annot <- xenoData$annot
  fisher.test(factor(annot$Cetuximab==3), factor(annot$Kras == "MUT" | annot$Braf=="MUT"))
  
  # QC eset
  modelNames <- gsub("([A-Z]{2}?-[A-Z]{2,3}?-[0-9A-Z]{4}).*","\\1",eset$Source.Name)
  xenoPassage <- as.numeric(gsub("[A-Z]{2}?-[A-Z]{2,3}?-[0-9A-Z]{4}-P([0-9]).*","\\1",eset$Source.Name))
  isXeno <- !is.na(xenoPassage)
  
  pc <- svd(exprs(eset) - rowMeans(exprs(eset)))
  
  # remove top PC
  adjEset <- eset
  exprs(adjEset) <- pc$u %*% diag(c(0,pc$d[-c(1)])) %*% t(pc$v)
  featureNames(adjEset) <- featureNames(eset)
  pcAdj <- svd(exprs(adjEset) - rowMeans(exprs(adjEset)))
  
  pdf("figs/xenoCorrection.pdf")
  par(mfrow=c(2,2))
  plot(pc$d^2 / sum(pc$d^2), ylab="%var", xlab="eigenrank", main="Pre adjust")
  plot(pc$v[, 1], pc$v[, 2], col=(as.numeric(isXeno)+1), xlab="PC1", ylab='PC2')
  
  plot(pc.adj$d^2 / sum(pc.adj$d^2), ylab="%var", xlab="eigenrank", main="Post adjust")
  plot(pc.adj$v[, 1], pc.adj$v[, 2], col=(as.numeric(is.xeno)+1), xlab="PC1", ylab='PC2')
  dev.off()
  
  yhat <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status == "MUT"), 
                              list(adjEset))$yhats[[1]]
  
  ## all data
  idxs <- match(modelNames, annot$Code)
  annotModel <- annot[na.omit(idxs), ]
  xenoModel <- xenoPassage[!is.na(idxs)]
  earlyModelMask <- xenoModel < 4 & !is.na(xenoModel)
  lateModelMask <- xenoModel >= 4 & !is.na(xenoModel)
  yhatModel <- yhat[!is.na(idxs)]
  
  colors = rep("azure4", nrow(annot.model))
  colors[annotModel$Kras == "MUT" & annotModel$Braf == "WT"] <- "red"
  colors[annotModel$Kras == "WT" & annotModel$Braf == "MUT"] <- "darkgoldenrod1"
  colors[annotModel$Kras == "MUT" & annotModel$Braf == "MUT"] <- "green"
  
  krasBrafFactor <- factor(annotModel$Kras == "MUT" | annotModel$Braf=="MUT")
  responseFactor <- factor(annotModel$Cetuximab==3)
  fisher.test(factor(annot$Cetuximab==3), krasBrafFactor)
  
  pdf("figs/xenomodelPrimaryandmouseCetux.pdf", width=8, height=6)
  par(mfrow=c(2, 2), oma=c(0, 0, 0, 8))
  model <- yhatModel[is.na(xenoModel)] ~ factor(annotModel$Cetuximab==3)[is.na(xenoModel)]
  mutPval <- fisher.test(responseFactor[is.na(xenoModel)], krasBrafFactor[is.na(xenoModel)])$p.value
  boxplot(model, ylab="RIS", names=c("Non-response", "Response"),
          pars=list(outpch=NA, medlwd=.5), main=paste("tumor, n=", sum(is.na(xenoModel)), sep=""))
  points(model, pch=19, cex=.9, col=colors[!is.na(xenoModel)])
  ptext <- ifelse(mutPval > .5,  " (p>.5)", paste(" (p=", format(mutPval, digits=2),")", sep=""))
  mtext(side=3, paste("p=", format(kruskal.test(model)$p.value, digits=2), ptext,sep=""), cex=.8)
  
  
  model <- yhatModel[earlyModelMask] ~ factor(annotModel$Cetuximab==3)[earlyModelMask]
  mutPval <- fisher.test(responseFactor[earlyModelMask], krasBrafFactor[earlyModelMask])$p.value
  boxplot(model, ylab="RIS", names=c("Non-response", "Response"),
          pars=list(outpch=NA, medlwd=.5), main=paste("xeno early, n=", sum(earlyModelMask), sep=""))
  points(model, pch=19, cex=.9, col=colors[earlyModelMask])
  ptext <- ifelse(mutPval > .5,  " (p>.5)", paste(" (p=", format(mutPval, digits=2),")", sep=""))
  mtext(side=3, paste("p=", format(kruskal.test(model)$p.value, digits=2), ptext, sep=""), cex=.8)
  
  
  model <- yhatModel[lateModelMask] ~ factor(annotModel$Cetuximab==3)[lateModelMask]
  mutPval <- fisher.test(responseFactor[lateModelMask], krasBrafFactor[lateModelMask])$p.value
  boxplot(model, ylab="RIS", names=c("Non-response", "Response"),
          pars=list(outpch=NA, medlwd=.5), main=paste("xeno late, n=", sum(lateModelMask), sep=""))
  points(model, pch=19, cex=.9, col=colors[lateModelMask])
  ptext <- ifelse(mutPval > .5,  " (p>.5)", paste(" (p=", format(mutPval, digits=2),")", sep=""))
  mtext(side=3, paste("p=", format(kruskal.test(model)$p.value, digits=2), ptext, sep=""), cex=.8)
  
  
  model <- yhatModel ~ factor(annotModel$Cetuximab==3)
  mutPval <- fisher.test(responseFactor, krasBrafFactor)$p.value
  boxplot(model, ylab="RIS", names=c("Non-response", "Response"),
          pars=list(outpch=NA, medlwd=.5), ylim=c(0, 1), main=paste("tumor + xeno, n=", length(yhatModel), sep=""))
  points(model, pch=19, cex=.9, col=colors)
  ptext <- ifelse(mutPval > .5, " (p>.5)", paste(" (p=", format(mutPval, digits=2),")", sep=""))
  mtext(side=3, paste("p=", format(kruskal.test(model)$p.value, digits=2), ptext, sep=""), cex=.8)
  
  
  
  par(xpd=NA)
  legend(par("usr")[2] + .2, mean(par("usr")[3:4]), legend=c("kras","braf","kras+braf","wt"),
         fill=c("red","darkgoldenrod1","green","azure4"))
  par(xpd=FALSE)
  dev.off()
  
  
  
  boxplot(yhatModel ~ factor(annotModel$Cetuximab), ylab="RIS", main="Cetuximab response", names=c(">42%", "10-42%", "-10-10%", "<-10%"))
  mtext(side=3, paste("p=", format(kruskal.test(yhatModel ~ factor(annotModel$Cetuximab))$p.value, digits=2), sep=""))
  
  summary(glm(factor(annotModel$Cetuximab==3) ~ factor(annotModel$Kras == "MUT" | annotModel$Braf), family="binomial"))
  summary(glm(factor(annotModel$Cetuximab==3) ~ factor(annotModel$Kras == "MUT" | annotModel$Braf=="MUT" | annotModel$Pik3ca=="MUT") + (yhatModel > .35), family="binomial"))
  
}

mekInhibition <- function(){
  
  tmp <- getEMEXP3557()
  kfsysccEset <- getKFSYSCC()
  
  # 1 array is major outlier
  pc <- svd(exprs(tmp) - rowMeans(exprs(tmp)))
  badIdxs <- which(pc$v[,1] > .8)
  mekEset <- tmp[,-badIdxs]
  
  yhat <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$krasStatus == "MUT"), 
                              list(mekEset))$yhats[[1]]
  
  mutMask <- mekEset$Factor.Value.GENOTYPE == "Kras mutant"
  cntrlMask <- mekEset$Factor.Value.COMPOUND. == "control"
  plot(yhat, type="p", pch=as.numeric(factor(mekEset$Factor.Value.GENOTYPE.)),
       col=as.numeric(factor(mekEset$Factor.Value.COMPOUND.)))
  
  boxplot(yhat ~ factor(mekEset$Factor.Value.COMPOUND.),main="Kras mut",ylab="RIS")
  points(factor(mekEset$Factor.Value.COMPOUND.)[mutMask], yhat, col=c(as.numeric(factor(mekEset$Factor.Value.COMPOUND.)) + 1))
  mtext(side=3, paste("p=", format(kruskal.test(yhat ~ factor(mekEset$Factor.Value.COMPOUND.))$p.value,digits=2), sep=""))
  
  pdf("figs/mekMouseXenograft.pdf", width=6, height=4)
  par(mfrow=c(1, 2))
  boxplot(yhat[mutMask] ~ factor(mekEset$Factor.Value.COMPOUND.)[mutMask], main="Kras mut", ylab="RIS", ylim=c(.25, .55))
  points(factor(mekEset$Factor.Value.COMPOUND.)[mutMask], 
         yhat[mutMask], col=c(as.numeric(factor(mekEset$Factor.Value.COMPOUND.)[mutMask]) + 1), pch=19)
  mtext(side=3, paste("p=", format(kruskal.test(yhat[mutMask] ~ factor(mekEset$Factor.Value.COMPOUND.)[mutMask])$p.value, digits=2), sep=""))
  boxplot(yhat[!mutMask] ~ factor(mekEset$Factor.Value.COMPOUND.)[!mutMask], main="Kras wt", ylab="RIS", ylim=c(.25, .55))
  points(factor(mekEset$Factor.Value.COMPOUND.)[!mutMask], 
         yhat[!mutMask], col=c(as.numeric(factor(mekEset$Factor.Value.COMPOUND.)[!mutMask]) + 1), pch=19)
  mtext(side=3, paste("p=", format(kruskal.test(yhat[!mutMask] ~ factor(mekEset$Factor.Value.COMPOUND.)[!mutMask])$p.value, digits=2), sep=""))
  dev.off()
}

getRISandMut <- function(){
  ## FIX THIS
  rasness <- read.table("tcga_crc_ras_scores.txt",sep="\t",header=F)
  
  tcgaAaMat <- getTCGARasMuts()
  
  idxs <- match(extractTcgaPatientIds(rasness[,1]), rownames(tcgaAaMat))
  rasScoreM <- rasness[!is.na(idxs),2]
  names(rasScoreM) <- rasness[!is.na(idxs),1]
  tcgaAaMatM <- tcgaAaMat[na.omit(idxs),]
  
  canonical <- c(canonicalKras, canonicalBraf, canonicalNras)
  idxs <- pmatch(canonical, colnames(tcgaAaMatM))
  hasMut <- rowSums(tcgaAaMatM[,idxs]) > 0
  return(list(ris=rasScoreM, mut=hasMut))
}

testRasnessInLuad <- function(){
  kfsysccEset <- getKFSYSCC()
  env <- new.env()
  load("../AZLung/tcga_luad_eset.rda", env)
  luadEset <- env$eset
  
  hgncSymbols <- sapply(featureNames(luadEset), function(x){ strsplit(x,"|",fixed=TRUE)[[1]][1] } )
  uniqSymbols <- setdiff(unique(hgncSymbols), "?")
  idxs <- match(uniqSymbols, hgncSymbols)
  luadEset <- luadEset[idxs,]
  featureNames(luadEset) <- hgncSymbols[idxs]
  yhats <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status=="MUT"), 
                                list(luadEset))$yhats
  
  pred <- prediction(yhats[[1]], list(as.logical(luadEset$kras)))
  aucKras <- performance(pred, 'auc')@y.values
  
  pred <- prediction(yhats[[1]], list(as.logical(luadEset$ras)))
  aucRas <- performance(pred, 'auc')@y.values
  
  pdf("figs/raspredictLuad.pdf", width=5, height=5)
  plot(1:4, c(NULL, aucKras, aucRas, NULL), xaxt="n", col="red", pch=19, main="Ras prediction in LUAD", ylab="AUC", xlab="Platform/Normalization")
  axis(1, 1:4, labels=c(NULL, "KRAS", "KRAS/BRAF/NRAS", NULL))
  dev.off()
}

computeRnaseqRasness <- function(){
  kfsysccEset <- getKFSYSCC()
  tcgaEset <- getTCGACRC()
  tcgaAaMat <- getTCGARasMuts()
  canonical <- c(canonicalKras, canonicalBraf, canonicalNras)
  
  idxs <- pmatch(canonical, colnames(tcgaAaMat))
  hasMutMask <- rowSums(tcgaAaMat[,idxs]) > 0
  
  coadEnv <- new.env()
  load("data~/tcga_coad_rnaseq.rda",env=coad.env)
  readEnv <- new.env()
  load("data~/tcga_read_rnaseq.rda",env=read.env)
  stopifnot(rownames(coadEnv$RPKM) == rownames(readEnv$RPKM))
  
  RPKM <- cbind(coadEnv$RPKM, readEnv$RPKM)
  Counts <- cbind(coadEnv$Counts, readEnv$Counts)
  
  hgncSymbols <- sapply(rownames(RPKM), function(x){ strsplit(x,"|",fixed=TRUE)[[1]][1] } )
  uniqSymbols <- setdiff(unique(hgncSymbols), "?")
  idxs <- match(uniqSymbols, hgncSymbols)
  RPKM <- RPKM[idxs,]
  Counts <- Counts[idxs,]
  rownames(RPKM) <- uniqSymbols
  rownames(Counts) <- uniqSymbols
  
  RPKMquant <- normalizeBetweenArrays(RPKM, method="quantile")
  
  # cqn normalization
  geneAnnot <- read.table("data~/gene.annot.txt",sep="\t",header=TRUE)
  idxs <- match(rownames(Counts), geneAnnot$hgnc_symbol)
  CountsM <- Counts[!is.na(idxs),]
  geneAnnotM <- geneAnnot[na.omit(idxs),]
  
  cqn <- cqn(CountsM, lengths=geneAnnotM$length, x=geneAnnotM$gcpercent, sizeFactors=colSums(Counts), verbose=TRUE)
  RPKMcqn <- cqn$y + cqn$offset	
  
  makeRasFactor <- function(matrix){
    if(class(matrix)=="ExpressionSet"){
      eset <- matrix
    }else{
      # remove duplicate rows
      idxs <- match(unique(extractTcgaPatientIds(colnames(matrix))),
                    extractTcgaPatientIds(colnames(matrix)))
      eset <- new("ExpressionSet",exprs=matrix[, idxs])
    }
    
    idxs <- match(extractTcgaPatientIds(rownames(tcgaAaMat)), extractTcgaPatientIds(sampleNames(eset)))
    tcgaAaMatM <- tcgaAaMat[!is.na(idxs),]
    esetM <- eset[, na.omit(idxs)]
    idxs <- pmatch(canonical, colnames(tcgaAaMatM))
    y <- rowSums(tcgaAaMatM[, idxs]) > 0
    return(list(eset=esetM, y=y))
  }
  
  rpkm <- makeRasFactor(RPKM)
  rpkmQuant <- makeRasFactor(RPKMquant)
  rpkmCqn <- makeRasFactor(RPKMcqn)
  agilent <- makeRasFactor(tcgaEset)
  
  factors <- list(rpkm=factor(!rpkm$y), rpkmQuant=factor(!rpkmQuant$y), rpkmCqn=factor(!rpkmCqn$y), agilent=factor(!agilent$y))
  aucs <- mclapply(1:10, function(idx){
    idx <- idx
    N <- dim(kfsysccEset)[2]
    mask <- sample(N, round(.66 * N), replace=TRUE)
    yhats <- binomialPredictEN(kfsysccEset[, mask], factor(kfsysccEset$kras_status)[mask],
                                  testEsets=list(rpkm$eset, rpkmQuant$eset, rpkmCqn$eset, agilent$eset), seed=(idx+1), alpha=.1)$yhats
    aucs <- NULL
    for(i in 1:length(yhats)){
      pred <- prediction(yhats[[i]], factors[[i]])	
      auc <- performance(pred, 'auc')@y.values[[1]]
      aucs <- c(aucs, auc)
    }
    aucs
  }, mc.cores=10, mc.set.seed=TRUE, mc.preschedule =FALSE)
  aucs <- t(sapply(aucs, function(x) x))
  aucList <- lapply(1:4, function(x) { aucs[,x] } )		
  names(aucList) <- c("rpkm", "rpkmQuant", "rpkmCqn", "agilent")
  
  pdf("figs/raspredictRnaseqVsAgilent.pdf",width=6,height=6)
  boxplot(aucList, main="Platform AUCs", ylab="AUC", xlab="Platform/Normalization",
          col=c(rep("red", 3), "cornflowerblue"),
          ylim=c(.81, .90))
  legend("topright", col=c("red", "cornflowerblue"), legend=c("RNAseq", "microarray"), lty=1, lwd=3)
  dev.off()
  
  return(list(pvalDiff=kruskal.test(aucList)$p.value))
}


# compute enrichment of RAS signatures in colon cancer data sets
compute_ras_signature_enrichment <- function(){
  require(GSVA)
  
  rasGsets <- getRASSigs()
  kfsysccEset <- getKFSYSCC()
  khambataEset <- getKhambata()
  gaedckeEset <- getGaedcke()
  tcgaEset <- getTCGACRC()
  tcgaEset <- tcgaEset[,!is.na(tcgaEset$kras)]
  
  kfsysccEsRas <- gsva(kfsysccEset, rasGsets,min.sz=10,max.sz=500)$es
  kfsysccEs <- t(scale(t(rbind(exprs(kfsysccEsRas), khambataEnrichmentMethod(kfsysccEset)))))
  khambataEsRas <- gsva(khambataEset,rasGsets,min.sz=10,max.sz=500)$es
  khambataEs <- t(scale(t(rbind(exprs(khambataEsRas), khambataEnrichmentMethod(khambataEset)))))
  gaedckeEsRas <- gsva(gaedckeEset,rasGsets,min.sz=10,max.sz=500)$es
  gaedckeEs <- t(scale(t(rbind(exprs(gaedckeEsRas), khambataEnrichmentMethod(gaedckeEset)))))
  tcgaEsRas <- gsva(tcgaEset, rasGsets, min.sz=10,max.sz=500)$es
  tcgaEs <- t(scale(t(rbind(exprs(tcgaEsRas), khambataEnrichmentMethod(tcgaEset)))))
  
  
  makeSigPlot <- function(sigIdx, title=FALSE,sig=""){
    names <- c("MU","WT")
    f1 <- formula(as.vector(kfsysccEs[sigIdx,]) ~ factor(pData(kfsysccEset)$kras_status))
    f2 <- formula(as.vector(tcgaEs[sigIdx,]) ~ factor(pData(tcgaEset)$kras==0))
    f3 <- formula(as.vector(khambataEs[sigIdx,]) ~ factor(pData(khambataEset)$kras_status))
    f4 <- formula(as.vector(gaedckeEs[sigIdx,]) ~ factor(pData(gaedckeEset)$kras_status))
    #par(mfrow=c(1,4))
    par(mar=c(0, 0, 0, 0))
    plot(0:10, 0:10, type="n", axes=FALSE)
    text(5, 5, sig, cex=1.5)
    par(mar=c(3, 2, 2, 2))
    par(cex.axis=.7)
    boxplot(f1, names=names, yaxt='n')
    if(title){ mtext("KFSYSCC", side=3, cex=.8)}
    mtext(paste("p=", format(kruskal.test(f1)$p.value, digits=2), sep=""), side=1, line=3, cex=.7)
    boxplot(f2, names=names, yaxt='n')
    if(title){ mtext("TCGA-CRC", side=3, cex=.8)}
    mtext(paste("p=", format(kruskal.test(f2)$p.value, digits=2), sep=""), side=1, line=3, cex=.7)
    boxplot(f3, names=names, yaxt='n')
    if(title){ mtext("Khambata-Ford", side=3, cex=.8)}
    mtext(paste("p=", format(kruskal.test(f3)$p.value, digits=2), sep=""), side=1, line=3, cex=.7)
    boxplot(f4, names=names, yaxt='n')
    if(title){ mtext("GAEDCKE", side=3, cex=.8)}
    mtext(paste("p=", format(kruskal.test(f4)$p.value, digits=2), sep=""), side=1, line=3, cex=.7)
  }
  
  
  
  pdf("figs/rasSigEval.pdf", width=5, 5, useDingbats=FALSE)
  par(omi=c(.5, 0, 0, 0))
  layout(matrix(1:20, 4, 5, byrow=TRUE))
  makeSigPlot(5, TRUE, sig="Loboda")
  makeSigPlot(3, sig="Meyerson")
  makeSigPlot(7, sig="Bild")
  makeSigPlot(8, sig="Baker")
  dev.off()
  
  
  makeRocPlot <- function(sigIdx,sig){
    kfsysccPerf <- performance(prediction(as.vector(kfsysccEs[sigIdx,]), 
                                           factor(kfsysccEset$kras_status == "MUT")),'tpr','fpr')
    tcgaPerf <- performance(prediction(tcgaEs[sigIdx,],
                                        factor(tcgaEset$kras==1)),'tpr','fpr')
    khambataPerf <- performance(prediction(khambataEs[sigIdx,],
                                            factor(khambataEset$kras_status == "MUT")),'tpr','fpr')
    gaedckePerf <- performance(prediction(gaedckeEs[sigIdx,],
                                           factor(gaedckeEset$kras_status=="MUT")),'tpr','fpr')
    
    plot(tcgaPerf, lty=1, lwd=1.5, col="red", main=sig)
    plot(khambataPerf, lty=1,lwd=1.5, col="cornflowerblue", add=TRUE)
    plot(gaedckePerf, lty=1, lwd=1.5, col="darkgoldenrod1", add=TRUE)
    plot(kfsysccPerf, lty=1, lwd=1.5, col="darkgreen", add=TRUE)
  }
  
  pdf("figs/rasSigEvalROC.pdf",width=12,height=3,useDingbats=FALSE)
  par(mfrow=c(1, 5))
  makeRocPlot(5, sig="Loboda")
  makeRocPlot(3, sig="Meyerson")
  makeRocPlot(7, sig="Bild")
  makeRocPlot(8, sig="Baker")
  plot(1:5, 1:, type="n", axes=F, xlab="", ylab="")
  legend("topleft",
         legend=c("TCGA CRC",
                  "Khambata-Ford", 
                  "Gaedcke",
                  "KFSYSCC"),
         col=c("red", "cornflowerblue", "darkgoldenrod1", "darkgreen"),
         lty=1, lwd=1.5)
  dev.off()
  
}

computeKfsysccBootstrappedAucs <- function(nbootstraps=100, numProcessors=5){
  require(ROCR)
  
  set.seed(2012)
  eset <- getKFSYSCC()
  N <- ncol(eset)
  
  splitAndTest <- function(i){
    idxs <- sample(N, N/2)
    kras <- factor(eset$kras_status)
    r <- binomialPredictEN(eset[,idxs], kras[idxs], list(eset[,-idxs]))
    yhat <- r$yhats[[1]]
    selectedFeatures <- r$featureVec[as.logical(coefficients(r$model) != 0)]
    auc <- performance(prediction(yhat, factor(kras[-idxs])), 'auc')@y.values[[1]]
    list(auc=auc, sf=selectedFeatures)
  }
  rlist <- mclapply(1:nbootstraps, splitAndTest, mc.cores=numProcessors, mc.set.seed=TRUE, mc.preschedule =FALSE)
  aucs <- unlist(lapply(rlist, function(x){ x$auc }))
  modelSizes <- unlist(lapply(rlist, function(x){ length(x$sf) }))
  
  df <- data.frame(gene=featureNames(eset), freq=rep(0, nrow(eset)))
  
  for(r in rlist){
    mask <- df$gene %in% r$sf
    df$freq[mask] <- df$freq[mask] + 1
  }
  df$freq <- df$freq / nbootstraps
  df <- df[order(df$freq, decreasing=T), ]
  rownames(df) <- NULL
  
  pdf("figs/kfsysccBootstrappedAucs.pdf", width=8, height=4, useDingbats=FALSE)
  par(mfrow=c(1, 2))
  boxplot(aucs, ylab="AUC", notch=TRUE, col="gold", main="Bootstrapped AUC")
  boxplot(modelSizes, ylab="# Genes Selected", col="red", main="Bootstrapped Model Size")
  dev.off()
}

computeKfsysccValidation <- function(){
  
  kfsysccEset <- getKFSYSCC()
  khambataEset <- getKhambata()
  gaedckeEset <- getGaedcke()
  tcgaEset <- getTCGACRC()
  tcgaEset <- tcgaEset[,!is.na(tcgaEset$kras)]
  
  yhats <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status=="MUT"), 
                                list(tcgaEset, khambataEset, gaedckeEset))$yhats
  
  aucAndPerf <- function(yhat, actual){
    pred <- prediction(yhat, actual)
    perf <- performance(pred, 'tpr', 'fpr')
    auc <- performance(pred, 'auc')@y.values[[1]]
    
    idx <- order(performance(pred,"acc")@y.values[[1]],decreasing=T)[1]
    sens <- pred@tp[[1]][idx] / (pred@tp[[1]][idx] + pred@fn[[1]][idx])
    spec <- pred@tn[[1]][idx] / (pred@tn[[1]][idx] + pred@fp[[1]][idx])
    
    list(pred=pred, perf=perf, auc=auc, sens=sens, spec=spec)
  }
  
  tcgaPerf <- aucAndPerf(yhats[[1]], factor(tcgaEset$kras==1))
  khambataPerf <- aucAndPerf(yhats[[2]], factor(khambataEset$kras_status=="MUT"))
  gaedckePerf <- aucAndPerf(yhats[[3]], factor(gaedckeEset$kras_status=="MUT"))
  
  pdf("figs/kfsysccKrasPredictionValidtion.pdf", width=5, height=5, useDingbats=FALSE) 
  plot(tcgaPerf$perf, lty=1, lwd=1.5, col="red")
  plot(khambataPerf$perf, lty=1, lwd=1.5, col="cornflowerblue", add=TRUE)
  plot(gaedckePerf$perf, lty=1, lwd=1.5, col="darkgoldenrod1", add=TRUE)
  legend(.4, .6,
         legend=c("TCGA CRC",
                  "Khambata-Ford", 
                  "Gaedcke"),
         col=c("red", "cornflowerblue", "darkgoldenrod1"),
         lty=1, lwd=1.5)
  
  dev.off()
 
  tmp <- list(tcgaPerf, khambataPerf, gaedckePerf)
  
  df <- data.frame(auc=unlist(sapply(tmp, function(x) format(x$auc, digits=3))), 
                   sens=unlist(sapply(tmp, function(x) format(x$sens, digits=3))),
                   spec=unlist(sapply(tmp, function(x) format(x$spec, digits=3))),
                   row.names=c("TCGA CRC", "Khambata-Ford", "Gaedcke"))
  
  return(list(tcga=tcgaPerf$auc, KF=khambataPerf$auc, gaedcke=gaedckePerf$auc))
}

computeTcgaRasactDifferencesForRasAminoacids <- function(){
  kfsysccEset <- getKFSYSCC()
  tcgaEset <- getTCGACRC()
  maf <- getCleanTcgaMaf()
  rasScore <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status=="MUT"), 
                                  list(tcgaEset))$yhats[[1]]
  
  
  overlapPats <- intersect(extractTcgaPatientIds(sampleNames(tcgaEset)), 
                           unique(extractTcgaPatientIds(unique(maf$Tumor_Sample_Barcode))))
  
  tcgaEsetM <- tcgaEset[,extractTcgaPatientIds(sampleNames(tcgaEset)) %in% overlapPats]
  rasScoreM <- rasScore[extractTcgaPatientIds(sampleNames(tcgaEset)) %in% overlapPats, 1]
  mafM <- maf[extractTcgaPatientIds(maf$Tumor_Sample_Barcode) %in% overlapPats, ]
  
  krasPats <- extractTcgaPatientIds(mafM[mafM$Hugo_Symbol == "KRAS",]$Tumor_Sample_Barcode)
  nrasPats <- extractTcgaPatientIds(mafM[mafM$Hugo_Symbol == "NRAS",]$Tumor_Sample_Barcode)
  brafPats <- extractTcgaPatientIds(mafM[mafM$Hugo_Symbol == "BRAF",]$Tumor_Sample_Barcode)
  rasPats <- c(krasPats, nrasPats, brafPats)
  brafKrasWt <- setdiff(brafPats, krasPats)
  pik3caPats <- extractTcgaPatientIds(mafM[mafM$Hugo_Symbol == "PIK3CA",]$Tumor_Sample_Barcode)
  pik3caRasWt <- setdiff(pik3caPats, rasPats)
  pik3caRasMut <- intersect(pik3caPats, rasPats)
  
  ##### KRAS,NRAS,BRAF
  
  krasMaf <- mafM[mafM$Hugo_Symbol == "KRAS" ,]
  krasMaf$AAChange <- paste("kras.", krasMaf$AAChange, sep="")
  brafMaf <- mafM[mafM$Hugo_Symbol == "BRAF", ]
  brafMaf$AAChange <- paste("braf.", brafMaf$AAChange, sep="")
  nrasMaf <- mafM[mafM$Hugo_Symbol == "NRAS", ]
  nrasMaf$AAChange <- paste("nras.", nrasMaf$AAChange, sep="")
  
  krasAaChanges <- sort(table(krasMaf$AAChange), decreasing=TRUE)
  krasAaChanges <- krasAaChanges[names(krasAaChanges) != "kras.."]
  brafAaChanges <- sort(table(brafMaf$AAChange), decreasing=TRUE)
  nrasAaChanges <- sort(table(nrasMaf$AAChange), decreasing=TRUE)
  
  bothMaf <- rbind(krasMaf, brafMaf, nrasMaf)
  aaChanges <- c(krasAaChanges, brafAaChanges, nrasAaChanges)
  
  aaFactor <- rep("WT", ncol(tcgaEsetM))
  names(aaFactor) <- extractTcgaPatientIds(sampleNames(tcgaEsetM))
  for(aa in names(aaChanges)){
    pats <- extractTcgaPatientIds(bothMaf[bothMaf$AAChange == aa,]$Tumor_Sample_Barcode)
    aaFactor[pats] <- aa
  }
  
  aaFactor <- factor(aaFactor, levels=c(names(aaChanges), "WT"))
  
  colors <- rep("azure4", length(aaFactor))
  colors[as.character(aaFactor) %in% names(krasAaChanges)] = "red"
  colors[as.character(aaFactor) %in% names(brafAaChanges)] = "cornflowerblue"
  colors[as.character(aaFactor) %in% names(nrasAaChanges)] = "darkgoldenrod1"
  
  pdf("figs/tcgaAminoacidRisKrasBraf.pdf", width=10, height=5, useDingbats=FALSE)
  par(mar=c(8, 4, 2, 8))
  boxplot(rasScoreM ~ aaFactor, las=2, pars=list(outpch=NA, medlwd=.5), ylim=c(0, 1), ylab="RIS")
  points(rasScoreM ~ aaFactor, pch=19, cex=.7, col=colors)
  abline(h= quantile(rasScoreM[aaFactor=="WT"])[4], lty=2)
  par(xpd=TRUE)
  legend(length(levels(aaFactor))+3, 1, legend=c("kras", "braf", "nras", "wt"),
         fill=c("red", "cornflowerblue", "darkgoldenrod1", "azure4"), xjust=0)
  par(xpd=FALSE)
  dev.off()
  
  ## PIK3CA
  pik3caMaf <- mafM[mafM$Hugo_Symbol == "PIK3CA" ,]
  pik3caExonFactor <- rep("WT", ncol(tcgaEsetM))
  pik3caExonChanges <- sort(table(pik3caMaf$Exon),decreasing=TRUE)
  names(pik3caExonFactor) <- extractTcgaPatientIds(sampleNames(tcgaEsetM))
  for(exon in names(pik3caExonChanges)){
    pats <- extractTcgaPatientIds(pik3caMaf[pik3caMaf$Exon == exon,]$Tumor_Sample_Barcode)
    pik3caExonFactor[pats] <- exon
  }
  
  pik3caExonFactor <- factor(pik3caExonFactor, levels=c(names(pik3caExonChanges), "WT"))
  pik3caExonFactor[pik3caExonFactor == "WT" & aaFactor != "WT"] <- NA
  rasWtMask <- pik3caExonFactor != "WT" & aaFactor == "WT"
  rasMutMask <- pik3caExonFactor != "WT" & aaFactor != "WT"
  names(pik3caExonFactor[rasMutMask]) <- paste("rasMut", names(pik3caExonFactor[rasMutMask]), sep="")
  names(pik3caExonFactor[rasWtMask]) <- paste("rasWt",names(pik3caExonFactor[rasWtMask]), sep="")
  
  model <- rasScoreM[pik3caExonFactor != "WT"] ~ (aaFactor != "WT")[pik3caExonFactor != "WT"]
  boxplot(model)
  kruskal.test(model)
  
  colors = rep("azure4", length(pik3caExonFactor))
  colors[pik3caExonFactor != "WT"] <- "cornflowerblue"
  colors[pik3caExonFactor != "WT" & aaFactor != "WT"] <- "red"
  pdf("figs/tcgaPik3caExonRis.pdf", width=8, height=5, useDingbats=FALSE)
  boxplot(rasScoreM ~ factor(pik3caExonFactor), las=2, pars=list(outpch=NA, medlwd=.5), ylim=c(0, 1), ylab="RIS")
  points(rasScoreM ~ factor(pik3caExonFactor), pch=19, cex=.6, col=colors)
  abline(h= quantile(na.omit(rasScoreM[pik3caExonFactor=="WT"]))[4], lty=2)
  legend(6, .95, legend=c("kras/braf wt", "kras/braf mut", "all wt"),
         fill=c("cornflowerblue", "red", "azure4"))
  dev.off()
  
  tmpFactor <- rep("WT", ncol(tcgaEsetM))
  names(tmpFactor) <- extractTcgaPatientIds(sampleNames(tcgaEsetM))
  tmpFactor[names(tmpFactor) %in% krasPats] <- "kras"
  tmpFactor[names(tmpFactor) %in% brafPats] <- "braf"
  tmpFactor[names(tmpFactor) %in% nrasPats] <- "nras"
  tmpFactor[names(tmpFactor) %in% pik3caRasMut] <- "pik3ca/ras+"
  tmpFactor[names(tmpFactor) %in% pik3caRasWt] <- "pik3ca/ras-"
  colors = rep("azure4",length(tmpFactor))
  colors[tmpFactor != "WT"] <- "azure4"
  colors[tmpFactor == "kras"] <- "red"
  colors[tmpFactor == "braf"] <- "cornflowerblue"
  colors[tmpFactor == "nras"] <- "darkgoldenrod1"
  colors[tmpFactor == "pik3ca/ras+"] <- "aquamarine3"
  colors[tmpFactor == "pik3ca/ras-"] <- "aquamarine4"
  
  factorColor <- c("cornflowerblue", "red", "darkgoldenrod1", "aquamarine3", "aquamarine4", "azure4")
  pdf("figs/tcgaKrasBrafHrasPik3ca.pdf", width=10, height=5, useDingbats=FALSE)
  par(mar=c(8, 4, 2, 10))
  boxplot(rasScoreM ~ factor(tmpFactor), las=2, pars=list(outpch=NA, medlwd=.5), ylim=c(0, 1), ylab="RIS")
  stripchart(rasScoreM ~ factor(tmpFactor), pch=19, cex=.8, col=factor_color, method="jitter", vertical=TRUE, add=TRUE)
  abline(h= quantile(na.omit(rasScoreM[tmpFactor=="WT"]))[4], lty=2)
  par(xpd=TRUE)
  legend(6.8, 1, legend=c("braf", "kras", "nras", "pik3ca/ras-", "pik3ca/ras+", "wt"),
         fill=c("cornflowerblue", "red", "darkgoldenrod1", "aquamarine3", "aquamarine4", "azure4"))
  par(xpd=FALSE)
  dev.off()  
  
  lvls <- levels(factor(tmpFactor))
  M <- matrix(NA, nrow=5,ncol=5)
  for(i in 1:4){
    for(j in (i+1):5){
      mask <- tmpFactor %in% lvls[c(i,j)]
      M[i,j] <- kruskal.test(rasScoreM[mask] ~ factor(tmpFactor[mask]))$p.value
    }
  }
  
  #################################
  # look for novel mutations associated in Ras WT samples
  wtPats <- names(aaFactor)[aaFactor == "WT"]
  mafWt <- mafM[extractTcgaPatientIds(mafM$Tumor_Sample_Barcode) %in% wtPats 
                  & mafM$Variant_Classification %in% c("MissenseMutation", "NonsenseMutation"),]
  mafPatids <- extractTcgaPatientIds(mafWt$Tumor_Sample_Barcode)
  rasScoreWt <- rasScoreM[aaFactor == "WT"]
  rasScoreWtPatids <- extractTcgaPatientIds(names(rasScoreWt))
  
  gene <- "CTNNB1"
  tmp <- unique(mafPatids[mafWt$Hugo_Symbol == gene])
  f <- factor(rasScoreWtPatids %in% tmp)
  pdf("figs/ctnnb1RisAssocInWt.pdf", width=6, height=6, useDingbats=FALSE)
  boxplot(rasScoreWt ~ f,las=2, pars=list(outpch=NA, medlwd=.5), ylim=c(0, 1), ylab="RIS", names=c("WT", "MUT"), main='CTNNB1')
  stripchart(rasScoreWt ~ f, pch=19, cex=.8, col=c("darkgoldenrod1", "cornflowerblue"), method="jitter", vertical=TRUE, add=TRUE)
  dev.off()
}

computeCcleDrugResponse <- function(sampleOrigin="LARGE_INTESTINE"){
  kfsysccEset <- getKFSYSCC()
  ccle <- getCCLE()
  
  ccleEset <- ccle[[1]]
  ccleDrug <- ccle[[2]]
  idxs <- match(sampleNames(ccleEset), sampleNames(ccleDrug))
  
  ccleEset <- ccleEset[,!is.na(idxs)]
  ccleResponse <- ccleDrug[na.omit(idxs)]
  
  if(!is.null(sampleOrigin)){
    subEset <- ccleEset[, grep(sampleOrigin, sampleNames(ccleEset))]
    subResponse <- pData(ccleResponse[grep(sampleOrigin, sampleNames(ccleResponse))])
  }else{
    subEset <- ccleEset
    subResponse <- pData(ccleResponse)
  }
  RASact <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status=="MUT"), 
                                list(subEset))$yhats[[1]]
  
  stats <- apply(subResponse[mask,], 2, function(x){ cor.test(RASact[mask], x, method="spearman")$estimate })
  pvals <- apply(subResponse[mask,], 2, function(x){ cor.test(RASact[mask], x, method="spearman")$p.value })
  hasRas <- length(levels(factor(subEset$KRAS))) > 1
  hasRaf <- length(levels(factor(subEset$BRAF))) > 1
  
  pvals <- apply(subResponse[mask,], 2, function(x){ cor.test(RASact[mask], x, method="spearman")$p.value })
  pvals <- apply(subResponse[mask,], 2, function(x){ 
    coef(summary(lm(RASact[mask] ~ x + has.mut[mask])))[2,4] })
  
  if(hasRas){
    mutKrasPval=apply(subResponse, 2, function(x) { kruskal.test(x ~ factor(subEset$KRAS))$p.value } )
  }else{
    mutKrasPval=rep(NA, length(pvals))
  }
  
  if(hasRaf){
    mutBrafPval=apply(subResponse, 2, function(x) { kruskal.test(x ~ factor(subEset$BRAF))$p.value } )
  }else{
    mutBrafPval <- rep(NA, length(pvals))
  }
  if(hasRas | hasRaf){
    mutKrasBrafPval=apply(subResponse, 2, function(x) { kruskal.test(x ~ factor(subEset$KRAS | subEset$BRAF))$p.value } )
  }else{
    mutKrasBrafPval <- rep(NA, length(pvals))
  }
  df <- data.frame(modelRho=stats,
                   modelPval=pvals,
                   mutKrasPval=mutKrasPval,
                   mutBrafPval=mutBrafPval,
                   mutKrasBrafPval=mutKrasBrafPval)
  
  idxs <- order(stats)
  pdf(paste("figs/ccleIc50RasactAssociation-", sampleOrigin, ".pdf", sep=""), width=8, height=8)
  par(mar=c(6, 8, 2, 2))
  col=rep("black", length(pvals))
  col[pvals[idxs] < .05] <- "red"
  bp <- barplot(stats[idxs], horiz=T, las=2, xlim=c(-1, 1), col=col, xlab="Spearman rho")
  text(.7, bp, labels=format(pvals[idxs], digits=2))
  dev.off()
  
  sigDrugs <- pvals < .01
  glyphs <- rep(19, length(pvals))
  glyphs[sigDrugs] <- 1:sum(sigDrugs)
  pdf(paste("figs/ccleIc50RasactAssociation-", sampleOrigin, "-volcano.pdf", sep=""), width=8, height=8)
  plot(stats, -log10(pvals), xlab="Spearman rho", pch=glyphs)
  abline(h=-log10(.01), lty=2, col="red")
  dev.off()
  
  return(df)
}

computeCcleDrugResponseMetagenomics <- function(sampleOrigin="LARGE_INTESTINE"){
  kfsysccEset <- getKFSYSCC()
  ccle <- getCCLEmetaGenomics()
  
  if(!is.null(sampleOrigin)){
    subEset <- ccle[, grep(sampleOrigin, sampleNames(ccle))]
  }else{
    subEset <- ccle
  }
  RASact <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status=="MUT"), 
                                list(subEset))$yhats[[1]]
  
  stats <- apply(pData(subEset), 2, function(x){ cor.test(RASact, x, method="spearman")$estimate })
  pvals <- apply(pData(subEset), 2, function(x){ cor.test(RASact, x, method="spearman")$p.value })
  
  ##TODO .....
}

predictCetuximabResponse <- function(){
  kfsysccEset <- getKFSYSCC()
  khambataEset <- getKhambata()
  
  yhat <- binomialPredictEN(kfsysccEset, factor(kfsysccEset$kras_status=="MUT"), 
                               list(khambataEset))$yhats[[1]]
  notUtdMask <- khambataEset$Best != "UTD"
  wtMask <- khambataEset$kras_status == "WT"
  
  grps <- factor(khambataEset$Best[notUtdMask]=="PD")
  byRasScore <- kruskal.test(yhat[notUtdMask] ~ grps)$p.value
  byKrasStatus <- fisher.test(factor(khambataEset$krasStatus[notUtdMask]), grps)$p.value
  
  byRasScoreWt <- kruskal.test(yhat[notUtdMask][wtMask] ~ grps[wtMask])$p.value
  byRasScoreMut <- kruskal.test(yhat[notUtdMask][!wtMask] ~ grps[!wtMask])$p.value
  
  pdf("figs/khambataFordCetuximabResponse.pdf",width=10, height=4)
  par(mfrow=c(1,3))
  model <- yhat[notUtdMask] ~ factor(khambataEset$Best[notUtdMask])
  boxplot(model,
          pars=list(outpch=NA, medlwd=.5),
          ylab="RIS",
          main="All Samples")
  stripchart(model, 
             vertical=TRUE,
             col=c("gold", "cornflowerblue", "cornflowerblue"),
             pch=19, cex=1,
             method = "jitter", add=TRUE)
  mtext(paste("p=", format(byRasScore, digits=3), sep=""), side=1, line=3)
  
  model <- yhat[notUtdMask][wtMask] ~ factor(khambataEset$Best[notUtdMask][wtMask])
  boxplot(model,
          pars=list(outpch=NA, medlwd=.5),
          ylab="RIS",
          main="Kras wild-type")
  stripchart(model, 
             vertical=TRUE,
             col=c("gold", "cornflowerblue", "cornflowerblue"),
             pch=19, cex=1,
             method = "jitter", add=TRUE)
  mtext(paste("p=", format(byRasScoreWt, digits=2), sep=""), side=1, line=3)
  
  model <- yhat[notUtdMask][!wtMask] ~ factor(khambataEset$Best[notUtdMask][!wtMask])
  boxplot(model,
          pars=list(outpch=NA, medlwd=.5),
          ylab="RIS",
          main="Kras mutant")
  stripchart(model, 
             vertical=TRUE,
             col=c("gold", "cornflowerblue", "cornflowerblue"),
             pch=19, cex=1,
             method = "jitter", add=TRUE)
  mtext(paste("p=", format(byRasScoreMut, digits=2), sep=""), side=1, line=3)
  dev.off()
  
  return (list(pvalByRasScore=byRasScore, 
               pvalByKrasStatus=byKrasStatus, 
               pvalByRasScoreWt=byRasScoreWt))
}

checkLobodaCCLEPredictions <- function(sampleOrigin=NULL){
  ccle <- getCCLE()
  
  ccleEset <- ccle[[1]]
  ccleDrug <- ccle[[2]]
  idxs <- match(sampleNames(ccleEset), sampleNames(ccleDrug))
  
  ccleEset <- ccleEset[,!is.na(idxs)]
  ccleResponse <- ccleDrug[na.omit(idxs)]
  
  if(!is.null(sampleOrigin)){
    subEset <- ccleEset[, grep(sampleOrigin, sampleNames(ccleEset))]
    subResponse <- pData(ccleResponse[grep(sampleOrigin, sampleNames(ccleResponse))])
  }else{
    subEset <- ccleEset
    subResponse <- pData(ccleResponse)
  }
  
  RASact <- lobodaEnrichmentMethod(subEset)
  
  stats <- apply(subResponse, 2, function(x){ cor.test(RASact, x, method="spearman")$estimate })
  pvals <- apply(subResponse, 2, function(x){ cor.test(RASact, x, method="spearman")$p.value })
  
  df <- data.frame(modelRho=stats,
                   modelPval=pvals)
  return(df)
}
