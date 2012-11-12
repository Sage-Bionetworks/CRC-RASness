source("code/data_functions.R")
source("code/util_functions.R")

getAllLung <- function(){
  tcga.luad <- get.luad.RNAseq(tumor.only=TRUE)
  tmp <- get.luad.exome(tcga.luad)
  tcga.luad.eset <- tmp$eset
  exome <- tmp$exome
  
  kras.mut.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode[exome$Hugo_Symbol == "KRAS" & exome$Variant_Classification=="Missense_Mutation"]))
  patids <- extract.tcga.patientIds(sampleNames(tcga.luad.eset))
  tcga.luad.eset$kras <- patids %in% kras.mut.patients
  
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
  
  return (list(tcga.luad=tcga.luad.eset,battle=battle,gse26939=lungA,chemores=chemores))
}

getCHEMORES <- function(){
	load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_EXP.RData")
	load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/CHEMORES_CLIN.RData")
	idxs <- match(colnames(CHEMORES_EXP), rownames(CHEMORES_CLIN))
	eset <- new("ExpressionSet",exprs=CHEMORES_EXP[,!is.na(idxs)])
	pData(eset) <- CHEMORES_CLIN[na.omit(idxs),]
	
	eset <- eset[,!is.na(eset$KRAS)]
	eset
}

getBattle <- function(){
	load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/BATTLE_EXP.RData")
	load("/home/cferte/FELLOW/cferte/KRAS_Project/OBJECTS/BATTLE_CLIN.RData")
	
	
	#### select the BATTLE dataset and get rid of the patients with KRAS status = na
	table(BATTLE_CLIN$KRAS)
	identical(colnames(BATTLE_EXP),rownames(BATTLE_CLIN))
	tmp <- which(!is.na(BATTLE_CLIN$KRAS))
	BATTLE_CLIN <- BATTLE_CLIN[tmp,]
	BATTLE_EXP <- BATTLE_EXP[,tmp]
	
	tmp <- combine_probes_2_gene(BATTLE_EXP, rownames(BATTLE_EXP))
	
	eset <- new("ExpressionSet",exprs=tmp)
	pData(eset) <- BATTLE_CLIN
	eset
}

get.luad.RNAseq <- function(tumor.only=TRUE){
	
	getOrFetch("__tcga.luad.rnaseq",{
		env <- new.env()
		load("data~/luad/luad_rnaseq_v3.1.3.rda",env)
		M <- env$RPKM.cqn
		
		if(tumor.only){
			mask <- unlist(sapply(colnames(M), function(x) { grepl("01",strsplit(x,"-",fixed=TRUE)[[1]][4]) }))
			M <- M[,mask]
		}
		eset <- new("ExpressionSet", exprs=M)
		clinical <- read.table("data~/luad/clinical_patient_public_luad.txt",sep="\t", stringsAsFactor=FALSE, 
                           na.strings=c("NA","[Not Applicable]","[Not Available]"),header=TRUE)
		idxs <- match(extract.tcga.patientIds(sampleNames(eset)), clinical$bcr_patient_barcode)
    clinical.m <- clinical[idxs,]
    rownames(clinical.m) <- sampleNames(eset)
    pData(eset) <- clinical.m
    eset
	})
	
	
}

get.tcga.survival <- function(eset){
   surv <- Surv(as.numeric(eset$days_to_death), eset$vital_status == "DECEASED")
   surv[surv[,2] == 0,1] <- eset$days_to_last_followup[surv[,2] == 0]
   surv
}

get.luad.exome <- function(eset=NULL){
	exome <- getOrFetch("__tcga.luad.exome", {
		exome <- read.table("data~/luad/LUAD.exome.cleaned.somatic.maf",sep="\t",quote="",as.is=TRUE,header=TRUE)
		exome
	})

	if(!is.null(eset)){
		exome.patients <- unique(extract.tcga.patientIds(exome$Tumor_Sample_Barcode))
		eset.patients <- unique(extract.tcga.patientIds(sampleNames(eset)))
		joint <- intersect(exome.patients, eset.patients)
		exome.m <- exome[extract.tcga.patientIds(exome$Tumor_Sample_Barcode) %in% joint,]
		eset.m <- eset[, extract.tcga.patientIds(sampleNames(eset)) %in% joint]
		return (list(exome=exome.m, eset=eset.m))
	}else{
		return(exome)
	}
}