## DATA PREPROCESSING FUNCTIONS
## WHICH STORE INTERMEDIATE DATA IN SYNAPSE
##
## ORIGINATING ANALYST: JUSTIN GUINNEY
## SUPPORTING ANALYST: BRIAN BOT
#####

require(synapseClient)
require(Biobase)


#####
## TCGA AGILENT LEVEL 3 DATA
mergeAgilentLevel3 <- function(folderId){
  
  ## GET ILLUMINA AND SOLID PLATFORMS
  agilentFolder <- synapseQuery(paste("SELECT id, name FROM entity WHERE entity.parentId=='", folderId, "'", sep=""))
  agilentFiles <- sapply(as.list(agilentFolder$entity.id), function(synId){
    tmp <- downloadEntity(synId)
    return(file.path(tmp$cacheDir, grep("XXXXXXXXXXXXX", tmp$files, fixed=T)))
  })
  
  allGene <- lapply(as.list(agilentFiles), function(gene){
    read.delim(gene, quote="", as.is=T, header=T, comment.char="", skip=1, fill=T)
  })
  masterGene <- do.call(cbind, allGene)
  
  return(masterGene)
}


#####
## TCGA RNA-SEQ LEVEL 3 DATA
mergeRNAseqLevel3 <- function(folderId){
  
  ## GET ILLUMINA AND SOLID PLATFORMS
  seqFolder <- synapseQuery(paste("SELECT id, name FROM entity WHERE entity.parentId=='", folderId, "'", sep=""))
  geneFiles <- sapply(as.list(seqFolder$entity.id), function(synId){
    tmp <- downloadEntity(synId)
    return(file.path(tmp$cacheDir, grep("XXXXXXXXXXXXX", tmp$files, fixed=T)))
  })
  
  allGene <- lapply(as.list(geneFiles), function(gene){
    read.delim(gene, quote="", as.is=T, header=T, comment.char="", skip=1, fill=T)
  })
  masterGene <- do.call(cbind, allGene)
  
  return(masterGene)
}


#####
## TCGA DNA-SEQ LEVEL 2 DATA (MAF FILES)
mergeMafLevel2 <- function(folderId){
  
  ## GET ILLUMINA AND SOLID PLATFORMS
  seqFolder <- synapseQuery(paste("SELECT id, name FROM entity WHERE entity.parentId=='", folderId, "'", sep=""))
  mafFiles <- sapply(as.list(seqFolder$entity.id), function(synId){
    tmp <- downloadEntity(synId)
    return(file.path(tmp$cacheDir, grep(".maf", tmp$files, fixed=T)))
  })
  
  allMafs <- lapply(as.list(mafFiles), function(maf){
    read.delim(maf, quote="", as.is=T, header=T, comment.char="", skip=1, fill=T)
  })
  masterMaf <- do.call(rbind, allMafs)
  
  return(masterMaf)
}

