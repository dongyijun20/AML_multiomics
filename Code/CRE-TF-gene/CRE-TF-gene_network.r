library(ArchR)
library(Signac)
library(Seurat)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(SingleCellExperiment)
library(cicero)
library(GenomicRanges)
library(tidyr)
library(stringr)
library(chromVARmotifs)
library(reshape2)
library(motifmatchr)
library(SummarizedExperiment)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

addArchRThreads(threads = 20)
addArchRGenome("hg38") #### choose your own species
h5disableFileLocking()
path <- "./project/ATAC/ArchRproj"
proj <- loadArchRProject(path=path)

### Add GeneScore
proj <- addImputeWeights(proj)
seGS <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
counts_RNA <- imputeMatrix(assay(seGS), getImputeWeights(proj))
rownames(counts_RNA) <- rowData(seGS)$name
obj_RNA <- CreateSeuratObject(
  counts = counts_RNA,
  assay = "RNA"
)
seRNA<-as.SingleCellExperiment(obj_RNA)

### Add peak-to-peak correlation
proj <- addCoAccessibility(
  ArchRProj = proj,
  maxDist = 5000000,
  reducedDims = "IterativeLSI"
)
cA <- getCoAccessibility(
  ArchRProj = proj,
  corCutOff = -1,
  resolution = 1,
  returnLoops = TRUE
)
p2p_all <- as.data.frame(cA[[1]],row.names = NULL)
p2p_all <- p2p_all[p2p_all$FDR<0.01,] #### then use it to calculate the insulation score using fanc under the ./insulation folder

### Add peak-to-gene correlation
proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA)
proj <- addIterativeLSI(
ArchRProj = proj, 
clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
),
saveIterations = FALSE,
useMatrix = "GeneExpressionMatrix", 
depthCol = "Gex_nUMI",
varFeatures = 2500,
firstSelection = "variable",
binarize = FALSE,
name = "LSI_RNA"
)
proj <- addPeak2GeneLinks(
ArchRProj = proj,
reducedDims = "IterativeLSI",
useMatrix = "GeneExpressionMatrix"
)

### Using GLASSO after TAD estimation
grs <- generate_windows(window, genomic_coords)
TAD <- read.table("insulation/100kb_all.insulation_boundaries_500kb.bed")
TAD <- TAD[,c(1:3,5)]
colnames(TAD) <- c("seqnames","start","end","insulation_score")
TAD <- nonOverlappingGR(GRanges(TAD[,1:3],insulation_score=TAD$insulation_score),by = "insulation_score")
TAD$TAD <- c(1:length(TAD))

knnIteration = 500
k = 100
overlapCutoff = 0.8
log2Norm = TRUE
scaleTo = 10^4

rD <- getReducedDims(proj, reducedDims = "IterativeLSI", corCutOff = 0.75, dimsToUse = 1:30)
idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
knnObj <- computeKNN(data = rD, query = rD[idx,], k = k)
keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
knnObj <- knnObj[keepKnn==0,]
knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
  rownames(rD)[knnObj[x, ]]
}) %>% SimpleList

peakSet <- getPeakSet(proj)
chri <- gtools::mixedsort(availableChr(getArrowFiles(proj), subGroup = "PeakMatrix"))
cS <- .getColSums(getArrowFiles(proj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
cS <- cS[proj$cellNames]
gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))

groupMat <- list()
for(x in chri){
  featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == x),]
  featureDF$seqnames <- x
  peak_chr <- peakSet[BiocGenerics::which(seqnames(peakSet) == x),]
  names(peak_chr) <- NULL
  peak_chr <- as.data.frame(peak_chr)
  peak_chr$mean_bp <-
    (as.numeric(as.character(peak_chr$start)) +
       as.numeric(as.character(peak_chr$end)))/2
  peak_chr <- unite(as.data.frame(peak_chr)[,c("seqnames","mean_bp")],col = "peaks",sep = "_")
  groupMat[[x]] <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(proj), 
    featureDF = featureDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    verbose = FALSE
  )
  groupMat[[x]] <- t(t(groupMat[[x]]) / gS) * scaleTo
  if(log2Norm){
    groupMat[[x]] <- log2(groupMat[[x]] + 1)
  }
  rownames(groupMat[[x]]) <- peak_chr$peaks
}

peaks <- lapply(chri,function(x){
  peaks <- proj@peakSet[proj@peakSet@seqnames==x]
  names(peaks) <- NULL
  peaks <- as.data.frame(peaks)[,1:5]
  colnames(peaks) <- c("chr","bp1","bp2","width","strand")
  peaks
})
names(peaks) <- chri

boundaries <- lapply(chri,function(x){
  boundaries <- TAD[TAD@seqnames==x]
  boundaries <- as.data.frame(boundaries)
  colnames(boundaries) <- c("chr","bp1","bp2","width","strand","insulation_score","TAD")
  boundaries
})
names(boundaries) <- chri

outlist <- lapply(seq_len(length(grs)),function(win){
  if(win%%100==0){
    message(paste("Compute", win, "windows",sep=" "))
  }

  if(as.character(grs[win]@seqnames) %in% chri == F){
    return("Unmatched chromosome")
    next()
  }
  
  peak <- peaks[[as.character(grs[win]@seqnames)]]
  boundary <- boundaries[[as.character(grs[win]@seqnames)]]
  win_range <- get_genomic_range(grs, peak, win)
  win_summits <- tidyr::unite(win_range[,c("chr","mean_bp")],col="peak",sep="_")
  boundaries_range <- get_genomic_range(grs, boundary, win)
  
  if(nrow(win_range) <= 1){
    return("Zero or one element in range")
    next()
  }
  
  if(nrow(boundaries_range) == 0){
    rho_mat <- matrix(0,nrow=nrow(win_range),ncol=nrow(win_range))
    expr <- groupMat[[as.character(grs[win]@seqnames)]]
    vals <- expr[win_summits$peak,]
    cov_mat <- cov(t(vals))
    diag(cov_mat) <- diag(cov_mat) + 1e-4
    
    GL <- glasso::glasso(cov_mat, rho_mat)
    
    colnames(GL$w) <- row.names(GL$w) <- win_summits$peak
    colnames(GL$wi) <- row.names(GL$wi) <- win_summits$peak

    return(GL)
    next()
  } else {
    dist_mat <- calc_dist_matrix(win_range,boundaries_range)
    rho_mat <- get_rho_mat(dist_mat, distance_parameter = 1)
    
    expr <- groupMat[[as.character(grs[win]@seqnames)]]
    vals <- expr[win_summits$peak,]
    
    cov_mat <- matrix(0,nrow=nrow(dist_mat),ncol=ncol(dist_mat))
    cov_mat[-which(rownames(dist_mat) == "b"),-which(rownames(dist_mat) == "b")] <- cov(t(vals))
    cov_mat[which(rownames(dist_mat) == "b"),] <- 0
    cov_mat[,which(rownames(dist_mat) == "b")] <- 0
    diag(cov_mat) <- diag(cov_mat) + 1e-4
    
    GL <- glasso::glasso(cov_mat, rho_mat)
    
    GL$w <- GL$w[-which(rownames(dist_mat) == "b"),-which(rownames(dist_mat) == "b")]
    GL$wi <- GL$wi[-which(rownames(dist_mat) == "b"),-which(rownames(dist_mat) == "b")]
    colnames(GL$w) <- row.names(GL$w) <- win_summits$peak
    colnames(GL$wi) <- row.names(GL$wi) <- win_summits$peak

    return(GL)
    next()
  }
})
results <- assemble_connections(outlist,groupMat,knnObj,0.01)

### Calculate and annotate peak-to-peak links
results<-results[results$FDR<0.01,]
p2p_links<-results[,1:3]
colnames(p2p_links)[3]<-"value"

peak<-as.data.frame(GRanges(unique(c(sub("_",":",p2p_links$Peak1),sub("_",":",p2p_links$Peak2)))))[,1:3]
write.table(peak,file = "./tmp/p2p_peaks.bed",sep="\t",quote = F,row.names = F,col.names = F)
promoter_distance = 1000
x = annotatePeak("./tmp/p2p_peaks.bed", tssRegion=c(-promoter_distance, promoter_distance), 
                 genomicAnnotationPriority =c("Promoter", "5UTR", "3UTR", "Exon","Intron","Intergenic"),
                 TxDb=txdb, annoDb=annodb)
annotation<-data.frame("peak"=paste(x@anno@seqnames,x@anno@ranges@start-1,sep="_"),
                       "ann"=x@anno$annotation,
                       "annotation"=x@anno$SYMBOL,stringsAsFactors = F)
p2p_links<-mutate(p2p_links,"annotation1"=annotation$ann[match(as.character(p2p_links$Peak1),annotation$peak)],
             "annotation2"=annotation$ann[match(as.character(p2p_links$Peak2),annotation$peak)],
             "gene1"=annotation$annotation[match(as.character(p2p_links$Peak1),annotation$peak)],
             "gene2"=annotation$annotation[match(as.character(p2p_links$Peak2),annotation$peak)])

p2p_links<-p2p_links[str_detect(p2p_links$annotation1,"Promoter")|str_detect(p2p_links$annotation2,"Promoter"),]
p2p_links<-p2p_links[!(str_detect(p2p_links$annotation1,"Promoter")&str_detect(p2p_links$annotation2,"Promoter")),]
p2p_links_final<-data.frame("Peak"=c(as.character(p2p_links$Peak1)[grep("Promoter",p2p_links$annotation2)],as.character(p2p_links$Peak2)[grep("Promoter",p2p_links$annotation1)]),
                     "Gene"=c(p2p_links$gene2[grep("Promoter",p2p_links$annotation2)],p2p_links$gene1[grep("Promoter",p2p_links$annotation1)]),
                     "value"=p2p_links$value[c(grep("Promoter",p2p_links$annotation2),grep("Promoter",p2p_links$annotation1))],stringsAsFactors = F)

value<-tapply(p2p_links_final$value,INDEX=paste(p2p_links_final$Peak,"+",p2p_links_final$Gene,sep=""),FUN=function(x){
  x[which.max(abs(x))]
})
p2p_links_final <- data.frame("Peak"=sub("^(.*)\\+.*$","\\1",names(value)),"Gene"=sub("^.*\\+(.*)$","\\1",names(value)),
                              "value"=value,"type"="peak2peak",stringsAsFactors = F)
p2p_links_final<-p2p_links_final[order(p2p_links_final$Gene),]

### Filter peak-to-gene links
p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = -1,
  FDRCutOff = 0.01,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1,
  returnLoops = FALSE
)

peakSet <- metadata(p2g)$peakSet
peakSet$mean_bp <- (as.numeric(as.character(peakSet@ranges@start)) + as.numeric(as.character(peakSet@ranges@start+peakSet@ranges@width-1)))/2
p2g_links <- data.frame(
  Peak = paste(peakSet@seqnames,peakSet$mean_bp,sep="_")[p2g$idxATAC],
  Gene = metadata(p2g)$geneSet$name[p2g$idxRNA],
  value = p2g$Correlation,
  type = "peak2gene"
)
p2g_links_final <- p2g_links

### Calculate total links
linksMat<-rbind(p2p_links_final,p2g_links_final)
colnames(linksMat)<-c("Peak","Gene","value","type")
linksMat[,3]<-as.numeric(linksMat[,3])
dup<-intersect(paste(linksMat$Peak,linksMat$Gene,sep="+")[linksMat$type=="peak2peak"],
               paste(linksMat$Peak,linksMat$Gene,sep="+")[linksMat$type=="peak2gene"])
linksMat_dup<-rbind(p2p_links_final[match(dup,paste(p2p_links_final$Peak,p2p_links_final$Gene,sep="+")),],
                    p2g_links_final[match(dup,paste(p2g_links_final$Peak,p2g_links_final$Gene,sep="+")),])
link_pos<-plot_df[plot_df$p2p>=0&plot_df$p2g>=0,]
link_pos$mean_score<-rowMeans(link_pos)
link_neg<-plot_df[plot_df$p2p<0&plot_df$p2g<0,]
link_neg$mean_score<-rowMeans(link_neg)
links_all<-rbind(link_pos,link_neg)

### Calculate binding score
width <- 200
data("human_pwms_v1") #filtered collection of human motifs from cisBP database
peakSet<-GRanges(sub("_",":",sub("^(.*)\\+.*$","\\1",rownames(links_all))))
peakSet<-resize(peakSet,width,"center")
peak<-paste(peakSet@seqnames,1+as.numeric(resize(peakSet@ranges,1,"center")@start),sep="_")

motifs <- human_pwms_v1
motif_ix <- matchMotifs(motifs, peakSet, out =  c("matches"),
                        genome = BSgenome.Hsapiens.UCSC.hg38)
motif_score <- matchMotifs(motifs, peakSet, out =  c("scores"),
                           genome = BSgenome.Hsapiens.UCSC.hg38)
TF_motif_score <- as.matrix(motif_ix@assays@data$motifMatches)
colnames(TF_motif_score) <- paste(motif_ix@colData$name,1:length(motif_ix@colData$name),sep="_")
rownames(TF_motif_score) <- peak
binding_matrix<-matrix(data=0,nrow = length(peak),ncol = length(motifs))
rownames(binding_matrix)<-peak
colnames(binding_matrix)<-colnames(TF_motif_score)
binding_matrix<-TF_motif_score[match(peak,rownames(TF_motif_score)),]

TF_motif_bindscore <- as.matrix(motif_score@assays@data$motifScores)
colnames(TF_motif_bindscore) <- paste(motif_score@colData$name,1:length(motif_score@colData$name),sep="_")
rownames(TF_motif_bindscore) <- peak
binding_score_matrix<-matrix(data=0,nrow = length(peak),ncol = length(motifs))
rownames(binding_score_matrix)<-peak
colnames(binding_score_matrix)<-colnames(TF_motif_bindscore)
binding_score_matrix<-TF_motif_bindscore[match(peak,rownames(TF_motif_bindscore)),]

### Enriching regulated TF for each gene
peak<-sapply(str_split(rownames(links_all),pattern = "\\+"),function(x){x[1]})
gene<-sapply(str_split(rownames(links_all),pattern = "\\+"),function(x){x[2]})
gene_peak_dt<-data.table("gene"=as.character(gene),"peak"=as.character(peak),"value"=as.numeric(links_all$mean_score))

width <- 200
genome <- BSgenome.Hsapiens.UCSC.hg38
peakSet <- GRanges(sub("_",":",sub("^(.*)\\+.*$","\\1",rownames(links_all))))
peakSet <- resize(peakSet,width,"center")
peak <- paste(peakSet@seqnames,1+as.numeric(resize(peakSet@ranges,1,"center")@start),sep="_")
names(peakSet) <- peak

number <- sapply(unique(gene_peak_dt$gene),function(x){nrow(gene_peak_dt[gene_peak_dt$gene==x,])})

computeEnrichment <- function(matches = NULL, compare = NULL, background = NULL){
  matchCompare <- matches[compare, ,drop=FALSE]
  matchBackground <- matches[background, ,drop=FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)

  pOut <- data.frame(
    feature = colnames(matches),
    CompareFrequency = matchCompareTotal,
    nCompare = nrow(matchCompare),
    CompareProportion = matchCompareTotal/nrow(matchCompare),
    BackgroundFrequency = matchBackgroundTotal,
    nBackground = nrow(matchBackground),
    BackgroundProporition = matchBackgroundTotal/nrow(matchBackground)
  )

  pOut$Enrichment <- pOut$CompareProportion / pOut$BackgroundProporition
  
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x){
    p <- -phyper(pOut$CompareFrequency[x], # Number of Successes
                 pOut$BackgroundFrequency[x], # Number of all successes in background
                 pOut$nBackground[x] - pOut$BackgroundFrequency[x], # Number of non successes in background
                 pOut$nCompare[x], # Number that were drawn
                 lower.tail = FALSE, log.p = TRUE)# P[X > x] Returns LN must convert to log10
    return(p/log(10))
  }) %>% unlist %>% round(4)
  
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]
  
  pOut
}

EnrichByGene <- function(genome, peakSet, neighbors = 500, matches, scores, gene_peak_dt){
  peakSplit <- split(peakSet, as.character(seqnames(peakSet)))
  chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
  
  message("Compute GC score")
  peakNuc <- lapply(seq_along(peakSplit), function(x){
    message(sprintf("%s of %s", x, length(peakSplit)))
    chrSeq <- Biostrings::getSeq(genome,chromSizes[which(as.character(seqnames(chromSizes))==names(peakSplit)[x])])
    grx <- peakSplit[[x]]
    aFreq <- alphabetFrequency(Biostrings::Views(chrSeq[[1]], ranges(grx)))
    mcols(grx)$GC <- rowSums(aFreq[, c("G","C")]) / rowSums(aFreq)
    mcols(grx)$AT <- rowSums(aFreq[, c("A","T")]) / rowSums(aFreq)
    return(grx)
  }) %>% GenomicRangesList %>% unlist %>% sortSeqlevels
  peakNuc$N <- 1 - (peakNuc$GC + peakNuc$AT)

  message("Compute Enrichment")
  EnrichList <- lapply(seq_along(unique(gene_peak_dt$gene)), function(x){
    if(x %% 100 == 0){
      message(sprintf("%s of %s", x, length(unique(gene_peak_dt$gene))))
    }
    
    df <- gene_peak_dt[gene_peak_dt$gene == unique(gene_peak_dt$gene)[x],]
    peaks <- df$peak

    BgdPeaks <- sapply(peaks, function(i) {
      diff <- abs(peakNuc$GC - peakNuc[i]$GC)
      if(sum(diff==0) > neighbors){
        idxNN <- sample(which(diff==0), neighbors + 1)
      } else idxNN <- head(order(diff), neighbors + 1)
      idxNN <- idxNN[!(idxNN %in% i)]
      idxNN[1:neighbors]
    })
    background <- names(peakNuc)[unique(as.vector(BgdPeaks))]
    
    out <- computeEnrichment(matches, peaks, background)
    
    if(length(peaks)>1 & nrow(out)>1){
      out$score <- colMeans(abs(scores[peaks,as.character(out$feature)])*df$value[df$peak%in%peaks])
    } else if(length(peaks)==1) {
      out$score <- abs(scores[peaks,as.character(out$feature)])*df$value[df$peak%in%peaks]
    } else if(nrow(out)==1) {
      out$score <- mean(abs(scores[peaks,as.character(out$feature)])*df$value[df$peak%in%peaks])
    }
    
    out <- out[out$score!=0,]
    out[order(abs(as.numeric(out$score)),decreasing = T),]
    
  }) %>% SimpleList
  
  names(EnrichList) <- unique(gene_peak_dt$gene)
  EnrichList
}

### Calculate regulons used for later identification
regulons <- EnrichByGene(genome = genome, peakSet = peakSet, neighbors = 500, 
                       matches = binding_matrix, scores = binding_score_matrix, gene_peak_dt = gene_peak_dt)

pCutoff <- 0.05
regulons <- lapply(regulons,function(x){
  x[x$mlog10p>(-log10(pCutoff)),]
})
