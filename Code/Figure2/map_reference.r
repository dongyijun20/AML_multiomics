####step2: map reference----------------------------------
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
load("Peak.matrix.seurat_final3.rdata")
load("HYX_qc.rdata")
peak_range <- GRanges(seqnames = sub("(^.*)?-.*-.*$","\\1",rownames(Peak.matrix.seurat)), 
                      ranges = IRanges(start = as.numeric(sub("^.*\\-(.*)\\-.*$","\\1",rownames(Peak.matrix.seurat))),
                                       end = as.numeric(sub("^.*\\-(.*)$","\\1",rownames(Peak.matrix.seurat)))))
useproj<-HYX_qc
counts_useproj <- FeatureMatrix(
  fragments = useproj@assays$peaks@fragments,
  features = peak_range,
  cells = colnames(useproj)
)
useproj_assay <- CreateChromatinAssay(
  counts = counts_useproj,
  fragments = useproj@assays$peaks@fragments
)
useproj_map <- CreateSeuratObject(counts = useproj_assay, assay = "peaks")
useproj_map <- FindTopFeatures(useproj_map, min.cutoff = 10)
useproj_map <- RunTFIDF(useproj_map)
useproj_map <- RunSVD(useproj_map)
useproj_map <- RunUMAP(useproj_map,reduction = "lsi", dims = 1:30)
transfer.anchors <- FindTransferAnchors(
  reference = Peak.matrix.seurat,
  query = useproj_map,
  query.assay="peaks",
  reference.assay ="peaks",
  reference.reduction = "lsi",
  features = rownames(useproj_map),
  reduction = "lsiproject",
  dims = 1:30
)
map.atac <- MapQuery(
  anchorset = transfer.anchors,
  reference = Peak.matrix.seurat,
  query = useproj_map,
  refdata = Peak.matrix.seurat$nbt_cluster,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)
### add celltype label------------------
map.atac$predicted.celltype<-map.atac$predicted.id
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster1"]<-"1-HSC/MPP"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster2"]<-"2-MEP"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster3"]<-"3-CMP/BMP"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster4"]<-"4-LMPP"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster5"]<-"5-CLP"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster6"]<-"6-ProB"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster7"]<-"7-PreB"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster8"]<-"8-GMP"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster9"]<-"9-MDP"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster10"]<-"10-pDC"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster11"]<-"11-cDC"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster12"]<-"12-Mono1"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster13"]<-"13-Mono2"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster14"]<-"14-NaiveB"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster15"]<-"15-MemoryB"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster16"]<-"16-Plasma"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster17"]<-"17-Basophil"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster18"]<-"18-immature_NK"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster19"]<-"19-Mature_NK"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster20"]<-"20-immature_NK"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster21"]<-"21-Naive_CD4_1"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster22"]<-"22-Naive_CD4_2"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster23"]<-"23-naive_treg"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster24"]<-"24-memory_CD4"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster25"]<-"25-treg"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster26"]<-"26-naive_CD8_1"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster27"]<-"27-naive_CD8_2"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster28"]<-"28-naive_CD8_3"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster29"]<-"29-Central_MemoryCD8"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster30"]<-"30-Effector_MemoryCD8"
map.atac$predicted.celltype[map.atac$predicted.id=="Cluster31"]<-"31-Gamma_deltaT"
map.atac$celltype<-map.atac$predicted.id
map.atac$celltype[map.atac$predicted.id=="Cluster6"|
                    map.atac$predicted.id=="Cluster7"]<-"proB"
map.atac$celltype[map.atac$predicted.id=="Cluster14"|
                    map.atac$predicted.id=="Cluster15"|
                    map.atac$predicted.id=="Cluster16"]<-"B cell"

map.atac$celltype[map.atac$predicted.id=="Cluster18"|
                    map.atac$predicted.id=="Cluster19"|
                    map.atac$predicted.id=="Cluster20"]<-"NK"

map.atac$celltype[map.atac$predicted.id=="Cluster21"|
                    map.atac$predicted.id=="Cluster22"|
                    map.atac$predicted.id=="Cluster23"|
                    map.atac$predicted.id=="Cluster24"|
                    map.atac$predicted.id=="Cluster25"]<-"CD4"

map.atac$celltype[map.atac$predicted.id=="Cluster26"|
                    map.atac$predicted.id=="Cluster27"|
                    map.atac$predicted.id=="Cluster28"|
                    map.atac$predicted.id=="Cluster29"|
                    map.atac$predicted.id=="Cluster30"|
                    map.atac$predicted.id=="Cluster31"]<-"CD8"

map.atac$celltype[map.atac$predicted.id=="Cluster1"]<-"HSC/MPP"
map.atac$celltype[map.atac$predicted.id=="Cluster2"]<-"MEP"
map.atac$celltype[map.atac$predicted.id=="Cluster3"]<-"CMP/BMP"
map.atac$celltype[map.atac$predicted.id=="Cluster4"]<-"LMPP"
map.atac$celltype[map.atac$predicted.id=="Cluster5"]<-"CLP"
map.atac$celltype[map.atac$predicted.id=="Cluster8"]<-"GMP"
map.atac$celltype[map.atac$predicted.id=="Cluster9"]<-"MDP"
map.atac$celltype[map.atac$predicted.id=="Cluster10"]<-"pDC"
map.atac$celltype[map.atac$predicted.id=="Cluster11"]<-"cDC"
map.atac$celltype[map.atac$predicted.id=="Cluster12"]<-"Mono1"
map.atac$celltype[map.atac$predicted.id=="Cluster13"]<-"Mono2"
map.atac$celltype[map.atac$predicted.id=="Cluster17"]<-"Basophil"


#DimPlot(Peak.matrix.seurat, reduction = "umap", group.by = "nbt_cluster",cols = as.character(color_bioclass$color), label = TRUE, repel = TRUE) + 
#   ggtitle("Reference")
#DimPlot(map.atac, reduction = "ref.umap", group.by = "predicted.id", cols = as.character(color_bioclass$color),label = TRUE, repel = TRUE) + 
#   ggtitle("Query")
#DimPlot(map.atac, reduction = "umap", group.by = "predicted.id", cols = as.character(color_bioclass$color),label = TRUE, repel = TRUE) + 
#  ggtitle("Query")
save(map.atac,file="map.atac.rdata")
save(transfer.anchors,file="transfer.anchors.rdata")
save(useproj_map,file="useproj_map.rdata")
rm(list=ls())
