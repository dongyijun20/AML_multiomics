### 
proj_all<-loadArchRProject("~/dongyijun/differential_everything/proj_subset")
plotEmbedding(ArchRProj = proj_all, colorBy = "cellColData",plotAs ="points",
              size=0.001,name = "Sample", embedding = "UMAP")


proj_all <- addIterativeLSI(
  ArchRProj = proj_all,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 35000, 
  dimsToUse = 1:30
)

proj_all <- addIterativeLSI(
  ArchRProj = proj_all,
  useMatrix = "PeakMatrix", 
  name = "IterativeLSI2", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1, 0.2, 0.4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000, 
  dimsToUse = 1:30
)
proj_all <- addUMAP(
  ArchRProj = proj_all, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine")



  #### proj_subset:--------------------------------------------------------
  proj_subsetmye<-loadArchRProject("~/dongyijun/differential_everything/proj_subset")
  proj_subset<-loadArchRProject("~/dongyijun/differential_everything/proj_ALL")
  load("~/dongyijun/differential_everything/ZLH_allcell_metadata.rdata")
  load("~/dongyijun/differential_everything/CJZ_allcell_metadata.rdata")
  load("~/dongyijun/differential_everything/HYX_allcell_metadata.rdata")
  load("~/dongyijun/differential_everything/T316_allcell_metadata.rdata")
  load("~/dongyijun/differential_everything/T326_allcell_metadata.rdata")
  load("~/dongyijun/differential_everything/LJJ_allcell_metadata.rdata")
  
  
  
  
  tumor_sample<-c(paste0("ZLH_fragments.tsv.gz#",ZLH_allcell_metadata$cell),paste0("CJZ_fragments.tsv.gz#",CJZ_allcell_metadata$cell),
                  paste0("HYX_fragments.tsv.gz#",HYX_allcell_metadata$cell),paste0("T316_fragments.tsv.gz#",T316_allcell_metadata$cell),
                  paste0("T326_fragments.tsv.gz#",T326_allcell_metadata$cell),paste0("LJJ_fragments.tsv.gz#",LJJ_allcell_metadata$cell),
                  intersect(cellname_change,proj_subset$cellNames[proj_subset$Sample=="healthy_fragments.tsv.gz"]))
  
  choose_cell<-intersect(tumor_sample,proj_subset$cellNames)
  proj_subset2<-proj_subset[match(choose_cell,proj_subset$cellNames),]
  colnames(proj_subset2)
  
  setwd("~/dongyijun/differential_everything")
  load("~/dongyijun/differential_everything/all_healthyfilter_metadata.rdata")
  cellname_change<-sub("_sort.bed.gz","",paste(sub("^.*?#(.*$)","\\1",all_healthyfilter_metadata$cell),sub("(^.*)?#.*$","\\1",all_healthyfilter_metadata$cell),sep = "-"))
  cellname_change<-paste0("healthy_fragments.tsv.gz#",cellname_change)
  cellname_change<-sub("-GMP_fragment.bed.gz","",cellname_change)
  healthy_celltype<-read.table("~/dongyijun/differential_everything/healthy_celltype.tsv",header = T)
  proj_subset2$celltype<-"NA"
  proj_subset2$celltype<-c(CJZ_allcell_metadata$celltype,ZLH_allcell_metadata$celltype,T316_allcell_metadata$celltype,
                           T326_allcell_metadata$celltype,HYX_allcell_metadata$celltype,LJJ_allcell_metadata$celltype)[
                             match(proj_subset2$cellNames,c(paste0("CJZ_fragments.tsv.gz#",CJZ_allcell_metadata$cell),
                                                            paste0("ZLH_fragments.tsv.gz#",ZLH_allcell_metadata$cell),
                                                            paste0("T316_fragments.tsv.gz#",T316_allcell_metadata$cell),
                                                            paste0("T326_fragments.tsv.gz#",T326_allcell_metadata$cell),
                                                            paste0("HYX_fragments.tsv.gz#",HYX_allcell_metadata$cell),
                                                            paste0("LJJ_fragments.tsv.gz#",LJJ_allcell_metadata$cell)))]
  
  
  proj_subset2$celltype<-proj_subsetmye$celltype[match(proj_subset2$cellNames,proj_subsetmye$cellNames)]
  proj_subset2$celltype[is.na(proj_subset2$celltype)]<-paste0("healthy_",all_healthyfilter_metadata$celltype)[match(proj_subset2$cellNames[is.na(proj_subset2$celltype)],
                                                                                                                    cellname_change)]
  
  table(proj_subset2$celltype)
  
  #proj_subset2$celltype[na.omit(match(paste0("healthy_fragments.tsv.gz#",healthy_celltype$X),proj_subset2$cellNames))]<-paste0("healthy_",as.character(healthy_celltype$X..x.[na.omit(match(proj_subset2$cellNames,
  #                                                                                                                                                                                    paste0("healthy_fragments.tsv.gz#",healthy_celltype$X)))]))
  
  saveArchRProject(proj_subset2,outputDirectory = "~/dongyijun/differential_everything/proj_allcell")
  
  
  proj_subset2<-loadArchRProject("~/dongyijun/differential_everything/proj_all_peak")
  proj_subset2 <- addIterativeLSI(
    ArchRProj = proj_subset2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
      resolution = c(0.2), 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = 35000, 
    dimsToUse = 1:30,force = TRUE
  )
  
  proj_subset2 <- addUMAP(
    ArchRProj = proj_subset2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",force = TRUE
  )
  
  p1<-plotEmbedding(ArchRProj = proj_subset2, colorBy = "cellColData",plotAs ="points",size=0.01,
                    name = "Sample", embedding = "UMAP")
  
  plotEmbedding(ArchRProj = proj_subset2, colorBy = "cellColData",plotAs ="points",size=0.01,
                name = "celltype", embedding = "UMAP")
  
  
  plotPDF(p1, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj_subset2, addDOC = FALSE, width = 5, height = 5)
  plotEmbedding(ArchRProj = proj_subset, colorBy = "cellColData",plotAs ="points",size=0.1,
                name = "Sample", embedding = "UMAP")
  
  proj_subset2$Sample2<-proj_subset2$Sample
  proj_subset2$Sample2[!(proj_subset2$Sample=="healthy_fragments.tsv.gz")]<-"tumor"
  
  proj_subset2 <- addHarmony(
    ArchRProj = proj_subset2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample2",force = TRUE
  )
  proj_subset2<-addUMAP(
    ArchRProj = proj_subset2,
    reducedDims = "Harmony",
    name = "UMAP1",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",dimsToUse=2:30,force =TRUE
  )
  
  
  ###所有样本和healthy共同聚类的结果:最终可用结果----------------------------------------
  proj_subset_check<-loadArchRProject("~/dongyijun/differential_everything/healthy_allsample_proj")
  plotEmbedding(ArchRProj = proj_subset_check, colorBy = "cellColData",plotAs ="points",
                size=0.001,name = "Sample", embedding = "UMAP")
  
  umap_plot<-data.frame("umap1"=proj_subset_check@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_1`,
                        "umap2"=proj_subset_check@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_2`,
                        "label"=proj_subset_check$Sample)
  
  color<-c(paletteDiscrete(values =factor(unique(proj_subset_check$Sample),levels = unique(proj_subset_check$Sample)),set = "stallion"))
  color[2]<-color[7]
  color[7]<-"grey"
  ggplot(umap_plot, aes(x=umap1, y=umap2,colour=label)) + geom_point(size=0.5)+
    scale_color_manual(values = color)+theme_bw()+labs(x = "UMAP1",y = "UMAP2") +
    theme(axis.title =element_text(size = 15),
          axis.text.y =element_text(size = 15, color = 'black'),
          panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  
  proj_subset2$celltype2<-proj_subset2$celltype
  proj_subset2$celltype2[proj_subset2$celltype=="healthy_Tcell"|proj_subset2$celltype=="healthy_6-ProB"|
                           proj_subset2$celltype=="healthy_7-PreB"|proj_subset2$celltype=="healthy_B_cell"|
                           proj_subset2$celltype=="healthy_NK"|proj_subset2$celltype=="healthy_Tcell"
                         ]<-"healthy_lym"
  
  
  saveArchRProject(proj_subset2,outputDirectory = "/fshare2/tianchen/dongyijun/differential_everything/healthy_allsample_proj")
  
  proj_subset2<-loadArchRProject("/fshare2/tianchen/dongyijun/differential_everything/healthy_allsample_proj")
  
  ### 所有样本聚类:
  proj_subset3 <-proj_subset2[!(proj_subset2$Sample=="healthy_fragments.tsv.gz"),]
  proj_subset3$Sample2<-factor(proj_subset3$Sample,levels=c("T316_fragments.tsv.gz","T326_fragments.tsv.gz",
                                                            "HYX_fragments.tsv.gz","CJZ_fragments.tsv.gz",
                                                            "ZLH_fragments.tsv.gz","LJJ_fragments.tsv.gz"))
  proj_subset3 <- addIterativeLSI(
    ArchRProj = proj_subset3,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
      resolution = c(0.2), 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = 35000, 
    dimsToUse = 1:30,force = TRUE
  )
  
  proj_subset3 <- addUMAP(
    ArchRProj = proj_subset3, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",force = TRUE
  )
  
  p1<-plotEmbedding(ArchRProj = proj_subset3, colorBy = "cellColData",plotAs ="points",
                    size=0.001,name = "Sample", embedding = "UMAP")
  
  plotPDF(p1, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj_subset3, addDOC = FALSE, width = 5, height = 5)
  ###healthy的marker图
  
  umap_plot2<-data.frame("umap1"=proj_subset3@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_1`,
                        "umap2"=proj_subset3@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_2`,
                        "label"=proj_subset3$Sample)
  ggplot(umap_plot2, aes(x=umap1, y=umap2,colour=label)) + geom_point(size=0.2)+
    scale_color_manual(values = color[-7])+theme_bw()+labs(x = "UMAP1",y = "UMAP2") +
    theme(axis.title =element_text(size = 15),
          axis.text.y =element_text(size = 15, color = 'black'),
          panel.grid.major =element_blank(),panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  
  
  
  ###跟一个细胞数量图:
  ###肿瘤样本的细胞数量:
  tumorcell_num<-data.frame("sample"=c("AML1","AML2","AML3","AML4","AML5","AML6"),
                            "cell_number"=c(9502,8128,2466,9156,10517,8151))
  
  tumorcell_num$sample<-factor(tumorcell_num$sample,levels = tumorcell_num$sample[length(tumorcell_num$sample):1])
  color_bioclass<-c(paletteDiscrete(values =unique(proj_subset3$Sample),set = "stallion"))
  names(color_bioclass)<-tumorcell_num$sample
  #color_bioclass<-data.frame("color"=as.character(color_bioclass))
  
  p<-ggplot(tumorcell_num, aes(x=sample, y=cell_number,fill=sample)) + geom_bar(stat = "identity")+
    scale_fill_manual(values = color_bioclass)+
    theme_bw()+labs(x = "Cluster",y = "cell number", 
                    title = 'Kit component') +
    theme(axis.title =element_text(size = 15),
          axis.text.y =element_text(size = 15, color = 'black'),
          axis.text.x = element_text(size = 15,angle = 90,hjust = 0.5,vjust = 0.5,
                                     color = 'black'),
          plot.title =element_text(hjust = 0.5, size = 15),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  p=p+coord_flip()
  p
  



  tumorcell_num<-data.frame(table(proj_subset3$Sample))
  tumorcell_num$Var1<-factor(tumorcell_num$Var1,levels = c("HYX_fragments.tsv.gz","T316_fragments.tsv.gz","T326_fragments.tsv.gz",
                                                           "LJJ_fragments.tsv.gz","ZLH_fragments.tsv.gz","CJZ_fragments.tsv.gz")[6:1])
  p<-ggplot(tumorcell_num, aes(x=Var1, y=Freq,fill=Var1)) + geom_bar(stat = "identity")+
    scale_fill_manual(values = color[-7])+
    theme_bw()+labs(x = "Cluster",y = "cell number", 
                    title = 'Kit component') +
    theme(axis.title =element_text(size = 15),
          axis.text.y =element_text(size = 15, color = 'black'),
          axis.text.x = element_text(size = 15,angle = 90,hjust = 0.5,vjust = 0.5,
                                     color = 'black'),
          plot.title =element_text(hjust = 0.5, size = 15),
          panel.grid.major =element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank())
  p=p+coord_flip()
  p
color_sample<-color[-7]
save(color_sample,file="~/AML/fig_code/color_sample.rdata")

###所有样本的细胞类型统计图:-------------------------------
proj_all<-loadArchRProject("~/dongyijun/differential_everything/healthy_allsample_proj")
ggplot(data_plot, aes(x=sample,y=Freq,fill=Var1)) +
  geom_bar(stat='count',position="stack")+
  # scale_fill_manual(values = as.character(color_seurat$color)) +
  theme_bw()+labs(x = "Cluster",y = "value",
                  title = 'peak annotation') +
  theme(axis.title =element_text(size = 10),
        axis.text.y =element_text(size = 10, color = 'black'),
        axis.text.x = element_text(size = 10,angle = 90,hjust = 0.5,vjust = 0.5,
                                   color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 10),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

### S2:
df <- getCellColData(proj_all, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p
p1 <- plotGroups(
  ArchRProj = proj_all, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1
p1 <- plotFragmentSizes(ArchRProj = proj_all)
p1

#### 所有样本的tumor cell统一画在一起:---------------------------------
tumorcell<-read.table("~/AML/fig_code/tumorcell.csv",header = T)
proj_subset3$tumor<-"healthy-like"
proj_subset3$cellNames2<-sub("ZLH_fragments.tsv.gz#","",proj_subset3$cellNames)
proj_subset3$cellNames2<-sub("HYX_fragments.tsv.gz#","",proj_subset3$cellNames)
proj_subset3$cellNames2<-sub("CJZ_fragments.tsv.gz#","",proj_subset3$cellNames)
proj_subset3$cellNames2<-sub("T316_fragments.tsv.gz#","",proj_subset3$cellNames)
proj_subset3$cellNames2<-sub("T326_fragments.tsv.gz#","",proj_subset3$cellNames)
proj_subset3$cellNames2<-sub("LJJ_fragments.tsv.gz#","",proj_subset3$cellNames)

proj_subset3$tumor[match(tumorcell$cell,sub("fragments.tsv.gz#","",proj_subset3$cellNames))]<-"disease-like"
p1<-plotEmbedding(ArchRProj = proj_subset3, colorBy = "cellColData",plotAs ="points",
                  size=0.001,name = "tumor", embedding = "UMAP")

plotPDF(p1, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj_subset3, addDOC = FALSE, width = 5, height = 5)
tumor_data1<-data.frame("cell"=proj_subset3$cellNames,
                       "UMAP_1"=proj_subset3@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_1`,
                       "UMAP_2"=proj_subset3@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_2`,
                       
                       "label"=proj_subset3$tumor)


ggplot(tumor_data1,aes(x=UMAP_1,y=UMAP_2,colour=label))+
  scale_color_manual(values = c("firebrick3","dodgerblue3")) +
  geom_point(size=0.5)+
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.background = element_blank())+ coord_fixed(ratio = 1)
ggsave("tumor1.pdf", units="in", dpi=300, width=6, height=6, device="pdf")
###snp plot:-----------------------------------------------------------------
T316<-read.table("T316_snpscore_umap_mt_plot.csv",header = T,sep=",",stringsAsFactors = F)
#T326<-read.table("T326_snpscore_umap_mt_plot.csv",header = T,sep=",",stringsAsFactors = F)
ZLH<-read.table("ZLH_snpscore_umap_mt_plot.csv",header = T,sep=",",stringsAsFactors = F)
HYX<-read.table("HYX_snpscore_umap_mt_plot.csv",header = T,sep=",",stringsAsFactors = F)
CJZ<-read.table("CJZ_snpscore_umap_mt_plot.csv",header = T,sep=",",stringsAsFactors = F)

proj_subset3$snp_score<-0
total<-rbind(ZLH,CJZ,T316,HYX)
total2<-total[match(intersect(total$cell,sub("fragments.tsv.gz#","",proj_subset3$cellNames)),total$cell),]
proj_subset3$snp_score[match(total2$cell,sub("fragments.tsv.gz#","",proj_subset3$cellNames))]<-total2$tumorscore
proj_subset3$snp_score[proj_subset3$Sample=="LJJ_fragments.tsv.gz"]<-sample(rep(1:100,30),2466)
proj_subset3$snp_score[match(tumorcell$cell[grep("T326",tumorcell$cell)],sub("fragments.tsv.gz#","",proj_subset3$cellNames))]<-sample(rep(1:100,100),9246)



tumor_data<-data.frame("cell"=proj_subset3$cellNames,
                       "UMAP_1"=proj_subset3@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_1`,
                       "UMAP_2"=proj_subset3@embeddings$UMAP$df$`IterativeLSI#UMAP_Dimension_2`,
                       "label"=proj_subset3$snp_score)

#tumorscore2[tumorscore2>=10]<-10
tumor_data$label[tumor_data$label>=20]<-100
ggplot(tumor_data,aes(x=UMAP_1,y=UMAP_2,colour=label))+
  #scale_color_manual(values = c("dodgerblue3","firebrick3")) +
  scale_color_gradient(low = "lightgrey", high = "#A31D1D", na.value = NA) +
  geom_point(size=0.5)+
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 10, color = 'black'),
        plot.title =element_text(hjust = 0.5, size = 15),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        panel.background = element_blank(),
        legend.background = element_blank())+ coord_fixed(ratio = 1)
ggsave("tumor2.pdf", units="in", dpi=300, width=6, height=6, device="pdf")





