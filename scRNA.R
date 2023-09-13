################# scRNA  #################
### 0.准备
### 1.导入
### 2.质控
### 3.整合
### 4.降维
### 5.聚类
### 6.标志物
### #灵活操作
### 7.富集
### 8.monocle
### 9.PPI
### 10.cellchat
### 11.RNAvelocity
### 12.scenic

#################  0.准备  ##################
if(F){
options()$repos
options()$BioC_mirror
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options("repos" = c(CseuratN="https://mirrors.tuna.tsinghua.edu.cn/CseuratN/")) 
options()$repos 
options()$BioC_mirror

rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(hdf5r)
library(clustree)
library(SingleR)
library(sctransform)
library(tidyverse) 
library(RCurl)
library(BiocManager)
library(SeuratDisk)
library(pryr)
library(ulimit)
library(metap)
}
install.packages('pryr')
devtools::install_github("krlmlr/ulimit")
#################  1.导入  ##################
if(F){
# 批量导入
  folders=list.files('data/as3/gz')
  folders
  sceList = lapply(folders,function(folder){ 
    CreateSeuratObject(counts = Read10X(paste0('data/ad2dup/',folder)),
                       project = folder)
  })
# merge
  ad2_dup <- merge(sceList[[1]], y = c(sceList[[2]],sceList[[3]],sceList[[4]],sceList[[5]],
                                       sceList[[6]],sceList[[7]],sceList[[8]],sceList[[9]],sceList[[10]]), 
                          add.cell.ids = folders, 
                          project = "ad2dup")
# 观察元数据和表达矩阵
  head(seurat_merge@meta.data)
  table(seurat_merge@meta.data$orig.ident)
  seurat_merge@assays$RNA[1:3,1:30]
  dim(seurat_merge@assays$RNA@counts)
  dim(seurat_merge@assays$RNA@data)
  dim(seurat_merge@assays$RNA@scale.data)
  seurat_merge@assays$RNA@counts[1:20,1:20]
  seurat_merge@assays$RNA@data[1:20,1:20]
  rownames(seurat_merge)
}
if(F){
  as2<- CreateSeuratObject(raw_counts, project = 'as2',assay = 'RNA',meta.data = meta)
}

#################  2.质控  ##################
### 2.1 整理元数据
if(F){
# 将每个细胞每个UMI的基因数目添加到元数据中
  seurat_merge$log10GenesPerUMI <- log10(seurat_merge$nFeature_RNA) / log10(seurat_merge$nCount_RNA)
  head(seurat_merge@meta.data)
# 计算线粒体比率 
  seurat_merge$mitoRatio <- PercentageFeatureSet(object = seurat_merge, pattern = "^MT-") 
  seurat_merge$mitoRatio <- seurat_merge@meta.data$mitoRatio / 100  
# 计算红细胞比率
  HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")     
  HB_m <- match(HB.genes_total,rownames(seurat_merge@assays$RNA))     
  HB.genes <- rownames(seurat_merge@assays$RNA)[HB_m]     
  HB.genes <- HB.genes[!is.na(HB.genes)]     
  seurat_merge[["percent.HB"]]<-PercentageFeatureSet(seurat_merge, assay = 'RNA', features=HB.genes)
  seurat_merge$percent.HB <- seurat_merge@meta.data$percent.HB / 100
# 创建元数据数据框
  metadata <- seurat_merge@meta.data
# 为元数据添加细胞ID 
  metadata$cells <- rownames(metadata)
# 重命名列 
  # metadata <- metadata %>% 
  #   dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
# 创建样本列 
  metadata$sample <- NA
  metadata$sample[which(str_detect(metadata$cells, "^HC1_"))] <- "HC" 
  metadata$sample[which(str_detect(metadata$cells, "^HC2_"))] <- "HC"
  metadata$sample[which(str_detect(metadata$cells, "^L1_"))] <- "L"
  metadata$sample[which(str_detect(metadata$cells, "^L2_"))] <- "L" 
  metadata$sample[which(str_detect(metadata$cells, "^L3_"))] <- "L" 
  metadata$sample[which(str_detect(metadata$cells, "^N1_"))] <- "N"
  metadata$sample[which(str_detect(metadata$cells, "^N2_"))] <- "N"
  metadata$sample[which(str_detect(metadata$cells, "^N3_"))] <- "N" 

  metadata$group <- paste0(metadata$tissue,'',metadata$sample)
  table(metadata$group)
  metadata$group <- NA
  metadata$group[which(str_detect(metadata$group, "branchHC"))] <- "ASHC" 
  metadata$group[which(str_detect(metadata$group, "branchL"))] <- "ASL"
  metadata$group[which(str_detect(metadata$group, "nasalHC"))] <- "CRSHC"
  metadata$group[which(str_detect(metadata$group, "nasalL"))] <- "CRSL" 
  metadata$group[which(str_detect(metadata$group, "skinHC"))] <- "ADHC" 
  metadata$group[which(str_detect(metadata$group, "skinL"))] <- "ADL"

# 将元数据添加回Seurat对象中 
  seurat_merge@meta.data <- metadata
  head(seurat_merge@meta.data)
  tail(seurat_merge@meta.data)
  
# 另一种直接添加元数据的方法
  # if(F){
  # Idents(seurat_merge) <- 'orig.ident'
  # levels(Idents(seurat_merge))
  # clinic <- Idents(seurat_merge)
  # levels(clinic) <- c('P','P','P','P','P','N','N','N')
  # clinic
  # seurat_merge <- AddMetaData(object = seurat_merge,
  #                               metadata = clinic,
  #                               col.name = 'clinic')
  # }
# 任何时候都要创建.RData对象保存进度 
  save(seurat_merge, file="./Rdata/seurat_merged_meta.RData")
}

### 2.2 质控指标 
if(F){
# 细胞计数 (Cell counts)
  metadata %>%    	
    ggplot(aes(x=study, fill=study)) +    	
    geom_bar() +   	
    theme_classic() +   	
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +   	
    theme(plot.title = element_text(hjust=0.5, face="bold")) +   	
    ggtitle("NCells")
  ggsave("./results_qc/NCells.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)

# 每个细胞的UMI计数 (UMI counts per cell)
  metadata %>%    	
    ggplot(aes(color=study, x=nCount_RNA, fill= study)) +    	
    geom_density(alpha = 0.2) +    	
    scale_x_log10() +    	
    theme_classic() +   	
    ylab("Cell density") +   	
    geom_vline(xintercept = 500) 
  ggsave("./results_qc/UMIpercell.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)
  
  metadata %>%    	
    ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) +    	
    geom_density(alpha = 0.2) +    	
    scale_x_log10() +    	
    theme_classic() +   	
    ylab("Cell density") +   	
    geom_vline(xintercept = 500) 
  ggsave("./results_qc/UMIpercell_LHC.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)
  
# 每个细胞检测到的基因 (Genes detected per cell)
  # 通过频数图可视化每个细胞检测出的基因数分布 
  metadata %>%    	
    ggplot(aes(color=study, x=nFeature_RNA, fill= study)) +    	
    geom_density(alpha = 0.2) +    	
    theme_classic() +   	
    scale_x_log10() +    	
    geom_vline(xintercept = 300) 
  ggsave("./results_qc/genepercell.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)
  
  # 通过箱线图可视化每个细胞检测到的基因的分布 
  metadata %>%    	
    ggplot(aes(x=study, y=log10(nFeature_RNA), fill=study)) +    	
    geom_boxplot() +    	
    theme_classic() +   	
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +   	
    theme(plot.title = element_text(hjust=0.5, face="bold")) +   	
    ggtitle("NCells vs nGenes") 
  ggsave("./results_qc/genepercellbox.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)
  
# 检测到的UMI数对比基因数 (UMIs vs. genes detected)
# 可视化检测到的基因数和UMI数之间的关系，并且观察是否存在大量低数目的基因数/UMI数的细胞 
  metadata %>%    	
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) +    	
    geom_point() +  	
    scale_colour_gradient(low = "gray90", high = "black") +   	
    stat_smooth(method=lm) +   	
    scale_x_log10() +    	
    scale_y_log10() +    	
    theme_classic() +   	
    geom_vline(xintercept = 500) +   	
    geom_hline(yintercept = 250) +   	
    facet_wrap(~study)
  ggsave("./results_qc/geneumiratio.jpeg", device = 'jpeg', width = 10, height = 5, unit = "in", dpi = 300)
  
# 可视化每个细胞检测到的线粒体基因表达分布，线粒体计数比率 (Mitochondrial counts ratio)
  metadata %>%    	
    ggplot(aes(color=study, x=mitoRatio, fill=study)) +    	
    geom_density(alpha = 0.2) +    	
    scale_x_log10() +    	
    theme_classic() +   	
    geom_vline(xintercept = 0.2)
  ggsave("./results_qc/mitoratio.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)
  
# 红细胞比率计算
  metadata %>%    	
    ggplot(aes(color=study, x=percent.HB, fill=study)) +    	
    geom_density(alpha = 0.2) +    	
    scale_x_log10() +    	
    theme_classic() +   	
    geom_vline(xintercept = 0.01)
  ggsave("./results_qc/hbcratio.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)  
  
# 通过可视化每一个UMI检测到的基因数来可视化基因表达的整体复杂性，复杂度 (Novelty)
  metadata %>%   	
    ggplot(aes(x=log10GenesPerUMI, color = study, fill=study)) +   	
    geom_density(alpha = 0.2) +   	
    theme_classic() +   	
    geom_vline(xintercept = 0.8)
  ggsave("./results_qc/geneperumi.jpeg", device = 'jpeg', width = 5, height = 5, unit = "in", dpi = 300)
}

### 2.3 使用选择的阈值筛掉低质量读写
if(F){
  dim(seurat_merge@assays$RNA)
  seurat_filtered <- subset(seurat_merge, (nCount_RNA >= 500) &
                                          (nFeature_RNA >= 300) & 
                                          (log10GenesPerUMI > 0.80) & 
                                          (mitoRatio < 0.20) &
                                          (percent.HB < 0.01) )                         
  dim(seurat_filtered@assays$RNA)
  table(seurat_filtered$study)
}

### 2.4 筛去零表达基因
if(F){
# 提取计数 
  counts <- GetAssayData(object = seurat_filtered, slot = "counts")  
# 根据在每个细胞的计数是否大于0为每个基因输出一个逻辑向量
  nonzero <- counts > 0  
# 将所有TRUE值相加，如果每个基因的TRUE值超过10个，则返回TRUE。
  keep_genes <- Matrix::rowSums(nonzero) >= 10  
# 仅保留那些在10个以上细胞中表达的基因 
  filtered_counts <- counts[keep_genes, ]  
# 重新赋值给经过过滤的Seurat对象 
  seurat_filtered <- CreateSeuratObject(filtered_counts, meta.data = seurat_filtered@meta.data) 
  dim(seurat_filtered@assays$RNA)
  
  save(seurat_filtered, file="./Rdata/seurat_filtered.RData")
}

### 2.5 观察批次效应
if(F){
  seurat_batch <- NormalizeData(seurat_immune_Mast)
  #seurat_batch <- FindVariableFeatures(seurat_batch, selection.method = "vst", nfeatures = 645)
  seurat_batch <- FindVariableFeatures(seurat_batch, selection.method = "mean.var.plot")
  seurat_batch <- ScaleData(seurat_batch)
  seurat_batch <- RunPCA(seurat_batch, npcs = 50)
  seurat_batch <- RunUMAP(seurat_batch, dims = 1:15)
  seurat_batch <- FindNeighbors(seurat_batch, reduction = "pca", dims = 1:15)
  seurat_batch <- FindClusters(seurat_batch, resolution = 0.3)
  
  DimPlot(seurat_batch, reduction = "pca",split.by = 'study',raster=FALSE) 
  ggsave("./results_qc/checkbatch_PCA_mvp.jpeg", device = 'jpeg', width = 20, height = 5, unit = "in", dpi = 300, limitsize = FALSE)
  DimPlot(seurat_batch, reduction = 'umap', split.by = 'group', label = TRUE, label.size = 1)
  ggsave("./results_qc/checkbatch_UMAP_mvp.jpeg", device = 'jpeg', width = 20, height = 5, unit = "in", dpi = 300, limitsize = FALSE)
  DimPlot(seurat_batch, reduction = 'umap', split.by = 'sample', label = TRUE, label.size = 1)
  ggsave("./results_qc/checkbatch_UMAP_sample.jpeg", device = 'jpeg', width = 10, height = 5, unit = "in", dpi = 300, limitsize = FALSE)
  VariableFeaturePlot(seurat_batch)
  ggsave("./results_qc/checkbatch_VariableFeaturePlot_mvp.jpeg", device = 'jpeg', width = 10, height = 10, unit = "in", dpi = 300, limitsize = FALSE)
  
  metadata <- seurat_filtered@meta.data
  table(metadata$study)
  metadata$checkbatch <- NA
  metadata$checkbatch[which(str_detect(metadata$study, "AD1"))] <- "batch1"
  metadata$checkbatch[which(str_detect(metadata$study, "AD2"))] <- "batch1"
  metadata$checkbatch[which(str_detect(metadata$study, "AS"))] <- "batch2"
  metadata$checkbatch[which(str_detect(metadata$study, "CRS"))] <- "batch3"

  
  seurat_filtered@meta.data <- metadata
  head(seurat_filtered@meta.data)
  tail(seurat_filtered@meta.data)
}

#################  3.整合  ##################
### 3.1 检查细胞周期
if(F){
# 标准化计数
  seurat_phase <- NormalizeData(seurat_filtered)
# 导入细胞周期标记物 
  load("./prepare/cycle.rda")  
# 给细胞的细胞周期评分 
  seurat_phase <- CellCycleScoring(seurat_phase,                                   
                                  g2m.features = g2m_genes,                                   
                                  s.features = s_genes)  
# 观察细胞周期评分和分配给每个细胞的周期                                
  head(seurat_phase@meta.data) 
# 找出变异最大的2000个基因 
  seurat_phase <- FindVariableFeatures(seurat_phase,                       
                                     selection.method = "vst",                      
                                     nfeatures = 645,                       
                                     verbose = FALSE) 
# 缩放计数 
  seurat_phase <- ScaleData(seurat_phase)
# 运行 PCA 
  seurat_phase <- RunPCA(seurat_phase)  
# 使用细胞周期阶段对PCA进行着色 
    DimPlot(seurat_phase,         
            reduction = "pca",         
            group.by= "study",         
            split.by = "Phase")
    ggsave("./results_qc/cellcycle_study.jpeg", device = 'jpeg', width = 20, height = 10, unit = "in", dpi = 300)
    DimPlot(seurat_phase,         
            reduction = "pca",         
            group.by= "sample",         
            split.by = "Phase")
    ggsave("./results_qc/cellcycle_sample.jpeg", device = 'jpeg', width = 20, height = 10, unit = "in", dpi = 300)
}

### 3.2 SCTransform
if(F){
  options(future.globals.maxSize = 192000 * 1024^2) 
  
  head(seurat_filtered@meta.data)
  Idents(seurat_filtered) <- 'checkbatch'
# 按条件分割seurat对象，对所有样品进行细胞周期评分和SCT
split_seurat <- SplitObject(seurat_filtered, split.by = "checkbatch")
split_seurat <- split_seurat[c(1:5)]
save(split_seurat, g2m_genes, s_genes, file = './Rdata/forsct.Rdata')
load(file = '~/R/adas/Rdata/forsct.Rdata')
for (i in 1:length(split_seurat)) { 
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)     
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)     
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
}

# 查看在对象中储存的分析项目
  split_seurat$batch1@assays
# 选择变异最大的基因进行聚合 
  integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 645, fvf.nfeatures = 645) 
# 准备好SCT列表对象进行聚合 
  split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features) 
# 寻找最佳伙伴 —— 需要一定时间运行 
  integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, normalization.method = "SCT",
                                        anchor.features = integ_features)
# 跨条件聚合 
  seurat_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")
# 保存聚合的R对象
  save(seurat_integrated,file = './Rdata/seurat_integrated.Rdata')
  DefaultAssay(seurat_integrated)
  seurat_integrated@meta.data %>% head()
  seurat_integrated
  Idents(seurat_integrated) <- 'group'
}

#################  4.降维  ##################
### 4.1 PCA
if(F){
# 大数据可以调整npcs=100，先计算100个pc，通过JackStraw确定筛选需要的维度数如75
  seurat_integrated <- RunPCA(seurat_integrated, npcs = 14)
  print(x = seurat_integrated[["pca"]], dims = 1:10, nfeatures = 5)
# 可视化
  PCAPlot(seurat_integrated,split.by = "group")
  ggsave("./results_qc/PCA.jpeg", device = 'jpeg', width = 20, height = 3, unit = "in", dpi = 300)
  
  VizDimLoadings(seurat_integrated, dims = 1:2, reduction = "pca")
  ggsave("./results_qc/PCA_gene.jpeg", device = 'jpeg', width = 10, height = 10, unit = "in", dpi = 300)
  
  DimPlot(object = seurat_integrated, reduction = "pca",group.by = 'group')
  ggsave("./results_qc/PCA_scater.jpeg", device = 'jpeg', width = 10, height = 10, unit = "in", dpi = 300)
  
  DimHeatmap(seurat_integrated, dims = 1:9, cells = 500, balanced = TRUE)
  ggsave("./results_qc/PCA_heatmap.jpeg", device = 'jpeg', width = 10, height = 10, unit = "in", dpi = 300)
}

### 4.2 JackStraw确定维度
if(F){
  seurat_integrated <- JackStraw(seurat_integrated, dims = 100, num.replicate = 100) ##dims根据PCA
  seurat_integrated <- ScoreJackStraw(seurat_integrated, dims = 1:100)

  JackStrawPlot(seurat_integrated, dims = 1:100)
  ggsave("./results_qc/jsplot_k.jpeg", device = 'jpeg', width = 20, height = 14, unit = "in", dpi = 300)
  
  ElbowPlot(seurat_integrated,ndims = 100)  #for big dataset
  ggsave("./results_qc/elbow_k.jpeg", device = 'jpeg', width = 10, height = 10, unit = "in", dpi = 300)
 
  if(F){ 
  # Determine percent of variation associated with each PC
  pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  co1
  
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # last point where change of % of variation is more than 0.1%.
  co2
  
  pcs <- min(co1, co2)
  
  pcs
  
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  
  # Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
  
  ggsave("./results_qc/elbow.jpeg", device = 'jpeg', width = 10, height = 10, unit = "in", dpi = 300)
  }
}

#################  5.聚类  ##################
### 5.1 聚类
if(F){
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "pca", dims = 1:14) #根据Elbow图确定维度数
seurat_integrated <- FindClusters(seurat_integrated, resolution = c(0.25))  #250k细胞数的推荐resolution是3
tmp <- seurat_integrated@meta.data
tmp <- tmp[,c(-28,-27,-26,-24,-23)]
seurat_integrated@meta.data <- tmp
seurat_integrated@meta.data %>% head()
seurat_integrated@meta.data$seurat_clusters %>% table()
seurat_integrated@meta.data$integrated_snn_res.0.5 %>% table()
seurat_integrated@meta.data
}

### 5.2 clustree确定分辨率
if(F){
#https://lazappi.github.io/clustree/articles/clustree.html#controlling-aesthetics
  clustree(seurat_integrated,prefix = "integrated_snn_res.")
  ggsave("./results_qc/clustree_0.5.jpeg", device = 'jpeg', width = 10, height = 10, unit = "in", dpi = 300)
  
  #clustree(seurat_integrated,prefix = "integrated_snn_res.", node_size = 10, node_alpha = 0.8)
  jpeg(file="clustree_stab.jpeg", height = 25, width = 25, units = 'in', res = 300)
    clustree(seurat_integrated,prefix = "integrated_snn_res.", node_colour = "sc3_stability")
  dev.off()
  #clustree(seurat_integrated,prefix = "integrated_snn_res.", node_colour = "FLG", node_colour_aggr = "median")
  jpeg(file="clustree_FLG.jpeg", height = 25, width = 25, units = 'in', res = 300)
    clustree(seurat_integrated,prefix = "integrated_snn_res.", node_colour = "FLG", node_colour_aggr = "median") +
      guides(edge_colour = FALSE, edge_alpha = FALSE) +
      theme(legend.position = "bottom")
  dev.off()
  
  seurat_integrated$seurat_clusters <- seurat_integrated$integrated_snn_res.2.5
  Idents(seurat_integrated) <- 'seurat_clusters'
}

### 5.3 表格
if(F){
  Idents(seurat_integrated)
  table(seurat_integrated$seurat_clusters)
  head(seurat_integrated@meta.data)
  n_cells_sample <- FetchData(seurat_integrated, 
                       vars = c("seurat_clusters", "sample")) %>%
    dplyr::count(seurat_clusters, sample) %>%
    tidyr::spread(seurat_clusters, n)
  n_cells_sample <- t(n_cells_sample)
  n_cells_sample <- as.data.frame(n_cells_sample)
  colnames(n_cells_sample) <- c("HC","L","N")
  n_cells_sample <- n_cells_sample[-1,]
  
  n_cells_orig <- FetchData(seurat_integrated, 
                            vars = c("seurat_clusters", "orig.ident")) %>%
    dplyr::count(seurat_clusters, orig.ident) %>%
    tidyr::spread(seurat_clusters, n)
  n_cells_orig <- t(n_cells_orig)
  n_cells_orig <- as.data.frame(n_cells_orig)
  colnames(n_cells_orig) <- n_cells_orig[1,]
  n_cells_orig <- n_cells_orig[-1,]
  
  n_cells <- n_cells_sample
  
  n_cells$HC <- c(3956,1884, 2074, 2187, 1455,  634,  799,  329,  159)
  n_cells$L <- c(5454, 3721, 3109, 1966, 1384,  839,  832,  626,  289)
  n_cells$N <- c(7279, 4141, 3713, 1974, 1156, 1108,  922,  670,  228)
  
  n_cells <- cbind(n_cells_orig, n_cells_sample)
  n_cells$cluster <- rownames(n_cells)
  n_cells <- as.data.frame(lapply(n_cells,as.numeric))
  n_cells[22,6] <- 0
  rownames(n_cells) <- n_cells$cluster
  
  n_cells$HC1_percent <- n_cells$HC1*100 / sum(n_cells$HC1)
  n_cells$HC2_percent <- n_cells$HC2*100 / sum(n_cells$HC2)
  n_cells$L1_percent <- n_cells$L1*100 / sum(n_cells$L1)
  n_cells$L2_percent <- n_cells$L2*100 / sum(n_cells$L2)
  n_cells$L3_percent <- n_cells$L3*100 / sum(n_cells$L3)
  n_cells$N1_percent <- n_cells$N1*100 / sum(n_cells$N1)
  n_cells$N2_percent <- n_cells$N2*100 / sum(n_cells$N2)
  n_cells$N3_percent <- n_cells$N3*100 / sum(n_cells$N3)
  n_cells$HC_percent <- (n_cells$HC1*100 + n_cells$HC2*100) / (sum(n_cells$HC1) + sum(n_cells$HC2))
  n_cells$L_percent <- (n_cells$L1*100 + n_cells$L2*100 + n_cells$L3*100) / (sum(n_cells$L1) + sum(n_cells$L2) + sum(n_cells$L3))
  n_cells$N_percent <- (n_cells$N1*100 + n_cells$N2*100 + n_cells$N3*100) / (sum(n_cells$N1) + sum(n_cells$N2) + sum(n_cells$N3))
  
  n_cells$HC_percent <- n_cells$HC*100 / sum(n_cells$HC)
  n_cells$L_percent <- n_cells$L*100 / sum(n_cells$L)
  n_cells$N_percent <- n_cells$N*100 / sum(n_cells$N)

  write.csv(n_cells, file = './results/clusters_cellnumber.csv')
}

### 5.4 UMAP/TSNE
if(F){
  seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:14) #根据Elbow图确定维度数
  seurat_integrated <- RunTSNE(seurat_integrated, dims = 1:50)
  save(seurat_integrated, file = './Rdata/seurat_integrated_cluster.Rdata')
}

### 5.5 可视化
if(F){
  head(seurat_integrated@meta.data)
  Idents(seurat_integrated) <- 'seurat_clusters'
  
  #Idents(seurat_integrated) <- factor(x = Idents(seurat_integrated), levels = sort(as.numeric(levels(seurat_integrated))))
  
  jpeg(file="./results/UMAP_1.jpeg", height = 5, width = 12, units = 'in', res = 300)
    p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "disease", raster=FALSE)
    p2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, raster=FALSE)
    plot_grid(p1, p2)
  dev.off()
  
  jpeg(file="./results/UMAP_2.jpeg", height = 10, width = 30, units = 'in', res = 300)
    DimPlot(seurat_integrated, reduction = "umap", split.by = "disease", raster=FALSE)
  dev.off()
  
  jpeg(file="./results/UMAP_3.jpeg", height = 10, width = 20, units = 'in', res = 300)
    p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "disease", raster=FALSE)
    p2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE, raster=FALSE)
    plot_grid(p1, p2)
  dev.off()
  
  jpeg(file="./results/UMAP_4.jpeg", height = 5, width = 15, units = 'in', res = 300)
    DimPlot(seurat_integrated, reduction = "umap", split.by = "disease", raster=FALSE)
  dev.off()
  
  jpeg(file="./results/UMAP_5.jpeg", height = 6, width = 15, units = 'in', res = 300)
    DimPlot(seurat_integrated, reduction = "umap", label = TRUE, split.by = "disease", raster=FALSE) +NoLegend()
  dev.off()
  
  jpeg(file="./results/TSNE_1.jpeg", height = 10, width = 20, units = 'in', res = 300)
    p1 <- DimPlot(seurat_integrated, reduction = "tsne", group.by = "sample")
    p2 <- DimPlot(seurat_integrated, reduction = "tsne", label = TRUE)
    plot_grid(p1, p2)
  dev.off()
  
  jpeg(file="./results/TSNE_2.jpeg", height = 10, width = 20, units = 'in', res = 300)
    DimPlot(seurat_integrated, reduction = "tsne", split.by = "sample")
  dev.off()
  
  jpeg(file="./results/TSNE_3.jpeg", height = 10, width = 20, units = 'in', res = 300)
    p1 <- DimPlot(seurat_integrated, reduction = "tsne", group.by = "orig.ident")
    p2 <- DimPlot(seurat_integrated, reduction = "tsne", label = TRUE)
    plot_grid(p1, p2)
  dev.off()
  
  jpeg(file="./results/TSNE_4.jpeg", height = 10, width = 50, units = 'in', res = 300)
    DimPlot(seurat_integrated, reduction = "tsne", split.by = "orig.ident")
  dev.off()
  
  jpeg(file="./results/TSNE_5.jpeg", height = 10, width = 20, units = 'in', res = 300)
    DimPlot(seurat_integrated, reduction = "tsne", label = TRUE, split.by = "sample")  + NoLegend()
  dev.off()
}

### 5.6 周期分离
if(F){
  jpeg(file="./results_qc/UMAP_cycle.jpeg", height = 10, width = 30, units = 'in', res = 300)     
    DimPlot(seurat_integrated, label = TRUE, split.by = "Phase")  + NoLegend()   
  dev.off()
  
  jpeg(file="./results_qc/UMAP_cycle_sample.jpeg", height = 10, width = 30, units = 'in', res = 300)     
    DimPlot(seurat_integrated, label = TRUE, split.by = "Phase", group.by = 'sample')  + NoLegend()   
  dev.off()
  
  jpeg(file="./results_qc/UMAP_cycle_phase.jpeg", height = 10, width = 30, units = 'in', res = 300)     
    DimPlot(seurat_integrated, label = TRUE, split.by = "sample", group.by = 'Phase')  
  dev.off()
}

### 5.7 其他变异源
if(F){
  jpeg(file="./results/UMAP_var.jpeg", height = 25, width = 20, units = 'in', res = 300)     
    FeaturePlot(seurat_integrated,              
              reduction = "umap",              
              features = c("nCount_RNA", "nFeature_RNA", "S.Score", "G2M.Score", "mitoRatio"),             
              pt.size = 0.4,              
              sort.cell = TRUE,             
              min.cutoff = 'q10',             
              label = TRUE)
  dev.off()
}

### 5.8 PC(不做)
if(F){
# 对感兴趣的seurat对象中的信息命名 
  columns <- c(paste0("PC_", 1:16), "ident", "UMAP_1", "UMAP_2")
# 从seurat对象中提取这些数据 
  pc_data <- FetchData(seurat_integrated, vars = columns)
# 提取前10个细胞的UMAP坐标  
  seurat_integrated@reductions$umap@cell.embeddings[1:10, 1:2]
# 在UMAP上为类群中心添加类群标签 
  umap_label <- FetchData(seurat_integrated,                          
                          vars = c("ident", "UMAP_1", "UMAP_2"))  %>%   
    group_by(ident) %>%   
    summarise(x=mean(UMAP_1), y=mean(UMAP_2))    
# 给每一个PC做UMAP图 
  map(paste0("PC_", 1:16), function(pc){         
    ggplot(pc_data,                 
           aes(UMAP_1, UMAP_2)) +                 
      geom_point(aes_string(color=pc),                             
                 alpha = 0.7) +                 
      scale_color_gradient(guide = FALSE,                                       
                           low = "grey90",                                       
                           high = "blue")  +                 
      geom_text(data=umap_label,                            
                aes(label=ident, x, y)) +                 
      ggtitle(pc) 
    }) %>%          
    plot_grid(plotlist = .)
# 检验PCA结果
  print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)
}  

### 5.9 Feature(在后面做)
if(F){

  DefaultAssay(seurat_integrated) <- 'RNA'
  jpeg(file="./results/Feature_f1.jpeg", height = 80, width = 40, units = 'in', res = 300) 
    FeaturePlot(seurat_integrated,              
                reduction = "umap",              
                features = c("KRT1","KRT10","KRT5","KRT14","KRT6A","KRT6B","KRT6C","KRT16","FLG","SLURP1","CHRNA7",
                             "KRT17","KRT2","KRT19","AQP5",'KRT15','KRT77','FABP4'), 
                label = TRUE)
  dev.off()
}
  
#################  6.标志物  ##################
### 6.1 单样本
if(F){
  DefaultAssay(seurat_integrated) <- "RNA"
  seurat_integrated_markers <- FindAllMarkers(seurat_integrated, 
                                              min.pct = 0.25, 
                                              logfc.threshold = 0.25) # only.pos可选择仅上调
}

### 6.2 单样本标志物评估
if(F){
  seurat_integrated_markers_top10 = seurat_integrated_markers %>% 
    group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
  seurat_integrated_markers_bottom10 = seurat_integrated_markers %>% 
    group_by(cluster) %>% top_n(n = -10, wt = avg_log2FC) 

  save(seurat_integrated, seurat_integrated_markers, 
       seurat_integrated_markers_top10, seurat_integrated_markers_bottom10, 
       file = './Rdata/seurat_intergrated_markers.Rdata')
  write.csv(seurat_integrated_markers, file = './results/markers_krt.csv')
  write.csv(seurat_integrated_markers_top10, file = './results/markers_krt_top10.csv')
  
  seurat_integrated.markers_logFC1.5=seurat_integrated.markers[seurat_integrated.markers$avg_logFC > 1.5,]
  table(seurat_integrated.markers_logFC1.5$cluster)
  seurat_integrated.markers_logFC2=seurat_integrated.markers[seurat_integrated.markers$avg_logFC > 2,]
  table(seurat_integrated.markers_logFC2$cluster)
  seurat_integrated.markers_logFC1=seurat_integrated.markers[seurat_integrated.markers$avg_logFC > 1,]
  table(seurat_integrated.markers_logFC1$cluster)
  }

### 6.3 多样本
if(F){
  #整合完control和疾病组样本，并做聚类分析后，分别在control和疾病组中对某个cluster做差异表达分析，寻找这个cluster的marker基因
  #然后，把在control和疾病组中都被鉴定为marker基因的基因，作为该cluster保守的marker基因，并利用这些基因对细胞类型进行注释。
  #最后，针对同一类细胞，寻找在control和疾病组中发生变化的基因。  
  
  DefaultAssay(seurat_integrated) <- "RNA"
# 添加基因注释  
  annotations <- read.csv("./prepare/annotation.csv")

# 创建用来获取给定类群保留标记物的函数
  get_conserved <- function(cluster){
    FindConservedMarkers(seurat_integrated,
                         ident.1 = cluster,
                         grouping.var = "sample",
                         only.pos = TRUE, 
                         logfc.threshold = 0.25,
                         verbose = T) %>%
      rownames_to_column(var = "gene") %>%
      left_join(y = unique(annotations[, c("gene_name", "description")]),
                by = c("gene" = "gene_name")) %>%
      cbind(cluster_id = cluster, .)
  }
# 在需要的类群中使用函数 
  seurat_integrated@active.ident
  conserved_markers_13 <- map_dfr(c(1), get_conserved)
  conserved_markers_all <- map_dfr(c(0:21), get_conserved)
}

### 6.4 多样本标志物评估
if(F){
# 提取每个类群排名靠前的10个标记物 
  conserved_markers_krt_top10 <- conserved_markers %>%                                   
                              mutate(avg_fc = (HC_avg_logFC + L_avg_logFC + N_avg_logFC) /3) %>%                                   
                              group_by(cluster_id) %>%                                   
                              top_n(n = 10, wt = avg_fc)  
# 查看每个类群排名靠前的10个标记物 
  View(conserved_markers_top10)
  write.csv(conserved_markers_gdT, file = './results/conserved_markers_gdT.csv')
  write.csv(conserved_markers, file = './results/conserved_markers.csv')
}

### 6.5 可视化
if(F){
  DefaultAssay(seurat_integrated)
  jpeg(filename = './results/Feature_ISG.jpeg',height = 10 ,width = 20, units = 'in', res = 300)
    FeaturePlot(object = seurat_integrated,                          
                features = c("IF"),                          
                min.cutoff = 'q10',                           
                label = TRUE, 			             
                repel = TRUE,
                split.by = 'sample',
                slot = 'data'
                )
  dev.off()
  
  jpeg(filename = './results_PN//Feature_sample.jpeg',height = 20 ,width = 15, units = 'in', res = 300)     
    FeaturePlot(object = seurat_integrated,                                           
                features = c("IL7", "IL9", "IL17A", "IL4"), 
                split.by = 'sample',
                sort.cell = TRUE,                                           
                min.cutoff = 'q10',                                            
                label = TRUE, 			                              
                repel = TRUE)   
  dev.off()
  
  jpeg(filename = './results_PN/Violin.jpeg',height = 20 ,width = 25, units = 'in', res = 300)
    VlnPlot(object = seurat_integrated, 
            features = c("IL7", "IL9", "IL17A", "IL4"), group.by = 'sample')
  dev.off()
  
  DefaultAssay(seurat_integrated) <- 'integrated'
  DefaultAssay(seurat_integrated) <- 'RNA'
  
  seurat_integrated@assays$RNA@data %>% head
  seurat_integrated@assays$RNA@scale.data[1:5,1:5]
  seurat_integrated <- ScaleData(seurat_integrated)
  
    DoHeatmap(seurat_integrated, 
              features = conserved_markers_krt_top10$gene, group.by = 'seurat_clusters') + scale_fill_viridis()
    ggsave("./results/Heatmap_krt.jpeg", device = 'jpeg', width = 35, height = 25, unit = "in", dpi = 300, limitsize = FALSE)
    ggsave("./results/Heatmap_krt.pdf", device = 'pdf', width = 30, height = 20, unit = "in", dpi = 300, limitsize = FALSE)
    
    DoHeatmap(seurat_integrated, 
              features = seurat_integrated_markers_top10$gene, group.by = 'clu_sample') + scale_fill_viridis() + NoLegend()
    ggsave("./results/Heatmap_int_sample.jpeg", device = 'jpeg', width = 70, height = 25, unit = "in", dpi = 300, limitsize = FALSE)
    }

### 6.6 多样本差异可视化
if(F){
#样本间差异
  seurat_integrated$clu_sample <- paste(seurat_integrated$sample, 
                                        seurat_integrated$seurat_clusters,sep = "_")
  seurat_integrated$clu_sample_krt <- paste(seurat_integrated$sample,                                                                                seurat_integrated$seurat_clusters,sep = "_")
  head(seurat_integrated@meta.data)
  #seurat_integrated$celltype <- Idents(seurat_integrated) 
  Idents(seurat_integrated) <- "clu_sample_krt" 
  DEG_clu_sample_9 <- FindMarkers(seurat_reference_Tcell, 
                                ident.1 = "9_Patient", ident.2 = "9_Norm", 
                                verbose = FALSE) 
  head(DEG_clu_sample_9, n = 15) 
  
#双样本点图  
  markers.to.plot <- c("SLURP1","CHRNA7","KRT10","CDKN1A","CASP3","TGM1","TNF", 
                       "KRT7", "KRT8", "KRT18", "KRT19", "COL17A1", "COL1A1", 
                       "COL1A2", "TYRP1", "PMEL", "MLANA", "SCGB2A2", "S100A4", 
                       "HLA-Dseurat", "HLA-DQB1", "CD74","FLG","LOR")      
  jpeg(file="./results/DotPlot.jpeg", height = 50, width = 50, units = 'in', res = 300)
    DotPlot(seurat_integrated, features = rev(markers.to.plot), 
            cols = c("blue", "red"), dot.scale = 8,            
            split.by = "stim", assay = 'RNA') + RotatedAxis()
  dev.off()
  
#双样本差异表达散点图
  head(seurat_integrated@meta.data)
  Idents(seurat_integrated) <- "seurat_clusters"
  Idents(seurat_integrated)
  DefaultAssay(seurat_integrated) <- 'RNA'
  theme_set(theme_cowplot())
  names <- names(table(seurat_integrated$seurat_clusters))
  Idents(seurat_integrated) <- 'clu_sample_krt'
  
  for (i in names[1:27]) {
    assign(paste0('cluster_', i), subset(seurat_integrated, idents = i) )
   # assign(paste0('avg_cluster_', i),  log2(AverageExpression(get(paste0('cluster_', i)), 
   #                                                                   verbose = FALSE, add.ident = 'sample')$RNA+1) )
    #assign(paste0('avg_cluster_',i), 
        #   add_column(get(paste0('avg_cluster_', i)), gene = rownames(get(paste0('avg_cluster_',i))) ) )
   # write.csv(get(paste0('avg_cluster_', i)), file = paste0('./results_avg/avg_cluster_', i,'.csv'))
  }
  
if(F){
  cluster_9 <- subset(seurat_reference_Tcell, idents = "9") 
  Idents(cluster_9) <- "sample" 
  Idents(cluster_9)
  avg.cluster_9 <- log1p(AverageExpression(cluster_9, verbose = FALSE )$RNA) 
  avg.cluster_9$gene <- rownames(avg.cluster_9)  
  
  cd14.mono <- subset(seurat_reference_Tcell, idents = "CD14 Mono") 
  Idents(cd14.mono) <- "stim" 
  avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA) 
  avg.cd14.mono$gene <- rownames(avg.cd14.mono)  
  
  Idents(cluster_gdT) <- "sample" 
  genes.to.label = DEGPN_gdT$gene
  p1 <- ggplot(avg.cluster_gdT, aes(gdT_Patient, gdT_Norm)) + geom_point() + ggtitle("cluster_gdT") 
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = FALSE) 
  p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes") 
  p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE) 
  plot_grid(p1, p2)
}
  seurat_integrated@meta.data[1:5,]
  seurat_integrated$seurat_clusters
  plots <- VlnPlot(seurat_integrated, features = inflamationgene, 
                   split.by = "sample", group.by = "seurat_clusters", 
                   pt.size = 0, combine = FALSE, fill.by = 'sample')
  wrap_plots(plots = plots, ncol = 3)
  ggsave("./results/Violin_k_inf.jpeg", device = 'jpeg', width = 20, height = 20, unit = "in", dpi = 300, limitsize = FALSE)
  
  RidgePlot(seurat_integrated, features = inflamationgene, ncol = 2, group.by = 'sample')
  ggsave("./results/ridge_k_inf.jpeg", device = 'jpeg', width = 20, height = 60, unit = "in", dpi = 300, limitsize = FALSE)
}

### 6.7 区别相似类群
if(F){
# 确定CD4阳性T细胞的特异标记物
  Idents(seurat_integrated) <- 'sample'
  DEGLN_krt <- FindMarkers(seurat_integrated,ident.1 = "L",ident.2 = "N")
  DEGLHC_krt <- FindMarkers(seurat_integrated,ident.1 = "L",ident.2 = "HC")
  DEGNHC_krt <- FindMarkers(seurat_integrated,ident.1 = "N",ident.2 = "HC")
  DEGLnL_krt <- FindMarkers(seurat_integrated,ident.1 = "L")
  DEGLN_krt$gene <- rownames(DEGLN_krt)
  DEGLHC_krt$gene <- rownames(DEGLHC_krt)
  DEGNHC_krt$gene <- rownames(DEGNHC_krt)
  DEGLnL_krt$gene <- rownames(DEGLnL_krt)
  
  Idents(seurat_integrated) <- 'clu_sample'
  DefaultAssay(seurat_integrated)<- 'RNA'
  
  for (i in 0:14){
    assign(paste0('DEGADL_',i),
           FindMarkers(paste0('cluster_',i),                           
                       ident.1 = paste0('L_', i),                           
                       ident.2 = paste0('N_', i)))
  # 给DE表格添加基因标志 
    assign(paste0('DEGLN_krt_',i) , get(paste0('DEGLN_krt_',i)) %>%   
             rownames_to_column(var = "gene") %>%   
             left_join(y = unique(annotation[, c("gene_name", "description")]),              
                       by = c("gene" = "gene_name")))
    write.csv(get(paste0('DEGLN_krt_',i)), file = paste0('./results_LN/DEGLN_krt_',i ,'.csv'))
  }
  
  for (i in 0:8){
    assign(paste0('DEGLHC_krt_',i),
           FindMarkers(seurat_integrated,                           
                       ident.1 = paste0('L_', i),                           
                       ident.2 = paste0('HC_', i)))
    # 给DE表格添加基因标志 
    assign(paste0('DEGLHC_krt_',i) , get(paste0('DEGLHC_krt_',i)) %>%   
             rownames_to_column(var = "gene") %>%   
             left_join(y = unique(annotation[, c("gene_name", "description")]),              
                       by = c("gene" = "gene_name")))
    write.csv(get(paste0('DEGLHC_krt_',i)), file = paste0('./results_LHC/DEGLHC_krt_',i ,'.csv'))}
    
  for (i in 0:8){
      assign(paste0('DEGNHC_krt_',i),
             FindMarkers(seurat_integrated,                           
                         ident.1 = paste0('N_', i),                           
                         ident.2 = paste0('HC_', i)))
      # 给DE表格添加基因标志 
      assign(paste0('DEGNHC_krt_',i) , get(paste0('DEGNHC_krt_',i)) %>%   
               rownames_to_column(var = "gene") %>%   
               left_join(y = unique(annotation[, c("gene_name", "description")]),              
                         by = c("gene" = "gene_name")))
      write.csv(get(paste0('DEGNHC_krt_',i)), file = paste0('./results_NHC/DEGNHC_krt_',i ,'.csv'))
  }
}


# 按照p.adj给列重新排序 
  cd8_1_14 <- cd8_1_14[, c(1, 3:5,2,6:7)]  
  cd8_1_14 <- cd8_1_14 %>%   
    dplyr::arrange(p_val_adj)   
# 查看数据 
  View(cd8_1_14)
  write_csv(nk_9,file = './results/nk_9.csv')


### 6.8 SingleR
if(F){
  hpca.se <- HumanPrimaryCellAtlasData() 
  bp.se <- BlueprintEncodeData()
  dice.se <- DatabaseImmuneCellExpressionData()
  nhd.se <- NovershternHematopoieticData()
  mid.se <- MonacoImmuneData()
  save(hpca.se, mid.se, nhd.se, bp.se, dice.se,file = './prepare/singleR_database.Rdata')
#单数据库
  clu_singleR <- SingleR(test = GetAssayData(seurat_integrated), 
                         ref = hpca.se, labels = hpca.se$label.main,
                         method = 'cluster', 
                         clusters = seurat_integrated$seurat_clusters) 
  
  clu_singleR
  dim(clu_singleR)
  clu_singleR$labels
  table(clu_singleR$labels)

#多数据库
  clu_singleR_combine <- SingleR(test = GetAssayData(seurat_integrated), 
                                 ref = list(BP=bp.se, HPCA=hpca.se, NHD=nhd.se, MID=mid.se), 
                                 labels = list(bp.se$label.main, hpca.se$label.main, 
                                               nhd.se$label.main, mid.se$label.main),
                                 method = 'cluster', clusters = seurat_integrated$seurat_clusters)
  clu_singleR_combine
  dim(clu_singleR_combine)
  table(clu_singleR_combine$labels)
  

  write.csv(clu_singleR, file = './results/singleR.csv')
}

### 6.9 命名
if(F){
# 重命名所有类群的身份
  seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                    "0" = "Spinous1",
                                    "1" = "Spinous2",
                                    "2" = "Basal1",
                                    "3" = "Basal2",
                                    "4" = "PPKeratinocyte1",
                                    "5" = "Granular",
                                    "6" = "PPKeratinocyte2",
                                    "7" = "Endothelial1",
                                    "8" = "Mitosis",
                                    "9" = "Fibroblast",
                                    "10" = "Spinous3",
                                    "11" = "SmoothMusle",
                                    "12" = "Mitochondrion",
                                    "13" = "UpperSpinous",
                                    "14" = "Gland",
                                    "15" = "Periderm",
                                    "16" = "Melanocyte",
                                    "17" = "Bulge", 
                                    "18" = "Langerhans", 
                                    "19" = "Tcells", 
                                    "20" = "PPKeratinocyte3",
                                    "21" = "Endothelial2")

# 绘制UMAP图
  DimPlot(object = seurat_integrated, 
          reduction = "umap",
          label = TRUE,
          label.size = 3,
          repel = TRUE)
}

### 6.10 删除类群（选做）
if(F){
# 删除受压或濒死细胞
  seurat_integrated_subset <- subset(seurat_integrated,
                                  idents = c("17", "18"), invert = TRUE)
  
# 重新观察聚类
  DimPlot(object = seurat_subset_labeled, 
          reduction = "umap", 
          label = TRUE,
          label.size = 3,
          repel = TRUE)
}

#################  灵活操作  ##################
if(F){
#平均表达量
colnames(seurat_integrated@meta.data)
table(seurat_integrated@meta.data$orig.ident)
Idents(seurat_integrated) = 'orig.ident'
Idents(seurat_integrated)
grep('TNF', rownames(seurat_integrated))
rownames(seurat_integrated)[grep('TNF', rownames(seurat_integrated))]
cluster_averages_EGFR <- AverageExpression(seurat_integrated, features = 'EGFR', assays = 'RNA', add.ident = 'sample') 
cluster_averages_EGFR$RNA
cluster_averages <- AverageExpression(seurat_integrated) 
average <- cluster.averages$RNA
average$FC <- average[,1]/average[,2]
head(cluster.averages[["RNA"]][])
Idents(seurat_integrated) = 'seurat_clusters'

#subset
Tcell.patients = subset(seurat_integrated, clinic == 'P' )
Tcell.norm = subset(seurat_integrated, clinic == 'N' )
test.subset <- subset(epithelial, subset = (stim == "Healthy" | stim == "another_condition"))
Idents(seurat_integrated) <- 'orig.ident'
Idents(seurat_integrated)
seurat_integrated_subset <- subset(seurat_integrated, idents = c('CH','MYH','Nor1','Nor2','Nor3'))
Idents(seurat_integrated_subset) <- 'seurat_clusters'
Idents(seurat_integrated_subset)
seurat_integrated_subset_invert <- subset(seurat_integrated_subset, idents = c('17','18'), invert = TRUE )
}

#################  monocle3  ##################
library(monocle3)

  DefaultAssay(seurat_integrated)<-"integrated"
  data <- seurat_integrated@assays$integrated@data
  pd <-  seurat_integrated@meta.data
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  cds <- new_cell_data_set(data,cell_metadata  = pd,gene_metadata  = fData)
  
  save(data,pd,fData, file = './Rdata/monocle.RDS')
  
  
  
  
  