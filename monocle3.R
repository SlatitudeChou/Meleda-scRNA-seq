library(BiocManager)
library(devtools)
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install('units')
sf
spdep
#####预处理
AnnotatedDataFrame(pf)
cds <- newCellDataSet(data,pd,fData)
#scale
cds <- preprocess_cds(cds, num_dim = 30,norm_method = "none")
#umap
cds <- reduce_dimension(cds,umap.n_neighbors = 20L)
#color by seurat cluster
plot_cells(cds,label_groups_by_cluster=FALSE,color_cells_by = "cell.type")
#聚类
cds <- cluster_cells(cds,resolution = 0.5)
#color by monocle cluster
plot_cells(cds, color_cells_by = "partition",label_groups_by_cluster=FALSE)

#####伪时间构建
cds <- learn_graph(cds)
plot_cells(cds,color_cells_by = "Cell.type",label_groups_by_cluster=FALSE,label_leaves=TRUE,label_branch_points=TRUE)
#根节点函数
get_earliest_principal_node <- function(cds, time_bin="Dendritic.cells.activated"){
  cell_ids <- which(colData(cds)[, "Cell.type"] == time_bin)
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
#定义根节点
#order cell
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
#pseudotime analysis
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,graph_label_size=1.5)
#寻找动态基因
diff_gene <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
id <- row.names(subset(diff_gene, q_value < 0.05))
head(id)
#可视化
CASE_genes <- c("CD74", "CXCL14", "HLA-DRA")
CASE_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,colData(cds)$stim %in% c("CASE")]
#plot genes
plot_genes_in_pseudotime(AFD_lineage_cds,color_cells_by="Cell.type")




