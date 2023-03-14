#!/usr/bin/env Rscript


lapply(c("dplyr","Seurat","HGNChelper", "ggplot2", "ggthemes"), library, character.only = T)

dat <- Read10X(data.dir = ".")  # edit if data not in current directory

# Initialize the Seurat object with the raw (non-normalized data).
dat <- CreateSeuratObject(counts = dat$`Gene Expression`, project = "test", min.cells = 3, min.features = 200)

# normalize data
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")
dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)

# scale and run PCA
dat <- ScaleData(dat, features = rownames(dat))
dat <- RunPCA(dat, features = VariableFeatures(object = dat))

# Check number of PC components 
g1 <- ElbowPlot(dat) + theme_few() + geom_hline(yintercept=2.5, linetype='dashed', alpha=0.5)
ggsave(g1, file='elbow_plot.png')


n_pcs <- 10     # went with 10 PCs

# cluster and visualize
dat <- FindNeighbors(dat, dims = 1:n_pcs)
dat <- FindClusters(dat, resolution = 0.8)
dat <- RunUMAP(dat, dims = 1:n_pcs)
g2 <- DimPlot(dat, reduction = "umap")
ggsave(g2, file='umap_1.png', width=25, height=18, units='cm')


# load gene set preparation function
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("gene_sets_prepare.R")
# load cell type annotation function
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source("sctype_score_.R")


# DB file
# db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
db_ = "ScTypeDB_full.xlsx"
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = dat[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either dat[["RNA"]]@scale.data (default), dat[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or dat[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(dat@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(dat@meta.data[dat@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(dat@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

dat@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  dat@meta.data$customclassif[dat@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

g3 <- DimPlot(dat, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') +
guides(color=FALSE) +
labs(x='UMAP 1', y='UMAP 2', title='')

ggsave(g3, file='umap_2.png', width=25, height=18, units='cm')


lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]
edges = cL_resutls
edges$type = paste0(edges$type,"_",edges$cluster)
edges$cluster = paste0("cluster ", edges$cluster)
edges = edges[,c("cluster", "type")]
colnames(edges) = c("from", "to")
rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]
nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster)
nodes_lvl1$Colour = "#f1f1ef"
nodes_lvl1$ord = 1
nodes_lvl1$realname = nodes_lvl1$cluster
nodes_lvl1 = as.data.frame(nodes_lvl1)
nodes_lvl2 = c()

ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a",
            "#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", 
            "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c",
            "#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")

for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]
  nodes_lvl2 = rbind(nodes_lvl2, 
  data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}

nodes = rbind(nodes_lvl1, nodes_lvl2)

nodes$ncells[nodes$ncells<1] = 1

files_db = openxlsx::read.xlsx(db_)[,c("cellName","cellName")]
colnames(files_db) <- c('cellName','shortName')
files_db = unique(files_db)

nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]
nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) +
  geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_few() +
  geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5))) +
  geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted") +
  theme(axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      plot.background = element_blank())

ggsave(gggr, file='bubble.png', width=25, height=25, units='cm')

library(cowplot)
g4 <-  plot_grid(g3, gggr, labels = c('A', 'B'), label_size = 12, rel_widths=c(0.45, 0.55))
ggsave(g4, file='clusters_combined.png', width=45, height=25, units='cm')


if(FALSE){
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = dat, scaled = TRUE, assay = "RNA")
}

quit()