#!/usr/bin/env Rscript

# This requires two argumens
# @arg[1] a pathway to a folder containing an outs folder of cellranger output
# @arg[2] a sample_name for this run

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# MAx number of PCs to explore
MAX_PCS <- 60

# load libraries
# sctype_env_2
lapply(c("dplyr","Seurat","HGNChelper", "openxlsx", "ggplot2",
         "sctransform", "Signac", "EnsDb.Hsapiens.v86", "data.table"), 
         library, character.only = T)

soup_rate <- 0.20
reticulate::source_python('scripts/utility/scrublet_py.py')

# params <- list()
print(paste0("Sample name: ", sample_name))
# NOTE: IF MANUAL do base_dir <- dir
dir <- '/data/CARD_singlecell/Brain_atlas/NABEC_multiome/batch1/Multiome/SH9608-ARC'
base_dir <- paste0(args[1], '/outs/')
# test with sample <- 'SH9608-ARC' which is in batch1
#
# base_dir <- paste0('/data/CARD_singlecell/Brain_atlas/NABEC_multiome/', batch, '/Multiome/',sample_name, '/outs/')
sample_name <- args[2]

# make directory in figures
if (!dir.exists('figures/') {
        dir.create('figures')
}
figure_folder <- paste0('figures/', sample_name)
dir.create(figure_folder)

# Initialize the Seurat object with the raw (non-normalized data).
input_data <- Read10X(data.dir = paste0(base_dir, "filtered_feature_bc_matrix/"))

    adj.matrix <- suppressWarnings(SoupCorrect(raw.input.list[[dataset]], filtered.input.list[[dataset]], contamination_rate=soup_rate))
    object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=dataset)

    doublet_rate <- (ncol(object) / 1000) * 0.008
    object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
    object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')

    m <- copy(object@meta.data)
    setDT(m, keep.rownames='cells')

    m[, laneID := rep(dataset, nrow(m))]

    m[,
        `:=` (
            laneID=tstrsplit(laneID, '_')[[1]],
            condition=tstrsplit(laneID, '_')[[2]],
            timepoint=tstrsplit(laneID, '_')[[3]]
            )
    ]

    lane <- m[, laneID]
    condition <- m[, condition]
    timepoint <- m[, timepoint]

    names(lane) <- names(condition) <- names(timepoint) <- m[, cells]

    object %>% 
        
        AddMetaData(metadata=factor(lane), col.name='laneID') %>% 
        AddMetaData(metadata=factor(condition), col.name='condition') %>% 
        AddMetaData(metadata=factor(timepoint), col.name='timepoint') %>%
        
        scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate) 

}, simplify=FALSE)

library(GenomicFeatures)
tx <- makeTxDvFromGFF()
transcripts(tx)

# Whether the data are multiome or unimodal
# If they are multiome set multiome <- TRUE
if (length(input_data) > 1){
    rna_counts <- input_data[['Gene Expression']]
    atac_counts <- input_data[['Peaks']]
    data <- CreateSeuratObject(counts = rna_counts, project = sample_name);
    #message("Using RNA data only")
    grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use), ]
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    genome(annotations) <- "hg38"
    seqlevelsStyle(annotations) <- 'UCSC'
    #data.atac <- CreateAssayObject(atac_counts[, colnames(x = data)])
    #data  <- NormalizeData(data, assay = "Peaks", normalization.method = "CLR")
    frag.file <- paste0(base_dir, "atac_fragments.tsv.gz")
    chrom_assay <- CreateChromatinAssay(
        counts = atac_counts,
        sep = c(":", "-"),
        genome = 'hg38',
        fragments = frag.file,
        min.cells = 3,
        annotation = annotations,
        validate.fragments = TRUE
    )
    data[['ATAC']] <- chrom_assay
} else {
    data <- CreateSeuratObject(counts = input_data, project = sample_name, min.cells = 3, min.features = 200)
    data <- ScaleData(data, features = rownames(data))
    data <- RunPCA(data, npcs = MAX_PCS, features = VariableFeatures(object = data))

    # Check number of PC components (we selected 10 PCs for downstream analysis, based on Elbow plot)
    pdf(paste0(figure_folder,'/', sample_name, '_elbow', '.pdf'))
    ElbowPlot(data, ndims=MAX_PCS,reduction='pca')
    dev.off()
    
    data <- JackStraw(data, dims=MAX_PCS, num.replicate=80, maxit=900)
    data <- ScoreJackStraw(data, dims=1:MAX_PCS)
    upper_b <- min(which(JS(data[['pca']], 'overall')[,2] >= 0.05))
    print(paste0("total number of PCs used. 
                Make sense? : ", upper_b))
    if (upper_b == 'Inf'){
        upper_b <- MAX_PCS
    }
    pdf(paste0(figure_folder,'/', sample_name, '_jackstraw', '.pdf'))
    JackStrawPlot(data,dims = 1:upper_b)
    dev.off()
    
    
    # cluster and visualize
    data <- FindNeighbors(data, dims = 1:upper_b)
    data <- FindClusters(data, resolution = 0.2)
    data <- RunUMAP(data, dims = 1:upper_b)
    
    pdf(paste0(figure_folder,'/', sample_name, '_cluster', '.pdf'))
    DimPlot(data, reduction = "umap")
    dev.off()
}

# NOTE: Need to figure out when to take these out and when to simply display
# People seem divided on this.
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
pdf(paste0(figure_folder,'/', sample_name, '_qc', '.pdf'))
VlnPlot(data, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
  log = TRUE, pt.size = 0) + NoLegend()
dev.off()

# data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
# data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

# scale and run PCA



# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
#db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
db_ = "data/ScTypeDB_full.xlsx"
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# Celltype
es.max = sctype_score(scRNAseqData = data[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: FOR NOW JUST USE THE SAME NAMES BECAUSE WE HAVE NO A PRIORI REASON TO MODIFY
# As soon as we do, we can change the names. But that will be a next step
# See https://stuartlab.org/signac/articles/pbmc_vignette.html#create-a-gene-activity-matrix for information on how to get ATAC data to reasonable matching to rnaseq for assimilation into this measure.
es.max_ATAC = sctype_score(scRNAseqData = data[['ATAC']]@counts, scaled=FALSE,gs_list$gs_positive, gs2=gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either data[["RNA"]]@scale.data (default), data[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or data[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results = do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
markers <- FindAllMarkers(data, densify=TRUE, min.pct=0.18)
marker_scores <- merge(sctype_scores, markers, on="cluster")
marker_scores <- marker_scores[c('cluster','type','ncells','p_val','avg_log2FC','pct.1','pct.2','p_val_adj','gene')]
# sort by p_Adj_values
markers <- arrange(marker_scores, cluster, p_val_adj)
write.csv(markers, file=paste0(figure_folder, '/', sample_name, '_markers', '.csv'), row.names=FALSE)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

data@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  data@meta.data$customclassif[data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
g3 <- DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')
ggsave(g3, file=paste0(figure_folder,'/', sample_name, '_umap', '.pdf'), width=25, height=18, units='cm')

# output the number of each cell type
cluster_counts <- data@meta.data %>%
    count(customclassif)
write.csv(cluster_counts, file=paste0(figure_folder, '/', sample_name, '_type_counts', '.csv'), row.names=FALSE)

# load libraries
#
#
lapply(c("ggraph","igraph","tidyverse", "data.tree", "ggthemes"), library, character.only = T)

# prepare edges
cL_results=cL_results[order(cL_results$cluster),]
edges = cL_results
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

for (i in 1:length(unique(cL_results$cluster))){
  dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]
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
gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
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

ggsave(gggr, file=paste0(figure_folder, '/', sample_name, '-bubble.png'), width=25, height=25, units='cm')

library(cowplot)
g4 <-  plot_grid(g3, gggr, labels = c('A', 'B'), label_size = 12, rel_widths=c(0.45, 0.55))
ggsave(g4, file=paste0(figure_folder, '/', sample_name, '-clusters_combined.png'), width=45, height=25, units='cm')


if(FALSE){
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
  db_ = "data/ScTypeDB_full.xlsx";
  tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = data, scaled = TRUE, assay = "RNA")
}

quit()


working_dir <- '/data/ShernData/CS032989_Taylor_Kondo/'; setwd(working_dir)



# Need to figure out how to get the integratedrna, apparently both things are in there
# At least in looking at the rpca.R and rlsi.R files

DefaultAssay(object) <- 'integratedrna'
VariableFeatures(object) <- rownames(object)

object <- object %>% 
            ScaleData(verbose=FALSE) %>% RunPCA(verbose=FALSE) %>% 
            FindNeighbors(reduction='pca', dims=1:50) %>% FindClusters(algorithm=3, resolution=0.5) %>% StashIdent(save.name='rna_clusters') %>% 
            
            RunUMAP(reduction='pca', dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_') %>%
            RunUMAP(reduction='integrated_lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')



object <- object %>% 
            FindMultiModalNeighbors(reduction.list=list('pca', 'integrated_lsi'), dims.list=list(1:50, 2:50)) %>% 
            FindClusters(graph.name='wsnn', algorithm=3, resolution=0.5) %>% StashIdent(save.name='wnn_clusters') %>% 
            RunUMAP(nn.name='weighted.nn', reduction.name='wnn.umap', reduction.key='wnnUMAP_')

object[['seurat_clusters']] <- NULL


saveRDS(object, snakemake@output[['seurat_object']])



