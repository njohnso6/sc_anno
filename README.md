# Testing single cell annotation


## Get data

Retrieve data from (already configured) rclone, hosted on google drive

```bash
module load rclone

rclone copy --drive-shared-with-me  dtigoogledrive:/Xylena_Test_Data/filtered_feature_matrix/features.tsv.gz .
rclone copy --drive-shared-with-me  dtigoogledrive:/Xylena_Test_Data/filtered_feature_matrix/barcodes.tsv.gz .
rclone copy --drive-shared-with-me  dtigoogledrive:/Xylena_Test_Data/filtered_feature_matrix/matrix.mtx.gz .
```


# Run sc-type in `R/4.2`
```R
# Note that you may need to install.packages 
# lapply(c("dplyr","Seurat","HGNChelper", "ggplot2", "ggthemes"), install.packages, character.only=TRUE)

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
```

![](elbow_plot.png)


Continuing on...

```R

n_pcs <- 10     # went with 10 PCs

# cluster and visualize
dat <- FindNeighbors(dat, dims = 1:n_pcs)
dat <- FindClusters(dat, resolution = 0.8)
dat <- RunUMAP(dat, dims = 1:n_pcs)
g1 <- DimPlot(dat, reduction = "umap")
ggsave(g1, file='umap_1.png', width=25, height=18, units='cm')

```

![](umap_1.png)


```R
# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
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

g2 <- DimPlot(dat, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') + guides(color=FALSE)

ggsave(g2, file='umap_2.png', width=25, height=18, units='cm')
```
![](umap_2.png)


```R
# guess tissue type
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
# db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
# tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = dat, scaled = TRUE, assay = "RNA")
```