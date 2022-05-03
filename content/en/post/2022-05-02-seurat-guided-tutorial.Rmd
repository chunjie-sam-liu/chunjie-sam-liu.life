---
title: "PlayGround - Seurat - Guided Clustering Tutorial"
author: "Chun-Jie Liu"
date: "2022-05-02"
---

The tutorial is from Seurat v4.0.6 10X genomics PBMC data, [here](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html). There are 2700 single cells that were sequenced on the Illumina NextSeq 500. Data is available [here](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).


We start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column)

## Load library

```{r}
# library(dplyr)
library(Seurat)
library(patchwork)
```

## Load data

```{r}
pbmc.data <- Seurat::Read10X(data.dir = "~/tmp/singlecell/filtered_gene_bc_matrices/hg19")
pbmc <- Seurat::CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
```
```{r}
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
```
```{r}
dense <- object.size(as.matrix(pbmc.data))
```

```{r}
sparse <- object.size(pbmc.data)
```
```{r}
dense
sparse
dense/sparse
```

## Standard pre-processing workflow

### QC and selecting cells for further analysis

1. The number of unique genes detected in each cell.
  - Low-quality cells or empty droplets will often have very few genes
  - Cell doublets or multiplets  may exhibit an aberrantly  high gene count
2. Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
3. The percetage of reads that map to the mitochondrial genome
  - Low-quality/dying cells often exhibit extensive mitochondrial contamination
  - We calculate mitochondrial QC metrics with `PercentageFeatureSet()` function, which calculates the percentage of count originating from a set of features
  - We use the set of all genes starting with `MT-` as a set of mitochondrial genes

```{r}
pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data)
```

- We filter cells that have unique feature count over 25000 or less than 200
- We filter cells that have >5% mitochondrial counts
```{r}
Seurat::VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```{r}
plot1 <- Seurat::FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- Seurat::FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```
```{r}
pbmc <- subset(pbmc, nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

### Normalizing the data

```{r}
pbmc
```

After removing unwanted cells from the dataset, next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" that normalize the feature expression measurements for each cell by the total expression, multiplies this by a scale factor, and log-transform the result.
```{r}
pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

### Identification of highly variable features

We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (they are highly expressed in some cells, and lowly expressed in others)

Our procedure in Seurat is described in detail here, and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures() function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

```{r}
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(Seurat::VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- Seurat::VariableFeaturePlot(pbmc)
plot2 <- Seurat::LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```
### Scaling data

Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.


- Shifts the expression of each gene, so that the mean expression across cells is 0
- Scales the expression of each gene, so that the variance across cells is 1
  - This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
- The results of this are stored in pbmc[["RNA"]]@scale.data

```{r}
#Scaling is an essential step in the Seurat workflow, but only on genes that will be used as input to PCA. Therefore, the default in ScaleData() is only to perform scaling on the previously identified variable features (2,000 by default). To do this, omit the features argument in the previous function call,
all.genes <- rownames(pbmc)
pbmc <- Seurat::ScaleData(pbmc, features = all.genes)
# Your PCA and clustering results will be unaffected. However, Seurat heatmaps (produced as shown below with DoHeatmap()) require genes in the heatmap to be scaled, to make sure highly-expressed genes don’t dominate the heatmap. To make sure we don’t leave any genes out of the heatmap later, we are scaling all genes in this tutorial.
```

## Perform linear dimensional reduction

Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.

```{r}
pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(object = pbmc))
```

```{r}
Seurat::VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
```


```{r}
Seurat::DimPlot(pbmc, reduction = "pca")
```

```{r}
Seurat::DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```

## Determine the ‘dimensionality’ of the dataset

```{r}
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- Seurat::JackStraw(pbmc, num.replicate = 100)
pbmc <- Seurat::ScoreJackStraw(pbmc, dims = 1:20)
```


```{r}
Seurat::JackStrawPlot(pbmc, dims = 1:15)
```


```{r}
Seurat::ElbowPlot(pbmc)
```


## Cluster cells


```{r}
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5)
```

## Run non-linear dimensional reduction (UMAP/tSNE)


```{r}
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10)
```
```{r}
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
Seurat::DimPlot(pbmc, reduction = "umap")
```


## Finding differentially expressed features (cluster biomarkers)

Seurat can help you find markers that define clusters via differential expression. By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. `FindAllMarkers()` automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

```{r}
# find all markers of cluster 2
cluster2.markers <- Seurat::FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
```
```{r}
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- Seurat::FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
```

```{r}
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- Seurat::FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

```
```{r}
pbmc.markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n = 2, order_by = avg_log2FC)
```
```{r}
cluster0.markers <- Seurat::FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

```{r}
Seurat::VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
```

```{r}
# you can plot raw counts as well
Seurat::VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
```

```{r}
Seurat::FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
```

```{r}
pbmc.markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 10, wt = avg_log2FC) -> top10
Seurat::DoHeatmap(pbmc, features = top10$gene) + Seurat::NoLegend()
```

```{r}
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
new.cluster.ids
```

```{r}
pbmc <- Seurat::RenameIdents(pbmc, new.cluster.ids)

```
```{r}
Seurat::DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + Seurat::NoLegend()
```

```{r}
sessionInfo()
```

