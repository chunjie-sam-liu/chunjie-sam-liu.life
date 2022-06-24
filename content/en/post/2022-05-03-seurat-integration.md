---
title: "PlayGround - Seurat - scRNA-seq integration"
author: "Chun-Jie Liu"
date: "2022-05-03"
---


# Introduction to scRNA-seq integration

The joint analysis of two or more single-cell datasets poses unique challenges. In particular, identifying cell populations that are present across multiple datasets can be problematic under standard workflows. Seurat v4 includes a set of methods to match (or ‘align’) shared cell populations across datasets. These methods first identify cross-dataset pairs of cells that are in a matched biological state (‘anchors’), can be used both to correct for technical differences between datasets (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions.

## Integration goals

The following tutorial is designed to give you an overview of the kinds of comparative analyses on complex cell types that are possible using the Seurat integration procedure. Here, we address a few key goals:

- Create an ‘integrated’ data assay for downstream analysis
- Identify cell types that are present in both datasets
- Obtain cell type markers that are conserved in both control and stimulated cells
- Compare the datasets to find cell-type specific responses to stimulation

## Setup the Seurat objects
```{r}
library(Seurat)
library(SeuratData)
library(patchwork)
```

```{r}
SeuratData::InstallData("ifnb")
```

```{r}
# load dataset
SeuratData::LoadData("ifnb")

```

```{r}
# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- Seurat::SplitObject(
  object = ifnb,
  split.by = "stim"
)
```

```{r}
ifnb.list %>% names()
```




```{r}
ifnb.list
```

```{r}
identical(ifnb.list$CTRL@assays$RNA@counts, ifnb.list$CTRL@assays$RNA@data)
```



```{r}
# normalize and identify variable features for each dataset independently

ifnb.list <- lapply(
  X = ifnb.list,
  FUN = function(x) {
    x <- Seurat::NormalizeData(x)

    x <- Seurat::FindVariableFeatures(
      x,
      selection.method="vst",
      nfeatures = 2000
    )

  }
)
```
```{r}
ifnb.list
```

```{r}
# select features that are repeatedly variable across datasets for integration
features <- Seurat::SelectIntegrationFeatures(
  object.list = ifnb.list
  )
```


```{r}
length(features)
```


```{r}
head(features)
length(features)
```

```{r}
Seurat::VariableFeatures(ifnb.list$CTRL) %>% head()

Seurat::VariableFeatures(ifnb.list$STIM) %>% head()
```

```{r}
head(features)
```

```{r}
ggvenn::ggvenn(
  data = list(
    "CTRL" = Seurat::VariableFeatures(ifnb.list$CTRL),
    "STIM" = Seurat::VariableFeatures(ifnb.list$STIM),
    "InteFeatures" = features
  ),
  show_percentage  = FALSE
)
```



## Perform integration

We then identify anchors using the `FindIntegrationAnchors()` function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with `IntegrateData()`.

```{r}
immune.anchors <- Seurat::FindIntegrationAnchors(
  object.list = ifnb.list,
  anchor.features = features
)
```

```{r}
immune.anchors
```


```{r}
slotNames(immune.anchors)
```

```{r}
immune.anchors@object.list
```
```{r}
immune.anchors@reference.cells
```
```{r}
immune.anchors@reference.objects
```

```{r}
immune.anchors@query.cells
```

```{r}
immune.anchors@anchors %>% head()
```

```{r}
immune.anchors@offsets
```



```{r}
immune.anchors@anchor.features %>% head()
```

```{r}
ggvenn::ggvenn(
  data = list(
    "CTRL" = Seurat::VariableFeatures(ifnb.list$CTRL),
    "STIM" = Seurat::VariableFeatures(ifnb.list$STIM),
    "InteFeatures" = features,
    "AnchorFeatures" = immune.anchors@anchor.features
  ),
  show_percentage  = FALSE
)
```


```{r}
immune.anchors@neighbors
```

## Perform integration

```{r}
immune.combined <- Seurat::IntegrateData(
  anchorset = immune.anchors
  )
```
```{r}
immune.combined
```


```{r}
immune.combined@assays
```


## Perform an integrated analysis
```{r}
Seurat::DefaultAssay(immune.combined) <- "integrated"
```
```{r}
immune.combined
```
```{r}
immune.combined@assays
```


```{r}
immune.combined@active.assay
```
```{r}
immune.combined@assays$RNA@counts[1:4,1:4]
```

```{r}
immune.combined@assays$RNA@data[1:4,1:4]
```

```{r}
identical(immune.combined@assays$RNA@counts, immune.combined@assays$RNA@data)
```



```{r}
immune.combined@assays$integrated@counts %>% dim()
```
```{r}
immune.combined@assays$integrated@data %>% dim()
```
```{r}
immune.combined@assays$RNA@data[features, ] %>% dim()
```

```{r}
immune.combined@assays$integrated@data[features, ] %>% dim()
```
```{r}
identical(colnames(immune.combined@assays$RNA@data), colnames(immune.combined@assays$integrated@data))
```

```{r}
identical(
  immune.combined@assays$RNA@data[features, ],
  immune.combined@assays$integrated@data[features, ]
  )
```
```{r}
immune.combined@active.assay
```
```{r}
immune.combined@meta.data %>% head()
```



```{r}
DefaultAssay(immune.combined) <- "integrated"
```

```{r}
# Run the standard workflow for visualizaiton and clustering
immune.combined <- Seurat::ScaleData(immune.combined, verbose = FALSE)
```
```{r}
immune.combined@assays$RNA@data %>% dim()
```

```{r}
immune.combined@assays$RNA@data[c("HBB", "HBA2", "HBA1", "CCL4"),1:4]
```
```{r}
immune.combined@assays$RNA@counts[c("HBB", "HBA2", "HBA1", "CCL4"),1:4]
```
```{r}
immune.combined@assays$RNA@scale.data[c("HBB", "HBA2", "HBA1", "CCL4"),1:4]
```


```{r}
immune.combined@assays$integrated@data[c("HBB", "HBA2", "HBA1", "CCL4"),1:4]
```

```{r}
immune.combined@assays$integrated@scale.data[c("HBB", "HBA2", "HBA1", "CCL4"),1:4]
```

```{r}
immune.combined <- RunPCA(
  immune.combined,
  npcs = 30,
  verbose = FALSE
)
```
```{r}
Seurat::VizDimLoadings(
  immune.combined,
  dims = 1:2,
  reduction = "pca"
)
```


```{r}
Seurat::DimPlot(immune.combined, reduction = "pca")
```
```{r}
Seurat::DimHeatmap(
  immune.combined,
  dims = 1,
  cells = 500,
  balanced = TRUE
)
```

```{r}
# immune.combined <- Seurat::JackStraw(immune.combined, num.replicate = 100)
```

```{r}
# immune.combined <- Seurat::ScoreJackStraw(immune.combined, dims = 1:20)
```

```{r}
# Seurat::JackStrawPlot(immune.combined, dims = 1:20)
```

```{r}
Seurat::ElbowPlot(immune.combined,ndims = 50)
```
```{r}
immune.combined
```


```{r}
immune.combined <- Seurat::RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
```


```{r}
immune.combined <- Seurat::FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
```


```{r}
immune.combined <- Seurat::FindClusters(immune.combined, resolution = 0.5)
```

```{r}
immune.combined@meta.data %>% head()
```



```{r}
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
```
```{r}
DimPlot(immune.combined, reduction = "umap", split.by = "stim")
```

## Identify conserved cell type markers

To identify canonical cell type marker genes that are conserved across conditions, we provide the FindConservedMarkers() function. This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. For example, we can calculated the genes that are conserved markers irrespective of stimulation condition in cluster 6 (NK cells).

```
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')

packages = c("cowplot_1.1.1", "ggplot2_3.3.5", "patchwork_1.1.1", "thp1.eccite.SeuratData_3.1.5", "stxBrain.SeuratData_0.1.1", "ssHippo.SeuratData_3.1.4", "pbmcsca.SeuratData_3.0.0", "pbmcMultiome.SeuratData_0.1.1", "pbmc3k.SeuratData_3.1.4", "panc8.SeuratData_3.0.2", "ifnb.SeuratData_3.1.0", "hcabm40k.SeuratData_3.0.0", "bmcite.SeuratData_0.3.0", "SeuratData_0.2.1", "SeuratObject_4.0.4", "Seurat_4.0.6", "systemfonts_1.0.2", "sn_2.0.0", "plyr_1.8.6", "igraph_1.2.11", "lazyeval_0.2.2", "splines_4.1.0", "listenv_0.8.0", "scattermore_0.7", "TH.data_1.0-10", "digest_0.6.29", "htmltools_0.5.2", "fansi_1.0.0", "magrittr_2.0.1", "memoise_2.0.0", "tensor_1.5", "cluster_2.1.2", "ROCR_1.0-11", "limma_3.48.0", "globals_0.14.0", "matrixStats_0.61.0", "sandwich_3.0-1", "pkgdown_1.6.1", "spatstat.sparse_2.1-0", "colorspace_2.0-2", "rappdirs_0.3.3", "ggrepel_0.9.1", "rbibutils_2.2", "textshaping_0.3.5", "xfun_0.25", "dplyr_1.0.7", "crayon_1.4.2", "jsonlite_1.7.2", "spatstat.data_2.1-2", "survival_3.2-11", "zoo_1.8-9", "glue_1.6.0", "polyclip_1.10-0", "gtable_0.3.0", "leiden_0.3.9", "future.apply_1.8.1", "BiocGenerics_0.38.0", "abind_1.4-5", "scales_1.1.1", "mvtnorm_1.1-2", "DBI_1.1.1", "miniUI_0.1.1.1", "Rcpp_1.0.7", "plotrix_3.8-1", "metap_1.4", "viridisLite_0.4.0", "xtable_1.8-4", "tmvnsim_1.0-2", "reticulate_1.22", "spatstat.core_2.3-2", "stats4_4.1.0", "htmlwidgets_1.5.4", "httr_1.4.2", "RColorBrewer_1.1-2", "TFisher_0.2.0", "ellipsis_0.3.2", "ica_1.0-2", "farver_2.1.0", "pkgconfig_2.0.3", "sass_0.4.0", "uwot_0.1.11", "deldir_1.0-6", "utf8_1.2.2", "tidyselect_1.1.1", "labeling_0.4.2", "rlang_0.4.12", "reshape2_1.4.4", "later_1.3.0", "munsell_0.5.0", "tools_4.1.0", "cachem_1.0.6", "cli_3.1.0", "generics_0.1.1", "mathjaxr_1.4-0", "ggridges_0.5.3", "evaluate_0.14", "stringr_1.4.0", "fastmap_1.1.0", "yaml_2.2.1", "ragg_1.1.3", "goftest_1.2-3", "knitr_1.33", "fs_1.5.2", "fitdistrplus_1.1-6", "purrr_0.3.4", "RANN_2.6.1", "pbapply_1.5-0", "future_1.23.0", "nlme_3.1-152", "mime_0.12", "formatR_1.11", "compiler_4.1.0", "plotly_4.10.0", "png_0.1-7", "spatstat.utils_2.3-0", "tibble_3.1.6", "bslib_0.3.1", "stringi_1.7.6", "highr_0.9", "desc_1.4.0", "RSpectra_0.16-0", "lattice_0.20-44", "Matrix_1.3-3", "multtest_2.48.0", "vctrs_0.3.8", "mutoss_0.1-12", "pillar_1.6.4", "lifecycle_1.0.1", "Rdpack_2.1.2", "spatstat.geom_2.3-1", "lmtest_0.9-39", "jquerylib_0.1.4", "RcppAnnoy_0.0.19", "data.table_1.14.2", "irlba_2.3.5", "httpuv_1.6.5", "R6_2.5.1", "promises_1.2.0.1", "KernSmooth_2.23-20", "gridExtra_2.3", "parallelly_1.30.0", "codetools_0.2-18", "MASS_7.3-54", "assertthat_0.2.1", "rprojroot_2.0.2", "withr_2.4.3", "mnormt_2.0.2", "sctransform_0.3.2", "multcomp_1.4-17", "mgcv_1.8-35", "parallel_4.1.0", "grid_4.1.0", "rpart_4.1-15", "tidyr_1.1.4", "rmarkdown_2.10", "Rtsne_0.15", "Biobase_2.52.0", "numDeriv_2016.8-1.1", "shiny_1.7.1")
```

## Identify conserved cell type markers


To identify canonical cell type marker genes that are conserved across conditions, we provide the FindConservedMarkers() function. This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. For example, we can calculated the genes that are conserved markers irrespective of stimulation condition in cluster 6 (NK cells).

```{r}
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
```


```{r}
nk.markers <- Seurat::FindConservedMarkers(
  immune.combined,
  ident.1 = 6,
  grouping.var = "stim",
  verbose = FALSE
)

head(nk.markers)
```

```{r}
Seurat::FeaturePlot(
  immune.combined,
  features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
    "CCL2", "PPBP"),
  min.cutoff = "q9"
)
```
```{r}
immune.combined@meta.data$seurat_clusters %>% table()
```
```{r}
immune.combined@active.ident %>% head()
```
```{r}
immune.combined$seurat_clusters
```


```{r}
immune.combined <- SetIdent(immune.combined, value = "seurat_clusters")
```

```{r}
immune.combined$seurat_clusters %>% table()
```


```{r}
immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
    `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
    `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
```
```{r}
immune.combined@active.ident %>% table()
```


```{r}
Seurat::DimPlot(immune.combined, label = TRUE)
```
```{r}
Idents(immune.combined) %>% head()
```

```{r}
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
    "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
    "CD4 Naive T", "CD4 Memory T"))
```

```{r}
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
```

```{r}
Seurat::DotPlot(
  immune.combined,
  features = markers.to.plot,
  cols = c("blue", "red"),
  dot.scale = 8,
  split.by = "orig.ident"
) +
  Seurat::RotatedAxis()
```

## Identify differential expressed genes across conditions

Now that we’ve aligned the stimulated and control cells, we can start to do comparative analyses and look at the differences induced by stimulation. One way to look broadly at these changes is to plot the average expression of both the stimulated and control cells and look for genes that are visual outliers on a scatter plot. Here, we take the average expression of both the stimulated and control naive T cells and CD14 monocyte populations and generate the scatter plots, highlighting genes that exhibit dramatic responses to interferon stimulation.
```{r}
immune.combined@active.ident %>% table()
```

```{r}
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
```

```{r}
t.cells
```


```{r}
t.cells@active.ident %>% head()
```
```{r}
t.cells@active.ident %>% head()
```
```{r}
Idents(t.cells) %>% head()
```


```{r}
Idents(t.cells) <- "stim"
```



```{r}
avg.t.cells <- as.data.frame(
  log1p(
    AverageExpression(t.cells, verbose = FALSE)$RNA
  )
)
avg.t.cells$gene <- rownames(avg.t.cells)
head(avg.t.cells)
```

```{r}
cd14.mono <- subset(
  immune.combined,
  idents = "CD14 Mono"
)
```
```{r}
Idents(cd14.mono) %>% table()
```
```{r}
cd14.mono@meta.data %>% head()
```


```{r}
AverageExpression(cd14.mono, verbose = FALSE)$RNA %>% head()
```

```{r}
log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA) %>% head()
```


```{r}
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)
```

```{r}
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
```
```{r}
head(avg.t.cells)
```

```{r}
theme_set(cowplot::theme_cowplot())
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2
```

As you can see, many of the same genes are upregulated in both of these cell types and likely represent a conserved interferon response pathway.

Because we are confident in having identified common cell types across condition, we can ask what genes change in different conditions for cells of the same type. First, we create a column in the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column. Then we use FindMarkers() to find the genes that are different between stimulated and control B cells. Notice that many of the top genes that show up here are the same as the ones we plotted earlier as core interferon response genes. Additionally, genes like CXCL10 which we saw were specific to monocyte and B cell interferon response show up as highly significant in this list as well.
```{r}
Idents(immune.combined) %>% head()
```

```{r}
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
```
```{r}
immune.combined$celltype.stim %>% table()
```


```{r}
immune.combined$celltype <- Idents(immune.combined)
```
```{r}
immune.combined$celltype %>% table()
```


```{r}
Idents(immune.combined) <- "celltype.stim"
```


```{r}
b.interferon.response <- FindMarkers(
  immune.combined,
  ident.1 = "B_STIM",
  ident.2 = "B_CTRL",
  verbose = FALSE
)
```


```{r}
head(b.interferon.response, n = 15)
```
Another useful way to visualize these changes in gene expression is with the split.by option to the FeaturePlot() or VlnPlot() function. This will display FeaturePlots of the list of given genes, split by a grouping variable (stimulation condition here). Genes such as CD3D and GNLY are canonical cell type markers (for T cells and NK/CD8 T cells) that are virtually unaffected by interferon stimulation and display similar gene expression patterns in the control and stimulated group. IFI6 and ISG15, on the other hand, are core interferon response genes and are upregulated accordingly in all cell types. Finally, CD14 and CXCL10 are genes that show a cell type specific interferon response. CD14 expression decreases after stimulation in CD14 monocytes, which could lead to misclassification in a supervised analysis framework, underscoring the value of integrated analysis. CXCL10 shows a distinct upregulation in monocytes and B cells after interferon stimulation but not in other cell types.

```{r}
FeaturePlot(
  immune.combined,
  features =  c("CD3D", "GNLY", "IFI6"),
  split.by = "stim",
  max.cutoff = 3,
  cols = c("grey", "red")
)
```

```{r}
plots <- VlnPlot(
  immune.combined,
  features = c("LYZ", "ISG15", "CXCL10"),
  split.by = "stim",
  group.by = "celltype",
  pt.size = 0,
  combine = FALSE
)
wrap_plots(plots = plots, ncol = 1)
```

## Performing integration on datasets normalized with SCTransform

In Hafemeister and Satija, 2019, we introduced an improved method for the normalization of scRNA-seq, based on regularized negative binomial regression. The method is named ‘sctransform’, and avoids some of the pitfalls of standard normalization workflows, including the addition of a pseudocount, and log-transformation. You can read more about sctransform in the manuscript or our SCTransform vignette.

Below, we demonstrate how to modify the Seurat integration workflow for datasets that have been normalized with the sctransform workflow. The commands are largely similar, with a few key differences:

Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
As discussed further in our SCTransform vignette, we typically use 3,000 or more features for analysis downstream of sctransform.
Run the PrepSCTIntegration() function prior to identifying anchors
When running FindIntegrationAnchors(), and IntegrateData(), set the normalization.method parameter to the value SCT.
When running sctransform-based workflows, including integration, do not run the ScaleData() function

```{r}
LoadData("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")
```

```{r}
ifnb.list$CTRL
```



```{r}
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
```


```{r}
ifnb.list$CTRL
```

```{r}
ifnb.list$CTRL@assays$RNA@counts %>% dim()
```

```{r}
ifnb.list$CTRL@assays$RNA@data %>% dim()
```
```{r}
identical(ifnb.list$CTRL@assays$RNA@counts, ifnb.list$CTRL@assays$RNA@data)
```

```{r}
identical(ifnb.list$CTRL@assays$SCT@counts, ifnb.list$CTRL@assays$SCT@data)
```

```{r}
identical(ifnb.list$CTRL@assays$SCT@counts, ifnb.list$CTRL@assays$RNA@counts)
```

```{r}
ifnb.list$CTRL@assays$RNA@counts[c("ISG15", "CXCL10", "LYZ"),1:50]

ifnb.list$CTRL@assays$SCT@counts[c("ISG15", "CXCL10", "LYZ"),1:50]

```

```{r}
ifnb.list$CTRL@assays$RNA@data[c("ISG15", "CXCL10", "LYZ"),1:50]

ifnb.list$CTRL@assays$SCT@data[c("ISG15", "CXCL10", "LYZ"),1:50]

ifnb.list$CTRL@assays$SCT@scale.data[c("ISG15", "CXCL10", "LYZ"), 1:50]

```


```{r}
ifnb.list
```





```{r}
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
```


```{r}
ggvenn::ggvenn(
  data = list(
    "CTRL" = ifnb.list$CTRL@assays$SCT@var.features,
    "STIM" = ifnb.list$STIM@assays$SCT@var.features,
    "Features" = features

  ),
  show_percentage  = F
)
```



```{r}
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
```

```{r}
ifnb.list
```

```{r}
immune.anchors <- FindIntegrationAnchors(
  object.list = ifnb.list,
  normalization.method = "SCT",
  anchor.features = features
)
```


```{r}
immune.combined.sct <- IntegrateData(
  anchorset = immune.anchors,
  normalization.method = "SCT"
)
```
```{r}
immune.combined.sct@assays$SCT@data[c("ISG15", "CXCL10", "LYZ"),1:50]

immune.combined.sct@assays$integrated@data[c("ISG15", "CXCL10", "LYZ"),1:50]
```

```{r}
immune.combined.sct@assays$SCT@scale.data[c("ISG15", "CXCL10", "LYZ"),1:5]


```
```{r}
immune.combined.sct@assays$integrated@scale.data[c("ISG15", "CXCL10", "LYZ"),1:5]
```

```{r}
immune.combined.sct
```


```{r}
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
```

```{r}
ElbowPlot(immune.combined.sct, ndims = 50)
```

```{r}
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
```


```{r}
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
    repel = TRUE)
p1 + p2
```

```{r}
sessionInfo()
```

