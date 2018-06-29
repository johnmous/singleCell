---
title: "Seurat - Guided Clustering Tutorial"
output:
  html_document:
    keep_md: true
    smart: false
    theme: united
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
---
***

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  cache.lazy = FALSE,
  tidy = TRUE
)
```
## In this Rmd all samples will be combined

### Setup the Seurat Object for Ovary

This is an attempt to run a [guided clustering tutorial](http://satijalab.org/seurat/pbmc3k_tutorial.html). Changes in various parameters and introduction of additional functions will be explored.

We start by reading in the data. All features in Seurat have been configured to work with sparse matrices which results in significant memory and speed savings for Drop-seq/inDrop/10x data. A list of cell cycle markers is also loaded


```{r init, message=FALSE}
library(plotly)
library(Seurat)
library(dplyr)
library(Matrix)

outputDir  = getwd()

# Also read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "/apthto/data/data/cellCycleMarkers/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S
# phase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

# Read the datasets one by one
# Read the datasets one by one
ovary.1.data <- Read10X(data.dir = "/path/to/data/sample1_B1_i12A/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.2.data <- Read10X(data.dir = "/path/to/data/sample2_B1_i12B/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.3.data <- Read10X(data.dir = "/path/to/data/sample3_B1_i12G/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.4.data <- Read10X(data.dir = "/path/to/data/sample4_B1_i12F/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.5.data <- Read10X(data.dir = "/path/to/data/sample5_B2_i10E/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.6a.data <- Read10X(data.dir = "/path/to/data/sample6a_B1_i12H/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.7.data <- Read10X(data.dir = "/path/to/data/sample7_B2_i10C/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.8a.data <- Read10X(data.dir = "/path/to/data/sample8a_B2_i10A/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.8b.data <- Read10X(data.dir = "/path/to/data/sample8b_B2_i10B/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.10.data <- Read10X(data.dir = "/path/to/data/sample10_B2_i10D/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.11.data <- Read10X(data.dir = "/path/to/data/sample11_B2_i10F/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.12.data <- Read10X(data.dir = "/path/to/data/sample12_B2_i10G/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.13.data <- Read10X(data.dir = "/path/to/data/sample13_B1_i12C/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.145.data <- Read10X(data.dir = "/path/to/data/sample145_B1_i12D/outs/filtered_gene_bc_matrices/GRCh38/")
ovary.C1.data <- Read10X(data.dir = "/path/to/data/sampleC1_B1_i12E/outs/filtered_gene_bc_matrices/GRCh38/")
```

```{r create.object, results='hide', message=FALSE}
# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells. Keep all cells with at
# least 100 detected genes
ovary.1 <- CreateSeuratObject(raw.data = ovary.1.data, min.cells = 3, min.genes = 100, project = "01")
ovary.2 <- CreateSeuratObject(raw.data = ovary.2.data, min.cells = 3, min.genes = 100, project = "02")
ovary.3 <- CreateSeuratObject(raw.data = ovary.3.data, min.cells = 3, min.genes = 100, project = "03")
ovary.4 <- CreateSeuratObject(raw.data = ovary.4.data, min.cells = 3, min.genes = 100, project = "04")
ovary.5 <- CreateSeuratObject(raw.data = ovary.5.data, min.cells = 3, min.genes = 100, project = "05")
ovary.6a <- CreateSeuratObject(raw.data = ovary.6a.data, min.cells = 3, min.genes = 100, project = "6a")
ovary.7 <- CreateSeuratObject(raw.data = ovary.7.data, min.cells = 3, min.genes = 100, project = "07")
ovary.8a <- CreateSeuratObject(raw.data = ovary.8a.data, min.cells = 3, min.genes = 100, project = "8a")
ovary.8b <- CreateSeuratObject(raw.data = ovary.8b.data, min.cells = 3, min.genes = 100, project = "8b")
ovary.10 <- CreateSeuratObject(raw.data = ovary.10.data, min.cells = 3, min.genes = 100, project = "10")
ovary.11 <- CreateSeuratObject(raw.data = ovary.11.data, min.cells = 3, min.genes = 100, project = "11")
ovary.12 <- CreateSeuratObject(raw.data = ovary.12.data, min.cells = 3, min.genes = 100, project = "12")
ovary.13 <- CreateSeuratObject(raw.data = ovary.13.data, min.cells = 3, min.genes = 100, project = "13")
ovary.145 <- CreateSeuratObject(raw.data = ovary.145.data, min.cells = 3, min.genes = 100, project = "145")
ovary.C1 <- CreateSeuratObject(raw.data = ovary.C1.data, min.cells = 3, min.genes = 100, project = "C1")

## Merge objects
ovary <- MergeSeurat(ovary.6a, ovary.8a, add.cell.id1 = "6a", add.cell.id2 = "8a", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.8b, add.cell.id2 = "8b", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.145, add.cell.id2 = "145", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.C1, add.cell.id2 = "C1", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.1, add.cell.id2 = "01", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.2, add.cell.id2 = "02", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.3, add.cell.id2 = "03", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.4, add.cell.id2 = "04", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.5, add.cell.id2 = "05", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.7, add.cell.id2 = "07", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.10, add.cell.id2 = "10", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.11, add.cell.id2 = "11", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.12, add.cell.id2 = "12", do.normalize = FALSE)
ovary <- MergeSeurat(ovary, ovary.13, add.cell.id2 = "13", do.normalize = FALSE)

table(ovary@meta.data$orig.ident)
```

***

### Standard pre-processing workflow
The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the creation of a Seurat object, the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable genes. In previous versions, we grouped many of these steps together in the `Setup` function, but in v2, we separate these steps into a clear and sequential workflow.


### Useful Definitions
* UMIs: Unique Molecular Identifiers (UMIs) are random oligonucleotide barcodes. Through a UMI, identical copies arising from distinct molecules can be distinguished from those arising through PCR-amplification of the same molecule.
* Multiplet: Lodsing and sequencing more than one cells (doublets, triplets)

### QC and selecting cells for further analysis
While the `CreateSeuratObject` imposes a basic minimum gene-cutoff, you may want to filter out cells at this stage based on technical or biological parameters. Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. In the example below, we visualize gene and molecule counts, plot their relationship, and exclude cells with a clear outlier number of genes detected as potential multiplets. Of course this is not a guaranteed method to exclude cell doublets, but we include this as an example of filtering user-defined outlier cells. We also filter cells based on the percentage of mitochondrial genes present.

```{r qc, results='hide', fig.height=20,fig.width=10,  }
# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = ovary@data), value = TRUE)
percent.mito <- Matrix::colSums(ovary@raw.data[mito.genes, ])/Matrix::colSums(ovary@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
ovary <- AddMetaData(object = ovary, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = ovary, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 1, point.size.use = 0.5)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = ovary, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = ovary, gene1 = "nUMI", gene2 = "nGene")


# We filter out cells that have unique gene counts over 2500 or less than
# 300 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
ovary <- FilterCells(object = ovary, subset.names = c("nGene", "percent.mito", "nUMI"), 
                     low.thresholds = c(200, -Inf, 300), high.thresholds = c(2500, 0.07, Inf))

```

```{r cell_number }
cat("### Number of cells after applying filtering:\n")
print(ovary@data@Dim[2])
table(ovary@meta.data$orig.ident)
```

***


***
### Normalizing the data

After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" that normalizes the gene expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
### This  is the c++ code implementing the Normalization calculations
`log1p(double(it.value()) / colSums[k] * scale_factor);`
`log1p( double arg )` Computes the natural (base e) logarithm of *1+arg*. This function is more precise than the expression `std::log(1+arg)` if arg is close to zero.
`colSums` is the sum of raw counts for one cell
The per cell sums of normalized counts are not identical between cells. 

```{r normalize, results='hide'}
ovary <- NormalizeData(object = ovary, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

```


## Dissociation genes on normalized data
```{r removeDissociatioAffecgedCellsNormalized, fig.height=6, fig.width=6}

## Gene list is from Van Den Brink et al. (2017) suplement. data R code

genesChrom <- c("Actg1__chr11","Ankrd1__chr19","Arid5a__chr1","Atf3__chr1","Atf4__chr15","Bag3__chr7","Bhlhe40__chr6",
"Brd2__chr17","Btg1__chr10","Btg2__chr1","Ccnl1__chr3","Ccrn4l__chr3","Cebpb__chr2","Cebpd__chr16",
"Cebpg__chr7","Csrnp1__chr9","Cxcl1__chr5","Cyr61__chr3","Dcn__chr10","Ddx3x__chrX","Ddx5__chr11",
"Des__chr1","Dnaja1__chr4","Dnajb1__chr8","Dnajb4__chr3","Dusp1__chr17","Dusp8__chr7",
"Egr1__chr18","Egr2__chr10","Eif1__chr11","Eif5__chr12","Erf__chr7","Errfi1__chr4","Fam132b__chr1",
"Fos__chr12","Fosb__chr7","Fosl2__chr5","Gadd45a__chr6","Gcc1__chr6","Gem__chr4","H3f3b__chr11",
"Hipk3__chr2","Hsp90aa1__chr12","Hsp90ab1__chr17","Hspa1a__chr17","Hspa1b__chr17","Hspa5__chr2",
"Hspa8__chr9","Hspb1__chr5","Hsph1__chr5","Id3__chr4","Idi1__chr13","Ier2__chr8","Ier3__chr17",
"Ifrd1__chr12","Il6__chr5","Irf1__chr11","Irf8__chr8","Itpkc__chr7","Jun__chr4","Junb__chr8",
"Jund__chr8","Klf2__chr8","Klf4__chr4","Klf6__chr13","Klf9__chr19","Litaf__chr16","Lmna__chr3",
"Maff__chr15","Mafk__chr5","Mcl1__chr3","Midn__chr10","Mir22hg__chr11","Mt1__chr8","Mt2__chr8",
"Myadm__chr7","Myc__chr15","Myd88__chr9","Nckap5l__chr15","Ncoa7__chr10","Nfkbia__chr12","Nfkbiz__chr16",
"Nop58__chr1","Nppc__chr1","Nr4a1__chr15","Odc1__chr12","Osgin1__chr8","Oxnad1__chr14","Pcf11__chr7",
"Pde4b__chr4","Per1__chr11","Phlda1__chr10","Pnp__chr14","Pnrc1__chr4","Ppp1cc__chr5","Ppp1r15a__chr7",
"Pxdc1__chr13","Rap1b__chr10","Rassf1__chr9","Rhob__chr12","Rhoh__chr5","Ripk1__chr13","Sat1__chrX",
"Sbno2__chr10","Sdc4__chr2","Serpine1__chr5","Skil__chr3","Slc10a6__chr5","Slc38a2__chr15",
"Slc41a1__chr1","Socs3__chr11","Sqstm1__chr11","Srf__chr17","Srsf5__chr12","Srsf7__chr17",
"Stat3__chr11","Tagln2__chr1","Tiparp__chr3","Tnfaip3__chr10","Tnfaip6__chr2","Tpm3__chr3",
"Tppp3__chr8","Tra2a__chr6","Tra2b__chr16","Trib1__chr15","Tubb4b__chr2","Tubb6__chr18",
"Ubc__chr5","Usp2__chr9","Wac__chr18","Zc3h12a__chr4","Zfand5__chr19","Zfp36__chr7","Zfp36l1__chr12",
"Zfp36l2__chr17","Zyx__chr6","Gadd45g__chr13","Hspe1__chr1","Ier5__chr1","Kcne4__chr1")

genes <- sapply(genesChrom, function(x){
  toupper( strsplit(x, "__")[[1]][1])
})

Data <- as.data.frame(as.matrix(ovary@data))
cat("All genes:\n")
print(unname(genes))
write.table(genes, paste0(outputDir, "/mouseDissocGenes.tsv"), sep ="\t", quote=FALSE, row.names = FALSE)

## Remove mouse only genes and put the corresponding human
genes <- genes[!genes %in% c("CCRN4L", "MT1", "MT2")]
genes <- c(genes, "NOCT", "MT1A", "MT2A")
cat("Genes from mouse we miss in human:\n")
unname(genes[!genes %in% row.names(Data)])

## Calculate the percentage of UMIs maping on dissociation genes
totalSum <- colSums(ovary@data)
selection <- Data[genes, ]
selection[is.na(selection)] <- 0
dissociationSums <- colSums(selection)  
countSums <- merge(totalSum, dissociationSums, by="row.names", all=TRUE, sort= FALSE)
rownames(countSums) <- countSums$Row.names
countSums <- countSums[-1]
colnames(countSums) <- c("totalCount", "dissociationCounts")
countSums$percentage <- countSums$dissociationCounts/countSums$totalCount
## Save in meta.data of object
ovary@meta.data$percent.dissoc <- countSums$percentage

plotDissociation <- function(sample){
  percentages <- ovary@meta.data$percent.dissoc[ovary@meta.data$orig.ident == sample]
  # Make a histogram showing the distribution of the metric "Percentage of dissociation-affected reads per cell":
  hist(percentages, breaks = 100, col = "lightgrey", main = paste("Expression dissociation-affected genes for sample: ", sample), 
       xlab = "Ratio of dissociation-affected genes to total gene count", ylab = "Number of cells", xlim = c(0, 0.20))
}

for (sample in levels(ovary@ident)){
  plotDissociation(sample)
}

## Draw histogram for all samples
percentages <- ovary@meta.data$percent.dissoc
hist(percentages, breaks = 100, col = "lightgrey", main = paste("Expression dissociation-affected genes for sample: ", sample), 
       xlab = "Ratio of dissociation-affected genes to total gene count", ylab = "Number of cells", xlim = c(0, 0.20))

```

## Keep cells with dissociation percentages below the threshold of 0.06 or 6%
```{r filterDissoc}
ovary <-FilterCells(ovary, subset.names = "percent.dissoc", low.thresholds = -Inf, high.thresholds = 0.06)
cat("### Number of cells after applying Dissociation-indiced genes filtering:\n")
print(ovary@data@Dim[2])
table(ovary@meta.data$orig.ident)
```

### Detection of variable genes across the single cells

Seurat calculates highly variable genes and focuses on these for downstream analysis. **`FindVariableGenes`** calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression. This function is unchanged from (Macosko *et al*.), but new methods for variable gene expression identification are coming soon. We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. The parameters here identify ~2,000 variable genes, and represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules.

*Exact parameter settings may vary empirically from dataset to dataset, and based on visual inspection of the plot. Setting the y.cutoff parameter to 2 identifies genes that are more than two standard deviations away from the average dispersion within a bin. The default X-axis function is the mean expression level, and for Y-axis it is the log(Variance/mean). log is natural logarithm (e=2.71). All mean/variance calculations are not performed in log-space, but the results are reported in log-space - see relevant functions for exact details.*

```{r var_genes, fig.height=7, fig.width=11, results='hide'}
ovary <- FindVariableGenes(object = ovary, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = 0.25)

## Write the genes that have an average log expression>3 to a table and save it.
topAveExpr = ovary@hvg.info[ovary@hvg.info[,1]>3, ]
topAveExprPath = paste0(outputDir, "/topAveExpr.tsv")
write.table(x = topAveExpr, file = topAveExprPath, sep = "\t")
```

```{r len_var_genes}
length(x = ovary@var.genes)
```

***

### Scaling the data and removing unwanted sources of variation

Your single cell dataset likely contains 'uninteresting' sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in Buettner *et al*, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on user-defined variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and  are used for dimensionality reduction and clustering. 

We can regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data),  the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a 'cell-cycle' score (see example [HERE]) and regress this out as well. In this simple example here for post-mitotic blood cells, we regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content.  

Seurat v2.0 implements this regression as part of the data scaling process. Therefore, the `RegressOut` function has been deprecated, and replaced with the vars.to.regress argument in `ScaleData`.

*Scales and centers genes in the dataset. If variables are provided in vars.to.regress, they are individually regressed against each gene, and the resulting residuals are then scaled and centered.*


```{r regress, fig.height=7, fig.width=11, results='hide'}
ovary <- ScaleData(object = ovary, vars.to.regress = c("nUMI", "percent.mito"))
```
***

```{r pca_viz_cell_cycle}
## Cell cycling scoring
ovary <- CellCycleScoring(object = ovary, s.genes = s.genes, g2m.genes = g2m.genes, 
                           set.ident = TRUE)
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells
# separate entirely by phase
ovary <- RunPCA(object = ovary, pc.genes = c(s.genes, g2m.genes), do.print = FALSE)
PCAPlot(object = ovary)
```

### Regress out cell cycle scores during data scaling
We now attempt to subtract (‘regress out’) this source of heterogeneity from the data. For users of Seurat v1.4, this was implemented in RegressOut. However, as the results of this procedure are stored in object@scale.data (therefore overwriting the output of ScaleData), we now merge this functionality into the ScaleData function itself.

For each gene, Seurat models the relationship between gene expression and the S and G2M cell cycle scores. The scaled residuals of this model represent a ‘corrected’ expression matrix, that can be used downstream for dimensional reduction.

```{r regress out cell cycle }
ovary <- ScaleData(object = ovary, vars.to.regress = c("S.Score", "G2M.Score"), 
    display.progress = FALSE)
```


### Perform linear dimensional reduction

Next we perform PCA on the scaled data. By default, the genes in `object@var.genes` are used as input, but can be defined using pc.genes. We have typically found that running dimensionality reduction on highly variable genes can improve performance. However, with UMI data - particularly after regressing out technical variables, we often see that PCA returns similar (albeit slower) results when run on much larger subsets of genes, including the whole transcriptome.

```{r pca}
ovary <- RunPCA(object = ovary, pc.genes = ovary@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
```

Seurat provides several useful ways of visualizing both cells and genes that define the PCA, including `PrintPCA`, `VizPCA`, `PCAPlot`, and `PCHeatmap`


```{r pca_viz}
# Examine and visualize PCA results a few different ways
## Print PCAs
PrintPCA(object = ovary, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
## Save the gene loading in a 
write.table( ovary@dr$pca@gene.loadings, file = paste0(outputDir, "/geneLoadings"), sep = "\t", quote = FALSE)

## 
VizPCA(object = ovary, pcs.use = 1:2)
PCAPlot(object = ovary, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation 
# with the calculated components. Though we don't use this further here, it can be used to identify markers that 
# are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. 
# The results of the projected PCA can be explored by setting use.full=T in the functions above
ovary <- ProjectPCA(object = ovary, do.print = FALSE)
```

In particular `PCHeatmap` allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and genes are ordered according to their PCA scores. Setting cells.use to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated gene sets.

```{r single-heatmap, warning=FALSE}
PCHeatmap(object = ovary, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
```

```{r multi-heatmap, fig.height=12, fig.width=9, warning=FALSE}
PCHeatmap(object = ovary, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
```

***

### Determine statistically significant principal components

To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metagene' that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

In Macosko *et al*, we implemented a resampling test inspired by the jackStraw procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of gene scores, and repeat this procedure. We identify 'significant' PCs as those who have a strong enrichment of low p-value genes.

```{r jackstraw, fig.height=6, fig.width=10, warning=FALSE, cache=TRUE, results = 'hide'}
# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in PCElbowPlot() can be used to reduce computation time
ovary <- JackStraw(object = ovary, num.replicate = 100 )
```

The `JackStrawPlot` function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 'Significant' PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line). In this case it appears that PCs 1-10 are significant.

```{r jsplots, fig.height=6, fig.width=10, warning=FALSE}
JackStrawPlot(object = ovary, PCs = 1:12)
```

A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. This can be done with `PCElbowPlot`. In this example, it looks like the elbow would fall around PC 9.

```{r elbow_plot, fig.height=6, fig.width=10, warning=FALSE}
PCElbowPlot(object = ovary)
```

PC selection -- identifying the true dimensionality of a dataset -- is an important step for Seurat, but can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-10 as a cutoff. We followed the jackStraw  here, admittedly buoyed by seeing the PCHeatmap returning interpretable signals (including canonical dendritic cell markers) throughout these PCs. Though the results are only subtly affected by small shifts in this cutoff (you can test below), we strongly suggest always explore the PCs they choose to include downstream.

***

### Cluster the cells

Seurat now includes an graph-based clustering approach compared to (Macosko *et al*.). Importantly, the *distance metric* which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data [[SNN-Cliq, Xu and Su, Bioinformatics, 2015]](http://bioinformatics.oxfordjournals.org/content/early/2015/02/10/bioinformatics.btv088.abstract) and CyTOF data [[PhenoGraph, Levine *et al*., Cell, 2015]](http://www.ncbi.nlm.nih.gov/pubmed/26095251). Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar gene expression patterns, and then attempt to partition this graph into highly interconnected 'quasi-cliques' or 'communities'. As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard distance). To cluster the cells, we apply modularity optimization techniques [[SLM, Blondel *et al*., Journal of Statistical Mechanics]](http://dx.doi.org/10.1088/1742-5468/2008/10/P10008), to iteratively group cells together, with the goal of optimizing the standard modularity function.

The `FindClusters` function implements the procedure, and contains a resolution parameter that sets the 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters are saved in the `object@ident` slot.


```{r cluster, fig.height=5, fig.width=7}

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
# but with a different resolution value (see docs for full details)
ovary <- FindClusters(object = ovary, reduction.type = "pca", dims.use = 1:15, 
                     resolution = 0.9, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
```

A useful feature in Seurat v2.0 is the ability to recall the parameters that were used in the latest function calls for commonly used functions. For FindClusters, we provide the function `PrintFindClustersParams` to print a nicely formatted formatted summary of the parameters that were chosen. 

```{r cluster.params}
PrintFindClustersParams(object = ovary)
# While we do provide function-specific printing functions, the more general function to 
# print calculation parameters is PrintCalcParams(). 
```

***

### Run Non-linear dimensional reduction (tSNE)

Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets. While we no longer advise clustering directly on tSNE components, cells within the graph-based clusters determined above should co-localize on the tSNE plot. This is because the tSNE aims to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. As input to the tSNE, we suggest using the same PCs as input to the clustering analysis, although computing the tSNE based on scaled gene expression is also supported using the genes.use argument.


## Draw t-SNE and PCA  
```{r clusterPCs1-15, fig.height=8, fig.width=9}

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
# but with a different resolution value (see docs for full details)

## Plot tSNEs and PCAs 
ovary <- RunTSNE(object = ovary, dims.use = 1:15, do.fast = TRUE )
pdfPath <- paste0(outputDir, "/pdfPlots/")
dir.create(pdfPath)
TSNEPlot(object = ovary, do.label = T)
ggsave(paste0(pdfPath, "tSNE.pdf"), width = 10, height = 7)
TSNEPlot(object = ovary, do.label = T, group.by = "orig.ident")
ggsave(paste0(pdfPath, "tSNEOrigIdent.pdf"), width = 10, height = 7)
PCAPlot(object = ovary, dim.1 = 1, dim.2 = 2, pt.size=1)
ggsave(paste0(pdfPath, "PCA.pdf"), width = 10, height = 7)
PCAPlot(object = ovary, dim.1 = 1, dim.2 = 2, pt.size=1, group.by = 'orig.ident')
ggsave(paste0(pdfPath, "PCAOrigIdent.pdf"), width = 10, height = 7)

## 

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
ovary.markers <- FindAllMarkers(object = ovary, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
topMarkers <- ovary.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
topGenesPath = paste0(outputDir, "/topMarkerGenes.tsv")
write.table(x = topMarkers, file = topGenesPath, sep = "\t", quote = FALSE, row.names = FALSE)

## Write ovary to rds
ovaryRds = paste0(outputDir, "/ovary.rds")
saveRDS(ovary, file =  ovaryRds)
```


### As an extra, make plot interactive
```{r interactivetSNEs, fig.height=8, fig.width=9}
p <- TSNEPlot(object = ovary, do.label = T, do.return = TRUE) #, do.hover = TRUE, do.identify = TRUE)
ggplotly(p)
p <- TSNEPlot(object = ovary, do.label = T, group.by = "orig.ident", do.return = TRUE) #, do.hover = TRUE, do.identify = TRUE)
ggplotly(p)
p <- PCAPlot(object = ovary, dim.1 = 1, dim.2 = 2, do.return = TRUE) # , do.hover = TRUE, do.identify = TRUE, do.label=TRUE, no.legend = FALSE)
ggplotly(p)
p <- PCAPlot(object = ovary, dim.1 = 1, dim.2 = 2, group.by = 'orig.ident', do.return = TRUE) #, do.hover = TRUE, do.identify = TRUE, no.legend = FALSE)
ggplotly(p)
```

