library(Seurat)
library(ggplot2)
library(stringr)

tabula_muris_counts <- read.table("/scratch/data-for_fast_access/pub-others/tabula_muris_180920/tabula_muris.umi.csv.gz",sep=",",header=TRUE,row.names=1)
# # tabula_muris_meta <- read.csv('/scratch/data-for_fast_access/pub-others/tabula_muris_180920/tabula_muris.cell_metadata.csv', row.names = 1)
# tabula_muris_meta <- read.csv('~/alegbe-project/tabula_muris_meta_with_height.csv', row.names = 1)
tabula_muris_meta <- read_csv('/scratch/data-for_fast_access/pub-others/tabula_muris_180920/tabula_muris.cell_metadata.csv')

tabula_muris_meta <- tabula_muris_meta %>%
  mutate(tissue_celltype = stringr::str_replace_all(tissue_celltype, pattern="_", " "))

# Processing CELLECT-pvals
EA <- read_tsv('/projects/alegbe/sc-genetics/out/CELLECT-LDSC-out/tabula-muris/tabula_muris.ESmu.continuous.benchmark190722__EA3_Lee2018.cell_type_results.txt') %>%
  separate(Name, c('benchmark', 'tissue_celltype'), '__') %>%
  select(tissue_celltype, Coefficient_P_value) %>%
  mutate(tissue_celltype = stringr::str_replace_all(tissue_celltype, pattern="_", " "))


tm_meta <- left_join(tabula_muris_meta, EA, "tissue_celltype") %>%
  mutate(log_pval = -log10(Coefficient_P_value))

tabula_muris_meta <- as.data.frame(tm_meta)
rownames(tabula_muris_meta) <- tabula_muris_meta$cell_id

tm <- CreateSeuratObject(counts = tabula_muris_counts, meta.data = tabula_muris_meta, project = "TabulaMuris")
# Only keep annotated cells
tm <- subset(tm, cells = names(which(!is.na(tm$tissue_celltype))))

tm

#Normalise
tm <- NormalizeData(tm, normalization.method = "LogNormalize", scale.factor = 10000)
tm <- FindVariableFeatures(tm)

# Regress out mitochondrial expression
tm[["percent.mt"]] <- PercentageFeatureSet(tm, pattern = "^mt-")
tm <- ScaleData(tm, vars.to.regress = "percent.mt")

# Make UMAP
tm <- RunPCA(tm, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)
tm <- RunUMAP(tm, dims = 1:75, min.dist = 0.75)


# Plot UMAP
p1 <- DimPlot(tm, reduction = "umap", pt.size = 0.1) + ggtitle(label = "UMAP")
p1

# Plot feature UMAP
p3 <- FeaturePlot(tm, features = c("log_pval"), reduction = 'umap', pt.size = 0.1) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
p3

### Save
file.out <- sprintf("figs/fig-EA-tabula-muris-umap.png")
ggsave(p3, filename=file.out, width=10, height=10)
