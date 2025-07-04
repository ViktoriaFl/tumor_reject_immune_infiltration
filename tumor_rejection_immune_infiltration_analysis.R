# SETUP ------------------------------------------------------------------------

# load libraries
library(stringr)
library(Seurat)
library(ggpubr)
library(viridis)
library(PICtR)
library(clustree)
library(ggrastr)
library(rstatix)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 2e9)

# directories
input.dir <- "./flowjo_data/"
data.dir <- "./data/"
plot.dir <- "./plots/"

# create directory structure if not already present 
for (d in c(input.dir, data.dir, plot.dir)) {
  ifelse(!dir.exists(d), dir.create(d), FALSE)
}

# LOAD DATA AND SKETCH ---------------------------------------------------------
# find channel value csv files 
files <- dir(paste0(input.dir), pattern="*.csv", recursive = T)

message("reading data...")
# read channel values and save in list
file_list <- lapply(files, function(csv) {
  # read csv
  tab <- read.table(paste0(input.dir, csv), header = T, sep = ",")
  # extract meta data 
  tab$sample <- csv
  return(tab)
})

# combine into one dataframe
data <- bind_rows(file_list)

# format 
data <- data %>% 
  mutate(sample = str_replace(sample, "\\s", "_")) %>% 
  mutate(ID = str_extract(sample, "samples_.._(.*)_202", group = 1)) %>% 
  mutate(group = case_when(
    ID %in% c("3203", "3745", "3746", "3747", "3948") ~ 1, 
    ID %in% c("3153", "3204", "3318", "3321", "3455") ~ 2, 
    ID %in% c("3152", "3205", "3206", "32x9", "3320") ~ 3)) %>% 
    mutate(group_name = case_when(
      group == 1 ~ "Dox", 
      group == 2 ~ "Dox_removed", 
      group == 3 ~ "Dox_reintroduced")) %>% 
  mutate(tissue = "tumor") %>% 
  mutate(date = str_extract(sample, "_(202._.._..)_", group = 1))

# remove sample 3206
# large batch effects probably due to the continuous clogging during sample prep 
data_qc <- data %>% dplyr::filter(ID != "3206")

# CREATE SEURAT OBJECT --------------------------------------------------------
data_qc$ratio <- scales::rescale(as.numeric(data_qc$FSC.A / data_qc$FSC.H), to = c(0, 1023))
obj <- CreateSeuratObject(counts = t(data_qc[,c(1,2,4,8, 11:44, 53)]), meta.data = data_qc, assay = "FACS")
obj@assays$FACS$data <- obj@assays$FACS$counts

# SKETCH PER GROUP ------------------------------------------------------------
split <- obj
split[["FACS"]] <- split(split[["FACS"]], f = split$group) # split layers 
DefaultAssay(split) <- "FACS"

# sketch 50 000 cells per group 
split <- FindVariableFeatures(split)
split <- SketchData(split, assay = "FACS", ncells = 50000, method = "LeverageScore", 
                    sketched.assay = "sketch", seed = 42)
DefaultAssay(split) <- "sketch"

# join layers 
split <- JoinLayers(split)

obj <- split
rm(split)

# remove uPAR due to unreliable signal
VariableFeatures(obj) <- VariableFeatures(obj) %>% str_subset("UPAR", negate = T)
obj <- ScaleData(obj) 
obj <- RunPCA(obj, npcs = 38, approx = F) 
obj <- FindNeighbors(obj, dims = 1:38)
obj <- FindClusters(obj, resolution = c(0.3, 0.5, 1, 2))
obj <- RunUMAP(obj, dims = 1:38, return.model = TRUE, seed.use = 42)

# save 
saveRDS(obj, file = paste0(data.dir, "obj_per_group_sketched_non_projected.rds"))

# QUALITY CONTROL -------------------------------------------------------------
ifelse(!dir.exists(paste0(plot.dir, "01_doublets/")), dir.create(paste0(plot.dir, "01_doublets/")), FALSE)
ifelse(!dir.exists(paste0(plot.dir, "02_qc/")), dir.create(paste0(plot.dir, "02_qc/")), FALSE)
ifelse(!dir.exists(paste0(plot.dir, "03_annotation/")), dir.create(paste0(plot.dir, "03_annotation/")), FALSE)

# per sample
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_ID_sketched_per_group.pdf"), height = 7, width = 9, onefile = T)
DimPlot(obj, group.by = "ID") 
dev.off() 

# per group
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_group_sketched_per_group.pdf"), height = 7, width = 9, onefile = T)
DimPlot(obj, group.by = "group") 
dev.off()

# per acquisition date
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_exp_date_sketched_per_group.pdf"), height = 7, width = 9, onefile = T)
DimPlot(obj, group.by = "date") 
dev.off()  # still some batch effects, but less pronounced without sample 3206, check features 

# feature plots
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group.pdf"), height = 7, width = 14, onefile = T)
FeaturePlot(obj, features = Features(obj)[c(1:12)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj, features = Features(obj)[13:24], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj, features = Features(obj)[c(25:36)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj, features = Features(obj)[c(37:39)], slot = "counts") & scale_color_viridis(option = 'magma')
dev.off() # predominantly CD73

# check area/height ratio distribution and cutoff
cutoff <- calculateThreshold(obj$ratio, method = "otsu")
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_ratio_histogram.pdf"), height = 6, width = 8, onefile = T)
hist(obj$ratio, breaks = 2000)
abline(v=cutoff)
dev.off() # check independently for CD45+ and CD45- cells

# CD45+
obj@meta.data %>% dplyr::filter(CD45 > 250) %>% ggplot(aes(x = ratio)) + 
  geom_histogram(bins = 1000) + 
  geom_vline(xintercept = calculateThreshold((obj@meta.data %>% dplyr::filter(CD45 > 250))$ratio))

# CD45-
obj@meta.data %>% dplyr::filter(CD45 < 250) %>% ggplot(aes(x = ratio)) + 
  geom_histogram(bins = 1000) + 
  geom_vline(xintercept = calculateThreshold((obj@meta.data %>% dplyr::filter(CD45 < 250))$ratio))

# default method is not ideal, try other thresholding methods 
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_ratio_histogram_evaluate_methods.pdf"), height = 6, width = 8, onefile = T)
for (m in c("otsu", "triangle", "kmeans", "IJDefault", "Huang", "Huang2", 
            "Intermodes", "IsoData", "Li", "Mean", "MinErrorI", 
            "Minimum", "Moments", "Percentile", "RenyiEntropy", "Shanbhag")) {
  obj@meta.data[[m]] <- if_else(obj$ratio >= calculateThreshold(obj$ratio, method = m), "high", "low")
  print(ggplot(obj@meta.data, aes(x = ratio, fill = .data[[m]])) +
    geom_histogram(bins = 2000) +
    labs(title = m))
  print(ggplot(obj@meta.data, aes(x = FSC.A, y = FSC.H, color = .data[[m]])) +
    geom_point() +
    labs(title = m))
}
dev.off()

p1 <- ggplot(obj@meta.data, aes(x = FSC.A, y = FSC.H, color = otsu)) + geom_point() + labs(title = "otsu")
p2 <- ggplot(obj@meta.data, aes(x = FSC.A, y = FSC.H, color = triangle)) + geom_point() + labs(title = "triangle")
p3 <- ggplot(obj@meta.data, aes(x = FSC.A, y = FSC.H, color = RenyiEntropy)) + geom_point() + labs(title = "RenyiEntropy")
p4 <- ggplot(obj@meta.data, aes(x = FSC.A, y = FSC.H, color = Moments)) + geom_point() + labs(title = "Moments")
ggarrange(p1, p2, p3, p4) 

# move forward with triangle
obj$ratio_anno <- obj$triangle
ratio_cluster_plot(obj, "sketch_snn_res.0.5", assay = "sketch") 

# ANNOTATION ------------------------------------------------------------------
# clustering resolutions
resolution <- c("sketch_snn_res.0.3", "sketch_snn_res.0.5", "sketch_snn_res.1", "sketch_snn_res.2")
ifelse(!dir.exists(paste0(plot.dir, "03_annotation/")), dir.create(paste0(plot.dir, "03_annotation/")), FALSE)

# dimplots grouped by cluster at each resolution 
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_dimplot_clusters_sketched_per_group.pdf"), height = 7, width = 9, onefile = T)
for (res in resolution) {
  print(DimPlot(obj, reduction = 'umap', group.by = res, label = T) + 
          labs(title = paste0(res)))
}
dev.off()

# dimplots grouped by cluster at each resolution 
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_dimplot_clusters.pdf"), height = 7, width = 9, onefile = T)
  print(DimPlot(obj, reduction = 'umap', group.by = "sub.cluster", label = T) + labs(title = "clusters"))
dev.off()

# clustree 
subset <- obj@meta.data %>% dplyr::filter(!is.na(seurat_clusters))
sketched <- subset(obj, cells = rownames(subset))

pdf(paste0(plot.dir, "01_doublets/", Sys.Date(), "_clustree_all_events_sketched_per_group.pdf"), height = 7, width = 14)
clustree(sketched, prefix = "sketch_snn_res.", exprs = 'data', assay = 'sketch')
dev.off()

# check resolution = 0.5 
Idents(obj) <- "sketch_snn_res.0.5"
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_dimplot_highlighted_clusters_res0.5_sketched_per_group.pdf"), height = 7, width = 9, onefile = T)
for (i in 0:(length(unique(obj$sketch_snn_res.0.5))-2)) {
  print(DimPlot(obj, cells.highlight = WhichCells(obj, idents = as.character(i))))
}
dev.off()

############################
##                        ##
##        choose          ##
##   sketch_snn_res.0.5   ##
##                        ##
############################

# ratio cluster plots and preliminary selection of doublet clusters
pdf(paste0(plot.dir, "01_doublets/", Sys.Date(), "_ratio_cluster_plots_all_events.pdf"), height = 7, width = 14, onefile = T)
for (res in resolution) {
  obj <- select_dbt(obj, clusters = res, selected_clusters = res)
  
  # plots
  print(ratio_cluster_plot(obj, clusters = res))
}
dev.off()

# subcluster 17, 19, 22 and 29 
Idents(obj) <- "sketch_snn_res.0.5"
obj <- FindSubCluster(obj, cluster = "17", resolution = 0.1, subcluster.name = "sub.cluster", graph.name = "sketch_snn")
Idents(obj) <- "sub.cluster"
obj <- FindSubCluster(obj, cluster = "19", resolution = 0.05, subcluster.name = "sub.cluster", graph.name = "sketch_snn")
Idents(obj) <- "sub.cluster"
obj <- FindSubCluster(obj, cluster = "29", resolution = 0.1, subcluster.name = "sub.cluster", graph.name = "sketch_snn")
Idents(obj) <- "sub.cluster"
obj <- FindSubCluster(obj, cluster = "22", resolution = 0.1, subcluster.name = "sub.cluster", graph.name = "sketch_snn") # also subcluster 22 to check Ly6G+ populations
Idents(obj) <- "sub.cluster"

# plot
DimPlot(obj, group.by = "sub.cluster", label =T)

# RidgePlot
Idents(obj) <- "sketch_snn_res.0.5"
features <- Features(obj)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_ridgeplots_resolution_0.5", ".pdf"), height = 8, width = 5, onefile = T)
for (f in features) {
  print(RidgePlot(obj, features = f, slot = 'data'))
}
dev.off()

# RidgePlot subclustered
Idents(obj) <- "sub.cluster"
DefaultAssay(obj) <- "sketch"
features <- Features(obj)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_ridgeplots_resolution_0.5_subclustered_all", ".pdf"), height = 8, width = 5, onefile = T)
for (f in features) {
  print(RidgePlot(obj, features = f, slot = 'data', assay = "sketch"))
}
dev.off()

# plot features 
features <- Features(obj) %>% str_replace("-", "_")
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_features_overall_landscape_rastered.pdf"), width = 12, height = 9)
plots = vector('list', length(features))
for(i in seq_along(features)){
  plots[[i]] = ggplot(obj@meta.data, aes(x = umap_1, y = umap_2, color = .data[[features[i]]])) + 
    ggrastr::geom_point_rast(size = 2.5, color = "black", raster.dpi = 700) + 
    ggrastr::geom_point_rast(aes(color = .data[[features[i]]]), size = 1, raster.dpi = 700) + 
    theme_void() +
    labs(title = features[i]) + 
    scale_color_viridis(option = "magma") +
    theme(legend.position = "none")
}
print(ggarrange(plotlist = plots))
dev.off()

# project chosen clustering to all cells
obj@meta.data <- obj@meta.data %>% 
  mutate(sub.cluster = str_replace(sub.cluster, "_", "."))

# join FACS layers
DefaultAssay(obj) <- "FACS"
obj <- JoinLayers(obj)

obj <- predict_data(obj = obj,
                    data_query = obj, 
                    ref_clusters = "sub.cluster", 
                    FSC.A="FSC.A",
                    FSC.H="FSC.H",
                    pred_name="pred_snn_res.0.5",
                    chunk_size=1000000, 
                    assay_ref = "sketch",
                    assay_query = "FACS")


# save projected obj
saveRDS(obj, file = paste0(data.dir, "obj_per_group_sketched_projected.rds"))

# annotation 
obj@meta.data <- mutate(obj@meta.data, celltype = case_when(
  pred_snn_res.0.5 == "0"~"M1 Macrophages",
  pred_snn_res.0.5 == "1"~"cDCs",
  pred_snn_res.0.5 == "2"~"Fibroblasts", 
  pred_snn_res.0.5 == "3"~"cDCs",
  pred_snn_res.0.5 == "4"~"NK cells", 
  pred_snn_res.0.5 == "5"~"Tumor cells", 
  pred_snn_res.0.5 == "6"~"M2 Macrophages",
  pred_snn_res.0.5 == "7"~"Fibroblasts", 
  pred_snn_res.0.5 == "8"~"Mesenchymal stromal cells?",
  pred_snn_res.0.5 == "9"~"Endothelial cells", 
  pred_snn_res.0.5 == "10"~"Tumor cells", 
  pred_snn_res.0.5 == "11"~"CD45+", 
  pred_snn_res.0.5 == "12"~"M1 Macrophages",
  pred_snn_res.0.5 == "13"~"Eosinophils",
  pred_snn_res.0.5 == "14"~"Neutrophils",
  pred_snn_res.0.5 == "15"~"PICs",
  pred_snn_res.0.5 == "16"~"Endothelial cells",
  pred_snn_res.0.5 == "17"~"Monocytes",
  pred_snn_res.0.5 == "17.1"~"M2 Macrophages",
  pred_snn_res.0.5 == "17.2"~"M1 Macrophages*Ery?", # potentially Erythrocytes stuck to Macrophages 
  pred_snn_res.0.5 == "18"~"unclear", 
  pred_snn_res.0.5 == "19"~"unclear", 
  pred_snn_res.0.5 == "19.1"~"unclear", 
  pred_snn_res.0.5 == "19.2"~"unclear", 
  pred_snn_res.0.5 == "20"~"regulatory Macrophages", 
  pred_snn_res.0.5 == "21"~"Tumor cells", 
  pred_snn_res.0.5 == "22"~"EMT-like tumor cells?", 
  pred_snn_res.0.5 == "22.1"~"unclear",
  pred_snn_res.0.5 == "23"~"Platelets",
  pred_snn_res.0.5 == "24"~"PICs",
  pred_snn_res.0.5 == "25"~"unclear", 
  pred_snn_res.0.5 == "26"~"unclear", 
  pred_snn_res.0.5 == "27"~"Mesenchymal stromal cells?",
  pred_snn_res.0.5 == "28"~"TAM Macrophages?",
  pred_snn_res.0.5 == "29"~"unclear",
  pred_snn_res.0.5 == "29.1"~"lowqual", # SiglecF in CD45neg, low quality 
  pred_snn_res.0.5 == "29.2"~"unclear", 
  pred_snn_res.0.5 == "30"~"unclear",
  pred_snn_res.0.5 == "31"~"lowqual" # low quality
))

## Macrophage subsets: 
### M1-like: CD86+, MHCII+, CD64+
### M2-like: CD206+, low CD86 and MHCII
### monocyte-derived: Ly6Chigh, intermediate F4/80
### regulatory: CD73+, CD11b+, CD64+ 
### relapse-specific: CD62P+, F4/80+, CD11b+, CD206+ 

# save annotated obj
saveRDS(obj, file = paste0(data.dir, "obj_per_group_sketched_projected_annotated.rds"))

# CLEANED DATASET -------------------------------------------------------------
# remove unclear populations and doublets, re-run analysis 
subset <- obj@meta.data %>% dplyr::filter(!celltype %in% c("unclear", "lowqual", "PICs"))
obj_clean <- subset(obj, cells = rownames(subset))

DefaultAssay(obj_clean) <- "sketch"
VariableFeatures(obj_clean) <- VariableFeatures(obj_clean) %>% str_subset("UPAR", negate = T)
obj_clean <- obj_clean %>%
  ScaleData() %>%
  RunPCA(npcs=38, approx=F, reduction.name = "pca") %>%
  FindNeighbors(dims = 1:38) %>%
  FindClusters(resolution = c(0.3, 0.5, 1), algorithm=1) %>% 
  RunUMAP(dims = 1:38, return.model = TRUE, reduction = "pca", reduction.name = "umap", seed.use = 42)

# save 
saveRDS(obj_clean, file = paste0(data.dir, "obj_per_group_sketched_clean.rds"))

obj_clean$umap_1 <- obj_clean@reductions$umap@cell.embeddings[,1]
obj_clean$umap_2 <- obj_clean@reductions$umap@cell.embeddings[,2]

DimPlot(obj_clean, group.by = "celltype", label = T)
DimPlot(obj_clean, group.by = "sketch_snn_res.0.5", label = T)
DimPlot(obj_clean, group.by = "sketch_snn_res.1", label = T)

# features 
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group_clean.pdf"), height = 7, width = 14, onefile = T)
FeaturePlot(obj_clean, features = Features(obj_clean)[c(1:12)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj_clean, features = Features(obj_clean)[13:24], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj_clean, features = Features(obj_clean)[c(25:36)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj_clean, features = Features(obj_clean)[c(37:39)], slot = "counts") & scale_color_viridis(option = 'magma')
dev.off() 

# features split by group
sketched <- obj_clean@meta.data %>% dplyr::filter(!is.na(sketch_snn_res.0.3))
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group_clean_by_group.pdf"), height = 7, width = 14, onefile = T)
features <- Features(obj_clean) %>% str_replace("-", "_")
for(i in seq_along(features)){
  print(ggplot(sketched, aes(x = umap_1, y = umap_2, color = .data[[features[i]]])) + 
    geom_point() + 
    facet_wrap(~group) + 
    theme_void() +
    labs(title = features[i]) + 
    scale_color_viridis(option = "magma") +
    theme(legend.position = "none"))
}
dev.off()

# check resolution = 0.5
Idents(obj_clean) <- "sketch_snn_res.0.5"
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_dimplot_highlighted_clusters_res0.5_sketched_per_group_clean.pdf"), height = 7, width = 9, onefile = T)
for (i in 0:(length(unique(obj_clean$sketch_snn_res.0.5))-2)) {
  print(DimPlot(obj_clean, cells.highlight = WhichCells(obj_clean, idents = as.character(i))))
}
dev.off()

# subcluster cluster 0 
obj_clean <- FindSubCluster(obj_clean, cluster = "0", resolution = 0.05, subcluster.name = "sub.cluster", graph.name = "sketch_snn") # refine Macrophage/Monocyte distinction
Idents(obj) <- "sub.cluster"

# plot
DimPlot(obj_clean, group.by = "sub.cluster", label =T)

# RidgePlot 
DefaultAssay(obj_clean) <- "sketch"
features <- Features(obj_clean)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_ridgeplots_resolution_0.5_all_clean", ".pdf"), height = 8, width = 5, onefile = T)
for (f in features) {
  print(RidgePlot(obj_clean, features = f, slot = 'data', assay = "sketch"))
}
dev.off()

# predict clustering --
obj_clean@meta.data <- obj_clean@meta.data %>% 
  mutate(sub.cluster = str_replace(sub.cluster, "_", "."))

DefaultAssay(obj_clean) <- "FACS"
obj_clean <- predict_data(obj = obj_clean,
                    data_query = obj_clean, 
                    ref_clusters = "sub.cluster", 
                    FSC.A="FSC.A",
                    FSC.H="FSC.H",
                    pred_name="pred_snn_res.0.5",
                    chunk_size=1000000, 
                    assay_ref = "sketch",
                    assay_query = "FACS")

# clean annotation ----
obj_clean@meta.data <- mutate(obj_clean@meta.data, celltype = case_when(
  pred_snn_res.0.5 == "0"~"Macrophages",
  pred_snn_res.0.5 == "0.1"~"Monocytes",
  pred_snn_res.0.5 == "1"~"cDCs",
  pred_snn_res.0.5 == "2"~"cDCs", 
  pred_snn_res.0.5 == "3"~"NK cells",
  pred_snn_res.0.5 == "4"~"Fibroblasts", 
  pred_snn_res.0.5 == "5"~"Tumor cells", 
  pred_snn_res.0.5 == "6"~"Macrophages",
  pred_snn_res.0.5 == "7"~"Fibroblasts", 
  pred_snn_res.0.5 == "8"~"CD73+ cells", # Mesenchymal stromal cells? but only hinges on CD73 signal
  pred_snn_res.0.5 == "9"~"Endothelial cells", 
  pred_snn_res.0.5 == "10"~"Tumor cells", 
  pred_snn_res.0.5 == "11"~"CD45+", 
  pred_snn_res.0.5 == "12"~"Macrophages",
  pred_snn_res.0.5 == "13"~"Eosinophils",
  pred_snn_res.0.5 == "14"~"Neutrophils",
  pred_snn_res.0.5 == "15"~"Endothelial cells",
  pred_snn_res.0.5 == "16"~"Macrophages",
  pred_snn_res.0.5 == "17"~"Platelets",
  pred_snn_res.0.5 == "18"~"Tumor cells", 
  pred_snn_res.0.5 == "19"~"Monocytes", 
  pred_snn_res.0.5 == "20"~"CD73+ cells", # Mesenchymal stromal cells? but only hinges on CD73 signal
  pred_snn_res.0.5 == "21"~"Ly6G+CD73+ cells", # EMT-like tumor cells? cannot be sure, but upregulation of CD73 (and Ly6G)
  pred_snn_res.0.5 == "22"~"remove", # very low CD45, low quality
  pred_snn_res.0.5 == "23"~"remove", # maybe Platelet/Erythrocyte aggregates 
  pred_snn_res.0.5 == "24"~"remove", # Erythrocytes?
))

# remove populations that cannot be annotated, re-calculate UMAP 
subset <- obj_clean@meta.data %>% dplyr::filter(!celltype %in% c("remove"))
obj_final <- subset(obj_clean, cells = rownames(subset))

DefaultAssay(obj_final) <- "sketch"
VariableFeatures(obj_final) <- VariableFeatures(obj_final) %>% str_subset("UPAR", negate = T)
obj_final <- obj_final%>%
  ScaleData() %>%
  RunPCA(npcs=38, approx=F, reduction.name = "pca") %>%
  RunUMAP(dims = 1:38, return.model = TRUE, reduction = "pca", reduction.name = "umap", seed.use = 42)

obj_final$umap_1 <- obj_final@reductions$umap@cell.embeddings[,1]
obj_final$umap_2 <- obj_final@reductions$umap@cell.embeddings[,2]

# save 
saveRDS(obj_final, file = paste0(data.dir, "obj__per_group_sketched_final.rds"))

# color palette
colors <- c(
  "cDCs" = "#CE87DE",                           
  "Tumor cells" = "magenta4",                    
  "Fibroblasts" = "green4",                     
  "CD45+" = "royalblue4",
  "Macrophages" = "goldenrod",  
  "Endothelial cells" = "#52A6B1",              
  "Eosinophils" = "#f2dc4b",                    
  "Neutrophils" = "pink2",                     
  "NK cells" = "#672C49",
  "Platelets" = "lightblue2",                     
  "Monocytes" = "gold4",                        
  "PICs" = "orange",                  
  "CD73+ cells" = "dodgerblue4", 
  "Ly6G+CD73+ cells" = "dodgerblue3"
)

# features
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group_final.pdf"), height = 7, width = 14, onefile = T)
FeaturePlot(obj_final, features = Features(obj_final)[c(1:12)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj_final, features = Features(obj_final)[13:24], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj_final, features = Features(obj_final)[c(25:36)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(obj_final, features = Features(obj_final)[c(37:39)], slot = "counts") & scale_color_viridis(option = 'magma')
dev.off() 

# supplement features figure layout  (Ext Figure 4B)
feat <- c("Ly6G", "Sca1", "CD45", "FAP", "CD105", "CD48", "CD41", "CD43", "CD86", 
          "CD16-32", "MHCII", "CD127", "CD62P", "CD71", "NK1", "CD64", "CD23", 
          "CD117", "CD80", "CD150", "CD73", "CD11c", "PDGFR", "CD135", "B220", 
          "CD206", "Ly6C", "CD34", "CD11b", "CD31", "CD21-35", "F4-80", "SiglecF")
Feature_Plot <- FeaturePlot(obj_final, features=feat, alpha=1, combine=F, raster=T, reduction="umap", pt.size=1.5, slot = "counts")
for(i in 1:length(Feature_Plot)) suppressMessages({
  Feature_Plot[[i]] <- Feature_Plot[[i]] + 
    scale_colour_gradientn(colours = pals::parula(n = 100)) +
    theme_classic() +
    NoLegend() + 
    NoAxes()
})
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group_final_rastered.pdf"), height = 8, width = 6)
print(cowplot::plot_grid(plotlist = Feature_Plot, nrow = 6))
dev.off()

# plot annotation (Figure 4D)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_dimplot_annotated_clean_overall_landscape_final.pdf"), height = 7, width = 10)
print(ggplot(obj_final@meta.data, aes(x = umap_1, y = umap_2, color = celltype)) + 
        ggrastr::geom_point_rast(size = 0.2, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(celltype)), size = 0.05, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = colors) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# plot annotation split by group (Ext Figure 4C)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_dimplot_annotated_clean_overall_landscape_split_by_group_landscape_final.pdf"), height = 3, width = 10)
print(ggplot(obj_final@meta.data, aes(x = umap_1, y = umap_2, color = celltype)) + 
        ggrastr::geom_point_rast(size = 0.2, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(celltype)), size = 0.05, raster.dpi = 700) + 
        facet_wrap(~group, nrow = 1) + 
        theme_classic() + 
        scale_color_manual(values = colors) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# IMMUNE COMPARTMENT -----------------------------------------------------------
immun <- subset(obj_clean, cells = rownames(obj_clean@meta.data %>% dplyr::filter(celltype %in% c("Macrophages", "Eosinophils", "Neutrophils", "cDCs", "CD45+", "Monocytes", "NK cells")))) 

DefaultAssay(immun) <- "sketch"
VariableFeatures(immun) <- VariableFeatures(immun) %>% str_subset("UPAR", negate = T)
immun <- immun %>%
  ScaleData() %>%
  RunPCA(npcs=38, approx=F, reduction.name = "pca") %>%
  FindNeighbors(dims = 1:38) %>%
  FindClusters(resolution = c(0.3, 0.5, 1), algorithm=1) %>% 
  RunUMAP(dims = 1:38, return.model = TRUE, reduction = "pca", reduction.name = "umap", seed.use = 42)

# save 
saveRDS(immun, file = paste0(data.dir, "obj__per_group_sketched_immune_only.rds"))

immun$umap_1 <- immun@reductions$umap@cell.embeddings[,1]
immun$umap_2 <- immun@reductions$umap@cell.embeddings[,2]

DimPlot(immun, group.by = "celltype", label = T)
DimPlot(immun, group.by = "sketch_snn_res.0.5", label = T)


# check resolution = 0.5
Idents(immun) <- "sketch_snn_res.0.5"
pdf(paste0(plot.dir, "02_qc/", Sys.Date(), "_dimplot_highlighted_clusters_res0.5_sketched_per_group_immun.pdf"), height = 7, width = 9, onefile = T)
for (i in 0:(length(unique(immun$sketch_snn_res.0.5))-2)) {
  print(DimPlot(immun, cells.highlight = WhichCells(immun, idents = as.character(i))))
}
dev.off()

# features
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group_immun.pdf"), height = 7, width = 14, onefile = T)
FeaturePlot(immun, features = Features(immun)[c(1:12)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(immun, features = Features(immun)[13:24], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(immun, features = Features(immun)[c(25:36)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(immun, features = Features(immun)[c(37:39)], slot = "counts") & scale_color_viridis(option = 'magma')
dev.off() 

# subcluster 3 (Macrophage subpopulations), 19
Idents(immun) <- "sketch_snn_res.0.5"
immun <- FindSubCluster(immun, cluster = "3", resolution = 0.1, subcluster.name = "sub.cluster", graph.name = "sketch_snn")
Idents(immun) <- "sub.cluster"
immun <- FindSubCluster(immun, cluster = "19", resolution = 0.1, subcluster.name = "sub.cluster", graph.name = "sketch_snn")
Idents(immun) <- "sub.cluster"

DimPlot(immun, group.by = "sub.cluster", label = T)

# RidgePlot 
DefaultAssay(immun) <- "sketch"
features <- Features(immun)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_ridgeplots_resolution_0.5_all_immun", ".pdf"), height = 8, width = 5, onefile = T)
for (f in features) {
  print(RidgePlot(immun, features = f, slot = 'data', assay = "sketch"))
}
dev.off()

# predict clustering --
immun@meta.data <- immun@meta.data %>% 
  mutate(sub.cluster = str_replace(sub.cluster, "_", "."))

immun <- predict_data(obj = immun,
                          data_query = immun, 
                          ref_clusters = "sub.cluster", 
                          FSC.A="FSC.A",
                          FSC.H="FSC.H",
                          pred_name="pred_snn_res.0.5",
                          chunk_size=1000000, 
                          assay_ref = "sketch",
                          assay_query = "FACS")


# refined annotation immune compartment ---
immun@meta.data <- mutate(immun@meta.data, celltype = case_when(
  pred_snn_res.0.5 == "0"~"M1-like Macrophages",
  pred_snn_res.0.5 == "1"~"activated cDCs",
  pred_snn_res.0.5 == "2"~"cDCs", 
  pred_snn_res.0.5 == "3"~"M1-like Macrophages",
  pred_snn_res.0.5 == "3.1"~"M2-like Macrophages",
  pred_snn_res.0.5 == "3.2"~"M1-like Macrophages",
  pred_snn_res.0.5 == "4"~"NK cells", 
  pred_snn_res.0.5 == "5"~"CD45+", 
  pred_snn_res.0.5 == "6"~"M2-like Macrophages",
  pred_snn_res.0.5 == "7"~"Eosinophils", 
  pred_snn_res.0.5 == "8"~"lowqual", # CD45 low 
  pred_snn_res.0.5 == "9"~"class. Monocytes", 
  pred_snn_res.0.5 == "10"~"cDCs", 
  pred_snn_res.0.5 == "11"~"Neutrophils", 
  pred_snn_res.0.5 == "12"~"NK cells",
  pred_snn_res.0.5 == "13"~"Regulatory Macrophages", # M2-like, CD73 upregulation 
  pred_snn_res.0.5 == "14"~"non-class. Monocytes",
  pred_snn_res.0.5 == "15"~"Undifferentiated Macrophages", # no M1/M2 polarization
  pred_snn_res.0.5 == "16"~"NK cells",
  pred_snn_res.0.5 == "17"~"lowqual", # CD45 low
  pred_snn_res.0.5 == "18"~"Eosinophils", 
  pred_snn_res.0.5 == "19"~"PICs", # remove
  pred_snn_res.0.5 == "19.1"~"PICs", # remove 
))

DimPlot(immun, group.by = "celltype", label = T)

# remove low quality populations and PICs, recalculate UMAP
subset <- immun@meta.data %>% dplyr::filter(!celltype %in% c("lowqual", "PICs"))
immun_final <- subset(immun, cells = rownames(subset))

DefaultAssay(immun_final) <- "sketch"
VariableFeatures(immun_final) <- VariableFeatures(immun_final) %>% str_subset("UPAR", negate = T)
immun_final <- immun_final %>%
  ScaleData() %>%
  RunPCA(npcs=38, approx=F, reduction.name = "pca") %>%
  RunUMAP(dims = 1:38, return.model = TRUE, reduction = "pca", reduction.name = "umap", seed.use = 42)

immun_final$umap_1 <- immun_final@reductions$umap@cell.embeddings[,1]
immun_final$umap_2 <- immun_final@reductions$umap@cell.embeddings[,2]

colors_immun <- c(
  "cDCs" = "#CE87DE", 
  "activated cDCs" = "orchid3", 
  "CD45+" = "royalblue4", 
  "class. Monocytes" = "yellowgreen", 
  "M1-like Macrophages" = "gold1", 
  "Eosinophils" = "#f2dc4b", 
  "Neutrophils" = "pink2", 
  "M2-like Macrophages" = "orange", 
  "Regulatory Macrophages" = "tomato2",
  "Undifferentiated Macrophages" = "sienna", 
  "NK cells" = "#672C49", 
  "non-class. Monocytes" = "gold4") 

# features
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group_immun_final.pdf"), height = 7, width = 14, onefile = T)
FeaturePlot(immun_final, features = Features(immun_final)[c(1:12)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(immun_final, features = Features(immun_final)[13:24], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(immun_final, features = Features(immun_final)[c(25:36)], slot = "counts") & scale_color_viridis(option = 'magma')
FeaturePlot(immun_final, features = Features(immun_final)[c(37:39)], slot = "counts") & scale_color_viridis(option = 'magma')
dev.off() 

# supplement features figure layout (Ext Figure 4E)
Feature_Plot <- FeaturePlot(immun_final, features=feat, alpha=1, combine=F, raster=T, reduction="umap", pt.size=1.5, slot = "counts")
for(i in 1:length(Feature_Plot)) suppressMessages({
  Feature_Plot[[i]] <- Feature_Plot[[i]] + 
    scale_colour_gradientn(colours = pals::parula(n = 100)) +
    theme_classic() +
    NoLegend() + 
    NoAxes()
})
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_feature_plots_sketched_per_group_immun_final_rastered.pdf"), height = 8, width = 6)
print(cowplot::plot_grid(plotlist = Feature_Plot, nrow = 6))
dev.off()

# plot (Figure 4E)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_dimplot_annotated_immun_landscape_final.pdf"), height = 7, width = 9)
print(ggplot(immun_final@meta.data, aes(x = umap_1, y = umap_2, color = celltype)) + 
        ggrastr::geom_point_rast(size = 0.2, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(celltype)), size = 0.05, raster.dpi = 700) + 
        theme_classic() + 
        scale_color_manual(values = colors_immun) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# split by group (Ext Figure 4F)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_dimplot_annotated_immun_landscape_final_by_group.pdf"), height = 3, width = 10)
print(ggplot(immun_final@meta.data, aes(x = umap_1, y = umap_2, color = celltype)) + 
        ggrastr::geom_point_rast(size = 0.2, color = "black", raster.dpi = 700) + 
        ggrastr::geom_point_rast(aes(color = as.factor(celltype)), size = 0.05, raster.dpi = 700) + 
        theme_classic() + 
        facet_wrap(~group) + 
        scale_color_manual(values = colors_immun) + 
        guides(colour = guide_legend(override.aes = list(size = 5))))
dev.off()

# QUANTIFICATION --------------------------------------------------------------
# cell type frequencies ----
# live, high quality event counts 
live <- obj@meta.data %>% 
  dplyr::filter(!celltype %in% c("lowqual")) %>% 
  group_by(sample) %>%
  count(name = "live")

# substract additional low quality events from the cleaned object
sub <- obj_clean@meta.data %>% 
  dplyr::filter(pred_snn_res.0.5 == "22") %>% 
  group_by(sample) %>%
  count(name = "low_quality")

live <- left_join(live, sub, by = "sample") %>% 
  mutate(low_quality = if_else(is.na(low_quality), 0, low_quality)) %>% 
  mutate(live = live - low_quality) %>% 
  select(-low_quality)

# quantification per group 
res <- obj_final@meta.data %>%
  # make factors for counting variables
  mutate(celltype = factor(celltype)) %>%
  mutate(sample = factor(sample)) %>% 
  # count while including zero counts
  group_by(celltype, sample, ID, group, .drop = F) %>%
  count() 

# add the meta.data 
res <- left_join(res, live, by = "sample") %>% 
  mutate(ID = str_extract(sample, "samples_.._(.*)_202", group = 1)) %>% 
  mutate(group = case_when(
    ID %in% c("3203", "3745", "3746", "3747", "3948") ~ 1, 
    ID %in% c("3153", "3204", "3318", "3321", "3455") ~ 2, 
    ID %in% c("3152", "3205", "3206", "32x9", "3320") ~ 3)) 

# calculate frequency
res <- res %>% 
  mutate(freq = n / live)

# plot 
plots = vector('list', length(unique(res$celltype)))
for(i in seq_along(unique(res$celltype))){
  plots[[i]] = ggstripchart(res %>% dplyr::filter(celltype == unique(res$celltype)[i]), 
                            x = "group", y = "freq", fill = "celltype", jitter = 0.2, size = 4, add = "mean_sd",
                            add.params = list(color = "black"), shape = 21, error.plot = "errorbar") +
    stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "errorbar", color = "black", width = 0.3, size = 0.5) + 
    labs(title = unique(res$celltype)[i]) + 
    theme(legend.position = "none") + 
    scale_fill_manual(values = colors)
}
print(ggarrange(plotlist = plots))

# stats 
stats <- res %>% 
  ungroup() %>% 
  group_by(celltype) %>%
  t_test(freq ~ group) %>% 
  adjust_pvalue(method = "holm")

# tumor vasculature ----
# CAVE: on day 1 (group 1&2) all cells were acquired, acquisition at day 2 was capped 
# -> all statements are relative 

# stats celltype frequency 
stats_vasc <- res %>% 
  dplyr::filter(celltype %in% c("Endothelial cells")) %>% 
  group_by(celltype) %>% 
  t_test(freq ~ group) %>% 
  adjust_pvalue(method = "holm") %>%
  add_xy_position()

# plot (Ext Figure 4D)
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_quantification_endothelial_cells.pdf"), height = 6, width = 8)
print(ggstripchart(res %>% dplyr::filter(celltype %in% c("Endothelial cells")), 
             x = "group", y = "freq", fill = "celltype", jitter = 0.2, 
             size = 6, add = "mean_sd", add.params = list(color = "black"), 
             shape = 21, error.plot = "errorbar", facet.by = "celltype") + 
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "errorbar", color = "black", width = 0.3, size = 0.5) + 
  theme_classic() +
  theme(legend.position = "none") + 
  stat_pvalue_manual(stats_vasc) + 
  scale_fill_manual(values = c("darkred")))
dev.off()

# QUANTIFICATION OF THE IMMUNE COMPARTMENT ONLY -------------------------------
# cell type frequencies ----
# live event counts for immune cells 
live_immun <- immun@meta.data %>% 
  # filter out low quality populations 
  dplyr::filter(!celltype %in% c("lowqual")) %>% 
  group_by(sample) %>%
  count(name = "live")

# quantification per group 
res <- immun_final@meta.data %>%
  # make factors for counting variables
  mutate(celltype = factor(celltype)) %>%
  mutate(sample = factor(sample)) %>% 
  # count while including zero counts
  group_by(celltype, sample, ID, group, .drop = F) %>%
  count() 

# add the meta.data 
res <- left_join(res, live_immun, by = "sample") %>% 
  mutate(ID = str_extract(sample, "samples_.._(.*)_202", group = 1)) %>% 
  mutate(group = case_when(
    ID %in% c("3203", "3745", "3746", "3747", "3948") ~ 1, # 3747 instead of 3947?
    ID %in% c("3153", "3204", "3318", "3321", "3455") ~ 2, 
    ID %in% c("3152", "3205", "3206", "32x9", "3320") ~ 3)) 

# calculate frequency
res <- res %>% 
  mutate(freq = n / live) %>% 
  # format 
  mutate(celltype = factor(celltype, levels = c("M1-like Macrophages", "M2-like Macrophages", "Regulatory Macrophages", 
                                              "Undifferentiated Macrophages", "activated cDCs", "cDCs", "class. Monocytes", 
                                              "non-class. Monocytes", "Eosinophils", "Neutrophils", "NK cells", "CD45+")))

# plot stripchart (Figure 4F) ----
plots = vector('list', length(unique(res$celltype)))
scaleFunc <- function(x) sprintf("%.3f", x)
celltypes <- c("M1-like Macrophages", "M2-like Macrophages", "Regulatory Macrophages", 
               "Undifferentiated Macrophages", "activated cDCs", "cDCs", "class. Monocytes", 
               "non-class. Monocytes", "Eosinophils", "Neutrophils", "NK cells", "CD45+")

for(i in seq_along(celltypes)){
  plots[[i]] = ggstripchart(res %>% dplyr::filter(celltype == celltypes[i]), 
                            x = "group", y = "freq", fill = "celltype", jitter = 0.2, size = 4, add = "mean_sd",
                            add.params = list(color = "black"), shape = 21, error.plot = "errorbar") +
    stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "errorbar", color = "black", width = 0.3, size = 0.5) + 
    labs(title = celltypes[i]) + 
    scale_y_continuous(n.breaks = 4, labels = scaleFunc) + 
    theme(legend.position = "none", axis.text.y = element_text()) + 
    scale_fill_manual(values = colors_immun) 
}
pdf(paste0(plot.dir, "03_annotation/", Sys.Date(), "_quantification_immune_cells.pdf"), height = 8, width = 7)
print(ggarrange(plotlist = plots, ncol = 3, nrow = 4))
dev.off()

# stats 
stats <- res %>% 
  ungroup() %>% 
  group_by(celltype) %>%
  t_test(freq ~ group) %>% 
  adjust_pvalue(method = "BH")