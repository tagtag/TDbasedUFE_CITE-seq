y0 <- read.csv("GSE297097_All.IGTs.Enhanced.Summary.Table.2026-02-12.csv.gz")
y <- read.csv("GSE297097_annotation_table_20260206_IGT1_104_cleaned.csv.gz")
source("compute_embedding_metrics.R")
source("compute_graph_knn_metrics.R")
require(Seurat)
#GSE301271
#---- TD ---
load("fit_GSE301271")
td_embedding <- fit$U3
#--- TD -----
#--- scMoMaT ---
load("cell_factors_py_GSE301271")
td_embedding <- cell_factors_py[[1]][,1:10] 
#--- scMoMaT ----
#--- Seurat ---
obj <- readRDS("results_GSE301271/GSE301271_seurat_wnn.rds")#RNA_only, ADT_only, Seurat Wnn
#-- Seurat ------


y1 <- read.csv("GSE301271_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT36-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT36/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc
col <- factor(y0$organ[match(y1[,1],y0$IGT.cellID)])

#GSE301960
#--- TD ---
load("fit_GSE301960")
td_embedding <- fit$U3
#---- TD ----
#---- scMoMaT ---
load("cell_factors_py_GSE301960")
td_embedding <- cell_factors_py[[1]][,1:10] 
#---- scMoMaT ---
#--- Seurat -----
obj <- readRDS("results_GSE301960/GSE301960_seurat_wnn.rds") #RNA_only, ADT_only, Seurat_Wnn
#--- Seurat ----

y1 <- read.csv("GSE301960_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT38-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT38/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc

#GSE301961
#---- TD ----
load("fit_GSE301961")
td_embedding <- fit$U3
#--- TD ----
#-- scMoMaT ----
load("cell_factors_py_GSE301961")
td_embedding <- cell_factors_py[[1]][,1:10] 
#---- scMoMaT ----
#--- Seurat ----
obj <- readRDS("results_GSE301961/GSE301961_seurat_wnn.rds") #RNA_only, ADT_only, Seurat Wnn
#--- Seurat ---

y1 <- read.csv("GSE301961_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT38-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT38/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc

#GSE281719
#---- TD ----
load("fit_GSE281719")
td_embedding <- fit$U3
#---- TD ---
#--- scMoMaT ----
load("cell_factors_py_GSE281719")
td_embedding <- cell_factors_py[[1]][,1:10]
#--- scMoMaT ----
#--- Seurat ---
obj <- readRDS("results_GSE281719/GSE281719_seurat_wnn.rds") #RNA_only, ADT_only, Seurat Wnn
#--- Seurat ----

y1 <- read.csv("GSE281719_cells.tsv.gz",sep="\t",header=F)
y2 <- read.csv("IGT15-cells-20260406.csv") #taken from https://rosetta.immgen.org/IGT15/tab/integration#?
col <- factor(y$level1[match(y1[,1],y$IGT.cellID)]) #CD4, CD8 etc

#---- Execution ---
labels <- col
names(labels) <- y1[,1]

#--- TD, scMoMaT ----
td_metrics <- compute_embedding_metrics(
  embedding = td_embedding[, 1:10],
  labels = labels,
  k = 20,
  scale_embedding = TRUE,
  max_silhouette_n = 5000
)
td_metrics

#---- RNA only ----
rna_embedding <- Embeddings(obj, "pca")[, 1:10]

rna_metrics <- compute_embedding_metrics(
  embedding = rna_embedding,
  labels = labels,
  k = 20,
  scale_embedding = TRUE,
  max_silhouette_n = 5000
)

#--- ADT only ----
adt_embedding <- Embeddings(obj, "apca")

adt_metrics <- compute_embedding_metrics(
  embedding = adt_embedding,
  labels = labels,
  k = 20,
  scale_embedding = TRUE,
  max_silhouette_n = 5000
)

#---Seurat Wnn ----
wnn_nn <- obj@neighbors$weighted.nn

nn_index <- wnn_nn@nn.idx

wnn_metrics <- compute_graph_knn_metrics_from_nn(
  nn_index = nn_index,
  labels = labels,
  cell_names = colnames(obj),
  k = 20
)
