run_seurat_wnn <- function(gse_id, data_dir = ".", out_dir = ".") {
  message("Running Seurat WNN for ", gse_id)

  # out_dir が存在しないと saveRDS/write.csv で落ちるので作る
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  rna_file  <- file.path(data_dir, paste0(gse_id, "_matrix_rna.mtx.gz"))
  adt_file  <- file.path(data_dir, paste0(gse_id, "_matrix_protein.mtx.gz"))
  cell_file <- file.path(data_dir, paste0(gse_id, "_cells.tsv.gz"))
  gene_file <- file.path(data_dir, paste0(gse_id, "_genes.tsv.gz"))

  # protein / protain / ADT などを含む tsv.gz を探す
  protein_candidates <- list.files(
    data_dir,
    pattern = paste0(
      "^", gse_id,
      ".*(protein|proteins|protain|protains|adt|ADT).*\\.tsv(\\.gz)?$"
    ),
    full.names = TRUE
  )

  if (length(protein_candidates) == 0) {
    stop(
      "Protein/ADT feature file was not found. Existing files for this GSE are:\n",
      paste(list.files(data_dir, pattern = gse_id), collapse = "\n")
    )
  }

  protein_file <- protein_candidates[1]
  message("Using protein feature file: ", protein_file)

  rna <- readMM(rna_file)
  adt <- readMM(adt_file)

  cells <- read.table(
    cell_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )[, 1]

  genes <- read.table(
    gene_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )[, 1]

  proteins <- read.table(
    protein_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )[, 1]

  rownames(rna) <- make.unique(genes)
  colnames(rna) <- cells

  rownames(adt) <- make.unique(proteins)
  colnames(adt) <- cells

  # 念のため細胞順を揃える
  common_cells <- intersect(colnames(rna), colnames(adt))
  rna <- rna[, common_cells]
  adt <- adt[, common_cells]

  obj <- CreateSeuratObject(
    counts = rna,
    assay = "RNA",
    project = gse_id
  )

  obj[["ADT"]] <- CreateAssayObject(counts = adt)

  # RNA側
  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, nfeatures = 4000)
  obj <- ScaleData(obj)
  obj <- RunPCA(
    obj,
    npcs = 50,
    reduction.name = "pca",
    reduction.key = "PC_"
  )

  # ADT側
  DefaultAssay(obj) <- "ADT"

  VariableFeatures(obj) <- rownames(obj[["ADT"]])

  obj <- NormalizeData(
    obj,
    normalization.method = "CLR",
    margin = 2
  )

  obj <- ScaleData(obj)

  n_adt_pc <- min(30, nrow(obj[["ADT"]]) - 1)

  obj <- RunPCA(
    obj,
    reduction.name = "apca",
    reduction.key = "apca_",
    npcs = n_adt_pc
  )

  # WNN
  obj <- FindMultiModalNeighbors(
    object = obj,
    reduction.list = list("pca", "apca"),
    dims.list = list(1:30, 1:n_adt_pc),
    modality.weight.name = c("RNA.weight", "ADT.weight")
  )

  obj <- RunUMAP(
    obj,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    seed.use = 1
  )

  obj <- FindClusters(
    obj,
    graph.name = "wsnn",
    resolution = 0.5,
    algorithm = 3,
    random.seed = 1
  )

  saveRDS(
    obj,
    file.path(out_dir, paste0(gse_id, "_seurat_wnn.rds"))
  )

  write.csv(
    Embeddings(obj, "wnn.umap"),
    file.path(out_dir, paste0(gse_id, "_seurat_wnn_umap.csv"))
  )

  write.csv(
    obj@meta.data,
    file.path(out_dir, paste0(gse_id, "_seurat_wnn_metadata.csv"))
  )

  return(obj)
}

obj <- run_seurat_wnn(
  "GSE301271",
  data_dir = "./",
  out_dir = "results_GSE301271"
)

obj <- run_seurat_wnn(
  "GSE301960",
  data_dir = "./",
  out_dir = "results_GSE301960"
)

obj <- run_seurat_wnn(
  "GSE301961",
  data_dir = "./",
  out_dir = "results_GSE301961"
)


obj <- run_seurat_wnn(
  "GSE281719",
  data_dir = "./",
  out_dir = "results_GSE281719"
)
