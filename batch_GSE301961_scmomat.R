# 必要なパッケージの読み込み
library(Matrix)
library(reticulate)
require(umap)

# ---------------------------------------------------------
# 1. Python環境のセットアップ
# ---------------------------------------------------------
# scmomatがインストールされているPython環境を指定します
# 例としてConda環境名 "scmomat_env" を使用する場合：
use_condaenv("scmomat_env", required = TRUE)
# use_python("/path/to/your/python") # パスで直接指定する場合

# Pythonパッケージのインポート
scmomat <- import("scmomat")

# ---------------------------------------------------------
# 2. データの読み込み (Matrixパッケージ)
# ---------------------------------------------------------
# Sparse Matrix (mtx形式) を読み込みます
rna_matrix <- readMM("GSE301961_matrix_rna.mtx.gz")
protein_matrix <- readMM("GSE301961_matrix_protein.mtx.gz")

# 【重要】データ次元の向きの調整
# R (Seuratなど) では「行=遺伝子/タンパク質, 列=細胞 (Features x Cells)」が標準ですが、
# Python (scMoMaTなど) では「行=細胞, 列=遺伝子/タンパク質 (Cells x Features)」が標準です。
# そのため、読み込んだ行列を転置(t)します。
rna_counts <- t(rna_matrix)
protein_counts <- t(protein_matrix)

# ---------------------------------------------------------
# 3. Pythonオブジェクトへの変換とデータ構造の準備（再々修正）
# ---------------------------------------------------------
# 一度Pythonの疎行列に変換します（ここはそのまま）
rna_py <- r_to_py(rna_counts)
protein_py <- r_to_py(protein_counts)

# 【修正】
# PyTorchが読み込めるように、疎行列 (Sparse) を密行列 (Dense/NumPy配列) に変換します。
rna_dense <- rna_py$toarray()
protein_dense <- protein_py$toarray()

# 密行列にしたデータを辞書に格納します
counts_dict <- reticulate::dict(
  nbatches = 1L,
  RNA = list(rna_dense),
  Protein = list(protein_dense)
)

# ---------------------------------------------------------
# 4. scMoMaTモデルの構築と学習
# ---------------------------------------------------------
# これでPyTorchが理解できる普通の配列として渡されます
model <- scmomat$scmomat_model(
  counts = counts_dict,
  K = 30L, 
  device = "cpu" 
)

cat("モデルの学習を開始します...\n")
losses <- model$train_func(
  T = 4000L 
)

save(file="model_GSE301961",model)
# ---------------------------------------------------------
# 5. 解析結果（潜在表現）の取得とRでの利用（修正版）
# ---------------------------------------------------------
# 統合された細胞の低次元エンベディング（Cell factors）を取得します
# 【修正】関数名: extract_cell_factors を使用します
cell_factors_py <- model$extract_cell_factors()
save(file="cell_factors_py_GSE301961",cell_factors_py)

UMAP <- umap(cell_factors_py[[1]])
save(file="UMAP_scmomat_GSE301961",UMAP)



