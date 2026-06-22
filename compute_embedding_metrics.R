library(FNN)
library(cluster)

compute_embedding_metrics <- function(
    embedding,
    labels,
    k = 20,
    scale_embedding = TRUE,
    max_silhouette_n = 5000,
    seed = 1
) {
  embedding <- as.matrix(embedding)
  labels <- as.character(labels)

  # cell IDで揃える
  if (!is.null(rownames(embedding)) && !is.null(names(labels))) {
    common <- intersect(rownames(embedding), names(labels))
    embedding <- embedding[common, , drop = FALSE]
    labels <- labels[common]
  }

  # NA labelを除く
  keep <- !is.na(labels) & labels != "" & labels != "NA"
  embedding <- embedding[keep, , drop = FALSE]
  labels <- labels[keep]

  # 分散ゼロの次元を除く
  vars <- apply(embedding, 2, var)
  embedding <- embedding[, vars > 0, drop = FALSE]

  if (scale_embedding) {
    embedding <- scale(embedding)
  }

  n <- nrow(embedding)

  if (n <= k + 1) {
    stop("Number of cells must be larger than k + 1.")
  }

  # ---------- kNN ----------
  # get.knn は自分自身を含まない近傍を返す
  knn <- FNN::get.knn(embedding, k = k)
  nn_index <- knn$nn.index

  purity_each <- numeric(n)
  pred <- character(n)

  for (i in seq_len(n)) {
    nn_labels <- labels[nn_index[i, ]]

    purity_each[i] <- mean(nn_labels == labels[i])

    tab <- sort(table(nn_labels), decreasing = TRUE)
    pred[i] <- names(tab)[1]
  }

  correct_each <- pred == labels

  knn_purity <- mean(purity_each)
  knn_classification_accuracy <- mean(correct_each)

  # ---------- per-cell-type metrics ----------
  label_levels <- sort(unique(labels))

  per_cell_type <- do.call(
    rbind,
    lapply(label_levels, function(x) {
      idx <- labels == x

      data.frame(
        cell_type = x,
        n_cells = sum(idx),
        knn_purity = mean(purity_each[idx]),
        knn_classification_accuracy = mean(correct_each[idx]),
        stringsAsFactors = FALSE
      )
    })
  )

  macro_knn_accuracy <- mean(per_cell_type$knn_classification_accuracy)
  macro_knn_purity <- mean(per_cell_type$knn_purity)

  # ---------- silhouette ----------
  set.seed(seed)

  if (n > max_silhouette_n) {
    sampled <- sample(seq_len(n), max_silhouette_n)
  } else {
    sampled <- seq_len(n)
  }

  emb_s <- embedding[sampled, , drop = FALSE]
  lab_s <- labels[sampled]

  # 1細胞しかないlabelはsilhouetteで問題になるので除く
  tab_s <- table(lab_s)
  keep_s <- lab_s %in% names(tab_s)[tab_s >= 2]

  emb_s <- emb_s[keep_s, , drop = FALSE]
  lab_s <- lab_s[keep_s]

  if (length(unique(lab_s)) < 2) {
    silhouette_mean <- NA_real_
  } else {
    d <- dist(emb_s)
    sil <- cluster::silhouette(as.integer(factor(lab_s)), d)
    silhouette_mean <- mean(sil[, "sil_width"])
  }

  # ---------- confusion matrix ----------
  confusion_matrix <- table(
    true = labels,
    predicted = pred
  )

  # ---------- cell-level prediction table ----------
  predictions <- data.frame(
    cell_id = if (!is.null(rownames(embedding))) rownames(embedding) else seq_len(n),
    true_label = labels,
    predicted_label = pred,
    correct = correct_each,
    knn_purity = purity_each,
    stringsAsFactors = FALSE
  )

  # ---------- summary ----------
  summary <- data.frame(
    n_cells = n,
    n_dimensions = ncol(embedding),
    k = k,
    knn_purity = knn_purity,
    macro_knn_purity = macro_knn_purity,
    knn_classification_accuracy = knn_classification_accuracy,
    macro_knn_accuracy = macro_knn_accuracy,
    silhouette = silhouette_mean,
    silhouette_n = length(lab_s)
  )

  return(
    list(
      summary = summary,
      per_cell_type = per_cell_type,
      confusion_matrix = confusion_matrix,
      predictions = predictions
    )
  )
}