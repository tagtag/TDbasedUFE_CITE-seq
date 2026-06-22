compute_graph_knn_metrics_from_nn <- function(
    nn_index,
    labels,
    cell_names = NULL,
    k = 20
) {
  labels <- as.character(labels)

  if (is.null(cell_names)) {
    if (!is.null(names(labels))) {
      cell_names <- names(labels)
    } else {
      cell_names <- seq_along(labels)
    }
  }

  if (is.null(names(labels))) {
    names(labels) <- cell_names
  }

  # nn_index は cells x neighbors の整数行列を想定
  nn_index <- as.matrix(nn_index)

  n <- nrow(nn_index)

  if (length(labels) != n) {
    stop("length(labels) must match nrow(nn_index).")
  }

  if (length(cell_names) != n) {
    stop("length(cell_names) must match nrow(nn_index).")
  }

  names(labels) <- cell_names

  keep <- !is.na(labels) & labels != "" & labels != "NA"
  labels <- labels[keep]
  cell_names <- cell_names[keep]
  nn_index <- nn_index[keep, , drop = FALSE]

  # 元のindexから、keep後のindexへ変換するための対応表
  old_to_new <- rep(NA_integer_, n)
  old_to_new[which(keep)] <- seq_len(sum(keep))

  n2 <- length(labels)

  purity_each <- numeric(n2)
  pred <- character(n2)

  for (i in seq_len(n2)) {
    # Seuratのnn.idxは元のcell indexを返すので、keep後のindexへ変換
    neigh_old <- nn_index[i, ]
    neigh_new <- old_to_new[neigh_old]

    # NAと自分自身を除く
    neigh_new <- neigh_new[!is.na(neigh_new)]
    neigh_new <- neigh_new[neigh_new != i]

    if (length(neigh_new) == 0) {
      purity_each[i] <- NA_real_
      pred[i] <- NA_character_
      next
    }

    neigh_new <- neigh_new[seq_len(min(k, length(neigh_new)))]
    nn_labels <- labels[neigh_new]

    purity_each[i] <- mean(nn_labels == labels[i])

    tab <- sort(table(nn_labels), decreasing = TRUE)
    pred[i] <- names(tab)[1]
  }

  correct_each <- pred == labels

  knn_purity <- mean(purity_each, na.rm = TRUE)
  knn_classification_accuracy <- mean(correct_each, na.rm = TRUE)

  label_levels <- sort(unique(labels))

  per_cell_type <- do.call(
    rbind,
    lapply(label_levels, function(x) {
      idx <- labels == x

      data.frame(
        cell_type = x,
        n_cells = sum(idx),
        graph_knn_purity = mean(purity_each[idx], na.rm = TRUE),
        graph_knn_classification_accuracy = mean(correct_each[idx], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )

  summary <- data.frame(
    n_cells = n2,
    k = k,
    graph_knn_purity = knn_purity,
    macro_graph_knn_purity = mean(per_cell_type$graph_knn_purity, na.rm = TRUE),
    graph_knn_classification_accuracy = knn_classification_accuracy,
    macro_graph_knn_accuracy = mean(
      per_cell_type$graph_knn_classification_accuracy,
      na.rm = TRUE
    )
  )

  confusion_matrix <- table(
    true = labels,
    predicted = pred
  )

  predictions <- data.frame(
    cell_id = cell_names,
    true_label = labels,
    predicted_label = pred,
    correct = correct_each,
    graph_knn_purity = purity_each,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    per_cell_type = per_cell_type,
    confusion_matrix = confusion_matrix,
    predictions = predictions
  )
}