#' @include utils.R
#' @include objects.R
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Multi Classifier ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Train Ensemble Classifier for Cell Type Prediction
#'
#' Trains multiple classifiers (SVM by default) on balanced subsets of data
#' to create robust ensemble model for cell type prediction.
#'
#' @import mlr3verse
#' @importFrom mlr3 as_task_classif lrn
#' @importFrom parallel detectCores
#' @importFrom pbapply pblapply
#'
#' @param feature.mat Signature score matrix (cells x features)
#' @param cell.types data.frame of cell type labels (one column per granularity)
#' @param do.norm Normalization method: 'L1', 'L2', or NULL. Default: 'L2'
#' @param mlr3.model MLR3 learner name (e.g. "classif.svm"). Default: "classif.svm"
#' @param batch.size Cells per model training. Default: 5000
#' @param n.models Number of models in ensemble. Default: 100
#' @param balance.cell.type Whether to balance cell types in sampling. Default: TRUE
#' @param cores Number of CPU cores (-1 = all available). Default: -1
#'
#' @return List of trained learners with class "Classif"
#'
#' @examples
#' \dontrun{
#' features <- matrix(rnorm(1000), nrow=100)
#' types <- data.frame(celltype=sample(LETTERS[1:3], 100, replace=TRUE))
#' model <- FitEnsemblMultiClassif(feature.mat=features, cell.types=types, n.models=10)
#' }
#'
#' @concept training_model
#' @export
FitEnsemblMultiClassif <- function(feature.mat, cell.types, do.norm = 'L2', 
                                  mlr3.model = "classif.svm", batch.size=5000,
                                  n.models=100, balance.cell.type = TRUE, cores=-1){
  ## Input validation
  if (!is.data.frame(feature.mat))) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  
  if (!all(rownames(feature.mat) == rownames(cell.types))) {
    stop("The rownames of 'feature.mat' and 'cell.types' are not matched.")
  }
  
  if (batch.size >= nrow(feature.mat))) {
    stop("Batch size is larger than the given cells. A smaller 'batch.size' is expected.")
  }
  
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  
  ## Clean column names
  bad.cols <- grepl("-", colnames(feature.mat))
  if (any(bad.cols))) {
    warning("Bad colnames in 'feature.mat'. Change the '-' to '_'.")
    colnames(feature.mat) <- gsub("-", "_", colnames(feature.mat))
  }
  
  bad.cols <- grepl("-", colnames(cell.types))
  if (any(bad.cols))) {
    warning("Bad colnames in 'cell.types'. Change the '-' to '_'.")
    colnames(cell.types) <- gsub("-", "_", colnames(cell.types))
  }
  
  ## Convert cell types to factors
  for (i in seq_along(colnames(cell.types))) {
    cell.types[[i]] %<>% as.factor()
  }
  
  ## Normalization
  if (do.norm == "L1")) {
    message("L1 normalization ...")
    feature.mat <- L1norm(feature.mat)
  } else if (do.norm == "L2")) {
    message("L2 normalization ...")
    feature.mat <- L2norm(feature.mat)
  } else if (is.null(do.norm))) {
    message("Skip normalization ...")
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  
  ## Create training tasks
  message(sprintf("Creating tasks: size=%s, n=%s", batch.size, n=n.models))
  
  # Sampling probability calculation for balanced sampling
  .getP <- function(labels) {
    labels.freq <- as.vector(table(labels))
    names(labels.freq) <- names(table(labels))
    1 / (labels.freq[labels] * length(labels.freq))
  }
  
  # Create tasks for each cell type granularity
  targets <- colnames(cell.types)
  tasks <- lapply(targets, function(ii) {
    if (balance.cell.type) {
      p <- .getP(cell.types[[ii]])
    } else {
      p <- NULL
    }
    
    pbapply::pblapply(1:n.models, function(jj) {
      row.ids <- sample(rownames(feature.mat), size = batch.size, replace = F, prob = p)
      data.use <- cbind(feature.mat[row.ids, , drop=FALSE], cell.types[row.ids, ii, drop=FALSE])
      as_task_classif(as.data.frame(data.use), target = ii)
    })
  })
  names(tasks) <- targets
  
  ## Training function
  .train <- function(dataset, batch.id) {
    task = dataset[[batch.id]]
    if (!mlr3.model %in% mlr3::mlr_learners$keys()) {
      stop(sprintf("Undefined 'mlr3.model': %s.", mlr3.model))
    }
    learner = lrn(mlr3.model)
    learner$train(task)
    learner
  }
  
  ## Train models in parallel
  message("Training ...")
  celltype_lrns = lapply(tasks, function(celltype_i) {
    parallel::mclapply(1:n.models, function(j) .train(celltype_i, j), mc.cores=cores)
  })
  names(celltype_lrns) = names(tasks)
  
  ## Return model object
  model <- list(
    model = celltype_lrns,
    features = colnames(feature.mat)
  )
  class(model) <- "Classif"
  return(model)
}

#' Predict Cell Types Using Trained Ensemble
#'
#' Predicts cell types for new data using trained ensemble classifier.
#' Handles feature matching and normalization to match training conditions.
#'
#' @param feature.mat Signature score matrix to predict (cells x features)
#' @param model Trained model from FitEnsemblMultiClassif()
#' @param do.norm Normalization method (must match training). Default: 'L2'
#' @param cores Number of CPU cores (-1 = all available). Default: -1
#'
#' @return data.frame of predicted cell types (one column per granularity)
#'
#' @examples
#' \dontrun{
#' new_data <- matrix(rnorm(500), nrow=50)
#' preds <- PredictNewdata(feature.mat=new_data, model=trained_model)
#' }
#'
#' @concept label_transfer
#' @export
PredictNewdata <- function(feature.mat, model, do.norm='L2', cores=-1){
  ## Input validation
  if (!is.data.frame(feature.mat))) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  
  # Check feature compatibility
  bad.features <- setdiff(colnames(feature.mat), model$features)
  if (length(bad.features) > 0) {
    stop(sprintf("These features are not in the model: %s", paste(bad.features, collapse = " ")))
  }
  
  # Handle missing features (fill with zeros)
  no.features <- setdiff(model$features, colnames(feature.mat))
  if (length(no.features) > 0) {
    warning(sprintf("These features are not in your inputs: %s\n. Filled with zeros.", 
                   paste(no.features, collapse = " ")))
    zeros <- matrix(rep(0, nrow(feature.mat)*length(no.features)), 
                  ncol = length(no.features)) %>% as.data.frame()
    colnames(zeros) <- no.features
    newdata <- cbind(feature.mat, zeros)[, model$features]
  } else {
    newdata <- feature.mat
  }
  
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  
  ## Normalization (must match training)
  if (do.norm == "L1")) {
    message("L1 normalization ...")
    newdata <- L1norm(newdata)
  } else if (do.norm == "L2")) {
    message("L2 normalization ...")
    newdata <- L2norm(newdata)
  } else if (is.null(do.norm))) {
    message("Skip normalization ...")
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  
  ## Parallel prediction
  n.models <- length(model$model[[1]])
  celltype_pred <- lapply(model$model, function(lrns) {
    celltype_i <- parallel::mclapply(
      X = 1:n.models,
      FUN = function(i) {
        lrns[[i]]$predict_newdata(newdata = as.data.frame(newdata))$data$response %>% as.character()
      }, mc.cores=cores)
    do.call(cbind, celltype_i)
  })
  names(celltype_pred) = names(model$model)
  
  ## Majority vote integration
  int.fun <- function(xx) {
    celltype.counts <- sort(table(xx), decreasing = T)
    names(celltype.counts)[1]
  }
  
  celltype_inter <- lapply(names(celltype_pred), function(xx) {
    apply(celltype_pred[[xx]], 1, int.fun)
  })
  names(celltype_inter) <- names(celltype_pred)
  celltype_inter <- as.data.frame(celltype_inter)
  rownames(celltype_inter) <- rownames(feature.mat)
  
  return(celltype_inter)
}

#' Perform K-means Over-clustering
#'
#' Helper function to cluster cells into k groups for majority voting.
#' Used to smooth noisy predictions at single-cell level.
#'
#' @param feature.mat Feature matrix for clustering (cells x features)
#' @param k Number of clusters
#'
#' @return Vector of cluster assignments
#'
#' @concept label_transfer
OverCluster <- function(feature.mat, k){
  clusters <- stats::kmeans(feature.mat, centers = k)$cluster
  return(clusters)
}

#' Majority Vote Cell Type Assignment
#'
#' Harmonizes cell type labels by majority voting within over-clustered groups.
#' Helps resolve ambiguous single-cell predictions.
#'
#' @param feature.mat Feature matrix for over-clustering (optional)
#' @param over.clusters Pre-computed cluster assignments (optional)
#' @param cell.types data.frame of original cell type labels
#' @param k Number of clusters if feature.mat provided. Default: 20
#' @param min.prop Minimum proportion to assign cluster to cell type. 
#'        Below this, labeled "Heterogeneous". Default: 0
#'
#' @return data.frame of consensus cell type assignments
#'
#' @examples
#' \dontrun{
#' features <- matrix(rnorm(1000), nrow=100)
#' types <- data.frame(celltype=sample(LETTERS[1:3], 100, replace=TRUE))
#' consensus <- MajorityVote(feature.mat=features, cell.types=types, k=10)
#' }
#'
#' @concept label_transfer
#' @export
MajorityVote <- function(feature.mat = NULL, over.clusters = NULL, cell.types, k = 20, min.prop = 0){
  ## Parameter validation
  if (is.null(feature.mat) && is.null(over.clusters))) {
    stop("Must provide one of the 'feature.mat' and 'over.clusters'.")
  }
  if (!is.null(over.clusters) && is.null(names(over.clusters))) {
    stop("'over.clusters' should be named.")
  }
  if (!is.null(feature.mat) && is.null(rownames(feature.mat))) {
    stop("No rownames of 'feature.mat'")
  }
  
  ## Perform over-clustering if needed
  if (is.null(over.clusters))) {
    common.cells <- intersect(rownames(feature.mat), rownames(cell.types))
    if (length(common.cells) != nrow(cell.types))) {
      stop("Inconsistent rownames between 'feature.mat' and 'cell.types'.")
    }
    over.clusters <- OverCluster(feature.mat, k = k)
  } else {
    common.cells <- intersect(names(over.clusters), rownames(cell.types))
    if (length(common.cells) != nrow(cell.types))) {
      stop("Inconsistent cell names between 'over.clusters' and 'cell.types'.")
    }
  }
  
  ## Perform majority voting per cluster
  all.vars <- colnames(cell.types)
  cluster.levels <- unique(over.clusters)
  
  major.votes <- lapply(cluster.levels, function(clu) {
    cells.use <- names(over.clusters[over.clusters == clu])
    row.mv <- sapply(seq_along(all.vars), function(ii) {
      freq <- table(cell.types[cells.use, all.vars[ii]]) / length(cells.use)
      top1 <- sort(freq, decreasing = TRUE)[1] %>% names()
      ifelse(top1 > min.prop, top1, "Heterogeneous")
    })
    names(row.mv) <- paste0(all.vars, ".major_votes")
    row.mv
  })
  
  major.votes <- do.call(rbind, major.votes)
  rownames(major.votes) <- cluster.levels
  
  ## Map consensus labels back to cells
  cell.types$over.clusters <- as.character(over.clusters[rownames(cell.types)])
  pred.df <- major.votes[cell.types$over.clusters, , drop=FALSE]
  rownames(pred.df) <- rownames(cell.types)
  
  return(as.data.frame(pred.df))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### KNN based label transfer ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' K-Nearest Neighbors Label Transfer
#'
#' Transfers labels from reference to query data using k-NN in embedding space.
#' Supports both categorical and continuous label transfer.
#'
#' @param query.emb Query embedding matrix (cells x dimensions)
#' @param ref.emb Reference embedding matrix (cells x dimensions)
#' @param ref.labels Named vector of reference labels (cell names = names)
#' @param k Number of nearest neighbors. Default: 100
#'
#' @return data.frame with:
#' \itemize{
#'   \item labels - Predicted labels
#'   \item votes - Number of supporting neighbors (for categorical)
#'   \item perc - Proportion of supporting neighbors
#' }
#'
#' @examples
#' \dontrun{
#' query_emb <- matrix(rnorm(200), nrow=50)
#' ref_emb <- matrix(rnorm(200), nrow=100)
#' ref_labs <- sample(LETTERS[1:3], 100, replace=TRUE)
#' names(ref_labs) <- paste0("cell_", 1:100)
#' knn_pred <- KnnLabelTransfer(query.emb=query_emb, ref.emb=ref_emb, 
#'                            ref.labels=ref_labs, k=30)
#' }
#'
#' @concept label_transfer
#' @export
KnnLabelTransfer <- function(query.emb, ref.emb, ref.labels, k=100) {
  ## Parameter validation
  if (is.null(names(ref.labels))) {
    stop("'ref.labels' should be named.")
  }
  
  common.cells <- intersect(rownames(ref.emb), names(ref.labels))
  if (length(common.cells) != nrow(ref.emb) || length(common.cells) != length(ref.labels))) {
    stop("The rownames of 'ref.emb' should be consistent with the names of 'ref.labels'.")
  }
  
  # Reorder labels to match reference embedding
  ref.labels <- ref.labels[rownames(ref.emb)]
  
  ## Find nearest neighbors
  nn.ranked <- NNHelper(data = ref.emb, query = query.emb, k = k+1, method = "rann")
  
  ## Transfer labels
  if (is.character(ref.labels) || is.factor(ref.labels))) {
    # Categorical labels - majority vote
    results <- pbapply::pblapply(1:nrow(nn.ranked$nn.idx), function(i) {
      xx <- nn.ranked$nn.idx[i, ]
      n.celltype <- sort(table(ref.labels[xx]), decreasing = TRUE)
      most.voted <- n.celltype[1]
      data.frame(
        labels = names(most.voted),
        votes = most.voted,
        perc = most.voted / length(xx)
      )
    })
  } else if (is.numeric(ref.labels))) {
    # Numeric labels - median value
    results <- pbapply::pblapply(1:nrow(nn.ranked$nn.idx), function(i) {
      xx <- nn.ranked$nn.idx[i, ]
      votes <- sum(is.na(ref.labels[xx]))
      median.val <- median(ref.labels[xx], na.rm = T)
      data.frame(
        labels = median.val,
        votes = votes,
        perc = votes / length(xx)
      )
    })
  } else {
    stop("ref.labels should be factor, character, or numeric.")
  }
  
  results <- do.call(rbind, results)
  rownames(results) <- nn.ranked$query.cell.names
  
  return(results)
}