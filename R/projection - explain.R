#' @include utils.R
#' @include objects.R
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Regression Model ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Training an ensemble SVM model for embedding coordinates regression
#' 
#' This function trains an ensemble of SVM models to predict embedding coordinates from feature matrices.
#' It performs data normalization, creates multiple regression tasks, and trains models in parallel.
#'
#' @import mlr3verse
#' @importFrom mlr3 as_task_regr lrn
#' @importFrom parallel detectCores
#' @importFrom pbapply pblapply
#' @importFrom stats setNames
#'
#' @param feature.mat A signature score matrix where rows are cells and columns are features.
#'        Should be coercible to data.frame.
#' @param emb.mat The embedding matrix to predict, with rows as cells and columns as dimensions.
#'        Must have matching row names with feature.mat.
#' @param cell.types Optional named vector recording cell types for balanced sampling.
#'        Names should match row names of feature.mat and emb.mat.
#' @param do.norm Normalization method for feature matrix: 'L1', 'L2', or NULL to skip.
#'        L1 normalization divides by row sums, L2 by row Euclidean norms. Default: 'L2'
#' @param batch.size Number of cells to sample for each model training. Default: 5000
#' @param n.models Number of SVM models to train in the ensemble. Default: 100
#' @param balance.cell.type Whether to perform balanced sampling by cell type. 
#'        Requires cell.types parameter. Default: FALSE
#' @param cores Number of CPU cores for parallel training. -1 uses all available. Default: -1
#'
#' @return A list with class "Regression" containing:
#' \itemize{
#'   \item model - List of trained learners (one list per coordinate)
#'   \item features - Vector of feature names used in training
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' feature_matrix <- as.data.frame(matrix(rnorm(1000), nrow=100)
#' emb_matrix <- as.data.frame(matrix(rnorm(200), nrow=100))
#' rownames(feature_matrix) <- rownames(emb_matrix) <- paste0("cell_", 1:100)
#' model <- FitEnsembleSVM(feature.mat = feature_matrix, 
#'                        emb.mat = emb_matrix,
#'                        n.models = 10)
#' }
#'
#' @concept training_model
#' @export
FitEnsembleSVM <- function(feature.mat, emb.mat, cell.types=NULL, do.norm='L2', batch.size=5000, n.models=100, balance.cell.type=FALSE, cores=-1) {
  ## Input validation and preprocessing
  if (!is.data.frame(feature.mat)) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  
  # Check for cell name consistency between feature and embedding matrices
  outliers.1 <- setdiff(rownames(feature.mat), rownames(emb.mat))
  outliers.2 <- setdiff(rownames(emb.mat), rownames(feature.mat))
  if (length(outliers.1)) {
    stop("Some cells in feature.mat are not in emb.mat.")
  }
  if (length(outliers.2)) {
    stop("Some cells in emb.mat are not in feature.mat.")
  }
  
  # Validate batch size
  if (batch.size >= nrow(feature.mat)) {
    stop("Batch size is larger than the samples. A smaller 'batch.size' is expected.")
  }
  
  # Set up parallel processing
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  
  # Clean column names (replace hyphens with underscores)
  bad.cols <- grepl("-", colnames(feature.mat))
  if (any(bad.cols)) {
    warning("Bad colnames in 'feature.mat'. Change the '-' to '_'.")
    colnames(feature.mat) <- gsub("-", "_", colnames(feature.mat))
  }
  bad.cols <- grepl("-", colnames(emb.mat))
  if (any(bad.cols)) {
    warning("Bad colnames in 'emb.mat'. Change the '-' to '_'.")
    colnames(emb.mat) <- gsub("-", "_", colnames(emb.mat))
  }
  
  ## Normalization
  if (is.null(do.norm)) {
    message("Skip normalization ...")
  } else if (do.norm == "L1")) {
    message("L1 normalization ...")
    feature.mat <- L1norm(feature.mat)
  } else if (do.norm == "L2")) {
    message("L2 normalization ...")
    feature.mat <- L2norm(feature.mat)
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  
  ## Sampling probability calculation for balanced sampling
  .getP <- function(labels) {
    labels.freq <- as.vector(table(labels))
    names(labels.freq) <- names(table(labels))
    1 / (labels.freq[labels] * length(labels.freq))
  }
  
  ## Create training datasets
  message(sprintf("Creating regression tasks: size=%s, n=%s", batch.size, n=n.models))
  if (balance.cell.type & !is.null(cell.types))) {
    message("Balanced sampling.")
    p <- .getP(cell.types)
  } else {
    p <- NULL
  }
  
  # Generate multiple training datasets by sampling
  train.datasets <- pbapply::pblapply(1:n.models, function(i) {
    row.ids <- sample(rownames(feature.mat), size = batch.size, replace = F, prob = p)
    as.data.frame(cbind(feature.mat[row.ids, , drop=FALSE], emb.mat[row.ids, , drop=FALSE]))
  })
  
  ## Create regression tasks for each coordinate
  vars <- colnames(feature.mat)
  targets <- colnames(emb.mat)
  tasks <- lapply(targets, function(coord_i) {
    pbapply::pblapply(train.datasets, function(data_j) {
      as_task_regr(data_j[, c(vars, coord_i)], target = coord_i)
    })
  })
  names(tasks) <- targets
  
  ## Training function for individual models
  .train <- function(dataset, batch.id) {
    task = dataset[[batch.id]]
    learner = lrn("regr.svm")
    learner$train(task)
    learner
  }
  
  ## Train models in parallel
  message("Regression ...")
  coords_lrns = lapply(tasks, function(coord_i) {
    parallel::mclapply(1:n.models, function(j) .train(coord_i, j), mc.cores=cores)
  })
  names(coords_lrns) = names(tasks)
  
  ## Return structured model object
  model <- list(
    model = coords_lrns,
    features = colnames(feature.mat)
  )
  class(model) <- "Regression"
  return(model)
}


#' Predict embedding coordinates using trained SVM ensemble
#'
#' Projects new data into the embedding space using a trained ensemble SVM model.
#' Handles feature matching, normalization, and coordinate prediction integration.
#'
#' @importFrom parallel detectCores
#' @importFrom pbapply pblapply
#' @importFrom stats median
#'
#' @param feature.mat A data.frame containing signature scores (rows=cells, cols=features).
#'        Will be coerced to data.frame if not already.
#' @param model Trained model object returned by FitEnsembleSVM().
#' @param do.norm Normalization method: 'L1', 'L2', or NULL. Must match training. Default: 'L2'
#' @param int.fun Function to integrate predictions across models: mean, median, or custom.
#'        Should accept a vector and return a scalar. Default: stats::median
#' @param cores Number of CPU cores for parallel prediction. -1 uses all. Default: -1
#'
#' @return A \code{CellProject} object containing:
#' \itemize{
#'   \item embeddings - Predicted coordinates matrix
#'   \item data - Original feature matrix (after processing)
#' }
#'
#' @examples
#' \dontrun{
#' # After training with FitEnsembleSVM():
#' new_features <- as.data.frame(matrix(rnorm(500), nrow=50))
#' rownames(new_features) <- paste0("new_cell_", 1:50)
#' projection <- ProjectNewdata(feature.mat = new_features, 
#'                            model = trained_model)
#' }
#'
#' @concept reference_mapping
#' @export
ProjectNewdata <- function(feature.mat, model, do.norm='L2', int.fun=stats::median, cores=-1) {
  ## Input validation and preprocessing
  if (!is.data.frame(feature.mat))) {
    warning("The 'feature.mat' should be a data.frame, enforce convert to data.frame.")
    feature.mat <- as.data.frame(feature.mat)
  }
  
  # Handle feature mismatch between input and model
  bad.features <- setdiff(colnames(feature.mat), model$features)
  if (length(bad.features) > 0) {
    warning(sprintf("These features are not in the model: %s. Dropping them (it).", 
                   paste(bad.features, collapse = " ")))
    feature.mat <- feature.mat[, !(colnames(feature.mat) %in% bad.features), drop = FALSE]
  }
  
  # Add missing features (filled with zeros)
  no.features <- setdiff(model$features, colnames(feature.mat))
  if (length(no.features) > 0) {
    warning(sprintf("These features are not in your inputs: %s\n. Filling with zeros.", 
                   paste(no.features, collapse = " ")))
    zeros <- matrix(rep(0, nrow(feature.mat)*length(no.features)), 
                  ncol = length(no.features)) %>% as.data.frame()
    colnames(zeros) <- no.features
    newdata <- cbind(feature.mat, zeros)[, model$features]
  } else {
    newdata <- feature.mat
  }
  
  # Set up parallel processing
  cores <- ifelse(cores > 0, cores, parallel::detectCores())
  
  ## Normalization (must match training normalization)
  if (is.null(do.norm))) {
    message("Skip normalization ...")
  } else if (do.norm == "L1")) {
    message("L1 normalization ...")
    newdata <- L1norm(newdata)
  } else if (do.norm == "L2")) {
    message("L2 normalization ...")
    newdata <- L2norm(newdata)
  } else {
    stop(sprintf("Undefined 'do.norm': %s.", do.norm))
  }
  
  ## Parallel prediction across all models
  n.models <- length(model$model[[1]])
  coords_pred <- lapply(model$model, function(lrns) {
    coords_i <- parallel::mclapply(
      X = 1:n.models,
      FUN = function(i) {
        lrns[[i]]$predict_newdata(newdata = as.data.frame(newdata))$data$response
      }, mc.cores=cores)
    do.call(cbind, coords_i)
  })
  names(coords_pred) = names(model$model)
  
  ## Integrate predictions across models
  coords_inter <- lapply(names(coords_pred), function(xx) {
    apply(coords_pred[[xx]], 1, int.fun)
  })
  names(coords_inter) <- names(coords_pred)
  coords_inter <- as.data.frame(coords_inter)
  rownames(coords_inter) <- rownames(feature.mat)
  
  ## Return CellProject object
  newCellProject(embeddings = coords_inter, data = feature.mat)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Mapping Quality ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Assess projection quality using k-NN consistency
#' 
#' Computes quality metrics for projected cells by comparing neighborhood consistency
#' between feature space and embedding space. Uses random background distribution
#' to calculate p-values for mapping quality.
#'
#' @importFrom pbapply pbsapply
#' @importFrom stats dist p.adjust
#'
#' @param object A \code{CellProject} object containing embeddings and original data
#' @param k Number of nearest neighbors to use for quality assessment. Default: 20
#' @param repeats Number of random samples for background distribution. Default: 1e4
#'
#' @return Modified \code{CellProject} object with quality metrics added to cellmeta:
#' \itemize{
#'   \item mean.knn.dist - Mean distance to k nearest neighbors in embedding space
#'   \item p.val - Empirical p-value for mapping quality
#'   \item p.adj - BH-adjusted p-values
#' }
#'
#' @references \url{https://www.nature.com/articles/s41467-021-25957-x}
#' @concept reference_mapping
#' @export
AddProjQual <- function(object, k=20, repeats=1e4) {
  pred.emb <- object@embeddings
  golden.emb <- object@data

  # Ensure consistent ordering
  if (!all(rownames(golden.emb) == rownames(pred.emb))) {
    pred.emb <- pred.emb[rownames(golden.emb), ]
  }

  # Build k-NN graph in original feature space
  message("Building k-NN graph in feature space ...")
  nn.golden = NNHelper(golden.emb, k = k, method = "sklearn", metric = "cosine")

  # Compute distances in predicted embedding space
  message("Calculating euclidean metric on pred.emb ...")
  pred.dist <- as.matrix(stats::dist(pred.emb))

  # Calculate quality metrics
  message("Calculating cell-based mapping quality metrics ...")
  
  # Mean distance to k-NNs in embedding space
  mean.knn.dist <- pbapply::pbsapply(1:nrow(pred.emb), function(i) {
    mean(pred.dist[i, nn.golden$nn.idx[i, ]])
  })
  names(mean.knn.dist) <- rownames(pred.emb)

  # Generate background distribution using random cells
  mean.krand.dist <- pbapply::pbsapply(1:repeats, function(i) {
    rand.cells <- sample(1:ncol(pred.dist), size = k+1, replace = F)
    mean(pred.dist[rand.cells[1], rand.cells[-1]])
  })
  
  # Calculate empirical p-values
  p.val <- pbapply::pbsapply(mean.knn.dist, function(q) {
    sum(q > mean.krand.dist) / length(mean.krand.dist) 
  })
  p.adj <- stats::p.adjust(p.val, method = "BH")
  
  # Store results
  results <- list(
    mean.knn.dist = mean.knn.dist,
    p.val = p.val,
    p.adj = p.adj
  )
  results <- as.data.frame(results)
  
  # Update object
  object@neighbors <- nn.golden
  cellmeta <- object@cellmeta
  
  # Handle potential column name conflicts
  res.names <- colnames(results)
  kept.names <- setdiff(colnames(cellmeta), res.names)
  cellmeta <- cellmeta[, kept.names, drop=FALSE]
  object@cellmeta <- cbind(cellmeta, results[rownames(cellmeta), ])
  
  return(object)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Refine Mapping Results ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Refine low-confidence mapping coordinates
#' 
#' Interpolates coordinates for poorly mapped cells based on their high-confidence
#' neighbors in feature space. Uses Gaussian kernel weighting for smooth interpolation.
#'
#' @param object A \code{CellProject} object with quality metrics (from AddProjQual)
#' @param p.val.cutoff Max p-value for cells to be considered high-quality.
#'        Cells above this will be refined. Default: 1 (refine all)
#' @param p.adj.cutoff Adjusted p-value cutoff for refinement. Default: 0.05
#' @param k Number of nearest neighbors to use for interpolation. Default: 10
#'
#' @return Modified \code{CellProject} object with refined embeddings stored in
#'         the refined.embeddings slot.
#'
#' @concept reference_mapping
#' @export
RefineProjection <- function(object, p.val.cutoff = 1, p.adj.cutoff = 0.05, k = 10){
  # Input validation
  if (is.null(object@neighbors))) {
    stop("No neighbors object provided. Run AddProjQual() first.")
  }
  
  # Extract data from object
  data <- object@data
  nn.dists <- object@neighbors$nn.dists
  query.cell.names <- object@neighbors$query.cell.names
  rownames(nn.dists) <- query.cell.names
  embeddings <- object@embeddings
  cellmeta <- object@cellmeta
  refined.embeddings <- embeddings
  
  # Validate k parameter
  if (k > ncol(nn.dists))) {
    warning("K =", k, "was too big, reset to", ncol(nn.dists)), "according to the built knn graph.")
    k <- ncol(nn.dists)
  }
  
  # Identify low-quality mappings
  sel.rows <- cellmeta$p.val > p.val.cutoff | cellmeta$p.adj > p.adj.cutoff
  low.qual.map <- cellmeta[sel.rows, ]
  
  if (nrow(low.qual.map))) {
    cat("There are", nrow(low.qual.map)), "cell(s) with low quality mapping coordinates under p <=", 
        p.val.cutoff, "and p.adj <=", p.adj.cutoff, "\n")
    
    ## Refinement process
    low.qual.cells <- rownames(low.qual.map)
    other.cells <- setdiff(rownames(cellmeta), low.qual.cells)
    
    # 1) Compute cosine distance matrix between low-quality and other cells
    D <- cosine.dist(as.matrix(data))[low.qual.cells, other.cells]
    
    # 2) Calculate adaptive bandwidth (sigma) for each cell
    sigma <- nn.dists[low.qual.cells, k]
    sigma <- matrix(rep(sigma, each=length(other.cells)), 
                   ncol = length(other.cells), byrow = T)
    
    # 3) Compute Gaussian affinity matrix
    A <- exp(-(D / sigma)^2)
    
    # 4) Row-normalize to get interpolation weights
    A <- t(apply(A, 1, function(xx) xx/sum(xx) ))
    
    # 5) Perform weighted interpolation
    U <- A %*% embeddings[other.cells, ]
    refined.embeddings[rownames(U), ] <- U
  }
  
  # Store refined embeddings
  object@refined.embeddings <- refined.embeddings
  return(object)
}