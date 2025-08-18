#' The CellProject Class
#'
#' A class to store single-cell projection results, extending basic Bioconductor 
#' ExpressionSet functionality. Stores embedding coordinates, original data, 
#' and metadata for projected cells.
#'
#' @slot embeddings Matrix containing predicted embedding coordinates (cells x dimensions)
#' @slot refined.embeddings Matrix containing refined embedding coordinates after quality control
#' @slot data data.frame containing original feature matrix used for projection
#' @slot cellmeta data.frame containing cell metadata and quality metrics
#' @slot neighbors \code{Neighbor} object storing nearest neighbor information for quality assessment
#'
#' @name CellProject-class
#' @rdname CellProject-class
#' @concept data_structure
#' @exportClass CellProject
setClass(
  Class = "CellProject",
  slots = c(
    embeddings = "matrix",
    refined.embeddings = "matrix",
    data = "data.frame",
    cellmeta = "data.frame",
    neighbors = "ANY")
)

#' Create a New CellProject Object
#'
#' Constructor function for CellProject class. Validates inputs and initializes
#' a new projection object with embeddings and original data.
#'
#' @param embeddings Matrix of predicted embedding coordinates (cells x dimensions)
#' @param data Original feature matrix (cells x features)
#' @param cellmeta Optional data.frame of cell metadata (default NULL creates empty)
#'
#' @return A new CellProject object
#'
#' @examples
#' \dontrun{
#' embeddings <- matrix(rnorm(200), nrow=100)
#' data <- matrix(rnorm(1000), nrow=100)
#' rownames(embeddings) <- rownames(data) <- paste0("cell_", 1:100)
#' proj <- newCellProject(embeddings, data)
#' }
#'
#' @concept data_structure
#' @export
newCellProject <- function(embeddings, data, cellmeta = NULL){
  # Input validation
  if (missing(embeddings) || missing(data)) {
    stop("Must provide 'embeddings' and 'data'")
  }
  if (is.null(rownames(data))) {
    stop("No cell names (rownames) present in the input data")
  }
  if (is.null(colnames(data))) {
    stop("No features names (colnames) present in the input data")
  }
  
  # Initialize cell metadata if not provided
  if (is.null(cellmeta)) {
    cellmeta <- data.frame(row.names = rownames(data))
  }
  
  # Convert inputs to proper formats
  if (!is.data.frame(data))) {
    data <- as.data.frame(data)
  }
  if (!is.matrix(embeddings))) {
    embeddings <- as.matrix(embeddings)
  }
  
  # Check cell name consistency
  if (!all(rownames(embeddings) == rownames(data))) {
    stop("The rownames of 'embeddings' and 'data' are not matched!")
  }
  
  # Create and return new object
  object <- new(
    Class = "CellProject",
    embeddings = embeddings,
    refined.embeddings = new(Class = "matrix"),
    data = data,
    cellmeta = cellmeta,
    neighbors = NULL)
  return(object)
}

#' Merge Two CellProject Objects
#'
#' Combines embeddings, data, and metadata from two CellProject objects.
#' Handles duplicate column names in metadata by dropping duplicates.
#'
#' @param x First CellProject object
#' @param y Second CellProject object
#' @param ... Additional arguments (currently unused)
#'
#' @return A merged CellProject object
#'
#' @method merge CellProject
#' @concept data_manipulation
#' @export
merge.CellProject <- function(x, y, ...) {
  # Combine embeddings and data
  embeddings <- rbind(x@embeddings, y@embeddings)
  data <- rbind(x@data, y@data)
  cellmeta <- rbind(x@cellmeta, y@cellmeta)

  # Handle duplicate column names in metadata
  check.dup.cn <- duplicated(colnames(cellmeta))
  if (any(check.dup.cn))) {
    warning("Duplicated colnames in 'cellmeta' will auto dropped.")
    cellmeta <- cellmeta[!check.dup.cn]
  }
  
  # Create merged object
  merged.object <- newCellProject(
    embeddings = embeddings,
    data = data,
    cellmeta = cellmeta)
  
  # Combine refined embeddings if present
  if (nrow(x@refined.embeddings) > 0 &&
      nrow(y@refined.embeddings) > 0 &&
      ncol(x@refined.embeddings) == ncol(y@refined.embeddings))) {
    refined.embeddings <- rbind(x@refined.embeddings, y@refined.embeddings)
    merged.object@refined.embeddings <- refined.embeddings
  }
  
  return(merged.object)
}

#' Subset a CellProject Object
#'
#' Creates a new CellProject object containing only specified cells.
#'
#' @param x A CellProject object
#' @param cells Vector of cell names to keep
#'
#' @return A subsetted CellProject object
#'
#' @concept data_manipulation
subset2 <- function(x, cells) {
  if (missing(x) || missing(cells))) {
    stop("Must provide 'x' and 'cells'")
  }
  
  # Create new object with subset data
  object <- newCellProject(
    embeddings = x@embeddings[cells, ],
    data = x@data[cells, ],
    cellmeta = x@cellmeta[cells, , drop=FALSE]
  )
  
  # Subset refined embeddings if present
  if (nrow(x@refined.embeddings) > 0) {
    object@refined.embeddings <- x@refined.embeddings[cells, ]
  }
  
  return(object)
}

#' Split CellProject Object by Metadata
#'
#' Divides a CellProject object into a list of objects based on a metadata field.
#'
#' @param x CellProject object to split
#' @param split.by Metadata column name to split by
#'
#' @return A list of CellProject objects, one per level of split.by
#'
#' @examples
#' \dontrun{
#' # Assuming proj has "batch" column in cellmeta:
#' batch_list <- SplitCellProject(proj, split.by = "batch")
#' }
#'
#' @concept data_manipulation
#' @export
SplitCellProject <- function(x, split.by) {
  if (missing(x) || missing(split.by))) {
    stop("Must provide 'x' and 'split.by'")
  }
  
  cellmeta <- x@cellmeta
  if (!split.by %in% colnames(cellmeta))) {
    stop(paste(split.by, "is not found in 'cellmeta'"))
  }
  
  if (!is.factor(cellmeta[[split.by]]) && !is.character(cellmeta[[split.by]])) {
    stop(paste(split.by, "field is not a factor or character."))
  }
  
  # Split by each unique value in the metadata column
  factor.levels <- unique(cellmeta[[split.by]])
  object.list <- lapply(factor.levels, function(xx) {
    subset2(x, cells = rownames(subset(cellmeta, get(split.by) == xx)))
  })
  names(object.list) <- factor.levels
  
  return(object.list)
}

# Show method for CellProject objects
setMethod(
  f = "show",
  signature = "CellProject",
  definition = function(object) {
    cat("An object of class", class(x = object), "\n")
    cat("@data: ", ncol(object@data), "features across", nrow(object@data), "cells.\n")
    cat("@embeddings: ", "(", nrow(object@embeddings), ",", ncol(object@embeddings), ").\n")
    cat("@refined.embeddings: ", "(", nrow(object@refined.embeddings), ",", ncol(object@refined.embeddings), ").\n")
    cat("@cellmeta: ", "(", nrow(object@cellmeta), ",", ncol(object@cellmeta), ").\n")
    nn.dim <- dim(object@neighbors$nn.idx)
    if (is.null(nn.dim))) {
      cat("@neighbors: NULL\n")
    } else {
      cat("@neighbors: ", nn.dim[2], "nearest neighbors graph on", nn.dim[1], "cells.\n")
    }
  }
)