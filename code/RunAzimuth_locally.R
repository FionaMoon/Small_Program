library(Matrix)
library(Seurat)

NNTransform <- function(
    object,
    meta.data,
    neighbor.slot = "query_ref.nn",
    key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}


Oxford <- function(..., join = c('and', 'or')) {
  join <- match.arg(arg = join)
  args <- as.character(x = c(...))
  args <- Filter(f = nchar, x = args)
  if (length(x = args) == 1) {
    return(args)
  } else if (length(x = args) == 2) {
    return(paste(args, collapse = paste0(' ', join, ' ')))
  }
  return(paste0(
    paste(args[1:(length(x = args) - 1)], collapse = ', '),
    paste0(', ', join, ' '),
    args[length(x = args)]
  ))
}

LoadReference <- function(path, seconds = 10L) {
  op <- options(Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  ref.names <- list(
    map = 'ref.Rds',
    ann = 'idx.annoy'
  )
  mapref <- file.path(path, ref.names$map)
  annref <- file.path(path, ref.names$ann)
  exists <- file.exists(c(mapref, annref))
    if (!all(exists)) {
      stop(
        "Missing the following files from the directory provided: ",
        Oxford(unlist(x = ref.names)[!exists], join = 'and')
      )
    }
  # Load the map reference
  map <- readRDS(file = mapref)
  
  # handle new parameters in uwot models beginning in v0.1.13
  if (!"num_precomputed_nns" %in% names(Misc(map[["refUMAP"]])$model)) {
    Misc(map[["refUMAP"]], slot="model")$num_precomputed_nns <- 1
  }
  
  # Load the annoy index into the Neighbor object in the neighbors slot
  map[["refdr.annoy.neighbors"]] <- LoadAnnoyIndex(
    object = map[["refdr.annoy.neighbors"]],
    file = annref
  )
  # Validate that reference contains required dims

  # Fix colnames of 'feature.loadings' in reference 
  key.pattern = "^[^_]*_"
  new.colnames <- gsub(pattern = key.pattern, 
                       replacement = Key(map[["refDR"]]), 
                       x = colnames(Loadings(
                         object = map[["refDR"]],
                         projected = FALSE)))
  colnames(Loadings(object = map[["refDR"]], 
                    projected = FALSE)) <- new.colnames
  # Create plotref
  ad <- Tool(object = map, slot = "AzimuthReference")
  plotref.dr <- GetPlotRef(object = ad)
  cm <- sparseMatrix(
    i = 1, j = 1, x = 0, dims = c(1, nrow(x = plotref.dr)),
    dimnames = list("placeholder", Cells(x = plotref.dr))
  )
  op <- options(Seurat.object.assay.version = "v3")
  on.exit(expr = options(op), add = TRUE)
  plot <- CreateSeuratObject(counts = cm)
  plot[["refUMAP"]] <- plotref.dr
  DefaultAssay(plot[["refUMAP"]]) <- DefaultAssay(plot)
  plot <- AddMetaData(object = plot, metadata = Misc(object = plotref.dr, slot = "plot.metadata"))
  gc(verbose = FALSE)
  return(list(
    map = map,
    plot = plot
  ))
}


GetPlotRef <- function(object, ...) {
  plotref <- slot(object = object, name = "plotref")
  # temporary fix for Key of tonsil ref - Ideally we will update the SeuratData object 
  Key(plotref) <- gsub("\\.", "", Key(plotref))
  return(plotref)
}

RunAzimuth.Seurat <- function(
    query,
    reference,
    query.modality = "RNA",
    annotation.levels = NULL,
    umap.name = "ref.umap",
    do.adt = FALSE,
    verbose = TRUE,
    assay = NULL,
    k.weight = 50,
    n.trees = 20,
    mapping.score.k = 100
) {
  assay <- assay %||% DefaultAssay(query) ## 空合并运算符
  ## 如果不是，则空合并运算符返回其左侧操作数的值null；否则，它计算右侧操作数并返回其结果。
  
  ## 
  if (query.modality == "RNA"){
    ## 判断参考文件是否存在
    if (dir.exists(reference)) {
      reference <- LoadReference(reference)$map
    } else {
      stop("Can't find reference")
    }
    
    dims <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
    
    if (isTRUE(do.adt) && !("ADT" %in% Assays(reference))) {
      warning("Cannot impute an ADT assay because the reference does not have antibody data")
      do.adt = FALSE
    }
    
    meta.data <- names(slot(reference, "meta.data"))
    
    # is annotation levels are not specify, gather all levels of annotation
    if (is.null(annotation.levels)) {
      annotation.levels <- names(slot(object = reference, name = "meta.data"))
      annotation.levels <- annotation.levels[!grepl(pattern = "^nCount", x = annotation.levels)]
      annotation.levels <- annotation.levels[!grepl(pattern = "^nFeature", x = annotation.levels)]
      annotation.levels <- annotation.levels[!grepl(pattern = "^ori", x = annotation.levels)]
    }
    
    # Calculate nCount_RNA and nFeature_RNA if the query does not
    if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))){
      stop("No nCount_RNA and nFeature_RNA")
    } 
    else if (any(grepl(pattern = '^MT-', x = rownames(x = query)))){
      query <- PercentageFeatureSet(
        object = query,
        pattern = '^MT-',
        col.name = 'percent.mt',
        assay = assay
      )
    }
    
    # Find anchors between query and reference （Seurat Function）
    anchors <- FindTransferAnchors(
      reference = reference,
      query = query,
      k.filter = NA,
      reference.neighbors = "refdr.annoy.neighbors",
      reference.assay = "refAssay",
      query.assay = assay,
      reference.reduction = "refDR",
      normalization.method = "SCT",
      features = rownames(Loadings(reference[["refDR"]])),
      dims = 1:dims,
      n.trees = n.trees,
      mapping.score.k = mapping.score.k,
      verbose = verbose
    )
    # LogNormalize / SCT
    # Transferred labels are in metadata columns named "predicted.*"
    # The maximum prediction score is in a metadata column named "predicted.*.score"
    # The prediction scores for each class are in an assay named "prediction.score.*"
    # The imputed assay is named "impADT" if computed
    # annotation.levels is meta.data
    refdata <- lapply(X = annotation.levels, function(x) {
      reference[[x, drop = TRUE]] ## drop=TRUE 取某一列的时候返回的是vector，而不是dataframe
    })
    names(x = refdata) <- annotation.levels
    
    if (isTRUE(do.adt)) {
      refdata[["impADT"]] <- GetAssayData(
        object = reference[["ADT"]],
        slot = "data"
      )
    }
    
    ## Seurat function
    query <- TransferData(
      reference = reference,
      query = query,
      query.assay = assay,
      dims = 1:dims,
      anchorset = anchors,
      refdata = refdata,
      n.trees = 20,
      store.weights = TRUE,
      k.weight = k.weight,
      verbose = verbose
    )
    
    # Calculate the embeddings of the query data on the reference SPCA
    ## Seurat function
    query <- IntegrateEmbeddings(
      anchorset = anchors,
      reference = reference,
      query = query,
      query.assay = assay,
      reductions = "pcaproject",
      reuse.weights.matrix = TRUE,
      verbose = verbose
    )
    
    # Calculate the query neighbors in the reference
    # with respect to the integrated embeddings
    query[["query_ref.nn"]] <- FindNeighbors(
      object = Embeddings(reference[["refDR"]]),
      query = Embeddings(query[["integrated_dr"]]),
      return.neighbor = TRUE,
      l2.norm = TRUE,
      verbose = verbose
    )
    
    # The reference used in the app is downsampled compared to the reference on which
    # the UMAP model was computed. This step, using the helper function NNTransform,
    # corrects the Neighbors to account for the downsampling.
    ## functions in helpers.R
    query <- NNTransform(
      object = query,
      meta.data = reference[[]]
    )
    
    # Project the query to the reference UMAP.
    query[[umap.name]] <- RunUMAP(
      object = query[["query_ref.nn"]],
      reduction.model = reference[["refUMAP"]],
      reduction.key = 'UMAP_',
      verbose = verbose
    )
    
    # Calculate mapping score and add to metadata
    query <- AddMetaData(
      object = query,
      metadata = MappingScore(anchors = anchors, ndim = dims),
      col.name = "mapping.score"
    )
  }
  return(query)
}
