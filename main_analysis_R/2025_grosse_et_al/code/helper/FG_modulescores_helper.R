#' Compat version of FindMarkers replicating the output from Seurat V4
#' for our analysis (including LogFC estimates like in V4)
#'
#' In v5, Seurat changed the way FoldChange and FindMarkers work,
#' essentially applying pseudo-counts after aggregation, which
#' seems to be beneficial for single-cell datasets but not so much
#' for our low-res spatial data. FindMarkers has a compat-version when
#' setting use.test = "wilcox_limma", but there is currently no compatibility
#' setting for FoldChange, which is included in the FindMarker output.
#'
#' This function is not general, it skips most checks from the normal
#' FindMarkers function and only works for this specific case.
#'
#' @param object Seurat object to calculate markers from, Idents(object) should contain the classes you care about.
#' @param features (names) list of genesets. If it has names and 'names' is not explicitly set, the names from the list are used.
#' @param names defaults to `NULL`. Either a character vector of length matching 'features', or a single character prefix, which is used for all genesets with an ascending number added. Overwrites names of the 'features' argument if not `NULL`.
#' @param assay Name of assay to use, defaults to `DefaultAssay(object)`
#' @param layer layer to use, usually (and defaults to) 'data'
#' @param search defaults to `FALSE`. Whether or not `UpdateSymbolList` should be used if any input genes are not in the object.
#' @param noise.multiplier multiplier for noise added before the cut_number function. Defaults to `.Machine$double.eps/2`
#' @inheritDotParams Seurat::AddModuleScore -object -k -features -name -kmeans.obj -assay
#'
#' @returns df containing logFC estimates and marker p.vals, see FindMarkers
#'
 customV4_FindMarkers <- function (object, ident.1, ident.2 = NULL, assay = NULL, slot = "data",
                                   min.pct = 0.01, logfc.threshold = 0.1, min.diff.pct = -Inf,
                                   features = NULL, pseudocount.use = 1, fc.slot = "data",
                                   verbose = F, only.pos = F, test.use = "wilcox_limma", ...)
 {
   # this part is from FindMarkers.Seurat
   assay = assay %||% DefaultAssay(object)
   assay.use = object[[assay]]
   cellnames.use = colnames(assay.use)
   cells <- Seurat:::IdentsToCells(
     object = object,
     ident.1 = ident.1,
     ident.2 = ident.2,
     cellnames.use = cellnames.use
   )
   cells.1 = cells$cells.1
   cells.2 = cells$cells.2
   norm.method = "LogNormalize" # always LogNormalize for us

   # this part is from FindMarkers.Assay
   data.slot <- ifelse(
     test = test.use %in% Seurat:::DEmethods_counts(),
     yes = "counts",
     no = slot
   )
   data.use <-  SeuratObject::GetAssayData(object = assay.use, layer = data.slot)
   # this function is reverted back to Seurat V4 behaviour here
   fc.results <- V4_FoldChangeAssay(
     object = assay.use,
     slot = fc.slot,
     cells.1 = cells.1,
     cells.2 = cells.2,
     features = features,
     pseudocount.use = pseudocount.use,
     norm.method = norm.method
   )

   # after calculating the compat fc.results, the rest can be delegated back
   # to the seurat method using the compat version "wilcox_limma" for a test
   de.results <- Seurat::FindMarkers(
     object = data.use,
     cells.1 = cells.1,
     cells.2 = cells.2,
     features = features,
     test.use = test.use,
     fc.results = fc.results,
     ...
   )
   return(de.results)
 }

 # Seurat V4.4.0 FoldChange.Assay function (renamed, copied https://github.com/satijalab/seurat/blob/d9f09de15ddf05fe89a8b16eaa100e3720ee122b/R/differential_expression.R#L1090)
 V4_FoldChangeAssay <- function(
    object,
    cells.1,
    cells.2,
    features = NULL,
    slot = "data",
    pseudocount.use = 1,
    fc.name = NULL,
    mean.fxn = NULL,
    base = 2,
    norm.method = NULL,
    ...
 ) {
   pseudocount.use <- pseudocount.use %||% 1
   data <- GetAssayData(object = object, slot = slot)
   # By default run as if LogNormalize is done
   log1pdata.mean.fxn <- function(x) {
     return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
   }
   scaledata.mean.fxn <- rowMeans
   counts.mean.fxn <- function(x) {
     return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
   }
   if (!is.null(x = norm.method)) {
     # For anything apart from log normalization set to rowMeans
     if (norm.method!="LogNormalize") {
       new.mean.fxn <- counts.mean.fxn
     } else {
       new.mean.fxn <- counts.mean.fxn
       if (slot == "data") {
         new.mean.fxn <- log1pdata.mean.fxn
       }  else if (slot == "scale.data") {
         new.mean.fxn <- scaledata.mean.fxn
       }
     }
   } else {
     # If no normalization method is passed use slots to decide mean function
     new.mean.fxn <- switch(
       EXPR = slot,
       'data' = log1pdata.mean.fxn,
       'scale.data' = scaledata.mean.fxn,
       'counts' = counts.mean.fxn,
       log1pdata.mean.fxn
     )
   }
   mean.fxn <- mean.fxn %||% new.mean.fxn
   # Omit the decimal value of e from the column name if base == exp(1)
   base.text <- ifelse(
     test = base == exp(1),
     yes = "",
     no = base
   )
   fc.name <- fc.name %||% ifelse(
     test = slot == "scale.data",
     yes = "avg_diff",
     no = paste0("avg_log", base.text, "FC")
   )
   FoldChange(
     object = data,
     cells.1 = cells.1,
     cells.2 = cells.2,
     features = features,
     mean.fxn = mean.fxn,
     fc.name = fc.name
   )
 }


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# update geneset names if necessary, check against feature names in all samples
# + module score adding function
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# there is probably a smarter way to do this
update_genes = function(gs, known.genes){
  known.inds = gs %in% known.genes
  if (any(!known.inds)) {
    known.part = gs[known.inds]
    unknown.part = gs[!known.inds]
    updated.part = UpdateSymbolList(symbols = unknown.part, timeout = 30)
    return(unique(c(known.part, updated.part)))
  } else {
    return(gs)
  }
}


#' Seurat's AddModuleScore function with fixes for cut_number issues
#'
#' Fixes issues with the computation of bins for background genes. In certain
#' cases, one of two errors may be seen with the original function:
#' - binning fails entirely when cut_number() throws an "Insufficient data values to produce XX bins." error
#' - binning may succeed but one of the bins may end up with less than #ctrl number of genes such that the sampling from this group fails afterwards
#'
#' IMHO the problem is, that the Seurat authors correctly foresaw problems with
#' many genes having very small average expression values and thus (correctly)
#' decided to add a little bit of noise to make the cut_number function work
#' more reliably in this case. Unfortunately, they decided to make the noise too
#' small of magnitude (.../1e30), which is way below the machine epsilon of R at
#' the moment: (see .Machine$double.eps) 2.220446049250313080847e-16 This means
#' that, in bad cases, R is only accurate to about 16 digits. Thereby, the added
#' noise will get rounded away as soon as it is added, nullifying it's effect.
#'
#' @param object Seurat object to calculate module scores on
#' @param features (names) list of genesets. If it has names and 'names' is not explicitly set, the names from the list are used.
#' @param names defaults to `NULL`. Either a character vector of length matching 'features', or a single character prefix, which is used for all genesets with an ascending number added. Overwrites names of the 'features' argument if not `NULL`.
#' @param assay Name of assay to use, defaults to `DefaultAssay(object)`
#' @param layer layer to use, usually (and defaults to) 'data'
#' @param search defaults to `FALSE`. Whether or not `UpdateSymbolList` should be used if any input genes are not in the object.
#' @param noise.multiplier multiplier for noise added before the cut_number function. Defaults to `.Machine$double.eps/2`
#' @inheritDotParams Seurat::AddModuleScore -object -k -features -name -kmeans.obj -assay
#'
#' @returns 'object' with module scores added to metadata
#'
#' @examples
#' # For a demonstration, see: (22 is the maximum number of printable digits in base R)
#'  # Assume we have a lot of genes with 1 count across 10 million cells.
#'  # We can easily show that all those genes end up with exactly the same
#'  # value entering the cut_number function, despite the noise being added.
#' print( rep(1, 100)/1e7 + (rnorm(100)/1e30) , digits = 22)
#'
#' unique(rep(1, 100)/1e7 + (rnorm(100)/1e30))
#' # [1] 9.999999999999999547481e-08
#'
#' # which shows that there is only one unique value despite the added noise
#' # For a fix, we can either use something like .Machine$double.eps as a
#' # multiplier for the noise, or, if somebody is concerned about this being
#' # not portable or too generous, take the smallest value greater than zero and
#' # then divide by 100 or 1000 for a noise multiplier.
#' # I decided to go the easy way and use the machine precision cutoff as a minimum noise multiplier
fixed_AddModuleScore <- function (object, features, names = NULL, layer = "data",
                                  pool = NULL, nbin = 24, ctrl = 100, assay = NULL,
                                  seed = 1, search = FALSE, noise.multiplier = .Machine$double.eps/2,
                                  ...)
{
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- LayerData(object = object, layer = layer)
  features.old <- features

  if (is.null(x = features)) {
    stop("Missing input feature list")
  }
  if (!(is.list(x = features))) {
    stop("`features` should be a list")
  }
  # resolve names
  .featuresnamed = !(is.null(names(features)))
  .explicitlynamed = !(is.null(names))
  .names = if (.explicitlynamed) {
    if (length(names)==length(features)) {
      names
    } else {
      if (length(names)==1) {
        paste0(names, 1:length(features))
      } else {
        stop("Length of `names` must be 1 or equal to the length of `features`.")
      }
    }
  } else {
    if (.featuresnamed) { # replaces empty names with UnnamedFeatureX dummies to avoid problems
      ifelse(names(features)=="",
             paste0("UnnamedFeature", 1:length(features)),
             names(features))
    } else {
      stop("Either features should be a named list or names explicitly given.")
    }
  }

  # the rest is practically copied from Seurat::AddModuleScore
  features <- lapply(X = features, FUN = function(x) {
    missing.features <- setdiff(x = x, y = rownames(x = object))
    if (length(x = missing.features) > 0) {
      warning("The following features are not present in the object: ",
              paste(missing.features, collapse = ", "), ifelse(test = search,
                                                               yes = ", attempting to find updated synonyms",
                                                               no = ", not searching for symbol synonyms"),
              call. = FALSE, immediate. = TRUE)
      if (search) {
        tryCatch(expr = {
          updated.features <- UpdateSymbolList(symbols = missing.features,
                                               ...)
          names(x = updated.features) <- missing.features
          for (miss in names(x = updated.features)) {
            index <- which(x == miss)
            x[index] <- updated.features[miss]
          }
        }, error = function(...) {
          warning("Could not reach HGNC's gene names database",
                  call. = FALSE, immediate. = TRUE)
        })
        missing.features <- setdiff(x = x, y = rownames(x = object))
        if (length(x = missing.features) > 0) {
          warning("The following features are still not present in the object: ",
                  paste(missing.features, collapse = ", "),
                  call. = FALSE, immediate. = TRUE)
        }
      }
    }
    return(intersect(x = x, y = rownames(x = object)))
  })
  cluster.length <- length(x = features)

  if (!all(Seurat:::LengthCheck(values = features))) {
    warning(paste("Could not find enough features in the object from the following feature lists:",
                  paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
                  "Attempting to match case..."))
    features <- lapply(X = features.old, FUN = CaseMatch,
                       match = rownames(x = object))
  }
  if (!all(Seurat:::LengthCheck(values = features))) {
    stop(paste("The following feature lists do not have enough features present in the object:",
               paste(names(x = which(x = !Seurat:::LengthCheck(values = features)))),
               "exiting..."))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))*noise.multiplier,
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut ==
                                                                              data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use),
                        ncol = ncol(x = object))
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use,
    ])
  }
  features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length,
                            ncol = ncol(x = object))
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- .names # insert proper names
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}


# performs additional checks and does not append numbers to names, make sure that names are unique!
add_ms_with_checks = function(se, features, names, ...){
  # need to catch case that no genes of a score are in the sample, set score value for all spots to 0 then
  score.computable = unlist(lapply(features, function(x)length(intersect(x, rownames(se)))>0), use.names = F)
  if (any(!score.computable)) {
    not.computable.features = features[!score.computable]
    not.computable.names = names[!score.computable]
    for (nm in not.computable.names) {
      se = AddMetaData(se, metadata = rep(0., ncol(se)), col.name = nm)
    }
    features = features[score.computable]
    names = names[score.computable]
  }
  se = fixed_AddModuleScore(se, features = features, names = names, ...)
  se
}


# do_hypergeometric_ora_tests <- function(se, genesets, tops.pct = 0.05, assay = "Spatial"){
#   # logFCs with Means calculated in non-log space
#   # same thing that Seurat does for FoldChange but with base natural log
#   allgenes = rownames(se)
#   lfcs = sweep(se@assays[[assay]]@data, MARGIN = 1,
#                STATS = log1p(rowMeans(expm1(se@assays[[assay]]@data))))
#   ntops = floor(length(lfcs[,1])*tops.pct)
#
#   enrich.pvals = bplapply(Cells(se), function(c){
#     tops = names(sort(lfcs[,c], decreasing = T)[1:ntops])
#     eres = fgsea::fora(pathways = genesets, genes = tops, universe = allgenes)
#     return(setNames(object = eres$pval, nm = eres$pathway))
#   })
#   enrich.pvals = simplify2array(enrich.pvals)
# }

## old stuff for testing
##
# my_cut_number <- function(x, n = NULL, ...){
# probs <- seq(0, 1, length.out = n + 1)
# brk <- stats::quantile(x, probs, na.rm = TRUE)
# if (anyDuplicated(brk))
#   abort(glue("Insufficient data values to produce {n} bins."))
# cut(x, brk, include.lowest = TRUE, ...)
# }
#
# plot_stuff <- function(data.avg, data.cut = NULL, ctrl.use = NULL, features = NULL){
#   # Plot the bins that have been created to split genes based on their average expression
#   plot(data.avg, pch=16, ylab="Average expression across all cells", xlab="All genes, ranked")
#   if (!(is.null(data.cut))){
#     nbin = length(unique(data.cut))
#     for(i in unique(data.cut)){
#       cut_pos <- which(data.cut==i)
#       col = rainbow(nbin)[i]
#       points(cut_pos, data.avg[cut_pos], pch=20, col=col)
#       if(i%%2==0){
#         rect(xleft = cut_pos[1], ybottom = min(data.avg), xright = cut_pos[length(cut_pos)], ytop = max(data.avg), col=scales::alpha("grey", 0.3))
#       } else {
#         rect(xleft = cut_pos[1], ybottom = min(data.avg), xright = cut_pos[length(cut_pos)], ytop = max(data.avg), col=scales::alpha("white", 0.3))
#       }
#     }
#   }
#
#   if (!(is.null(ctrl.use))){
#     # Add red points for selected control genes
#     points(which(names(data.avg)%in%ctrl.use[[1]]), data.avg[which(names(data.avg)%in%ctrl.use[[1]])], pch=16, col="red")
#   }
#
#   if (!(is.null(features))){
#     # Add blue points for genes in the input gene list
#     points(which(names(data.avg)%in%features[[1]]), data.avg[which(names(data.avg)%in%features[[1]])], pch=16, col="blue")
#   }
#
#   # Add a legend
#   legend(x = "topleft",
#          legend = c("gene", "selected control gene", "gene in geneset"),
#          col = c("black", "red", "blue"),
#          pch = 16)
# }
