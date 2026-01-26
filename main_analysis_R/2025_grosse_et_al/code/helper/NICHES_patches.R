########################################################################################################################################
# PATCHED VERSION OF NICHES FUNCTION
# because the transformation to Seurat objects inside of NICHES mangles the names of complexes
# by replacing all _ with - and this interferes with genes like HLA-DRA, which contain -
# COMPLEX GENES ARE NOW SEPARATED BY +
# Addition: more memory-friendly computations for large HD samples
# Addition: edgelist_fun can be set to change the function that computes
# the edgelist to something more efficient, for large datasets
########################################################################################################################################
library(NICHES)
# patches as of NICHES v1.1.0

fixed_RunCellToCellSpatial <- function(filtered.obj,
                                 ground.truth,
                                 assay,
                                 meta.data.to.map,
                                 edgelist,
                                 output_format
){


  # Make ligand matrix

  # needed to make this more efficient, do everything without casting to dense matrix if possible
  # maximum number of subunits
  nsubunits = ncol(ground.truth$source.subunits)
  source_subunits_subset = ground.truth$source.subunits
  lig.dats = list()
  source_indices = 1:nrow(source_subunits_subset)
  final_indices = c()
  # for each "size class" of complexes (1 unit, 2 subunits, 3, ...), extract only those columns,
  # aggregate and rbind everything at the end, this should avoid the dense matrix filled with 1s issue
  for (s in nsubunits:1) {
    applicable_complexes = !is.na(source_subunits_subset[,s])
    source_complexes_filt = source_subunits_subset[applicable_complexes,]
    source_subunits_subset = source_subunits_subset[!applicable_complexes,] # remove from list for lower size classes
    subset_indices = source_indices[applicable_complexes]
    source_indices= source_indices[!applicable_complexes] # same for indices for restoring order at the end

    if (s == 1) {
      lig.dat = getSeuratAssay(filtered.obj,assay,"data")[source_complexes_filt[,1],edgelist$from]
    } else {
      subunit.list <- lapply(1:s, function(si){
        getSeuratAssay(filtered.obj,assay,"data")[source_complexes_filt[,si],edgelist$from]
      })
      lig.dat.names = apply(FUN = paste, MARGIN = 1, sapply(subunit.list, rownames), collapse = "_")
      lig.dat = Reduce('*',subunit.list)
      rownames(lig.dat) = lig.dat.names
    }
    lig.dats[[s]] = lig.dat
    final_indices = c(subset_indices, final_indices)
  }
  lig.dat = do.call(rbind, lig.dats)[order(final_indices),] # restore order

  # Make receptor matrix

  nsubunits = ncol(ground.truth$target.subunits)
  target_subunits_subset = ground.truth$target.subunits
  rec.dats = list()
  target_indices = 1:nrow(target_subunits_subset)
  final_indices = c()
  # for each "size class" of complexes (1 unit, 2 subunits, 3, ...), extract only those columns,
  # aggregate and rbind everything at the end, this should avoid the dense matrix filled with 1s issue
  for (s in nsubunits:1) {
    applicable_complexes = !is.na(target_subunits_subset[,s])
    target_complexes_filt = target_subunits_subset[applicable_complexes,]
    target_subunits_subset = target_subunits_subset[!applicable_complexes,] # remove from list for lower size classes
    subset_indices = target_indices[applicable_complexes]
    target_indices = target_indices[!applicable_complexes] # same for indices for restoring order at the end

    if (s == 1) {
      rec.dat = getSeuratAssay(filtered.obj,assay,"data")[target_complexes_filt[,1],edgelist$to]
    } else {
      subunit.list <- lapply(1:s, function(si){
        getSeuratAssay(filtered.obj,assay,"data")[target_complexes_filt[,si],edgelist$to]
      })
      rec.dat.names = apply(FUN = paste, MARGIN = 1, sapply(subunit.list, rownames), collapse = "_")
      rec.dat = Reduce('*',subunit.list)
      rownames(rec.dat) = rec.dat.names
    }
    rec.dats[[s]] = rec.dat
    final_indices = c(subset_indices, final_indices)
  }
  rec.dat = do.call(rbind, rec.dats)[order(final_indices),] # restore order


  # Make SCC matrix
  scc <- lig.dat*rec.dat
  rownames(scc) <- paste(rownames(lig.dat),rownames(rec.dat),sep = '—')
  # patch here!
  rownames(scc) <- gsub(rownames(scc), pattern = "_", replacement = "+", fixed = T)
  colnames(scc) <- paste(colnames(lig.dat),colnames(rec.dat),sep = '—')
  sending.cell.idents <- as.character(Seurat::Idents(filtered.obj)[colnames(lig.dat)])
  receiving.cell.idents <- as.character(Seurat::Idents(filtered.obj)[colnames(rec.dat)])

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = scc, assay = 'CellToCellSpatial')

  # JC: Seurat V5 will not create data slot automatically, the following step is to manually add this slot
  if(SeuratObject::Version(demo) >= "5.0.0"){
    demo <- Seurat::NormalizeData(demo,assay = "CellToCellSpatial")  # Seura Object need to be >= 5.0.1
    demo@assays$CellToCellSpatial@layers$data <- demo@assays$CellToCellSpatial@layers$counts # Seura Object need to be >= 5.0.1

  }

  # Add key metadata

  meta.data.to.add <- data.frame(SendingType = sending.cell.idents,
                                 ReceivingType = receiving.cell.idents)

  rownames(meta.data.to.add) <- colnames(scc)
  meta.data.to.add$VectorType <- paste(meta.data.to.add$SendingType,
                                       meta.data.to.add$ReceivingType,
                                       sep = '—')

  #Add metadata to the Seurat object
  demo <- Seurat::AddMetaData(demo,metadata = meta.data.to.add)

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- colnames(lig.dat)
    receiving.barcodes <- colnames(rec.dat)
    # Pull and format sending and receiving metadata
    # jc: possible bug, change object to filtered.obj
    sending.metadata <- as.matrix(filtered.obj@meta.data[,meta.data.to.map,drop=FALSE][sending.barcodes,])
    receiving.metadata <- as.matrix(filtered.obj@meta.data[,meta.data.to.map,drop=FALSE][receiving.barcodes,])
    # Make joint metadata
    datArray <- abind::abind(sending.metadata,receiving.metadata,along=3)
    joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- cbind(sending.metadata,receiving.metadata,joint.metadata)
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,receiving.barcodes,sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$VectorType

  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$VectorType)),'distinct VectorTypes were computed, out of',length(table(Seurat::Idents(filtered.obj)))^2,'total possible'))


  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("CellToCellSpatialMatrix","metadata")
    output_list[["CellToCellSpatialMatrix"]] <- getSeuratAssay(demo,"CellToCellSpatial","counts")
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }

}
environment(fixed_RunCellToCellSpatial) <- asNamespace("NICHES")
assignInNamespace("RunCellToCellSpatial", fixed_RunCellToCellSpatial, asNamespace("NICHES"))

fixed_RunCellToNeighborhood <- function(filtered.obj,
                                  ground.truth,
                                  assay,
                                  meta.data.to.map,
                                  blend="mean",
                                  edgelist,
                                  output_format
){
  # Make ligand matrix

  # needed to make this more efficient, do everything without casting to dense matrix if possible
  # maximum number of subunits
  nsubunits = ncol(ground.truth$source.subunits)
  source_subunits_subset = ground.truth$source.subunits
  lig.dats = list()
  source_indices = 1:nrow(source_subunits_subset)
  final_indices = c()
  # for each "size class" of complexes (1 unit, 2 subunits, 3, ...), extract only those columns,
  # aggregate and rbind everything at the end, this should avoid the dense matrix filled with 1s issue
  for (s in nsubunits:1) {
    applicable_complexes = !is.na(source_subunits_subset[,s])
    source_complexes_filt = source_subunits_subset[applicable_complexes,]
    source_subunits_subset = source_subunits_subset[!applicable_complexes,] # remove from list for lower size classes
    subset_indices = source_indices[applicable_complexes]
    source_indices= source_indices[!applicable_complexes] # same for indices for restoring order at the end

    if (s == 1) {
      lig.dat = getSeuratAssay(filtered.obj,assay,"data")[source_complexes_filt[,1],edgelist$from]
    } else {
      subunit.list <- lapply(1:s, function(si){
        getSeuratAssay(filtered.obj,assay,"data")[source_complexes_filt[,si],edgelist$from]
      })
      lig.dat.names = apply(FUN = paste, MARGIN = 1, sapply(subunit.list, rownames), collapse = "_")
      lig.dat = Reduce('*',subunit.list)
      rownames(lig.dat) = lig.dat.names
    }
    lig.dats[[s]] = lig.dat
    final_indices = c(subset_indices, final_indices)
  }
  lig.dat = do.call(rbind, lig.dats)[order(final_indices),] # restore order

  # Make receptor matrix

  nsubunits = ncol(ground.truth$target.subunits)
  target_subunits_subset = ground.truth$target.subunits
  rec.dats = list()
  target_indices = 1:nrow(target_subunits_subset)
  final_indices = c()
  # for each "size class" of complexes (1 unit, 2 subunits, 3, ...), extract only those columns,
  # aggregate and rbind everything at the end, this should avoid the dense matrix filled with 1s issue
  for (s in nsubunits:1) {
    applicable_complexes = !is.na(target_subunits_subset[,s])
    target_complexes_filt = target_subunits_subset[applicable_complexes,]
    target_subunits_subset = target_subunits_subset[!applicable_complexes,] # remove from list for lower size classes
    subset_indices = target_indices[applicable_complexes]
    target_indices = target_indices[!applicable_complexes] # same for indices for restoring order at the end

    if (s == 1) {
      rec.dat = getSeuratAssay(filtered.obj,assay,"data")[target_complexes_filt[,1],edgelist$to]
    } else {
      subunit.list <- lapply(1:s, function(si){
        getSeuratAssay(filtered.obj,assay,"data")[target_complexes_filt[,si],edgelist$to]
      })
      rec.dat.names = apply(FUN = paste, MARGIN = 1, sapply(subunit.list, rownames), collapse = "_")
      rec.dat = Reduce('*',subunit.list)
      rownames(rec.dat) = rec.dat.names
    }
    rec.dats[[s]] = rec.dat
    final_indices = c(subset_indices, final_indices)
  }
  rec.dat = do.call(rbind, rec.dats)[order(final_indices),] # restore order


  # Make SCC matrix
  scc <- lig.dat*rec.dat
  rownames(scc) <- paste(rownames(lig.dat),rownames(rec.dat),sep = '—')
  # patch here!
  rownames(scc) <- gsub(rownames(scc), pattern = "_", replacement = "+", fixed = T)

  # Condense by column name
  colnames(scc) <- colnames(lig.dat) # Make colnames equal to sending cell
  # scc <- as.matrix(scc) # just no, sparse eval is much faster
  # if(blend == "sum") scc <- t(rowsum(t(scc), colnames(scc)))
  # else if(blend == "mean")  scc <- sapply(unique(colnames(scc)), function(sending_name)
  #   rowMeans(scc[,colnames(scc)== sending_name,drop=FALSE], na.rm=TRUE) )
  # more efficiently, do rowsum unconditionally and then look back at colnames(lig.dat)
  # to find out what to devide by
  scc <- Matrix::t(rowsum(Matrix::t(scc), colnames(scc)))
  if(blend == "mean") {
    scc <- sweep(scc, STATS=table(colnames(lig.dat))[colnames(scc)], MARGIN=2, FUN = "/")
  }

  # Label columns properly
  barcodes <- colnames(scc)
  colnames(scc) <- paste(barcodes,'Neighborhood',sep = '—')

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = scc,assay = 'CellToNeighborhood')
  # JC: Seurat V5 will not create data slot automatically, the following step is to manually add this slot
  if(SeuratObject::Version(demo) >= "5.0.0"){
    demo <- Seurat::NormalizeData(demo,assay = "CellToNeighborhood")  # Seura Object need to be >= 5.0.1
    demo@assays$CellToNeighborhood@layers$data <- demo@assays$CellToNeighborhood@layers$counts # Seura Object need to be >= 5.0.1

  }

  # Add metadata based on ident slot
  # bug fix: add the Neighborhood - prefix
  sending_type.meta <- data.frame(SendingCell = barcodes,
                                  SendingType = Seurat::Idents(filtered.obj)[barcodes],
                                  row.names = paste(barcodes,"Neighborhood",sep = '—'))

  demo <- Seurat::AddMetaData(demo,metadata = sending_type.meta,col.name = c("SendingCell","SendingType"))

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- barcodes # Only sending cell metadata applies for this function
    #receiving.barcodes <- colnames(rec.map)
    # Pull and format sending and receiving metadata
    # jc: possible bug, change object to filtered.obj
    sending.metadata <- as.matrix(filtered.obj@meta.data[,meta.data.to.map,drop=FALSE][sending.barcodes,])
    #receiving.metadata <- as.matrix(object@meta.data[,meta.data.to.map][receiving.barcodes,])
    # Make joint metadata
    #datArray <- abind(sending.metadata,receiving.metadata,along=3)
    #joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    #colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    #colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    #colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- sending.metadata
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,'Neighborhood',sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$SendingType
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$SendingCell)),'Cell-To-Neighborhood edges were computed, across',length(unique(demo$SendingType)),'cell types'))

  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("CellToNeighborhoodMatrix","metadata")
    output_list[["CellToNeighborhoodMatrix"]] <- getSeuratAssay(demo,"CellToNeighborhood","counts")
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }

}
environment(fixed_RunCellToNeighborhood) <- asNamespace("NICHES")
assignInNamespace("RunCellToNeighborhood", fixed_RunCellToNeighborhood, asNamespace("NICHES"))

fixed_RunNeighborhoodToCell <- function(filtered.obj,
                                  ground.truth,
                                  assay,
                                  meta.data.to.map,
                                  blend="mean",
                                  edgelist,
                                  output_format
){

  # Make ligand matrix

  # needed to make this more efficient, do everything without casting to dense matrix if possible
  # maximum number of subunits
  nsubunits = ncol(ground.truth$source.subunits)
  source_subunits_subset = ground.truth$source.subunits
  lig.dats = list()
  source_indices = 1:nrow(source_subunits_subset)
  final_indices = c()
  # for each "size class" of complexes (1 unit, 2 subunits, 3, ...), extract only those columns,
  # aggregate and rbind everything at the end, this should avoid the dense matrix filled with 1s issue
  for (s in nsubunits:1) {
    applicable_complexes = !is.na(source_subunits_subset[,s])
    source_complexes_filt = source_subunits_subset[applicable_complexes,]
    source_subunits_subset = source_subunits_subset[!applicable_complexes,] # remove from list for lower size classes
    subset_indices = source_indices[applicable_complexes]
    source_indices= source_indices[!applicable_complexes] # same for indices for restoring order at the end

    if (s == 1) {
      lig.dat = getSeuratAssay(filtered.obj,assay,"data")[source_complexes_filt[,1],edgelist$from]
    } else {
      subunit.list <- lapply(1:s, function(si){
        getSeuratAssay(filtered.obj,assay,"data")[source_complexes_filt[,si],edgelist$from]
      })
      lig.dat.names = apply(FUN = paste, MARGIN = 1, sapply(subunit.list, rownames), collapse = "_")
      lig.dat = Reduce('*',subunit.list)
      rownames(lig.dat) = lig.dat.names
    }
    lig.dats[[s]] = lig.dat
    final_indices = c(subset_indices, final_indices)
  }
  lig.dat = do.call(rbind, lig.dats)[order(final_indices),] # restore order

  # Make receptor matrix

  nsubunits = ncol(ground.truth$target.subunits)
  target_subunits_subset = ground.truth$target.subunits
  rec.dats = list()
  target_indices = 1:nrow(target_subunits_subset)
  final_indices = c()
  # for each "size class" of complexes (1 unit, 2 subunits, 3, ...), extract only those columns,
  # aggregate and rbind everything at the end, this should avoid the dense matrix filled with 1s issue
  for (s in nsubunits:1) {
    applicable_complexes = !is.na(target_subunits_subset[,s])
    target_complexes_filt = target_subunits_subset[applicable_complexes,]
    target_subunits_subset = target_subunits_subset[!applicable_complexes,] # remove from list for lower size classes
    subset_indices = target_indices[applicable_complexes]
    target_indices = target_indices[!applicable_complexes] # same for indices for restoring order at the end

    if (s == 1) {
      rec.dat = getSeuratAssay(filtered.obj,assay,"data")[target_complexes_filt[,1],edgelist$to]
    } else {
      subunit.list <- lapply(1:s, function(si){
        getSeuratAssay(filtered.obj,assay,"data")[target_complexes_filt[,si],edgelist$to]
      })
      rec.dat.names = apply(FUN = paste, MARGIN = 1, sapply(subunit.list, rownames), collapse = "_")
      rec.dat = Reduce('*',subunit.list)
      rownames(rec.dat) = rec.dat.names
    }
    rec.dats[[s]] = rec.dat
    final_indices = c(subset_indices, final_indices)
  }
  rec.dat = do.call(rbind, rec.dats)[order(final_indices),] # restore order


  # Make SCC matrix
  scc <- lig.dat*rec.dat
  rownames(scc) <- paste(rownames(lig.dat),rownames(rec.dat),sep = '—')
  # patch here!
  rownames(scc) <- gsub(rownames(scc), pattern = "_", replacement = "+", fixed = T)

  # Condense by column name
  colnames(scc) <- colnames(rec.dat) # Make colnames equal to receiving cell
  # scc <- as.matrix(scc)
  # if(blend == "sum") scc <- t(rowsum(t(scc), colnames(scc)))
  # else if(blend == "mean") scc <- sapply(unique(colnames(scc)), function(receiving_name)
  #   rowMeans(scc[,colnames(scc)== receiving_name,drop=FALSE], na.rm=TRUE) )
  scc <- Matrix::t(rowsum(Matrix::t(scc), colnames(scc)))
  if(blend == "mean") {
    scc <- sweep(scc, STATS=table(colnames(rec.dat))[colnames(scc)], MARGIN=2, FUN = "/")
  }

  # Label columns properly
  barcodes <- colnames(scc)
  colnames(scc) <- paste('Neighborhood',barcodes,sep = '—')

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = scc,assay = 'NeighborhoodToCell')
  # JC: Seurat V5 will not create data slot automatically, the following step is to manually add this slot
  if(SeuratObject::Version(demo) >= "5.0.0"){
    demo <- Seurat::NormalizeData(demo,assay = "NeighborhoodToCell")  # Seura Object need to be >= 5.0.1
    demo@assays$NeighborhoodToCell@layers$data <- demo@assays$NeighborhoodToCell@layers$counts # Seura Object need to be >= 5.0.1

  }

  # Add metadata based on ident slot
  # bug fix: add the Neighborhood - prefix
  receiving_type.meta <- data.frame(ReceivingCell = barcodes,
                                    ReceivingType = Seurat::Idents(filtered.obj)[barcodes],
                                    row.names = paste("Neighborhood",barcodes,sep = '—'))

  demo <- Seurat::AddMetaData(demo,metadata = receiving_type.meta,col.name = c("ReceivingCell","ReceivingType"))

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    #sending.barcodes <- barcodes
    receiving.barcodes <- barcodes # Only receiving cell metadata applies for this function
    # Pull and format sending and receiving metadata
    #sending.metadata <- as.matrix(object@meta.data[,meta.data.to.map][sending.barcodes,])
    # jc: possible bug, change object to filtered.obj
    receiving.metadata <- as.matrix(filtered.obj@meta.data[,meta.data.to.map,drop=FALSE][receiving.barcodes,])
    # Make joint metadata
    #datArray <- abind(sending.metadata,receiving.metadata,along=3)
    #joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    #colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    #colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    #colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- receiving.metadata
    rownames(meta.data.to.add.also) <- paste('Neighborhood',receiving.barcodes,sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$ReceivingType
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$ReceivingCell)),'Neighborhood-To-Cell edges were computed, across',length(unique(demo$ReceivingType)),'cell types'))

  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("NeighborhoodToCellMatrix","metadata")
    output_list[["NeighborhoodToCellMatrix"]] <- getSeuratAssay(demo,"NeighborhoodToCell","counts")
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }

}
environment(fixed_RunNeighborhoodToCell) <- asNamespace("NICHES")
assignInNamespace("RunNeighborhoodToCell", fixed_RunNeighborhoodToCell, asNamespace("NICHES"))

# add edgelist_fun
augmented_RunNICHES.default <- function(object,
                              assay="RNA",
                              LR.database="fantom5",
                              species,
                              min.cells.per.ident = NULL,
                              min.cells.per.gene = NULL,
                              meta.data.to.map = NULL,
                              position.x = NULL,
                              position.y = NULL,
                              custom_LR_database = NULL,
                              k = 4,
                              rad.set = NULL,
                              blend = 'mean',
                              CellToCell = T,
                              CellToSystem = F,
                              SystemToCell = F,
                              CellToCellSpatial = F,
                              CellToNeighborhood = F,
                              NeighborhoodToCell = F,
                              output_format = "seurat",
                              edgelist_fun = compute_edgelist,
                              ...
){
  # TODO: check the parameter validity here, then register the parameters
  # 1: check the data format of the required parameters and the optional parameters (if(!is.null())): int, character, etc.
  # 2. check some of the dependencies of the parameters
  #   (1) If LR.database is 'custom', then 'custom_LR_database' can't be NULL
  #   (2) If CellToCellSpatial,CellToNeighborhood, or NeighborhoodToCell is T, then rad.set, 'position.x' and 'position.y' can't be null
  # 3. check whether other functions have done similar checks
  # 4. think whether to put stop or warning


  # check lr database, species, custom_LR_database
  if(LR.database %in% c("fantom5","omnipath","custom")){
    if(LR.database == "fantom5"){
      if(!(species %in% c('human','mouse','rat','pig')))
        stop(paste0("Unsupported species ",species, " for fantom 5 database, only 'human','mouse','rat',or 'pig' supported."))
    }
    if(LR.database == "omnipath"){
      if(!(species %in% c('human','mouse','rat')))
        stop(paste0("Unsupported species ",species, " for omnipath database, only 'human','mouse',or 'rat' supported."))
    }
    if(LR.database == "custom"){
      # check the format of custom_LR_database
      if(is.null(custom_LR_database)) stop("custom_LR_database is NULL")
      # TODO: check gene names
      message("Custom Ligand Receptor database enabled...")
      message("Checking the format of the custom database...")
      if(!is.data.frame(custom_LR_database)){
        warning("Custom database provided is not in dataframe format.")
        warning("Converting to dataframe format...")
        custom_LR_database <- as.data.frame(custom_LR_database)
      }
      if(ncol(custom_LR_database) < 2) stop("Custom database provided contains less than 2 columns.")
    }
  }else stop('\n LR.receptor argument not recognized. Only accepts "omnipath","fantom5" or "custom".')

  if(!is.null(custom_LR_database) & LR.database!="custom"){
    warning("custom_LR_database is provided but LR.databse is not specified as 'custom'")
  }
  # Convert any non-integer inputs to integers. Still allows NULL as option.
  if (!is.null(min.cells.per.ident)){
    min.cells.per.ident <- as.integer(min.cells.per.ident)
  }
  if (!is.null(min.cells.per.gene)){
    min.cells.per.gene <- as.integer(min.cells.per.gene)
  }

  # check indicators
  # jc: Add organization names to the list
  org_names_indicator <- c(CellToCell,CellToSystem,SystemToCell,CellToCellSpatial,CellToNeighborhood,NeighborhoodToCell)
  org_names_indicator <- sapply(org_names_indicator,function(org){
    org_out <- as.logical(org)
    if(is.na(org_out)) warning("Organization indicator ",org," is not set to TRUE or FALSE.")
    return(org_out)
  })
  names(org_names_indicator) <- c("CellToCell","CellToSystem","SystemToCell","CellToCellSpatial","CellToNeighborhood","NeighborhoodToCell")


  # If requested, additionally calculate spatially-limited NICHES organizations
  if (org_names_indicator["CellToCellSpatial"] == T | org_names_indicator["CellToNeighborhood"] == T | org_names_indicator["NeighborhoodToCell"] == T){

    if (is.null(position.x) | is.null(position.y)){stop("\n Position information not provided. Please specify metadata columns containing x- and y-axis spatial coordinates.")}
    if(!is.null(k)) k <- as.integer(k)
    if(!is.null(rad.set)) rad.set <- as.numeric(rad.set)
  }

  if((!is.null(position.x) | !is.null(position.y)) & org_names_indicator["CellToCellSpatial"] == F & org_names_indicator["CellToNeighborhood"] == F & org_names_indicator["NeighborhoodToCell"] == F)
    warning("Spatial positions are provided but the spatial organization functions: 'CellToCellSpatial','CellToNeighborhood', and 'NeighborhoodToCell' are set to FALSE.")


  if(org_names_indicator["CellToSystem"] == T | org_names_indicator["SystemToCell"] == T)
    if(!blend %in% c("sum","mean","mean.adj")) stop("blend paramter is not recognized: need to be 'sum' or 'mean" )
  if(blend == "sum") warning("Operator `sum` will be deprecated in the later release.")
  if(blend == "mean.adj") warning("Operator `mean.adj` is still in experimental stage, use in caution")

  # Initialize output structure
  output <- list()


  # jc: move the shared preprocessing steps here to avoid redundancy and reduce the number of parameters to be passed to other functions

  # NOTE: relies on Idents(object) to be cell types to subset
  filtered.obj <- prepSeurat(object,assay,min.cells.per.ident,min.cells.per.gene)
  ground.truth <- lr_load(LR.database,custom_LR_database,species,rownames(filtered.obj@assays[[assay]]))
  if (org_names_indicator["CellToCellSpatial"] == T | org_names_indicator["CellToNeighborhood"] == T | org_names_indicator["NeighborhoodToCell"] == T){
    ## 1. Move the neighbor graph construction here
    ## 2. Enable a k-nearest-neighbor parameter as an alternative
    edgelist <- edgelist_fun(filtered.obj,position.x,position.y,k,rad.set)
  }

  # check the output format
  if(!output_format %in% c("seurat","raw"))
    stop(paste0("Unsupported output format: ",output_format,", Currently only 'seurat' and 'raw' are supported."))


  # Calculate NICHES organizations without spatial restrictions
  # jc: only pass the processed data to each function
  # NOTE: RunCellToCell relies on Idents(object) to be cell types to subset
  #       Also each RunXXX function needs Idents(object) to build VectorType meta data

  if (CellToCell == T){output[[length(output)+1]] <- RunCellToCell(filtered.obj=filtered.obj,
                                                                   ground.truth=ground.truth,
                                                                   assay = assay,
                                                                   meta.data.to.map = meta.data.to.map,
                                                                   output_format = output_format
  )}
  if (CellToSystem == T){output[[length(output)+1]] <- RunCellToSystem(filtered.obj=filtered.obj,
                                                                       ground.truth=ground.truth,
                                                                       assay = assay,
                                                                       meta.data.to.map = meta.data.to.map,
                                                                       blend = blend,
                                                                       output_format = output_format
  )}
  if (SystemToCell == T){output[[length(output)+1]] <- RunSystemToCell(filtered.obj=filtered.obj,
                                                                       ground.truth=ground.truth,
                                                                       assay = assay,
                                                                       meta.data.to.map = meta.data.to.map,
                                                                       blend = blend,
                                                                       output_format = output_format
  )}


  if (CellToCellSpatial == T){output[[length(output)+1]] <- RunCellToCellSpatial(filtered.obj=filtered.obj,
                                                                                 ground.truth=ground.truth,
                                                                                 assay = assay,
                                                                                 meta.data.to.map = meta.data.to.map,
                                                                                 edgelist = edgelist,
                                                                                 output_format = output_format
  )} #Spatially-limited Cell-Cell vectors
  if (CellToNeighborhood == T){output[[length(output)+1]] <- RunCellToNeighborhood(filtered.obj=filtered.obj,
                                                                                   ground.truth=ground.truth,
                                                                                   assay = assay,
                                                                                   meta.data.to.map = meta.data.to.map,
                                                                                   blend = blend,
                                                                                   edgelist = edgelist,
                                                                                   output_format = output_format
  )} #Spatially-limited Cell-Neighborhood vectors
  if (NeighborhoodToCell == T){output[[length(output)+1]] <- RunNeighborhoodToCell(filtered.obj=filtered.obj,
                                                                                   ground.truth=ground.truth,
                                                                                   assay = assay,
                                                                                   meta.data.to.map = meta.data.to.map,
                                                                                   blend = blend,
                                                                                   edgelist = edgelist,
                                                                                   output_format = output_format
  )} #Spatially-limited Neighborhood-Cell vectors (niches)

  # jc: Add organization names to the list
  names(output) <- names(org_names_indicator)[org_names_indicator]

  # Compile objects for output
  return(output)
}
environment(augmented_RunNICHES.default) <- asNamespace("NICHES")
assignInNamespace("RunNICHES.default", augmented_RunNICHES.default, asNamespace("NICHES"))
