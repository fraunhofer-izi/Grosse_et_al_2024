# function for plotting densities for score models:
# density of a mixture of gaussians
dnormmixt <- function(t, lam, mu, sig) {
  m <- length(lam); f <- 0
  for (j in 1:m) f <- f + lam[j]*dnorm(t,mean=mu[j],sd=sig[j])
  f
  }


# combine nichenet and liana curated LR-DBs
# use cached if possible
# tried to load all samples once and look for updates for genes not in our samples

get_combined_lr_dbs <- function(cached = T){
  if (cached && file.exists("./data/LRI_networks/cached_liana_nichenet_consensus_db.rds")) {
    combined_LR_resources = readRDS("./data/LRI_networks/cached_liana_nichenet_consensus_db.rds")
  } else {
    liana_LR_db = liana::select_resource("Consensus")$Consensus
    liana_LR_db = liana_LR_db %>%
      select(source_genesymbol, target_genesymbol) %>%
      rename(from = source_genesymbol, to = target_genesymbol) %>%
      mutate(source = "liana-consensus")

    # in the NicheNet framework, ligand-target links are predicted based on collected biological knowledge on ligand-receptor, signaling and gene regulatory interactions
    # The complete networks can be downloaded from Zenodo
    if (file.exists("./data/LRI_networks/NicheNet/lr_network.rds")) {
      nichenet_LR_db = readRDS("./data/LRI_networks/NicheNet/lr_network.rds")
    } else {
      nichenet_LR_db = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
      dir.create("./data/LRI_networks/NicheNet/", recursive = T)
      saveRDS(nichenet_LR_db, file = "./data/LRI_networks/NicheNet/lr_network.rds")
    }
    nichenet_LR_db = nichenet_LR_db %>%
      select(from, to) %>%
      mutate(source = "nichenet-curated")

    combined_LR_resources = rbind(liana_LR_db, nichenet_LR_db)
    combined_LR_resources = combined_LR_resources %>%
      mutate(mechanism = paste(from, to, sep = "—")) %>%
      group_by(mechanism) %>%
      summarize(to = data.table::first(to), from = data.table::first(from),
                source = paste(sort(unique(source)), collapse = ", "))

    ### translate, save cache
    # original entry
    combined_LR_resources$orig_entry = combined_LR_resources$mechanism
    known_genes = c()

    if (!exists("manifest")) {
      warning("Can't load Seurat objects to update LRI names")
    } else {
      print("Updating Liana/Nichenet consensus LRI gene names based on our samples:")
      require("Seurat")
      for (sname in manifest$samples) {
        se_tmp = readRDS(file.path(file.path(manifest$workdata_parent, manifest[[sname]]$workdata), "02_seurat.rds"))
        known_genes = unique(c(known_genes, rownames(se_tmp@assays$Spatial)))
      }
      rm(se_tmp)
      combined_LR_resources_is_complex = str_detect(combined_LR_resources$to, "_") | str_detect(combined_LR_resources$from, "_")
      combined_LR_resources_complexes = combined_LR_resources[combined_LR_resources_is_complex,]
      combined_LR_resources_easy = combined_LR_resources[!combined_LR_resources_is_complex,]
      LRI_genes = unique(c(combined_LR_resources_easy$to, combined_LR_resources_easy$from,
                           unlist(str_split(combined_LR_resources_complexes$to, pattern = "_")), unlist(str_split(combined_LR_resources_complexes$from, pattern = "_"))))
      LRI_genes_unknown = LRI_genes[!(LRI_genes %in% known_genes)]
      # since we mix annotations, gene names can be too old or too young, try to resolve this
      LRI_genes_tr = resolve_translations(LRI_genes_unknown, prefer = known_genes)
      # translate only what is in this list
      combined_LR_resources_easy_tr = combined_LR_resources_easy
      combined_LR_resources_easy_tr$to = setNames(coalesce(LRI_genes_tr[combined_LR_resources_easy_tr$to], combined_LR_resources_easy_tr$to), nm = NULL)
      combined_LR_resources_easy_tr$from = setNames(coalesce(LRI_genes_tr[combined_LR_resources_easy_tr$from], combined_LR_resources_easy_tr$from), nm = NULL)
      combined_LR_resources_easy_tr$mechanism = paste0(combined_LR_resources_easy_tr$from, "—", combined_LR_resources_easy_tr$to)

      combined_LR_resources_complexes_tr = combined_LR_resources_complexes
      combined_LR_resources_complexes_tr$to = sapply(str_split(combined_LR_resources_complexes_tr$to, pattern = "_"),
             function(gl){return(paste(setNames(coalesce(LRI_genes_tr[gl], gl), nm = NULL), collapse = "_"))})
      combined_LR_resources_complexes_tr$from = sapply(str_split(combined_LR_resources_complexes_tr$from, pattern = "_"),
                                                     function(gl){return(paste(setNames(coalesce(LRI_genes_tr[gl], gl), nm = NULL), collapse = "_"))})
      combined_LR_resources_complexes_tr$mechanism = paste0(combined_LR_resources_complexes_tr$from, "—", combined_LR_resources_complexes_tr$to)

      # update
      combined_LR_resources = rbind(combined_LR_resources_complexes_tr, combined_LR_resources_easy_tr)
      # cache
      saveRDS(combined_LR_resources, file = "./data/LRI_networks/cached_liana_nichenet_consensus_db.rds")
    }

  }
  return(combined_LR_resources)
}

# try to resolve translations in either direction
# prefer symbols in this order:
# - an updated symbol is in prefer
# - a previous symbol is in prefer
# - if not in prefer, use the most up-to-date symbol
#
# there may be complex failure modes with weird renamings, so tread carefully
resolve_translations <- function(LRI_genes_unknown, prefer){
  old_names = LRI_genes_unknown
  newer_names = setNames(Seurat::UpdateSymbolList(old_names), nm = old_names)
  older_names = ReverseGenesymbolUpdates(old_names, search.types = "prev_symbol")

  translations = lapply(old_names, function(nm){
    if (newer_names[nm] %in% prefer) {
      return(newer_names[nm])
    } else {
      if (any(older_names[[nm]] %in% prefer)) {
        older_sub = older_names[[nm]]
        return(older_sub[older_sub %in% prefer][1])
      } else {
        return(newer_names[nm])
      }
    }
  })

  return(setNames(unlist(translations), nm = old_names))
}


geomean = function(x){exp(sum(log(x))/length(x))}

## Return Seurat Object with an added assay with a binarized count matrix (1/0 present or absent)
add_binarized_counts <- function(se, layer = "counts", assay = NULL, new.assay.name = "binarycounts"){
  if (is.null(assay)){
    assay = DefaultAssay(se)
  }
  cnts = GetAssayData(object = se, layer = layer, assay = assay)
  bincounts = cnts>0
  se[[new.assay.name]] = CreateAssayObject(counts = bincounts, check.matrix = T)
  return(se)
}

### reverse lookup of old genesymbols from newer ones
### Modeled after the GeneSymbolThesaurus function from Seurat
ReverseGenesymbolUpdates <- function(symbols, timeout = 20,
                                     search.types = c("alias_symbol", "prev_symbol"),
                                     verbose = TRUE, ...){
  db.url <- "http://rest.genenames.org/fetch"
  search.types <- match.arg(arg = search.types, several.ok = TRUE)
  synonyms <- vector(mode = "list", length = length(x = symbols))
  not.found <- vector(mode = "logical", length = length(x = symbols))
  names(x = not.found) <- names(x = synonyms) <- symbols
  if (verbose) {
    pb <- txtProgressBar(max = length(x = symbols), style = 3,
                         file = stderr())
  }
  for (symbol in symbols) {
    sym.syn <- character()

    response <- httr::GET(url = paste(db.url, "symbol", symbol,
                                sep = "/"), config = c(httr::accept_json(), httr::timeout(seconds = timeout)),
                    ...)
    if (!identical(x = httr::status_code(x = response), y = 200L)){
      # querying error
    } else {
      response <- httr::content(x = response)
      if (response$response$numFound < 1) {
        # no entry found
      } else {
        # I believe this case is impossible for querying current symbols,
        # so this throws now
        if (response$response$numFound > 1) {
          stop("Multiple hits found for symbol ", symbol,
               ". This should not be possible, please check this HUGO entry manually.")
        }
        if (("alias_symbol" %in% search.types) &&
            ("alias_symbol" %in% names(response$response$docs[[1]]))) {
          sym.syn <- c(sym.syn, unlist(response$response$docs[[1]]$alias_symbol))
        }
        if (("prev_symbol" %in% search.types) &&
            ("prev_symbol" %in% names(response$response$docs[[1]]))) {
          sym.syn <- c(sym.syn, unlist(response$response$docs[[1]]$prev_symbol))
        }
      }
    }

    not.found[symbol] <- length(x = sym.syn) < 1
    if (length(x = sym.syn) >= 1) {
      synonyms[[symbol]] <- sym.syn
    }
    if (verbose) {
      setTxtProgressBar(pb = pb, value = pb$getVal() +
                          1)
    }
  }
  if (verbose) {
    close(con = pb)
  }
  if (sum(not.found) > 0) {
    warning("The following symbols had no synonyms: ", paste(names(x = which(x = not.found)),
                                                             collapse = ", "), call. = FALSE, immediate. = TRUE)
  }
  synonyms <- Filter(f = Negate(f = is.null), x = synonyms)
  return(synonyms)
}


## patched to separate complexes by + instead of mangling names to be - separated
## because the - sometimes conflicts with - in gene names

# returns output object of NICHES analysis
#
# pass the chosen clustering resolution to conserved.metadata to get them in the NICHES output
# automatically tries to conserve Kunz (histo annotations) and
# active.region.annotated (from thresholding) if those are in the metadata of se
#
# all ... arguments go to the RunNICHES call
do_NICHES_with_patches <- function(se, lr_db, conserved.metadata = NULL, assay = NULL){
  library("NICHES")
  if (is.null(assay)) {
    assay = DefaultAssay(se)
  }

  # common preparations
  for (conserveme in c("Kunz", "active.region.annotated")){
    if (conserveme %in% colnames(se@meta.data)) {
      conserved.metadata = c(conserved.metadata, conserveme)
    }
  }
  conserved.metadata = unique(conserved.metadata)

  DefaultAssay(se) = assay
  se@meta.data$ypos <- se@tools$Staffli@meta_data$x # for whatever reason,
  se@meta.data$xpos <- se@tools$Staffli@meta_data$y # these are swapped
  # because we use the array rows/cols as x and y coordinates, but the
  # physical distance is different between those directions.
  # leaving out the modification and setting rad.set=2 in NICHES
  # (-> neighborhoods of euclidean distance <= 2 in the provided coords)
  # leads to diamond-shaped neighborhoods, xpos_modified makes sure
  # the neighborhoods are hexagonal
  se[["xpos_modified"]] = 1.5 * se$xpos

  NICHES_output <- RunNICHES(object = se,
                             LR.database = "custom",
                             custom_LR_database = lr_db,
                             assay = assay,
                             position.x = 'xpos_modified',
                             position.y = 'ypos',
                             k = NULL,
                             rad.set = 2, # distance, Geometry dependent
                             min.cells.per.ident = 0,
                             min.cells.per.gene = NULL,
                             meta.data.to.map = conserved.metadata,
                             blend = "sum",
                             CellToCell = F, CellToSystem = F, SystemToCell = F,
                             CellToCellSpatial = T, CellToNeighborhood = T, NeighborhoodToCell = T)

  return(NICHES_output)
}


# adds CellToNeighborhood and NeighborhoodToCell from NICHES_output to se
# returns modified se
add_NICHES_to_se <- function(se, NICHES_output){
  # Add NICHES obj CellToNeighborhood
  niche = NICHES_output$CellToNeighborhood
  niches.data <- GetAssayData(object =  niche[['CellToNeighborhood']], layer = 'counts')
  colnames(niches.data) <- niche[['SendingCell']]$SendingCell
  se[["CellToNeighborhood"]] <- CreateAssayObject(data = niches.data)

  DefaultAssay(se) <- "CellToNeighborhood"
  se <- ScaleData(se, features = rownames(se))

  # Add NICHES obj NeighborhoodToCell
  niche = NICHES_output$NeighborhoodToCell
  niches.data <- GetAssayData(object =  niche[['NeighborhoodToCell']], layer = 'counts')
  colnames(niches.data) <- niche[['ReceivingCell']]$ReceivingCell
  se[["NeighborhoodToCell"]] <- CreateAssayObject(data = niches.data)

  DefaultAssay(se) <- "NeighborhoodToCell"
  se <- ScaleData(se, features = rownames(se))

  return(se)
}


lr_de_spots <- function(se, group.by,
                        active.groups, background.groups,
                        assay = NULL, mean_type = "arithmetic", R = 9999, ...){
  if (!(group.by %in% colnames(se@meta.data))) {
    stop("group.by must be in the metadata of se and contain at least one entry of each of active.groups and background.groups")
  }

  if (is.null(assay)) {
    assay = DefaultAssay(se)
  }
  DefaultAssay(se) <- assay
  Idents(se) = group.by

  # we must make sure that active.groups and background.groups labels are really in the data
  # or else diff_mean_test may just throw errors, so subset and throw helpful errors if necessary
  lvls = unique(se@meta.data[,group.by])
  active.not.in.se = setdiff(active.groups, lvls)
  if (length(active.not.in.se)>0) {warning(paste0("Active class label(s) not found in grouping column: [",paste(active.not.in.se, collapse = ", "),"]",
                                                  "\n  Checking if other labels are present."), immediate. = T)}
  active.groups.cleaned = intersect(active.groups, lvls)
  if (length(active.groups.cleaned)==0) {stop("None of the active group labels appear in the grouping column, aborting.")}

  background.not.in.se = setdiff(background.groups, lvls)
  if (length(background.not.in.se)>0) {warning(paste0("Background class label(s) not found in grouping column: [",paste(background.not.in.se, collapse = ", "),"]",
                                                      "\n  Checking if other labels are present."), immediate. = T)}
  background.groups.cleaned = intersect(background.groups, lvls)
  if (length(background.groups.cleaned)==0) {stop("None of the background group labels appear in the grouping column, aborting.")}


  # format according to what diff_mean_test expects
  # character-vector of length 2 = test 2 groups
  # length 2 list of character-vectors = make bigger groups out of the labels in each position of the list
  if (length(active.groups)>1 | length(background.groups)>1){
    compare = list(active = c(active.groups.cleaned), background = c(background.groups.cleaned))
  } else {
    compare = c(active.groups.cleaned, background.groups.cleaned)
  }
  border.markers = sctransform::diff_mean_test(y = se[[assay]]@data,                        # pulls from data slot without doing something stupid
                                               group_labels = se@meta.data[[group.by]],     # metadatum to split groups
                                               compare = compare,                           # groups to compare
                                               mean_type = mean_type,                       # how to calculate the mean
                                               R = R, ...)                                  # number of random permutations to derive p-values
  border.markers$tested = T

  # rerun with reduced R just for the mean+lfc stats, set all pvals to
  # 1 because we don't do any tests on those, add a column marking testes mechanisms
  border.markers.stats = sctransform::diff_mean_test(y = se[[assay]]@data,
                                                     group_labels = se@meta.data[[group.by]],
                                                     compare = compare,
                                                     mean_type = mean_type,
                                                     R = 99, verbosity = 0,
                                                     mean_th = 0., log2FC_th = 0., cells_th = 0.)
  border.markers.stats$tested = F
  border.markers.stats$emp_pval = 1.
  border.markers.stats$emp_pval_adj = 1.
  border.markers.stats$pval = 1.
  border.markers.stats$pval_adj = 1.
  border.markers.stats = border.markers.stats[setdiff(rownames(border.markers.stats), rownames(border.markers)),]

  # order by adjusted p value, add the stats of the other mechanisms
  border.markers = rbind(border.markers[order(border.markers$pval_adj),], border.markers.stats)

  return(border.markers)
}


lr_de_edges <- function(input,
                        groups.metadatum,
                        active.groups,
                        background.groups,
                        mean_type = "arithmetic", R = 999){
  # if the whole NICHES_output if given, extract CellToCellSpatial first
  if (length(input)==3) {input = input$CellToCellSpatial}
  if (!(names(input)=="CellToCellSpatial")) {stop("Input should be the NICHES output CellToCellSpatial")}

  # needed so diff_mean_test doesn't die with default setting
  # drop non-existant levels from active.groups and background.groups
  group_labels = input@meta.data[[groups.metadatum]]
  active.groups = active.groups[active.groups %in% group_labels]
  background.groups = background.groups[background.groups %in% group_labels]
  stopifnot("active.groups contains no entries present in the groups.metadatum column"= length(active.groups)>0,
            "background.groups contains no entries present in the groups.metadatum column"= length(background.groups)>0)

  border.markers = sctransform::diff_mean_test(y = input[["CellToCellSpatial"]]@data,
                                               group_labels = group_labels,                      # metadatum to split groups
                                               compare = list(active.groups, background.groups), # groups to compare
                                               mean_type = mean_type,                       # how to calculate the mean
                                               R = R)                                       # number of random permutations to derive p-values
  return(border.markers)
}




########################################################################################################################################
# PATCHED VERSION OF NICHES FUNCTION
# because the transformation to Seurat objects inside of NICHES mangles the names of complexes
# by replacing all _ with - and this interferes with genes like HLA-DRA, which contain -
# COMPLEX GENES ARE NOW SEPARATED BY +
# Addition: more memory-friendly computations for large HD samples
# Addition: edgelist_fun can be set to change the function that computes
# the edgelist to something more efficient, for large datasets
########################################################################################################################################

source("code/helper/NICHES_patches.R")
