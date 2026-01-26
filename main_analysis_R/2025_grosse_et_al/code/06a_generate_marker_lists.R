# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Download necessary resources, transform stuff and cache it.
# Prerequisite for the next step, but does not do much on it's own.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("tidyverse")
library("Seurat")
library("yaml")
library("readxl")
library("stringi")
library("curl")

# bioc
library("GSEABase")
library("BiocParallel")

library("BayesPrism")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

manifest = yaml.load_file("assets/manifest.yaml")
workdata = manifest$workdata_parent
all_samples = setdiff(manifest$samples, manifest$rejected_samples)
random_seed = 1234

if (!(dir.exists("./data/genesets/"))) {
  dir.create("./data/genesets/")
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get panglao DB
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Panglao DB Markers prep
if (!(dir.exists("./data/panglaoDB"))) {
  dir.create("./data/panglaoDB", recursive = T)
}
if (!(file.exists("./data/panglaoDB/PanglaoDB_markers_27_Mar_2020.tsv.gz"))) {
  library(curl)
  curl_download(url = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz",
                destfile = "./data/panglaoDB/PanglaoDB_markers_27_Mar_2020.tsv.gz")
} else {
  print("Panglao DB was already downloaded, skipping step.")
}

# Panglao DB Markers prep
panglao <- read_tsv("./data/panglaoDB/PanglaoDB_markers_27_Mar_2020.tsv.gz")
panglao = panglao %>% dplyr::filter(str_detect(species, "Hs"))

# convert dataframes to a named list of marker genes per celltype
panglao <- panglao %>%
  group_by(`cell type`) %>%
  summarise(geneset = list(`official gene symbol`))
all_gs <- setNames(panglao$geneset, panglao$`cell type`)

# select scores of interest
panglao_selection = c("Fibroblasts", "Epithelial cells", "Endothelial cells", "Basal cells",
                      "B cells", "NK cells", "T cells", "T helper cells", "T cytotoxic cells",
                      "Monocytes", "Macrophages", "Keratinocytes", "Melanocytes",
                      "Stromal cells", "Smooth muscle cells", "Sebocytes")
panglao_gs = all_gs[panglao_selection]
# add Panglao to names
cleaned_names = make.names(names = names(panglao_gs), unique = T)
names(panglao_gs) = cleaned_names
names(panglao_gs) = paste0(names(panglao_gs), "_Panglao")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get jerby-arnon data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

jerby_data_dir = paste0(workdata, "jerby_arnon_2018_GSE115978")
if (!(dir.exists(jerby_data_dir))) {dir.create(jerby_data_dir, recursive = T)}

if (!(file.exists("./data/genesets/jerby_se.rds"))){

    # annotations
    jerby.geo.annotations.archive = paste0(jerby_data_dir, "cell.annotations.csv.gz")
    jerby.geo.annotations.file = paste0(jerby_data_dir, "cell.annotations.csv")
    jerby.geo.annotations.link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115978&format=file&file=GSE115978%5Fcell%2Eannotations%2Ecsv%2Egz"
    if (!(file.exists(jerby.geo.annotations.file))) {
      curl_download(destfile = jerby.geo.annotations.archive, url = jerby.geo.annotations.link, quiet = F)
      # unpack
      R.utils::gunzip(jerby.geo.annotations.archive)
      if (!(file.exists(jerby.geo.annotations.file))) {stop("Something went wrong with extracting the metadata file.")}
    }

    # counts
    jerby.geo.counts.archive = paste0(jerby_data_dir, "cell.counts.csv.gz")
    jerby.geo.counts.file = paste0(jerby_data_dir, "cell.counts.csv")
    jerby.geo.counts.link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115978&format=file&file=GSE115978%5Fcounts%2Ecsv%2Egz"
    if (!(file.exists(jerby.geo.counts.file))) {
      curl_download(destfile = jerby.geo.counts.archive, url = jerby.geo.counts.link, quiet = F)
      # unpack
      R.utils::gunzip(jerby.geo.counts.archive)
      if (!(file.exists(jerby.geo.counts.file))) {stop("Something went wrong with extracting the counts file.")}
    }

    # TPM
    jerby.geo.tpm.archive = paste0(jerby_data_dir, "cell.tpm.csv.gz")
    jerby.geo.tpm.file = paste0(jerby_data_dir, "cell.tpm.csv")
    jerby.geo.tpm.link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE115978&format=file&file=GSE115978%5Ftpm%2Ecsv%2Egz"
    if (!(file.exists(jerby.geo.tpm.file))) {
      curl_download(destfile = jerby.geo.tpm.archive, url = jerby.geo.tpm.link, quiet = F)
      # unpack
      R.utils::gunzip(jerby.geo.tpm.archive)
      if (!(file.exists(jerby.geo.tpm.file))) {stop("Something went wrong with extracting the TPM file.")}
    }

  # piece together a seurat object

  jerby_annotations = read_csv(jerby.geo.annotations.file)
  jerby_annotations = as.data.frame(jerby_annotations)
  rownames(jerby_annotations) = jerby_annotations$cells
  jerby_annotations = jerby_annotations[,-1]

  jerby_counts = read_csv(jerby.geo.counts.file)
  jerby_counts.names = jerby_counts %>% pull(1)
  jerby_counts = as.sparse(as.matrix(jerby_counts[-1]))
  rownames(jerby_counts) = jerby_counts.names
  rm(jerby_counts.names)

  jerby_tpm = read_csv(jerby.geo.tpm.file)
  jerby_tpm.names = jerby_tpm %>% pull(1)
  jerby_tpm = as.sparse(as.matrix(jerby_tpm[-1]))
  rownames(jerby_tpm) = jerby_tpm.names
  rm(jerby_tpm.names)
  # the tpm come in log2(tpm+1), so undo for raw counts slot
  # probably unnecessary, but for the sake  of completeness I do it here
  jerby_tpm_raw = as.sparse( (2^(jerby_tpm) - 1) * 10)

  jerby_se = CreateSeuratObject(counts = jerby_counts, assay = "RNA",
                                project = "jerby_arnon_data", meta.data = jerby_annotations)
  jerby_tpm_assay = CreateAssayObject(counts = jerby_tpm_raw, check.matrix = T)
  jerby_tpm_assay = SetAssayData(object = jerby_tpm_assay, layer = "data",
                                 new.data = jerby_tpm)

  jerby_se[["TPM"]] = jerby_tpm_assay
  DefaultAssay(jerby_se) = "TPM"
  jerby_se = FindVariableFeatures(jerby_se)
  jerby_se = ScaleData(jerby_se)
  jerby_se = RunPCA(jerby_se)
  jerby_se = RunUMAP(jerby_se, dims = 1:25)

  DefaultAssay(jerby_se) = "RNA"
  jerby_se = NormalizeData(jerby_se)

  saveRDS(jerby_se, "./data/genesets/jerby_se.rds")
} else {
  jerby_se = readRDS("./data/genesets/jerby_se.rds")
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get Solé-Boldo data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

sole_boldo_data_dir = paste0(workdata, "sole_boldo_2020_GSE130973")
if (!(dir.exists(sole_boldo_data_dir))) {dir.create(sole_boldo_data_dir, recursive = T)}

if (!(file.exists("./data/genesets/sole_boldo_se.rds"))){

  # download preprocessed Seurat .rds file
  sole_boldo.geo.archive = paste0(sole_boldo_data_dir, "GSE130973_seurat_analysis_lyko.rds.gz")
  sole_boldo.geo.file = paste0(sole_boldo_data_dir, "GSE130973_seurat_analysis_lyko.rds")
  sole_boldo.geo.link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130973&format=file&file=GSE130973%5Fseurat%5Fanalysis%5Flyko%2Erds%2Egz"
  if (!(file.exists(sole_boldo.geo.file))) {
    curl_download(destfile = sole_boldo.geo.archive, url = sole_boldo.geo.link, quiet = F)
    # unpack
    R.utils::gunzip(sole_boldo.geo.archive)
    if (!(file.exists(sole_boldo.geo.file))) {stop("Something went wrong with extracting the metadata file.")}
  }

  # piece together a seurat object
  sole_boldo_se = readRDS(file = sole_boldo.geo.file)
  
  # prevents issues with using Seurat v5 on older objects
  sole_boldo_se = SeuratObject::UpdateSeuratObject(object = sole_boldo_se)
  
  sole_boldo_se$celltype = str_replace(sole_boldo_se$celltype.age, pattern = "_.*", replacement = "")
  # identify cluster numbers with cell types from the paper
  cell.type.translation = unlist(list("0" = "macrophages and dendritic cells_1",
                                      "1" = "secretory-reticular fibroblasts",
                                      "2" = "pro-inflammatory fibroblasts",
                                      "3" = "secretory-papillary fibroblasts",
                                      "4" = "vascular endothelial cells",
                                      "5" = "differentiated keratinocytes",
                                      "6" = "T cells",
                                      "7" = "epidermal stem cells and other undifferentiated progenitors_1",
                                      "8" = "pericytes_1",
                                      "9" = "mesenchymal fibroblasts",
                                      "10" = "pericytes_2",
                                      "11" = "erythrocytes",
                                      "12" = "lymphatic endothelial cells",
                                      "13" = "macrophages and dendritic cells_2",
                                      "14" = "melanocytes",
                                      "15" = "epidermal stem cells and other undifferentiated progenitors_2",
                                      "16" = "macrophages and dendritic cells_3"
  ))
  cell.type.translation.broad = unlist(list("0" = "macrophages and dendritic cells",
                                            "1" = "fibroblasts",
                                            "2" = "fibroblasts",
                                            "3" = "fibroblasts",
                                            "4" = "endothelial cells",
                                            "5" = "keratinocytes",
                                            "6" = "T cells",
                                            "7" = "keratinocytes",
                                            "8" = "pericytes",
                                            "9" = "fibroblasts",
                                            "10" = "pericytes",
                                            "11" = "erythrocytes",
                                            "12" = "endothelial cells",
                                            "13" = "macrophages and dendritic cells",
                                            "14" = "melanocytes/malignant cells",
                                            "15" = "keratinocytes",
                                            "16" = "macrophages and dendritic cells"
  ))
  sole_boldo_se$celltype.broad = setNames(cell.type.translation.broad[sole_boldo_se$celltype], nm = NULL)
  sole_boldo_se$celltype = setNames(cell.type.translation[sole_boldo_se$celltype], nm = NULL)

  saveRDS(sole_boldo_se, file = "./data/genesets/sole_boldo_se.rds")

} else {
  sole_boldo_se = readRDS("./data/genesets/sole_boldo_se.rds")
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get prepared Stubenvoll-Schmidt data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

stub_schmidt_se = readRDS("./data/sc_ref_Stubenvoll_Schmidt_2025.rds")
mela_ref = stub_schmidt_se$mela_ref
nevi_ref = stub_schmidt_se$nevi_ref

if (file.exists("./data/genesets/JA_SB_StubS_harmonized_for_markersearch.rds")) {
  obj = readRDS(file = "./data/genesets/JA_SB_StubS_harmonized_for_markersearch.rds")
  jerby_counts = obj$jerby_counts
  jerby_states = obj$jerby_states
  jerby_types = obj$jerby_types
  sole_boldo_counts = obj$sole_boldo_counts
  sole_boldo_states = obj$sole_boldo_states
  sole_boldo_types = obj$sole_boldo_types
  stubs_counts = obj$stubs_counts
  stubs_states = obj$stubs_states
  stubs_types = obj$stubs_types
} else {
  # combined JA, SB and StubS datasets
  # remove "?" cells and collapse t cell subtypes in a new annotation
  jerby_clean = subset(x = jerby_se, cells = WhichCells(jerby_se, expression = cell.types=="?", invert = T))
  jerby_counts = GetAssayData(jerby_clean, layer = "counts", assay = "RNA")
  jerby_states = jerby_clean$cell.types
  # harmonize group names with sole-boldo
  jerby_types = jerby_states
  jerby_types = ifelse(jerby_types %in% c("T.cell", "T.CD8", "T.CD4"), "T cells", jerby_types)
  jerby_types = ifelse(jerby_types == "NK", "NK cells", jerby_types)
  jerby_types = ifelse(jerby_types == "Mal", "melanocytes/malignant cells", jerby_types)
  jerby_types = ifelse(jerby_types == "B.cell", "B cells", jerby_types)
  jerby_types = ifelse(jerby_types == "CAF", "fibroblasts", jerby_types)
  jerby_types = ifelse(jerby_types == "Endo.", "endothelial cells", jerby_types)
  jerby_types = ifelse(jerby_types == "Macrophage", "macrophages and dendritic cells", jerby_types)
  jerby_clean$celltype.broad = jerby_types

  sole_boldo_clean = sole_boldo_se
  sole_boldo_counts = GetAssayData(sole_boldo_clean, assay = "RNA", layer = "counts")
  sole_boldo_states = sole_boldo_clean$celltype
  sole_boldo_types = sole_boldo_clean$celltype.broad

  # group cells form StubS, too, subset for less imbalance and faster processing
  Idents(stub_schmidt_se$mela_ref) = "celltype_medium"
  stubs_clean = subset(stub_schmidt_se$mela_ref, downsample = 5000, seed = random_seed)
  stubs_counts = GetAssayData(stubs_clean, layer = "counts", assay = "RNA")
  stubs_states = stubs_clean$celltype_medium
  stubs_types = stubs_states
  stubs_types = ifelse(stubs_types %in% c("Tcyt_activated", "Tdp", "Teff_mem", "Treg"), "T cells", stubs_types)
  stubs_types = ifelse(stubs_types == "NK", "NK cells", stubs_types)
  stubs_types = ifelse(stubs_types %in% c("Mel_nc-like", "Mel_trans", "Mel_trans_melan"), "melanocytes/malignant cells", stubs_types)
  stubs_types = ifelse(stubs_types == "B_cells", "B cells", stubs_types)
  stubs_types = ifelse(stubs_types == "Plasma_cells", "Plasma cells", stubs_types)
  stubs_types = ifelse(stubs_types == "Mast_cells", "Mast cells", stubs_types)
  stubs_types = ifelse(stubs_types %in% c("Fb_1", "Fb_2", "Fb_3"), "fibroblasts", stubs_types)
  stubs_types = ifelse(stubs_types %in% c("LE", "VE_1", "VE_2"), "endothelial cells", stubs_types)
  stubs_types = ifelse(stubs_types %in% c("Mac", "pDC", "cDC_1", "cDC_2"), "macrophages and dendritic cells", stubs_types)
  stubs_types = ifelse(stubs_types %in% c("Kc_postmit", "Kc_premit"), "keratinocytes", stubs_types)
  stubs_types = ifelse(stubs_types %in% c("Pericyte_1", "Pericyte_2"), "pericytes", stubs_types)
  stubs_types = ifelse(stubs_types %in% c("Sebaceous_gland", "Eccrine_sweat_gland"), "Glandular cells", stubs_types)
  stubs_clean$celltype.broad = stubs_types

  # bayesprism functions for cleaning and computing markers (uses tests from scran internally)
  # sanity checks for correct orientation
  stopifnot("jerby counts matrix seems broken or transposed" = "MALAT1" %in% rownames(jerby_counts),
            "sole-boldo counts matrix seems broken or transposed" = "MALAT1" %in% rownames(sole_boldo_counts),
            "stubenvoll&schmidt counts matrix seems broken or transposed" = "MALAT1" %in% rownames(stubs_counts))
  jerby_counts <- cleanup.genes( input=as.matrix(t(jerby_counts)),
                                 input.type="count.matrix",
                                 species="hs",
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                 exp.cells=5)
  jerby_genes = colnames(jerby_counts)
  sole_boldo_counts <- cleanup.genes(input=as.matrix(t(sole_boldo_counts)),
                                     input.type="count.matrix",
                                     species="hs",
                                     gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                     exp.cells=5)
  sole_boldo_genes = colnames(sole_boldo_counts)
  stubs_counts <- cleanup.genes(input=as.matrix(t(stubs_counts)),
                                     input.type="count.matrix",
                                     species="hs",
                                     gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                     exp.cells=5)
  stubs_genes = colnames(stubs_counts)
  # try to find updates for non-shared genes in case symbols changed
  shared_genes = unique(c(jerby_genes, sole_boldo_genes, stubs_genes))
  shared_genes = shared_genes[((shared_genes %in% jerby_genes) + (shared_genes %in% sole_boldo_genes) + (shared_genes %in% stubs_genes))>=2]
  non.shared_genes_jerby = setdiff(jerby_genes, shared_genes)
  non.shared_genes_sole_boldo = setdiff(sole_boldo_genes, shared_genes)
  non.shared_genes_stubs = setdiff(stubs_genes, shared_genes)

  updatedsymbols_jerby = UpdateSymbolList(symbols = non.shared_genes_jerby, timeout = 30, verbose = T)
  updatedsymbols_jerby = setNames(updatedsymbols_jerby, nm = non.shared_genes_jerby)
  updatedsymbols_sole_boldo = UpdateSymbolList(symbols = non.shared_genes_sole_boldo, timeout = 30, verbose = T)
  updatedsymbols_sole_boldo = setNames(updatedsymbols_sole_boldo, nm = non.shared_genes_sole_boldo)
  updatedsymbols_stubs = UpdateSymbolList(symbols = non.shared_genes_stubs, timeout = 30, verbose = T)
  updatedsymbols_stubs = setNames(updatedsymbols_stubs, nm = non.shared_genes_stubs)

  # discards updates that would collide with an existing symbol and also completely discard two or more updates
  # that would result in a new name clash. Sometimes nomenclature is just weird. Those are usually not important anyways.
  # Finally, merges list with
  fix_and_merge_updates <- function(updates, non.updates){
    clashes.with.existing = updates %in% non.updates
    clashing.old.names = names(updates[clashes.with.existing])
    updates = updates[!clashes.with.existing]
    new.clashes = duplicated(updates)
    newly.clashing.oldnames = c()
    if (any(new.clashes)) {
      newly.clashing.names = unique(updates[new.clashes])
      clashing.matches = updates %in% newly.clashing.names
      newly.clashing.oldnames = names(updates[clashing.matches])
      updates = updates[!clashing.matches]
      newly.clashing.oldnames = setNames(newly.clashing.oldnames, nm = newly.clashing.oldnames)
      # as weird as it is, there can be "chains" of update clashes, so check again
      if (any(updates %in% newly.clashing.oldnames)){
        new.new.clashes = updates %in% newly.clashing.oldnames
        new.new.clashing.oldnames = names(updates[new.new.clashes])
        updates = updates[!new.new.clashes]
        newly.clashing.oldnames = c(newly.clashing.oldnames, setNames(new.new.clashing.oldnames, nm = new.new.clashing.oldnames))
      }
    }
    clashing.old.names = setNames(clashing.old.names, nm = clashing.old.names)
    non.updates = setNames(non.updates, nm = non.updates)
    final.updates = c(updates, non.updates, clashing.old.names, newly.clashing.oldnames)
    return(final.updates)
  }

  # update and clean again, in case that any of the changed symbols are in the cleaning list
  # also, it's practically free
  jerby_updates = fix_and_merge_updates(updatedsymbols_jerby, shared_genes)
  colnames(jerby_counts) = setNames(jerby_updates[colnames(jerby_counts)], nm = NULL)
  jerby_counts <- cleanup.genes( input=jerby_counts,
                                 input.type="count.matrix",
                                 species="hs",
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                 exp.cells=5)


  sole_boldo_updates = fix_and_merge_updates(updatedsymbols_sole_boldo, shared_genes)
  colnames(sole_boldo_counts) = setNames(sole_boldo_updates[colnames(sole_boldo_counts)], nm = NULL)
  sole_boldo_counts <- cleanup.genes(input=sole_boldo_counts,
                                     input.type="count.matrix",
                                     species="hs",
                                     gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                     exp.cells=5)

  stubs_updates = fix_and_merge_updates(updatedsymbols_stubs, shared_genes)
  colnames(stubs_counts) = setNames(stubs_updates[colnames(stubs_counts)], nm = NULL)
  stubs_counts <- cleanup.genes(input=stubs_counts,
                                 input.type="count.matrix",
                                 species="hs",
                                 gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                 exp.cells=5)

  saveRDS(list(jerby_counts = jerby_counts,
               jerby_states = jerby_states,
               jerby_types = jerby_types,
               sole_boldo_counts = sole_boldo_counts,
               sole_boldo_states = sole_boldo_states,
               sole_boldo_types = sole_boldo_types,
               stubs_counts = stubs_counts,
               stubs_states = stubs_states,
               stubs_types = stubs_types),
          file = "./data/genesets/JA_SB_StubS_harmonized_for_markersearch.rds")
}

# after restoring, we are sure to have
# jerby_counts
# jerby_states
# jerby_types
# sole_boldo_counts
# sole_boldo_states
# sole_boldo_types
# stubs_counts
# stubs_states
# stubs_types

# this is exactly the BayesPrism function, so all credit goes there and to scran
# for some reason, this does not throw but BayesPrism::get.exp.stat sometimes does on this machine

get.exp.stat <- function (sc.dat, cell.type.labels, cell.state.labels, psuedo.count = 0.1,
                          cell.count.cutoff = 50, n.cores = 1)
{
  ct.to.cst <- unique(cbind(cell.type = cell.type.labels, cell.state = cell.state.labels))
  cst.count.table <- table(cell.state.labels)
  low.count.cst <- names(cst.count.table)[cst.count.table <
                                            cell.count.cutoff]
  lib.size <- rowSums(sc.dat)
  lib.size <- lib.size/median(lib.size)
  dat.tmp <- sc.dat/lib.size
  dat.tmp <- log2(dat.tmp + psuedo.count) - log2(psuedo.count)
  fit.up <- scran::pairwiseTTests(x = t(dat.tmp), groups = cell.state.labels,
                           direction = "up", BPPARAM = MulticoreParam(n.cores))
  pairs.celltype.first <- ct.to.cst[match(fit.up$pairs$first,
                                          ct.to.cst[, "cell.state"]), "cell.type"]
  pairs.celltype.second <- ct.to.cst[match(fit.up$pairs$second,
                                           ct.to.cst[, "cell.state"]), "cell.type"]
  filter.idx <- pairs.celltype.first != pairs.celltype.second &
    !fit.up$pairs$second %in% low.count.cst
  fit.up[[1]] <- fit.up[[1]][filter.idx]
  fit.up[[2]] <- fit.up[[2]][filter.idx, ]
  output.up <- scran::combineMarkers(fit.up$statistics, fit.up$pairs,
                              pval.type = "all", min.prop = NULL, log.p.in = F, log.p.out = F,
                              full.stats = F, pval.field = "p.value", effect.field = "logFC",
                              sorted = F)
  all.ct <- unique(ct.to.cst[, "cell.type"])
  ct.stat.list <- lapply(all.ct, function(ct.i) {
    cst.i <- ct.to.cst[ct.to.cst[, "cell.type"] == ct.i,
                       "cell.state"]
    output.up.i <- output.up[cst.i]
    pval.up.i <- do.call(cbind, lapply(output.up.i, "[",
                                       "p.value"))
    pval.up.min.i <- apply(pval.up.i, 1, min)
    lfc.i <- apply(do.call(cbind, lapply(output.up.i, function(output.up.i.j) {
      apply(output.up.i.j[, grepl("logFC", colnames(output.up.i.j)),
                          drop = F], 1, min)
    })), 1, max)
    data.frame(pval.up.min = pval.up.min.i, min.lfc = lfc.i)
  })
  names(ct.stat.list) <- all.ct
  return(ct.stat.list)
}

# find markers for groups
jerby_exp_stats = get.exp.stat(sc.dat = jerby_counts,
                               cell.type.labels = jerby_types,
                               cell.state.labels = jerby_states,
                               n.cores = 30)

sole_boldo_exp_stats = get.exp.stat(sc.dat = sole_boldo_counts,
                                    cell.type.labels = sole_boldo_types,
                                    cell.state.labels = sole_boldo_states,
                                    n.cores = 30)

stubs_exp_stats = get.exp.stat(sc.dat = stubs_counts,
                                    cell.type.labels = stubs_types,
                                    cell.state.labels = stubs_states,
                                    n.cores = 30)

# filter genes above a minimal lfc and with pval below threshold in both datsets if possible
# or one if these cells only appear in one dataset
lfc.threshold = 0.3
pval.threshold = 0.01
all.types = unique(c(names(jerby_exp_stats), names(sole_boldo_exp_stats), names(stubs_exp_stats)))

combine.2.dfs = function(df1, df2){
  if (is.null(df1)) {return(df2 %>% filter((pval.up.min<pval.threshold)&(min.lfc>lfc.threshold)) %>% arrange(desc(min.lfc)))}
  if (is.null(df2)) {return(df1 %>% filter((pval.up.min<pval.threshold)&(min.lfc>lfc.threshold)) %>% arrange(desc(min.lfc)))}
  # else combine
  # subset to shared entries and arrange in the same order
  shared.entries = intersect(rownames(df1), rownames(df2))
  df1 = df1[shared.entries, ]
  df2 = df2[shared.entries, ]
  # do filtering
  df3 = data.frame(pval.up.min = pmax(df1$pval.up.min, df2$pval.up.min),
                   min.lfc = pmin(df1$min.lfc, df2$min.lfc),
                   row.names = shared.entries)
  return(df3 %>% filter((pval.up.min<pval.threshold)&(min.lfc>lfc.threshold)) %>% arrange(desc(min.lfc)))
}

combine.3.dfs = function(df1, df2, df3){
  if (is.null(df1)) {return(combine.2.dfs(df2, df3))}
  if (is.null(df2)) {return(combine.2.dfs(df1, df3))}
  if (is.null(df3)) {return(combine.2.dfs(df1, df2))}
  # all filled, need a 3-way merge
  shared.entries = unique(c(rownames(df1), rownames(df2), rownames(df3)))
  shared.entries = shared.entries[((shared.entries %in% rownames(df1)) + (shared.entries %in% rownames(df2)) + (shared.entries %in% rownames(df3)))>=2]

  # subset to shared entries and arrange in the same order
  shared.entries = intersect(rownames(df1), rownames(df2))
  df1 = df1[shared.entries, ]
  df2 = df2[shared.entries, ]
  df3 = df3[shared.entries, ]
  # do filtering
  df4 = data.frame(pval.up.min = pmax(df1$pval.up.min, df2$pval.up.min, df3$pval.up.min, na.rm = T),
                   min.lfc = pmin(df1$min.lfc, df2$min.lfc, df3$min.lfc, na.rm = T),
                   row.names = shared.entries)
  return(df4 %>% filter((pval.up.min<pval.threshold)&(min.lfc>lfc.threshold)) %>% arrange(desc(min.lfc)))
}


combined.marker.stats = lapply(setNames(all.types, nm = all.types),
                               function(t)combine.3.dfs(jerby_exp_stats[[t]], sole_boldo_exp_stats[[t]], stubs_exp_stats[[t]]))

# for cell types present in only one dataset, it would (in theory) be possible that different cell types
# have shared markers. If any such cases are present, remove them.
combined.marker.stats.cleaned = lapply(seq_along(combined.marker.stats), function(i){
  df = combined.marker.stats[[i]]
  genes = rownames(df)
  others.genes = unique(unlist(lapply(combined.marker.stats[-i], rownames), use.names = F))
  ambiguous.markers = genes %in% others.genes
  return(df[!ambiguous.markers, ])
})
names(combined.marker.stats.cleaned) = names(combined.marker.stats)

# additionally add foldchange and percent_expressed in cell of that type and other cells as annotations
# computes over the all cells from sole-boldo AND jerby-arnon datasets
get.percent.expressed = function(genes, celltype, fill = NA){
  # non-existant genes are assumed to have 0 counts everywhere
  jerby.measured = intersect(genes, colnames(jerby_counts))
  jerby.missing = setdiff(genes, jerby.measured)
  jerby_counts_sub = cbind(jerby_counts[,jerby.measured],
                           matrix(data = fill, ncol = length(jerby.missing), nrow = nrow(jerby_counts),
                                  dimnames = list(rownames(jerby_counts), jerby.missing)))
  sole_boldo.measured = intersect(genes, colnames(sole_boldo_counts))
  sole_boldo.missing = setdiff(genes, sole_boldo.measured)
  sole_boldo_counts_sub = cbind(sole_boldo_counts[,sole_boldo.measured],
                                matrix(data = fill, ncol = length(sole_boldo.missing), nrow = nrow(sole_boldo_counts),
                                       dimnames = list(rownames(sole_boldo_counts), sole_boldo.missing)))
  stubs.measured = intersect(genes, colnames(stubs_counts))
  stubs.missing = setdiff(genes, stubs.measured)
  stubs_counts_sub = cbind(stubs_counts[,stubs.measured],
                                matrix(data = fill, ncol = length(stubs.missing), nrow = nrow(stubs_counts),
                                       dimnames = list(rownames(stubs_counts), stubs.missing)))
  counts = rbind(jerby_counts_sub[,genes], sole_boldo_counts_sub[,genes], stubs_counts_sub[,genes])
  cell.types = c(jerby_types, sole_boldo_types, stubs_types)
  in.group = cell.types == celltype
  out.group = !in.group
  return(data.frame(row.names = genes,
                    pct.expressed.celltype = colMeans(counts[in.group,]>0, na.rm = T),
                    pct.expressed.other = colMeans(counts[out.group,]>0, na.rm = T)))
}


combined.marker.stats.annotated = lapply(seq_along(combined.marker.stats.cleaned), function(i){
  nm = names(combined.marker.stats)[i]
  df = combined.marker.stats.cleaned[[i]]
  genes = rownames(df)
  return(cbind(df, get.percent.expressed(genes, nm)))
})
names(combined.marker.stats.annotated) = names(combined.marker.stats)

# filter detection rate > 0.2 in cell type of interest, < 0.2 in other cells
# and at least 4-fold higher in the cell type of interest than in background
combined.marker.stats.filtered = lapply(combined.marker.stats.annotated,
                                        function(df){df %>%
                                            filter((pct.expressed.celltype > 0.2) &
                                                     (pct.expressed.other < 0.2) &
                                                     (pct.expressed.celltype > 4*pct.expressed.other))})

saveRDS(combined.marker.stats.filtered, file = "./data/genesets/collected_markers_stats.rds")
# take 30 markers max
final.markers = lapply(combined.marker.stats.filtered, function(x)rownames(x) %>% head(n = 30))
print("Identified markers:")
print(final.markers)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# update all marker lists, save in one spot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

saveRDS(list(panglao_gs = panglao_gs, final.markers = final.markers), file = "./data/genesets/collected_markers.rds")

print("Done!")
