# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("httr")
library("tidyverse")
library("Seurat")
library("yaml")
library("naturalsort")
library("reshape2")

library("liana")

source("code/helper/FG_helper_NICHES.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

lr_db = get_combined_lr_dbs()
all_genes = unique(c(lr_db$to, lr_db$from))

# there is only a handful of recent gene symbols with underscores
# I hope this stays up-to-date for long
# we filter those out before separating complex components on underscores

genes_with_underscores = c("C4B_2", "GTF2H2C_2", "APOBEC3A_B")
gene_w_underscore_in_list = genes_with_underscores %in% all_genes

all_genes = setdiff(all_genes, genes_with_underscores)
all_complexes = unique(all_genes[str_detect(all_genes, pattern = "_")])

all_genes = unique(unlist(str_split(all_genes, pattern = "_")))
all_genes = c(all_genes, genes_with_underscores[gene_w_underscore_in_list])

###############################
### GTEx
###############################

gtex9 = readRDS("./assets/GTEx_v9_snRNA_seq_Eraslan_et_al_2022/GTEx_8_tissues_snRNAseq_atlas_071421.converted_to.Seurat.subset_to_tissue_is_skin.Rds")
gtex9 = subset(x = gtex9, subset = Broad.cell.type != "Unknown")
# table(droplevels(gtex9$Broad.cell.type))

present_genes = intersect(all_genes, rownames(gtex9))
missing_genes = setdiff(all_genes, present_genes)

missing_genes_lookup = resolve_translations(missing_genes, prefer = rownames(gtex9))
# if any could be translated, rename the gtex features to match
missing_genes_lookup = missing_genes_lookup[missing_genes_lookup!=names(missing_genes_lookup)]
updated.genenames = rownames(gtex9)
for (nm in names(missing_genes_lookup)) {
  gtexname = missing_genes_lookup[nm]
  db_name = nm
  # if db_name is also present, this is an error, skip
  if (db_name %in% rownames(gtex9)) {
    warning(paste("Changing symbol",gtexname,"to",db_name,"from LR DB for gtex9 would lead to a clash."))
  } else {
    updated.genenames[match(gtexname, updated.genenames)] = db_name
  }
}
rownames(gtex9@assays$originalexp@counts) = updated.genenames
rownames(gtex9@assays$originalexp@data) = updated.genenames
rownames(gtex9@assays$originalexp@meta.features) = updated.genenames

present_genes = intersect(all_genes, rownames(gtex9))
missing_genes = setdiff(all_genes, present_genes)

# compute, clean and add zero as expressions for genes not expressed in this dataset
log.avg.normalized.expression = AverageExpression(gtex9, features = present_genes, group.by = "Broad.cell.type", slot = "data")$originalexp
colnames(log.avg.normalized.expression) = setNames(colnames(log.avg.normalized.expression), nm = NULL)
log.avg.normalized.expression = rbind(log.avg.normalized.expression,
                                      matrix(data = 0., nrow = length(missing_genes), ncol = ncol(log.avg.normalized.expression),
                                             dimnames = list(missing_genes, colnames(log.avg.normalized.expression))))
# handle complexes
# log.avg.normalized.expression -> geomean
log.avg.normalized.expression_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = log.avg.normalized.expression[parts,]
  return(apply(dat, 2, geomean))
}))
log.avg.normalized.expression = rbind(log.avg.normalized.expression, log.avg.normalized.expression_complexes)


# percent of cells in a type with non-zero expression
gtex9 = add_binarized_counts(gtex9)
fraction.expressed = AverageExpression(gtex9, features = present_genes, group.by = "Broad.cell.type", slot = "counts")$binarycounts
colnames(fraction.expressed) = setNames(colnames(fraction.expressed), nm = NULL)
fraction.expressed = rbind(fraction.expressed,
                           matrix(data = 0., nrow = length(missing_genes), ncol = ncol(fraction.expressed),
                                  dimnames = list(missing_genes, colnames(fraction.expressed))))
# handle complexes
# fraction.expressed -> geomean
fraction.expressed_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = fraction.expressed[parts,]
  return(apply(dat, 2, geomean))
}))
fraction.expressed = rbind(fraction.expressed, fraction.expressed_complexes)

# empirical quantiles of expression per celltype
summed.counts.per.type = Seurat:::PseudobulkExpression(gtex9, layer = "counts", assay = "originalexp", method = "aggregate", group.by = "Broad.cell.type")$originalexp
per.celltype.ecdf.function = lapply(seq_along(colnames(summed.counts.per.type)), FUN = function(i){ecdf(summed.counts.per.type[,i])})
names(per.celltype.ecdf.function) = colnames(summed.counts.per.type)

per.celltype.empirical.quantile = sapply(seq_along(colnames(summed.counts.per.type)), FUN = function(i){
  x = summed.counts.per.type[,i]
  qfun = per.celltype.ecdf.function[[i]]
  return(qfun(x))
})
rownames(per.celltype.empirical.quantile) = rownames(summed.counts.per.type)
colnames(per.celltype.empirical.quantile) = colnames(summed.counts.per.type)
per.celltype.empirical.quantile = rbind(per.celltype.empirical.quantile[present_genes,],
                                        matrix(data = apply(per.celltype.empirical.quantile, 2, min), # we can safely assume that each celltype has at least one gene with zero counts -> min value of ecdf even for non-observed genes
                                               byrow = T, nrow = length(missing_genes), ncol = ncol(per.celltype.empirical.quantile),
                                               dimnames = list(missing_genes, colnames(per.celltype.empirical.quantile))))
# handle complexes
# per.celltype.empirical.quantile -> expm1(geomean(log(counts))) of complex components -> then recompute with ecdf calibrated on observed counts
# need to pad zeros to observed counts matrix for non-observed genes in the query
padded_counts = rbind(summed.counts.per.type,
                      matrix(data = 0, nrow = length(missing_genes), ncol = ncol(summed.counts.per.type),
                             dimnames = list(missing_genes, colnames(summed.counts.per.type))))
per.celltype.empirical.quantile_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = padded_counts[parts,]
  return(expm1(apply(log1p(dat), 2, geomean)))
}))
# apply ecdf of each celltype separately
per.celltype.empirical.quantile_complexes = sapply(seq_along(colnames(per.celltype.empirical.quantile_complexes)), FUN = function(i){
  x = per.celltype.empirical.quantile_complexes[,i]
  qfun = per.celltype.ecdf.function[[i]]
  return(qfun(x))
})
rownames(per.celltype.empirical.quantile_complexes) = all_complexes
colnames(per.celltype.empirical.quantile_complexes) = colnames(padded_counts)
per.celltype.empirical.quantile = rbind(per.celltype.empirical.quantile, per.celltype.empirical.quantile_complexes)

# also empirical quantiles, but over all celltypes after normalizing by total counts
across.celltype.normalized.counts = apply(summed.counts.per.type, 2, FUN = function(x)x/sum(x))
across.normalized.celltypes.ecdf = ecdf(as.vector(across.celltype.normalized.counts))
across.normalized.celltypes.empirical.quantiles = apply(across.celltype.normalized.counts, 2, FUN = function(x)across.normalized.celltypes.ecdf(x))
rownames(across.normalized.celltypes.empirical.quantiles) = rownames(summed.counts.per.type)
across.normalized.celltypes.empirical.quantiles = rbind(across.normalized.celltypes.empirical.quantiles[present_genes,],
                                        matrix(data = min(across.normalized.celltypes.empirical.quantiles), # same deal, there is at least one zero in the whole dataset which represents out minimum value for non-observed genes
                                               nrow = length(missing_genes), ncol = ncol(across.normalized.celltypes.empirical.quantiles),
                                               dimnames = list(missing_genes, colnames(across.normalized.celltypes.empirical.quantiles))))
# handle complexes
# across.normalized.celltypes.empirical.quantiles -> expm1(geomean(log(counts))) of complex components -> then recompute with ecdf calibrated on observed counts
# need to pad zeros to observed counts matrix for non-observed genes in the query
padded_counts = rbind(across.celltype.normalized.counts,
                      matrix(data = 0, nrow = length(missing_genes), ncol = ncol(across.celltype.normalized.counts),
                             dimnames = list(missing_genes, colnames(across.celltype.normalized.counts))))
across.normalized.celltypes.empirical.quantiles_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = padded_counts[parts,]
  geomeandat = expm1(apply(log1p(dat), 2, geomean))
  return(across.normalized.celltypes.ecdf(geomeandat))
}))
across.normalized.celltypes.empirical.quantiles = rbind(across.normalized.celltypes.empirical.quantiles, across.normalized.celltypes.empirical.quantiles_complexes)


# save
expression.annotations = list(log.avg.normalized.expression = log.avg.normalized.expression,
                              fraction.expressed = fraction.expressed,
                              per.celltype.empirical.quantile = per.celltype.empirical.quantile,
                              per.celltype.ecdf.function = per.celltype.ecdf.function,
                              across.normalized.celltypes.empirical.quantiles = across.normalized.celltypes.empirical.quantiles,
                              across.normalized.celltypes.ecdf = across.normalized.celltypes.ecdf)

manifest = yaml.load_file("./assets/manifest.yaml")
destination = manifest$workdata_parent
saveRDS(expression.annotations, file = file.path(destination, "GTEx9_LRI_gene_expression_annos.Rds"))

###############################
### Jerby-Arnon
###############################

# assume that the Jerby-Arnon et al. dataset has been properly prepared in step 05
jerby.se = readRDS("./data/genesets/jerby_se.rds")
jerby.se = subset(x = jerby.se, subset = cell.types != "?")
# table(droplevels(jerby.se$cell.types))

present_genes = intersect(all_genes, rownames(jerby.se))
missing_genes = setdiff(all_genes, present_genes)

missing_genes_lookup = resolve_translations(missing_genes, prefer = rownames(jerby.se))
# if any could be translated, rename the gtex features to match
missing_genes_lookup = missing_genes_lookup[missing_genes_lookup!=names(missing_genes_lookup)]
updated.genenames = rownames(jerby.se)
for (nm in names(missing_genes_lookup)) {
  jerbyname = missing_genes_lookup[nm]
  db_name = nm
  # if db_name is also present, this is an error, skip
  if (db_name %in% rownames(jerby.se)) {
    warning(paste("Changing symbol",jerbyname,"to",db_name,"from LR DB for Jerby-Arnon would lead to a clash."))
  } else {
    updated.genenames[match(jerbyname, updated.genenames)] = db_name
  }
}
rownames(jerby.se@assays$RNA) = updated.genenames
rownames(jerby.se@assays$TPM) = updated.genenames
jerby.se = NormalizeData(jerby.se, assay = "RNA")
jerby.se = ScaleData(jerby.se, assay = "RNA")

# recalc
present_genes = intersect(all_genes, rownames(jerby.se))
missing_genes = setdiff(all_genes, present_genes)

# compute, clean and add zero as expressions for genes not expressed in this dataset
log.avg.normalized.expression = AverageExpression(jerby.se, features = present_genes, group.by = "cell.types", layer = "data")$RNA
colnames(log.avg.normalized.expression) = setNames(colnames(log.avg.normalized.expression), nm = NULL)
log.avg.normalized.expression = rbind(log.avg.normalized.expression,
                                      matrix(data = 0., nrow = length(missing_genes), ncol = ncol(log.avg.normalized.expression),
                                             dimnames = list(missing_genes, colnames(log.avg.normalized.expression))))
# handle complexes
# log.avg.normalized.expression -> geomean
log.avg.normalized.expression_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = log.avg.normalized.expression[parts,]
  return(apply(dat, 2, geomean))
}))
log.avg.normalized.expression = rbind(log.avg.normalized.expression, log.avg.normalized.expression_complexes)


# percent of cells in a type with non-zero expression
jerby.se = add_binarized_counts(jerby.se, assay = "RNA")
fraction.expressed = AverageExpression(jerby.se, features = present_genes, group.by = "cell.types", layer = "counts")$binarycounts
colnames(fraction.expressed) = setNames(colnames(fraction.expressed), nm = NULL)
fraction.expressed = rbind(fraction.expressed,
                           matrix(data = 0., nrow = length(missing_genes), ncol = ncol(fraction.expressed),
                                  dimnames = list(missing_genes, colnames(fraction.expressed))))
# handle complexes
# fraction.expressed -> geomean
fraction.expressed_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = fraction.expressed[parts,]
  return(apply(dat, 2, geomean))
}))
fraction.expressed = rbind(fraction.expressed, fraction.expressed_complexes)

# empirical quantiles of expression per celltype
summed.counts.per.type = Seurat:::PseudobulkExpression(jerby.se, layer = "counts", assay = "RNA", method = "aggregate", group.by = "cell.types")$RNA
per.celltype.ecdf.function = lapply(seq_along(colnames(summed.counts.per.type)), FUN = function(i){ecdf(summed.counts.per.type[,i])})
names(per.celltype.ecdf.function) = colnames(summed.counts.per.type)

per.celltype.empirical.quantile = sapply(seq_along(colnames(summed.counts.per.type)), FUN = function(i){
  x = summed.counts.per.type[,i]
  qfun = per.celltype.ecdf.function[[i]]
  return(qfun(x))
})
rownames(per.celltype.empirical.quantile) = rownames(summed.counts.per.type)
colnames(per.celltype.empirical.quantile) = colnames(summed.counts.per.type)
per.celltype.empirical.quantile = rbind(per.celltype.empirical.quantile[present_genes,],
                                        matrix(data = apply(per.celltype.empirical.quantile, 2, min), # we can safely assume that each celltype has at least one gene with zero counts -> min value of ecdf even for non-observed genes
                                               byrow = T, nrow = length(missing_genes), ncol = ncol(per.celltype.empirical.quantile),
                                               dimnames = list(missing_genes, colnames(per.celltype.empirical.quantile))))
# handle complexes
# per.celltype.empirical.quantile -> expm1(geomean(log(counts))) of complex components -> then recompute with ecdf calibrated on observed counts
# need to pad zeros to observed counts matrix for non-observed genes in the query
padded_counts = rbind(summed.counts.per.type,
                      matrix(data = 0, nrow = length(missing_genes), ncol = ncol(summed.counts.per.type),
                             dimnames = list(missing_genes, colnames(summed.counts.per.type))))
per.celltype.empirical.quantile_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = padded_counts[parts,]
  return(expm1(apply(log1p(dat), 2, geomean)))
}))
# apply ecdf of each celltype separately
per.celltype.empirical.quantile_complexes = sapply(seq_along(colnames(per.celltype.empirical.quantile_complexes)), FUN = function(i){
  x = per.celltype.empirical.quantile_complexes[,i]
  qfun = per.celltype.ecdf.function[[i]]
  return(qfun(x))
})
rownames(per.celltype.empirical.quantile_complexes) = all_complexes
colnames(per.celltype.empirical.quantile_complexes) = colnames(padded_counts)
per.celltype.empirical.quantile = rbind(per.celltype.empirical.quantile, per.celltype.empirical.quantile_complexes)

# also empirical quantiles, but over all celltypes after normalizing by total counts
across.celltype.normalized.counts = apply(summed.counts.per.type, 2, FUN = function(x)x/sum(x))
across.normalized.celltypes.ecdf = ecdf(as.vector(across.celltype.normalized.counts))
across.normalized.celltypes.empirical.quantiles = apply(across.celltype.normalized.counts, 2, FUN = function(x)across.normalized.celltypes.ecdf(x))
rownames(across.normalized.celltypes.empirical.quantiles) = rownames(summed.counts.per.type)
across.normalized.celltypes.empirical.quantiles = rbind(across.normalized.celltypes.empirical.quantiles[present_genes,],
                                                        matrix(data = min(across.normalized.celltypes.empirical.quantiles), # same deal, there is at least one zero in the whole dataset which represents out minimum value for non-observed genes
                                                               nrow = length(missing_genes), ncol = ncol(across.normalized.celltypes.empirical.quantiles),
                                                               dimnames = list(missing_genes, colnames(across.normalized.celltypes.empirical.quantiles))))
# handle complexes
# across.normalized.celltypes.empirical.quantiles -> expm1(geomean(log(counts))) of complex components -> then recompute with ecdf calibrated on observed counts
# need to pad zeros to observed counts matrix for non-observed genes in the query
padded_counts = rbind(across.celltype.normalized.counts,
                      matrix(data = 0, nrow = length(missing_genes), ncol = ncol(across.celltype.normalized.counts),
                             dimnames = list(missing_genes, colnames(across.celltype.normalized.counts))))
across.normalized.celltypes.empirical.quantiles_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = padded_counts[parts,]
  geomeandat = expm1(apply(log1p(dat), 2, geomean))
  return(across.normalized.celltypes.ecdf(geomeandat))
}))
across.normalized.celltypes.empirical.quantiles = rbind(across.normalized.celltypes.empirical.quantiles, across.normalized.celltypes.empirical.quantiles_complexes)


# save
expression.annotations = list(log.avg.normalized.expression = log.avg.normalized.expression,
                              fraction.expressed = fraction.expressed,
                              per.celltype.empirical.quantile = per.celltype.empirical.quantile,
                              per.celltype.ecdf.function = per.celltype.ecdf.function,
                              across.normalized.celltypes.empirical.quantiles = across.normalized.celltypes.empirical.quantiles,
                              across.normalized.celltypes.ecdf = across.normalized.celltypes.ecdf)

manifest = yaml.load_file("./assets/manifest.yaml")
destination = manifest$workdata_parent
saveRDS(expression.annotations, file = file.path(destination, "jerby.se_LRI_gene_expression_annos.Rds"))

###############################
### Stubenvoll&Schmidt
###############################

# assume that the Stubenvoll, Schmidt et al. dataset has been properly prepared in step 05
stubs.se = readRDS("./data/sc_ref_Stubenvoll_Schmidt_2025.rds")
# take melanoma reference
stubs.se = stubs.se$mela_ref

present_genes = intersect(all_genes, rownames(stubs.se))
missing_genes = setdiff(all_genes, present_genes)

missing_genes_lookup = resolve_translations(missing_genes, prefer = rownames(stubs.se@assays$RNA))
# if any could be translated, rename the gtex features to match
missing_genes_lookup = missing_genes_lookup[missing_genes_lookup!=names(missing_genes_lookup)]
if (length(missing_genes_lookup)>0) {
  updated.genenames = rownames(stubs.se@assays$RNA)
  for (nm in names(missing_genes_lookup)) {
    stubsname = missing_genes_lookup[nm]
    db_name = nm
    # if db_name is also present, this is an error, skip
    if (db_name %in% rownames(stubs.se@assays$RNA)) {
      warning(paste("Changing symbol",stubsname,"to",db_name,"from LR DB for Stubenvoll-Schmidt would lead to a clash."))
    } else {
      updated.genenames[match(stubsname, updated.genenames)] = db_name
    }
  }
  rownames(stubs.se@assays$RNA) = updated.genenames
}
stubs.se = NormalizeData(stubs.se, assay = "RNA")
stubs.se = ScaleData(stubs.se, assay = "RNA")

# recalc
present_genes = intersect(all_genes, rownames(stubs.se))
missing_genes = setdiff(all_genes, present_genes)

# compute, clean and add zero as expressions for genes not expressed in this dataset
log.avg.normalized.expression = AverageExpression(stubs.se, features = present_genes, group.by = "celltype_heatmap_anno", layer = "data")$RNA
colnames(log.avg.normalized.expression) = setNames(colnames(log.avg.normalized.expression), nm = NULL)
log.avg.normalized.expression = rbind(log.avg.normalized.expression,
                                      matrix(data = 0., nrow = length(missing_genes), ncol = ncol(log.avg.normalized.expression),
                                             dimnames = list(missing_genes, colnames(log.avg.normalized.expression))))
# handle complexes
# log.avg.normalized.expression -> geomean
log.avg.normalized.expression_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = log.avg.normalized.expression[parts,]
  return(apply(dat, 2, geomean))
}))
log.avg.normalized.expression = rbind(log.avg.normalized.expression, log.avg.normalized.expression_complexes)


# percent of cells in a type with non-zero expression
stubs.se = add_binarized_counts(stubs.se, assay = "RNA")
fraction.expressed = AverageExpression(stubs.se, features = present_genes, group.by = "celltype_heatmap_anno", layer = "counts")$binarycounts
colnames(fraction.expressed) = setNames(colnames(fraction.expressed), nm = NULL)
fraction.expressed = rbind(fraction.expressed,
                           matrix(data = 0., nrow = length(missing_genes), ncol = ncol(fraction.expressed),
                                  dimnames = list(missing_genes, colnames(fraction.expressed))))
# handle complexes
# fraction.expressed -> geomean
fraction.expressed_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = fraction.expressed[parts,]
  return(apply(dat, 2, geomean))
}))
fraction.expressed = rbind(fraction.expressed, fraction.expressed_complexes)

# empirical quantiles of expression per celltype
summed.counts.per.type = Seurat:::PseudobulkExpression(stubs.se, layer = "counts", assay = "RNA", method = "aggregate", group.by = "celltype_heatmap_anno")$RNA
per.celltype.ecdf.function = lapply(seq_along(colnames(summed.counts.per.type)), FUN = function(i){ecdf(summed.counts.per.type[,i])})
names(per.celltype.ecdf.function) = colnames(summed.counts.per.type)

per.celltype.empirical.quantile = sapply(seq_along(colnames(summed.counts.per.type)), FUN = function(i){
  x = summed.counts.per.type[,i]
  qfun = per.celltype.ecdf.function[[i]]
  return(qfun(x))
})
rownames(per.celltype.empirical.quantile) = rownames(summed.counts.per.type)
colnames(per.celltype.empirical.quantile) = colnames(summed.counts.per.type)
per.celltype.empirical.quantile = rbind(per.celltype.empirical.quantile[present_genes,],
                                        matrix(data = apply(per.celltype.empirical.quantile, 2, min), # we can safely assume that each celltype has at least one gene with zero counts -> min value of ecdf even for non-observed genes
                                               byrow = T, nrow = length(missing_genes), ncol = ncol(per.celltype.empirical.quantile),
                                               dimnames = list(missing_genes, colnames(per.celltype.empirical.quantile))))
# handle complexes
# per.celltype.empirical.quantile -> expm1(geomean(log(counts))) of complex components -> then recompute with ecdf calibrated on observed counts
# need to pad zeros to observed counts matrix for non-observed genes in the query
padded_counts = rbind(summed.counts.per.type,
                      matrix(data = 0, nrow = length(missing_genes), ncol = ncol(summed.counts.per.type),
                             dimnames = list(missing_genes, colnames(summed.counts.per.type))))
per.celltype.empirical.quantile_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = padded_counts[parts,]
  return(expm1(apply(log1p(dat), 2, geomean)))
}))
# apply ecdf of each celltype separately
per.celltype.empirical.quantile_complexes = sapply(seq_along(colnames(per.celltype.empirical.quantile_complexes)), FUN = function(i){
  x = per.celltype.empirical.quantile_complexes[,i]
  qfun = per.celltype.ecdf.function[[i]]
  return(qfun(x))
})
rownames(per.celltype.empirical.quantile_complexes) = all_complexes
colnames(per.celltype.empirical.quantile_complexes) = colnames(padded_counts)
per.celltype.empirical.quantile = rbind(per.celltype.empirical.quantile, per.celltype.empirical.quantile_complexes)

# also empirical quantiles, but over all celltypes after normalizing by total counts
across.celltype.normalized.counts = apply(summed.counts.per.type, 2, FUN = function(x)x/sum(x))
across.normalized.celltypes.ecdf = ecdf(as.vector(across.celltype.normalized.counts))
across.normalized.celltypes.empirical.quantiles = apply(across.celltype.normalized.counts, 2, FUN = function(x)across.normalized.celltypes.ecdf(x))
rownames(across.normalized.celltypes.empirical.quantiles) = rownames(summed.counts.per.type)
across.normalized.celltypes.empirical.quantiles = rbind(across.normalized.celltypes.empirical.quantiles[present_genes,],
                                                        matrix(data = min(across.normalized.celltypes.empirical.quantiles), # same deal, there is at least one zero in the whole dataset which represents out minimum value for non-observed genes
                                                               nrow = length(missing_genes), ncol = ncol(across.normalized.celltypes.empirical.quantiles),
                                                               dimnames = list(missing_genes, colnames(across.normalized.celltypes.empirical.quantiles))))
# handle complexes
# across.normalized.celltypes.empirical.quantiles -> expm1(geomean(log(counts))) of complex components -> then recompute with ecdf calibrated on observed counts
# need to pad zeros to observed counts matrix for non-observed genes in the query
padded_counts = rbind(across.celltype.normalized.counts,
                      matrix(data = 0, nrow = length(missing_genes), ncol = ncol(across.celltype.normalized.counts),
                             dimnames = list(missing_genes, colnames(across.celltype.normalized.counts))))
across.normalized.celltypes.empirical.quantiles_complexes = t(sapply(all_complexes, function(comp){
  parts = str_split(comp, pattern = "_", simplify = T)
  dat = padded_counts[parts,]
  geomeandat = expm1(apply(log1p(dat), 2, geomean))
  return(across.normalized.celltypes.ecdf(geomeandat))
}))
across.normalized.celltypes.empirical.quantiles = rbind(across.normalized.celltypes.empirical.quantiles, across.normalized.celltypes.empirical.quantiles_complexes)


# save
expression.annotations = list(log.avg.normalized.expression = log.avg.normalized.expression,
                              fraction.expressed = fraction.expressed,
                              per.celltype.empirical.quantile = per.celltype.empirical.quantile,
                              per.celltype.ecdf.function = per.celltype.ecdf.function,
                              across.normalized.celltypes.empirical.quantiles = across.normalized.celltypes.empirical.quantiles,
                              across.normalized.celltypes.ecdf = across.normalized.celltypes.ecdf)

saveRDS(expression.annotations, file = file.path(destination, "stubs.se_LRI_gene_expression_annos.Rds"))

