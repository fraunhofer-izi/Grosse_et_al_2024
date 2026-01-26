# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Once qc cutoffs are final, apply them and store subsetted seurat objects
# Take:
#   02_seurat.rds
#   manifest$QC_cutoffs
# Store:
#   04_seurat.rds
#   04_extended_clustering.rds (the same as above with multiple possible clustering variants)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("tidyverse")
library("Seurat")
library("yaml")
library("naturalsort")
library("reshape2")
library("scico")
library("semla")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
all_samples = c(manifest$samples, manifest$external_samples)
random_seed = 1234

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Meaningful choices are in the manifest
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
genes_cutoff <- manifest$QC_cutoffs$genes
counts_cutoff <- manifest$QC_cutoffs$counts
mito_cutoff <- manifest$QC_cutoffs$mito

for (sample_id in all_samples) {

  print(paste("Active sample is now:", sample_id))
  sample_manifest = manifest[[sample_id]]
  work.dir = file.path(manifest$workdata_parent, sample_manifest$workdata)
  se = readRDS(file.path(work.dir, "02_seurat.rds"))

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Add selection criteria to the whole dataset
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  se[["genes_exceed_per_spot_cutoff"]] <- se$nFeature_Spatial >= genes_cutoff
  se[["counts_exceed_per_spot_cutoff"]] <- se$nCount_Spatial >= counts_cutoff
  se[["perc_of_mito_genes_below_cutoff"]] <- se$Perc_of_mito_spot <= mito_cutoff

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Make a spot subset
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  print(paste("Subsetting sample:", sample_id))
  se_qc_subset <- SubsetSTData(object = se,
                               expression =
                                 genes_exceed_per_spot_cutoff &
                                 counts_exceed_per_spot_cutoff &
                                 perc_of_mito_genes_below_cutoff
                               )

  # store number of genes and cells after QC cutoff
  Misc(se_qc_subset, "post_QC_dims_before_RiboMito_separation") <- dim(se_qc_subset[["Spatial"]])

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Regenerate scaling, SCT, PCA and UMAPs
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  print(paste("Recompute assays:", sample_id))

  # split off mitochondrial and ribosomal genes into a separate assay,
  # because we usually don't care about those
  mito.ribo.counts = GetAssayData(se_qc_subset, assay = "Spatial", layer = "counts")
  mito.ribo.counts = mito.ribo.counts[grep(rownames(se_qc_subset), pattern = "^M?RPL|^M?RPS|^MT-", value = T),]
  se_qc_subset[["MitoRibo"]] = CreateAssayObject(counts = mito.ribo.counts)

  countsmat = GetAssayData(se_qc_subset, assay = "Spatial", layer = "counts")
  countsmat = countsmat[setdiff(rownames(countsmat), rownames(mito.ribo.counts)),]
  se_qc_subset[["Spatial"]] = CreateAssayObject(counts = countsmat)

  # regenerate usual Seurat preprocessing for RNA assay
  DefaultAssay(se_qc_subset) = "Spatial"
  se_qc_subset = NormalizeData(se_qc_subset, verbose = F)
  se_qc_subset = FindVariableFeatures(se_qc_subset, selection.method = "vst")
  se_qc_subset = ScaleData(se_qc_subset, features = rownames(se_qc_subset))

  # additionally, introduce log-transformed counts without size-normalization as a separate assay
  se_qc_subset[["RNA_logcounts"]] = se_qc_subset@assays$Spatial
  DefaultAssay(se_qc_subset) = "RNA_logcounts"
  countsmat = as.matrix(GetAssayData(se_qc_subset, assay = "RNA_logcounts", layer = "counts"))
  countsmat = log1p(countsmat)
  se_qc_subset = SetAssayData(se_qc_subset, assay = "RNA_logcounts", layer = "data", new.data = countsmat)
  se_qc_subset = ScaleData(se_qc_subset, features = rownames(se_qc_subset), assay = "RNA_logcounts")
  Misc(se_qc_subset, "post_QC_dims_after_RiboMito_separation") <- dim(se_qc_subset)

  # Recompute SCT assay from scratch (also without ribo/mito genes)
  se_qc_subset[["SCT"]] = NULL
  se_qc_subset = SCTransform(se_qc_subset, vst.flavor = "v2",
                             seed = random_seed, assay = "Spatial",
                             variable.features.n = 4000)


  print(paste("Recomputing UMAPs and clustering:", sample_id))
  DefaultAssay(se_qc_subset) = "SCT"
  se_qc_subset = RunPCA(se_qc_subset)
  nbr.dims = manifest$pca_components_clustering_sct
  se_qc_subset = FindNeighbors(se_qc_subset, reduction = "pca", dims = 1:nbr.dims, verbose = F)
  se_qc_subset = FindClusters(object = se_qc_subset, random.seed = random_seed, resolution = manifest$clustree_resolutions, verbose=F)
  se_qc_subset = RunUMAP(se_qc_subset, reduction = "pca", dims = 1:nbr.dims, seed.use = random_seed, verbose = F)

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Save
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  saveRDS(se_qc_subset, file = file.path(work.dir, "04_seurat.rds"))
}
