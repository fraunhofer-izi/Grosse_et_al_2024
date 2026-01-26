# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Match data from the previous step with annotations, removing unannotated spots
# which can be used to remove sections that have been wrongly annotated as
# tissue by spaceranger
# Take:
#   annotation.csv
#   01_seurat.rds
# Store:
#   02_seurat.rds # independent of whether annotation worked or not
#                 # is_annotated in Misc(se) stores this
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("tidyverse")
library("Seurat")
library("yaml")
library("ggplot2")
library("cowplot")
library("clustree")
library("naturalsort")
library("reshape2")
library("patchwork")
library("semla")
library("hdf5r")

# bioc
library("BiocParallel")
library("glmGamPoi")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
all_samples = c(manifest$samples, manifest$external_samples)
random_seed = 1234
set.seed(random_seed)

for (sample_id in all_samples) {

  print(paste("Active sample is now:", sample_id))
  sample_manifest = manifest[[sample_id]]
  work.dir = file.path(manifest$workdata_parent, sample_manifest$workdata)

  annotations_dir = sample_manifest$annotations
  annotations_file = file.path(annotations_dir, "annotations.csv")

  if (!file.exists(annotations_file)) {
    print(paste0("No annotation for sample ", sample_id, " available."))
    # annotations are missing, so instead of subsetting we just cluster and save the file under a different name
    # an annotation indicates to later parts of the pipeline that there are no pictures and things for comparison to the clustering
    se = readRDS(file.path(work.dir, "01_seurat.rds"))
    DefaultAssay(se) = "SCT"
    se = FindClusters(object = se, random.seed = random_seed,
                      resolution = manifest$clustree_resolutions)
    se = RunUMAP(se, reduction = "pca", dims = 1:nbr.dims,
                 seed.use = random_seed)
    # store indication of missing annotations
    Misc(se, "is_annotated") <- F
    # store number of genes and cells after processing (should not have changed if no annotations were transferred)
    Misc(se, "annotation_transfer_dims") <- dim(se[["Spatial"]])
    # save back to file
    saveRDS(se, file = file.path(work.dir, "02_seurat.rds"))
  } else {

    # read in data
    se = readRDS(file.path(work.dir, "01_seurat.rds"))

    annotations_kunz <- read.csv(annotations_file, row.names = 1)

    # drop cells without annotation if any
    shared_cells = intersect(Cells(se), rownames(annotations_kunz))
    if (length(shared_cells)<length(Cells(se))) {
      se <- SubsetSTData(se, spots = shared_cells)
      DefaultAssay(se) = "Spatial"
      se[["SCT"]] = NULL
      # redo part of the pipeline now because we may have lost spots,
      # which affects scaling and clustering
      se = Seurat::SCTransform(se, verbose = FALSE, vst.flavor = "v2",
                               seed = random_seed, assay = "Spatial",
                               variable.features.n = 4000) # redo because cells changed
    }
    # add annotations
    se <- AddMetaData(se, annotations_kunz)

    DefaultAssay(se) = "SCT"
    nbr.dims = manifest$pca_components_clustering_sct
    se = RunPCA(se)
    se = FindNeighbors(se, reduction = "pca", dims = 1:nbr.dims)
    se = FindClusters(object = se, random.seed = random_seed,
                      resolution = manifest$clustree_resolutions)
    se = RunUMAP(se, reduction = "pca", dims = 1:nbr.dims,
                 seed.use = random_seed)

    # store indication of successful annotations
    Misc(se, "is_annotated") <- T
    # store number of genes and cells after processing (should not have changed if no annotations were transferred)
    Misc(se, "annotation_transfer_dims") <- dim(se[["Spatial"]])
    # save back to file
    saveRDS(se, file = file.path(work.dir, "02_seurat.rds"))
  }
}
