# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Read in data from the cluster (from secondary analysis) and do basic preprocessing
# Store:
#   01_seurat.rds (only spots under tissue according to spaceranger)
#   01_seurat_full.rds (all spots for comparison with background)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("hdf5r")
library("jsonlite")
library("Seurat")
library("yaml")
library("ggplot2")

# bioc
library("BiocParallel")

library("semla")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FOR EACH OF THE SAMPLES IN THE MANIFEST
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
all_samples = c(manifest$samples, manifest$external_samples)
# also, set a seed for some operations
random_seed = 1234

source("assets/cellCycleMarkers.R")
source("assets/houseKeepingMarkers.R")

for (sample_id in all_samples) {

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # make object
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  sample_manifest = manifest[[sample_id]]
  data.dir = sample_manifest$secondary_analysis
  work.dir = file.path(manifest$workdata_parent, sample_manifest$workdata)

  infoTable = data.frame(
    samples = file.path(data.dir,"raw_feature_bc_matrix.h5"),
    imgs = file.path(data.dir,"spatial/tissue_hires_image.png"),
    spotfiles = file.path(data.dir,"spatial/tissue_positions_list.csv"),
    json = file.path(data.dir,"spatial/scalefactors_json.json")
  )

  # For semla there is another S4 class object stored inside the Seurat object,
  # called Staffli. This object contains all the STutility/semla specific meta data, like
  # pixel cooridinates, sample IDs, platform types etc. (GetStaffli(se)).
  se = ReadVisiumData(infoTable)
  se = LoadImages(se)
  # exclude cells without any counts
  se = SubsetSTData(se, expression = nCount_Spatial != 0)

  # also save number of genes and cells before any of the processing steps occur
  Misc(se, "raw_input_dims") <- dim(se)

  transform_df = tibble(sampleID = 1,
                        mirror_x = as.logical(sample_manifest$mirror_x),
                        mirror_y = as.logical(sample_manifest$mirror_y),
                        angle = sample_manifest$rotate,
                        tr_x = 1/2000, # hack to avoid semla complaining if we don't need to apply any transform
                        tr_y = 0,
                        scalefactor = 1)

  se = RigidTransformImages(se, transforms = transform_df)

  # cropping is now handled during plotting

  # all spots for visualization of out of sample spots
  se.full = ReadVisiumData(infoTable, remove_spots_outside_tissue = F)
  se.full = LoadImages(se.full, image_height = 1000)
  # Samples in se.full contain all (including empty) spots and remain uncropped + has higher resolution
  se.full = RigidTransformImages(se.full, transforms = transform_df)

  # LogNormalize usual Assay
  se = NormalizeData(se, verbose = FALSE)
  se = FindVariableFeatures(se, assay = "Spatial", selection.method = "vst")
  se = ScaleData(se, assay = "Spatial", features = rownames(se))

  # Additionally create SCTransform Assay
  se = SCTransform(se, verbose = FALSE, vst.flavor = "v2",
                   seed = random_seed, assay = "Spatial",
                   variable.features.n = 4000) # Variance stabilization

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Add basic QC stuff to meta.data
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  DefaultAssay(se) = "Spatial"
  se[["Perc_of_mito_spot"]] = PercentageFeatureSet(se, pattern = "^MT-")
  se[["Perc_of_ribosomal_spot"]] = PercentageFeatureSet(se, pattern = "^RPL|^RPS")

  # Cell cycle states
  se = CellCycleScoring(se, s.features = s.genes, g2m.features = g2m.genes, search = TRUE)

  # Average expression for housekeeping genes
  hk.exprs = FetchData(se, housekeeping, layer = "data")
  se[["hk_exprs_ave"]] = apply(hk.exprs, 1, function(x){log(x = mean(x = exp(x = x) - 1) + 1)})

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # PCA
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  DefaultAssay(se) = "SCT"

  se = RunPCA(se, assay = "SCT")
  nbr.dims = manifest$pca_components_clustering_sct
  se = FindNeighbors(se, reduction = "pca", dims = 1:nbr.dims)
  # clustering and UMAP and so on happen after annotation and QC

  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Save
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if (!(dir.exists(work.dir))) dir.create(work.dir, recursive = T)
  saveRDS(se, file = file.path(work.dir, "01_seurat.rds"))
  saveRDS(se.full, file = file.path(work.dir, "01_seurat_full.rds"))

}
