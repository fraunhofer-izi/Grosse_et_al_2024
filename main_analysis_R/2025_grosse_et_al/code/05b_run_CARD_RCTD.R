# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Decon...TODO
# Take:
#   02_seurat.rds #TODO change this
# Store:
#   05_seurat.rds
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("Seurat")
library("yaml")
library("naturalsort")
library("reshape2")
library("semla")
library("hdf5r")
library("Matrix")
library("tictoc")
library("MuSiC")
library("CARD")
library("spacexr")
library("tidyverse")

# bioc
library("BiocParallel")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
all_samples = c(manifest$samples, manifest$external_samples)

# load and prep single-cell refs
refs = readRDS(file = "data/sc_ref_Stubenvoll_Schmidt_2025.rds")

random_seed = 1234
set.seed(random_seed)

for (sample_id in all_samples) {

  print(paste("Active sample is now:", sample_id))
  tic()
  sample_manifest = manifest[[sample_id]]
  work.dir = file.path(manifest$workdata_parent, sample_manifest$workdata)
  se = readRDS(file.path(work.dir, "02_seurat.rds"))

  # get appropriate reference
  ref = if (sample_id %in% manifest$naevus_samples) refs$nevi_ref else refs$mela_ref
  sc_counts = GetAssayData(ref, assay = "RNA", layer = "counts")
  sc_meta = ref@meta.data %>%
    select(sampleID, celltype, celltype_medium, celltype_broad)


  sp_counts = GetAssayData(se, assay = "Spatial", layer = "counts")
  sp_locations = data.frame(row.names = GetStaffli(se)@meta_data$barcode,
                            x = GetStaffli(se)@meta_data$pxl_row_in_fullres_transformed,
                            y = GetStaffli(se)@meta_data$pxl_col_in_fullres_transformed)

  # Run CARD, it's relatively fast

  # Original fine annotation
  CARD_obj_fine = createCARDObject(
    sc_count = sc_counts,
    sc_meta = sc_meta,
    spatial_count = sp_counts,
    spatial_location = sp_locations,
    ct.varname = "celltype",
    ct.select = unique(sc_meta$celltype),
    sample.varname = "sampleID",
    minCountGene = 50,
    minCountSpot = 5)

  CARD_obj_fine = CARD_deconvolution(CARD_object = CARD_obj_fine)

  # annotation with collapsed subtypes
  CARD_obj_medium = createCARDObject(
    sc_count = sc_counts,
    sc_meta = sc_meta,
    spatial_count = sp_counts,
    spatial_location = sp_locations,
    ct.varname = "celltype_medium",
    ct.select = unique(sc_meta$celltype_medium),
    sample.varname = "sampleID",
    minCountGene = 100,
    minCountSpot = 5)

  CARD_obj_medium = CARD_deconvolution(CARD_object = CARD_obj_medium)

  # and even broader annotation
  CARD_obj_broad = createCARDObject(
    sc_count = sc_counts,
    sc_meta = sc_meta,
    spatial_count = sp_counts,
    spatial_location = sp_locations,
    ct.varname = "celltype_broad",
    ct.select = unique(sc_meta$celltype_broad),
    sample.varname = "sampleID",
    minCountGene = 100,
    minCountSpot = 5)

  CARD_obj_broad = CARD_deconvolution(CARD_object = CARD_obj_broad)


  # RCTD

  sp = SpatialRNA(sp_locations, sp_counts)
  barcodes = colnames(sp@counts)
  # default is 25, set it down to 5. for some reason, we need to filter manually
  min_cells_per_type = 5
  if (!(all(table(sc_meta$celltype_medium) >= min_cells_per_type))) {
    # when some types in the reference have too few cells
    type_min_mask = sc_meta$celltype_medium %in% names(table(sc_meta$celltype_medium)[table(sc_meta$celltype_medium)>=5])

    sc_types_sub = sc_meta$celltype_medium[type_min_mask]
    sc_counts_sub = sc_counts[,type_min_mask]
    sc_cells_sub = rownames(sc_meta)[type_min_mask]
  } else {
    sc_types_sub = sc_meta$celltype_medium
    sc_counts_sub = sc_counts
    sc_cells_sub = rownames(sc_meta)
  }

  ref_types = setNames(factor(sc_types_sub), nm = sc_cells_sub)

  # run medium resolution reference
  reference_medium <- spacexr::Reference(counts = sc_counts_sub,
                                         cell_types = ref_types)

  myRCTD_medium <- create.RCTD(sp, reference_medium,
                               MAX_MULTI_TYPES = 4, # assume a maximum of 4 different cell types per spot
                               UMI_min = manifest$QC_cutoffs$counts, # for consistency with our QC
                               CELL_MIN_INSTANCE = min_cells_per_type,
                               max_cores = 20)
  # we may want multi or full mode, but as per the manual, multi mode output
  # should contain the full mode composition estimates,
  # so we can just extract them later if need be
  myRCTD_medium <- run.RCTD(myRCTD_medium, doublet_mode = 'multi')

  obj = list(
    se = se,
    CARD_obj_broad = CARD_obj_broad,
    CARD_obj_medium = CARD_obj_medium,
    CARD_obj_fine = CARD_obj_fine,
    RCTD_obj_medium = myRCTD_medium
  )

  toc() # time keeping
  saveRDS(obj, file = file.path(work.dir, "05_seurat_deconv_results.rds"))

}
