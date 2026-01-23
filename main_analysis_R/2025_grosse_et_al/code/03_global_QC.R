# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Take all valid spots and do QC metrics
# Take:
#   02_seurat.rds of all samples
# Store:
#   03_naively_merged_samples_se.rds and
#   03_seurat_integrated_samples_se.rds
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("tidyverse")
library("Seurat")
library("yaml")
library("naturalsort")
library("reshape2")
library("semla")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
manifest = yaml.load_file("assets/manifest.yaml")
all_samples = manifest$samples
random_seed = 1234

ses = lapply(all_samples, function(s){
  sample_manifest = manifest[[s]]
  work.dir = file.path(manifest$workdata_parent, sample_manifest$workdata)
  se = readRDS(file.path(work.dir, "02_seurat.rds"))
  DefaultAssay(se) = "Spatial"
  se = DietSeurat(se) # load all samples
  se[["Sample"]] = s
  se
})

# merge samples naively
prefixed_ses = lapply(ses, function(se){
  sample = unique(se[["Sample"]]$Sample)
  se = RenameCells(object = se, add.cell.id = sample)
  return(se)
})

# merging is broken at the time of writing,
# so do this manually
counts = do.call(what = cbind, lapply(prefixed_ses, function(se)GetAssayData(se, "Spatial", layer = "counts")))
meta = Reduce(full_join, lapply(prefixed_ses, function(se)se@meta.data %>%
                              select(!(contains("SCT") | ends_with("_Spatial") | orig.ident | seurat_clusters)) %>%
                              mutate(cellname = rownames(se@meta.data))))
merged_ses = CreateSeuratObject(counts = counts, assay = "Spatial", project = "NaivelyMergedMelanoma")
merged_ses = AddMetaData(merged_ses, meta %>% select(!cellname))
Idents(merged_ses) = "Sample"

# Assess mixing
merged_ses = NormalizeData(merged_ses)
VariableFeatures(merged_ses) = split(row.names(merged_ses@meta.data), merged_ses@meta.data$Sample) %>% lapply(function(cells_use) {
  merged_ses[,cells_use] %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 1500) %>%
    VariableFeatures()
}) %>% unlist %>% unique
merged_ses = ScaleData(merged_ses, features = VariableFeatures(merged_ses))
merged_ses = RunPCA(merged_ses, features = VariableFeatures(merged_ses))
merged_ses = RunUMAP(merged_ses, dims = 1:20)
#UMAPPlot(merged_ses)
max.k = 500
mm = MixingMetric(merged_ses, "Sample", reduction = "pca", dims = 1:10, max.k = max.k)
merged_ses$MixingMetric = max.k - mm
#VlnPlot(merged_ses, "MixingMetric") # higher is better, so this is badly mixed

# thus, we need to integrate the samples
# ses = lapply(ses, function(se){
#   DefaultAssay(se) = "SCT"
#   se
# })
# integration_features = SelectIntegrationFeatures(object.list = ses, nfeatures = 3000)
# ses = PrepSCTIntegration(object.list = ses, anchor.features = integration_features)
# integration_anchors = FindIntegrationAnchors(object.list = ses, normalization.method = "SCT", anchor.features = integration_features)
# integrated_ses = IntegrateData(anchorset = integration_anchors, normalization.method = "SCT")
# integrated_ses = RunPCA(integrated_ses)
# integrated_ses = RunUMAP(integrated_ses, dims = 1:30)
# # hopefully it's better now
# Idents(integrated_ses) = "Sample"
# #UMAPPlot(integrated_ses)
# max.k.integ = 500
# mm.integ = MixingMetric(integrated_ses, "Sample", reduction = "umap", max.k = max.k)
# integrated_ses[["MixingMetric"]] = max.k.integ - mm.integ
# integrated_ses = FindNeighbors(integrated_ses, reduction = "pca", dims = 1:20)
# resolutions = manifest$global_qc_resolutions
# integrated_ses = FindClusters(integrated_ses, random.seed = random_seed, resolution = resolutions, algorithm = 2)

# Note: shifted this away from the normal Seurat integration workflow towards Harmony in order to save time
integrated_ses = merged_ses
integrated_ses[["Spatial"]] <- split(integrated_ses[["Spatial"]], f = integrated_ses$Sample)
integrated_ses <- FindVariableFeatures(integrated_ses)
integrated_ses <- ScaleData(integrated_ses)
integrated_ses <- RunPCA(integrated_ses)
integrated_ses <- FindNeighbors(integrated_ses, dims = 1:20, reduction = "pca")

# Modifying Parameters, see HarmonyIntegration help
# We can also add arguments specific to Harmony such as theta, to give more diverse clusters
integrated_ses <- IntegrateLayers(object = integrated_ses, assay = "Spatial",
                                  method = HarmonyIntegration, orig.reduction = "pca",
                                  new.reduction = 'harmony', verbose = FALSE)

integrated_ses <- FindNeighbors(integrated_ses, reduction = "harmony", dims = 1:20)
integrated_ses <- FindClusters(integrated_ses, resolution = manifest$global_qc_resolutions)
integrated_ses <- RunUMAP(integrated_ses, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")

Idents(integrated_ses) = "Sample"
max.k.integ = 500
mm.integ = MixingMetric(integrated_ses, "Sample", reduction = "umap.harmony", max.k = max.k)
integrated_ses$MixingMetric = max.k.integ - mm.integ

# # Integrating SCTransformed data
# obj <- SCTransform(object = obj)
# obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
#                        orig.reduction = "pca", new.reduction = 'harmony',
#                        assay = "SCT", verbose = FALSE)



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
destination =  file.path(manifest$workdata_parent, manifest$merge_and_integration_dir)
if (!(dir.exists(destination))) {dir.create(destination, recursive = T)}
saveRDS(merged_ses, file = file.path(destination, "03_naively_merged_samples_se.rds"))
saveRDS(integrated_ses, file = file.path(destination, "03_seurat_integrated_samples_se.rds"))
