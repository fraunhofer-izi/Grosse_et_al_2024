# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Prepare single cell reference, load NCBI GEO dataset,
# parse annotation tsv, QC and save ref object
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("tidyverse")
library("Seurat")
library("yaml")
library("BayesPrism")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

manifest = yaml.load_file("assets/manifest.yaml")
random_seed = 1234

ref_data_dir = file.path(manifest$single_cell_reference_dir, "GSE277165_extracted")
ref_anno_file = "./assets/single_cell_reference/Stubenvoll_Schmidt_2025_cell_annotation_parsed_and_converted.tsv"
ref_collapsing_file = "./assets/single_cell_reference/cell_type_state_matching.tsv"

ref_anno = read_tsv(ref_anno_file)
celltype_collapsing = read_tsv(ref_collapsing_file)
ref_anno = ref_anno %>%
  left_join(celltype_collapsing,
            by = "celltype",
            relationship = "many-to-one")

sample_ids = ref_anno %>% select(sampleID, GEO_samplename) %>% distinct()

ref_data_files = grep(x = dir(ref_data_dir), "_ATAC_", invert = T, fixed = T, value = T)

barcode_files = grep(x = ref_data_files, "barcodes", fixed = T, value = T)
feature_files = grep(x = ref_data_files, "features", fixed = T, value = T)
matrix_files  = grep(x = ref_data_files, "matrix", fixed = T, value = T)

ses = mapply(sample_ids$GEO_samplename,
             sample_ids$sampleID,
             FUN = function(geo_sid, prefix) {
     barcodes = file.path(ref_data_dir, grep(geo_sid, barcode_files, value = T))
     features = file.path(ref_data_dir, grep(geo_sid, feature_files, value = T))
     matrix = file.path(ref_data_dir, grep(geo_sid, matrix_files, value = T))
     expr_mtx = ReadMtx(mtx = matrix, features = features, cells = barcodes)
     se = CreateSeuratObject(counts = expr_mtx)
     annos = ref_anno %>%
       filter(sampleID == prefix)
     annos = data.frame(celltype = annos$celltype,
                        celltype_medium = annos$celltype_medium,
                        celltype_broad = annos$celltype_broad,
                        celltype_heatmap_anno = annos$celltype_heatmap_anno,
                        sum_type = annos$sum_type,
                        sampleID = geo_sid,
                        row.names = annos$in_sample_barcode)
     se = AddMetaData(se, metadata = annos)
     # cells not in the list did not pass QC in the publication
     se = subset(se, cells = rownames(annos))
     se = RenameCells(se, add.cell.id = geo_sid)
     return(se)
  })

# merge and save back
# melanoma and nevi separately
melas = ses[sort(grep("MM", names(ses), value = T))]
nevi = ses[sort(grep("Nae", names(ses), value = T))]

# currently, there seems to be a problem with merging counts, so we need this workaround
mela_ref = CreateSeuratObject(
  counts = do.call(cbind, lapply(melas, function(se)GetAssayData(se, "RNA", "counts"))),
  meta.data = do.call(rbind, lapply(melas, function(se)se@meta.data[,c("celltype", "celltype_medium", "celltype_broad", "celltype_heatmap_anno", "sum_type", "sampleID")])))
nevi_ref = CreateSeuratObject(
  counts = do.call(cbind, lapply(nevi, function(se)GetAssayData(se, "RNA", "counts"))),
  meta.data = do.call(rbind, lapply(nevi, function(se)se@meta.data[,c("celltype", "celltype_medium", "celltype_broad", "celltype_heatmap_anno", "sum_type", "sampleID")])))

obj = list(
  mela_ref = mela_ref,
  nevi_ref = nevi_ref
)

saveRDS(obj, file = "data/sc_ref_Stubenvoll_Schmidt_2025.rds")

# additionally find markers for celltypes

# subset for speed-up and less imbalance
Idents(mela_ref) = "celltype_medium"
mela_ref_sub = subset(mela_ref, downsample = 2500, seed = random_seed)

mela_counts = GetAssayData(mela_ref_sub, assay = "RNA", layer = "counts")
mela_states = mela_ref_sub$celltype_medium
mela_types = mela_ref_sub$sum_type

mela_counts_clean = cleanup.genes(input=t(mela_counts),
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                  exp.cells=10)

mela_exp_stats = get.exp.stat(sc.dat = mela_counts_clean,
                              cell.type.labels = mela_types,
                              cell.state.labels = mela_states,
                              n.cores = 8)

# raw for backup, only lfc>0 as minimal filtering
mela_exp_stats_raw = lapply(mela_exp_stats, function(df)df %>%
                              filter(min.lfc>0) %>%
                              arrange(desc(min.lfc)))

lfc.threshold = 0.3
pval.threshold = 0.01

mela_exp_stats_filtered = lapply(mela_exp_stats, function(df)df %>%
                                  filter((pval.up.min<pval.threshold)&(min.lfc>lfc.threshold)) %>%
                                  arrange(desc(min.lfc)))

mela_exp_stats_filt_annot = lapply(names(mela_exp_stats_filtered), function(nm){
  df = mela_exp_stats_filtered[[nm]]
  in.group = mela_types == nm
  out.group = !in.group
  counts_use = mela_counts_clean[,rownames(df)]
  df.pct.expressed = data.frame(row.names = rownames(df),
                                pct.expressed.celltype = colMeans(counts_use[in.group,]>0),
                                pct.expressed.other = colMeans(counts_use[out.group,]>0))
  return(cbind(df, df.pct.expressed))
})

# take genes with consistently high expression & low off-target expression,
# subset to at most 30 genes with highest logFC
mela_chosen_markers = lapply(mela_exp_stats_filt_annot, function(df){
  df %>%
    filter((pct.expressed.celltype > 0.3) &
            (pct.expressed.other < 0.2) &
            (pct.expressed.celltype > 4*pct.expressed.other)) %>%
    rownames() %>%
    head(n = 30)
})
names(mela_chosen_markers) = names(mela_exp_stats_filtered)

# same for nevi samples

Idents(nevi_ref) = "celltype_medium"
nevi_ref_sub = subset(nevi_ref, downsample = 2500, seed = random_seed)

nevi_counts = GetAssayData(nevi_ref_sub, assay = "RNA", layer = "counts")
nevi_states = nevi_ref_sub$celltype_medium
nevi_types = nevi_ref_sub$sum_type

nevi_counts_clean = cleanup.genes(input=t(nevi_counts),
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY","act") ,
                                  exp.cells=3)

nevi_exp_stats = get.exp.stat(sc.dat = nevi_counts_clean,
                              cell.type.labels = nevi_types,
                              cell.state.labels = nevi_states,
                              n.cores = 8)

# raw for backup, only lfc>0 as minimal filtering
nevi_exp_stats_raw = lapply(nevi_exp_stats, function(df)df %>%
                              filter(min.lfc>0) %>%
                              arrange(desc(min.lfc)))

lfc.threshold = 0.3
pval.threshold = 0.01

nevi_exp_stats_filtered = lapply(nevi_exp_stats, function(df)df %>%
                                   filter((pval.up.min<pval.threshold)&(min.lfc>lfc.threshold)) %>%
                                   arrange(desc(min.lfc)))

nevi_exp_stats_filt_annot = lapply(names(nevi_exp_stats_filtered), function(nm){
  df = nevi_exp_stats_filtered[[nm]]
  in.group = nevi_types == nm
  out.group = !in.group
  counts_use = nevi_counts_clean[,rownames(df)]
  df.pct.expressed = data.frame(row.names = rownames(df),
                                pct.expressed.celltype = colMeans(counts_use[in.group,]>0),
                                pct.expressed.other = colMeans(counts_use[out.group,]>0))
  return(cbind(df, df.pct.expressed))
})

# take genes with consistently high expression & low off-target expression,
# subset to at most 30 genes with highest logFC
nevi_chosen_markers = lapply(nevi_exp_stats_filt_annot, function(df){
  df %>%
    filter((pct.expressed.celltype > 0.3) &
             (pct.expressed.other < 0.2) &
             (pct.expressed.celltype > 4*pct.expressed.other)) %>%
    rownames() %>%
    head(n = 30)
})
names(nevi_chosen_markers) = names(nevi_exp_stats_filtered)
# due to so few cells in the nevi ref, some types (glandular, B cells) end up with
# hardly any markers. We currently don't use those downstream but if we tried to,
# this would be a problem. So rescue those here by padding with the markers from
# the mela_ref, which should be consistent for B and glandular cells
nevi_chosen_markers = mapply(nevi_chosen_markers, mela_chosen_markers[names(nevi_chosen_markers)], FUN = function(nevi_m, mela_m){head(c(nevi_m, setdiff(mela_m, nevi_m)), n = 30)})


saveRDS(list(mela_exp_stats_raw = mela_exp_stats_raw,
             mela_exp_stats_filt_annot = mela_exp_stats_filt_annot,
             mela_chosen_markers = mela_chosen_markers,
             nevi_exp_stats_raw = nevi_exp_stats_raw,
             nevi_exp_stats_filt_annot = nevi_exp_stats_filt_annot,
             nevi_chosen_markers = nevi_chosen_markers), file = "data/sc_markers_Stubenvoll_Schmidt_2025.rds")

