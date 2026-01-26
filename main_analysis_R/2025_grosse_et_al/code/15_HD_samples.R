# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Analysis of preprocessed Visium HD data (QC, define regions, LIANA, NICHES)
# Take:
#   paths from manifest
# Store:
#   HD_preprocessed_{sample_id}.rds
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("conflicted")
library("Seurat")
library("yaml")
library("tidyverse")
library("reshape2")
library("robustbase")
library("semla")
library("anndata")
anndata::install_anndata()

library("NICHES")
library("liana")

library("distances")
library("SingleCellExperiment")

# collection of important functions
source("code/helper/FG_helper_NICHES.R")

conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::first,
  dplyr::rename,
  base::setdiff,
  matrixStats::rowSds
  )


min_dist_to <- function(distobj, distobj_names, spots_df, to_types_vector,
                        spotcol = "object_id", celltypecol = "celltype"){
  target_rows = match(spots_df %>% dplyr::pull(!!spotcol), table = distobj_names)
  matching_spots = spots_df[[celltypecol]] %in% to_types_vector
  target_cols = match(spots_df[matching_spots,] %>% pull(!!spotcol), table = distobj_names)
  ret = apply(distobj[target_rows, target_cols], FUN = "min", MARGIN = 1)
  return(ret)
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# First: prep LRI based on marker gene scores
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

manifest = yaml.load_file("assets/manifest.yaml")

# HD files from python processing pipeline
hd_samples = manifest$HD_data_ids_used
visiumHD_files_named = lapply(setNames(nm = hd_samples), function(x)manifest$HD_data_input_paths[[x]])

random_seed = 1234567

ref_collapsing_file = "assets/single_cell_reference/cell_type_state_matching.tsv"
celltype_collapsing = read_tsv(ref_collapsing_file)

# load once at the start, prepares if necessary (briefly loads seurat objects,
# looks for genes and searches updates for LRI genes not in any of out samples)
liana_nichenet_combined_db = get_combined_lr_dbs()
lr_db = liana_nichenet_combined_db %>% select(from, to) %>% distinct()

liana_format_combined_db = liana_nichenet_combined_db %>%
  transmute(source_genesymbol = from, target_genesymbol = to, db_source = source) %>%
  group_by(source_genesymbol, target_genesymbol) %>%
  summarize(db_source=paste(sort(db_source), collapse = ", "), .groups = "drop")

medium2broad = celltype_collapsing %>% select(celltype_medium, celltype_broad) %>% distinct()
medium2broad = setNames(medium2broad$celltype_broad, nm = medium2broad$celltype_medium)


# palette for all celltype annos, based on the paul tol palettes but manually altered to represent groups better
cellcolors = c("Plasma_cells" = "#114477", #
  "B_cells" = "#4477AA", #
  "pDC" = "#77AADD", #
  "cDC" =  "#228895", #
  "cDC_1" =  "#228895",
  "cDC_2" = "#55BBBB",
  "Mac" = "#33D9B1", #
  "NK" = "#B9E532", #
  "T" = "#EE6310", #
  "Tcyt_activated" = "#FF5013",
  "Teff_mem" = "#EE6310",
  "Treg" = "#F1831D",
  "Tdp" = "#FEA514",
  "VE" = "#AA4455", #
  "VE_1" = "#661122",
  "VE_2" = "#AA4455",
  "LE" = "#EE7788", #
  "Mast_cells" = "#FEEB88", #
  "Mel" = "#773311", #
  "Mel_nc-like" = "#331100", #
  "Mel_trans" = "#AA7744", #
  "Mel_trans_melan" = "#773311",
  "Pericyte" = "#BABABA",  #
  "Pericyte_1" = "#BABABA",
  "Pericyte_2" = "#808080",
  "Fb" = "#AAAA44", #
  "Fb_1" = "#777711",
  "Fb_2" = "#AAAA44",
  "Fb_3" = "#DDDD77",
  "Kc" = "#9A2E71", #
  "Kc_postmit" = "#9A2E71",
  "Kc_premit" = "#FDA7DF",
  "Eccrine_sweat_gland" = "#CC99BB", #
  "Sebaceous_gland" = "#BB1199") #



## processing steps: load, extract, filter, annotate border region, run LIANA, run NICHES
all_samples = hd_samples
random_seed = 1234

for (sample_id in all_samples) {
  print(paste("Loading sample:", sample_id))

  obj.segmented_hema.path = visiumHD_files_named[[sample_id]]
  obj.segmented_hema = read_h5ad(obj.segmented_hema.path)

  obs_for_alluvials = obj.segmented_hema$obs
  # pull transformed spatial coords (epidermis up)
  # mirror y axis (to account for ggplot having (0,0) bottom left instead of top left)
  # divide by stored mpp to get coordinates in µm with some arbitrary
  # offset -> makes cutoff distances are easy to calculate
  mpp = obj.segmented_hema$uns$spatial[[1]]$scalefactors$microns_per_pixel
  spatial_coords_transformed = data.frame(x_tr = obj.segmented_hema$obsm$spatial[,1],
                                          y_tr = - obj.segmented_hema$obsm$spatial[,2])
  spatial_coords_transformed = spatial_coords_transformed*mpp
  obs_for_alluvials = cbind(obs_for_alluvials, spatial_coords_transformed)
  gc()

  # sanity check, distances should be approximately the same if calculated from
  # array or image coordinates, irrespective of rotation or mirroring
  tmp = obs_for_alluvials[sample(seq_len(nrow(obs_for_alluvials)), size = 100), c("array_row", "array_col", "x_tr", "y_tr")]
  mean_factor_array_coords_to_image_coords = mean(as.matrix(distances::distances(tmp, dist_variables = c("array_row", "array_col"))*2) /
                                                    as.matrix(distances::distances(tmp, dist_variables = c("x_tr", "y_tr"))), na.rm=T)
  stopifnot(mean_factor_array_coords_to_image_coords>0.99)

  # hema picture segmentations, rotated such that dermis is generally up
  obs_for_alluvials$tacco_medium_translated = setNames(medium2broad[as.character(obs_for_alluvials$tacco_OT_celltype_medium)], nm = NULL)
  obs_for_alluvials$tacco_annos_consistent = obs_for_alluvials$tacco_medium_translated==as.character(obs_for_alluvials$tacco_OT_celltype_broad)
  tacco_counts_destriped = obj.segmented_hema$X
  obs_for_alluvials$total_counts = rowSums(tacco_counts_destriped)
  tacco_colnames = obj.segmented_hema$var_names # genes
  tacco_rownames = obj.segmented_hema$obs_names # cells

  # filter valid spots
  obs_for_alluvials$pass_QC = (obs_for_alluvials$total_counts >= manifest$HD_filter$min_counts) &
                              (obs_for_alluvials$bin_count >= manifest$HD_filter$min_bins) &
                              (obs_for_alluvials$bin_count < manifest$HD_filter$max_bins) &
                              (obs_for_alluvials$labels_joint_source %in% manifest$HD_filter$allowed_sources) &
                              obs_for_alluvials$tacco_annos_consistent

  print("Spots passing QC:")
  print(table(obs_for_alluvials$pass_QC))

  obs_filtered = obs_for_alluvials %>%
    filter(pass_QC)
  tacco_counts_destriped_filt = tacco_counts_destriped[obs_for_alluvials$pass_QC,]

  # calculate distances
  spots_filtered_dists = distances::distances(data.frame(spot = tacco_rownames[obs_for_alluvials$pass_QC],
                                                         x = obs_filtered$x_tr,
                                                         y = obs_filtered$y_tr),
                                              id_variable = "spot", dist_variables = c("x", "y"))
  distobj_names = names(spots_filtered_dists[,1])
  spots_filtered_border = data.frame(object_id = tacco_rownames[obs_for_alluvials$pass_QC],
                                     celltype = as.character(obs_filtered$tacco_OT_celltype_broad))
  spots_filtered_border = spots_filtered_border %>%
    mutate(dist_to_nearest_TNK = min_dist_to(spots_filtered_dists, distobj_names, spots_filtered_border,
                                             to_types_vector = c("NK", "T")),
           dist_to_nearest_Mel = min_dist_to(spots_filtered_dists, distobj_names, spots_filtered_border,
                                             to_types_vector = c("Mel"))) %>%
    mutate(close_to_TNK = dist_to_nearest_TNK <= manifest$HD_border_region_dist_cutoff_mym,
           close_to_Mel = dist_to_nearest_Mel <= manifest$HD_border_region_dist_cutoff_mym,
           bordering = ifelse(close_to_TNK, ifelse(close_to_Mel, "Tumor-Immune border", "T/NK neighborhood"),
                              ifelse(close_to_Mel, "Tumor neighborhood", "Background")),
           object_id = as.numeric(object_id))

  stopifnot(obs_filtered$object_id == spots_filtered_border$object_id)
  obs_filtered = obs_filtered %>%
    inner_join(spots_filtered_border %>% select(!c(celltype, close_to_TNK, close_to_Mel)), by = "object_id")
  gc()

  tacco_filt_se = CreateSeuratObject(counts = t(tacco_counts_destriped_filt), assay = "HD",
                                     meta.data = data.frame(row.names = tacco_rownames[obs_for_alluvials$pass_QC],
                                                            x = obs_filtered$x_tr,
                                                            y = obs_filtered$y_tr,
                                                            bincounts = obs_filtered$bin_count,
                                                            celltype = obs_filtered$tacco_OT_celltype_medium,
                                                            celltype_broad = obs_filtered$tacco_OT_celltype_broad,
                                                            bordering = obs_filtered$bordering))
  tacco_filt_se = NormalizeData(tacco_filt_se)
  tacco_filt_se = ScaleData(tacco_filt_se)

  #tacco_filt_se = subset(tacco_filt_se, bordering=="Tumor-Immune border")
  Idents(tacco_filt_se) = "celltype"
  gc()

  ## NICHES
  # more efficient edgelist computation for NICHES
  # it's O(n) or O(n log(n)) with a slight overhead instead of O(n^2) in the number of spots
  my_edgelist_2 = function(filtered.obj,position.x,position.y,k,rad.set, k.max=30){
    # kmax = maximum number of neighbors when using rad.set
    stopifnot("nn based edgelist not implemented in this function" = is.null(k))
    require(distances)
    require(dplyr)

    df <- data.frame(x = filtered.obj[[position.x]], y = filtered.obj[[position.y]])
    df$barcode <- rownames(df)
    # df$x <- as.character(df$x)
    # df$y <- as.character(df$y) # why all of this??
    df$x <- as.numeric(df$x)
    df$y <- as.numeric(df$y)
    # df <- df[,c('x','y')] #why?

    dist_obj = distances::distances(df, id_variable = "barcode", dist_variables = c("x", "y"))
    distobj_names = names(dist_obj[,1])

    edges = nearest_neighbor_search(dist_obj, k = k.max)
    edge_df = tibble(from = colnames(edges), as_tibble(t(edges)))
    edge_df = edge_df %>%
      pivot_longer(cols = !from, values_to = "to") %>%
      select(!name) %>%
      mutate(to = distobj_names[to])
    edge_df = edge_df %>%
      left_join(df %>% rename(from = barcode, x_from = x, y_from = y), by = "from") %>%
      left_join(df %>% rename(to = barcode, x_to = x, y_to = y), by = "to")
    edge_df = edge_df %>%
      transmute(from = from, to = to,
                dist = sqrt(abs(x_from - x_to)^2 + abs(y_from - y_to)^2)) %>%
      filter(dist <= rad.set)

    return(data.frame(from = edge_df$from, to = edge_df$to))
  }
  # run NICHES
  set.seed(random_seed)
  NICHES_output <- RunNICHES(object = tacco_filt_se,
                             LR.database = "custom",
                             custom_LR_database = lr_db,
                             assay = "HD",
                             position.x = 'x',
                             position.y = 'y',
                             k = NULL,
                             cell_types = "celltype",
                             rad.set = manifest$HD_border_region_dist_cutoff_mym,
                             min.cells.per.ident = 5,
                             min.cells.per.gene = 3,
                             meta.data.to.map = c('bordering','celltype_broad'),
                             CellToCell = F,CellToSystem = F,SystemToCell = F,
                             CellToCellSpatial = T,CellToNeighborhood = T,NeighborhoodToCell = F, # w/o NeighborhoodToCell to save memory
                             edgelist_fun = my_edgelist_2)

  gc()
  ## LIANA
  set.seed(random_seed)
  liana_res = liana_wrap(tacco_filt_se %>% subset(bordering == "Tumor-Immune border"),
                         resource = "custom",
                         external_resource = liana_format_combined_db,
                         expr_prop = 0.02)
  gc()

  res.obj = list(obs_for_alluvials = obs_for_alluvials,
                 NICHES_output = NICHES_output,
                 liana_res = liana_res)

  saveRDS(res.obj, file = paste0(manifest$workdata_parent,
                                 manifest$HD_intermediates,
                                 "HD_preprocessed_",sample_id,".rds"))
  print(paste("Finished sample:", sample_id))

}


## venn diagram data
# visium_venn = readr::read_csv(file.path(manifest$results_directory_parent, manifest$global_index_dir, "LRI_results_melanomas_with_annos.csv.gz")) %>%
#   select(mechanism, count_signif) %>%
#   filter(count_signif > 0) %>%
#   mutate(interaction_liana_format = gsub(pattern = "—", replacement = " -> ", mechanism)) %>%
#   transmute(interaction = interaction_liana_format, visium_signif_max = count_signif)
#
# hd_venn = do.call(rbind, lapply(seq_along(liana_aggs), function(i){
#   agg = liana_aggs[[i]][["primary_border_region_medium"]] %>%
#     filter(aggregate_rank<=0.01 #&
#            #source %in% c("Mel_nc-like", "Mel_trans", "Mel_trans_melan") &
#            #target %in% c("NK", "Treg", "Tcyt_activated", "Tdp", "Teff_mem")
#     )
#   agg$sample=names(liana_aggs)[i]
#   agg
# })) %>% group_by(source, target, ligand.complex, receptor.complex) %>%
#   mutate(n=n()) %>%
#   ungroup() %>%
#   filter(n >= 1) %>%
#   unite(source, target, sep = " --> ", col = "celltypes") %>%
#   select(celltypes, ligand.complex, receptor.complex, n) %>%
#   group_by(ligand.complex, receptor.complex) %>%
#   slice_max(n, n = 1, with_ties = T) %>%
#   summarize(n = first(n), celltypes = paste(unique(celltypes), collapse = ", ")) %>%
#   unite(ligand.complex, receptor.complex, sep = " -> ", col = "interaction") %>%
#   rename(hd_signif_max = n)
#
# hd_venn_TNK_Mel = do.call(rbind, lapply(seq_along(liana_aggs), function(i){
#   agg = liana_aggs[[i]][["primary_border_region_medium"]] %>%
#     filter(aggregate_rank<=0.01 & (
#       (source %in% c("Mel_nc-like", "Mel_trans", "Mel_trans_melan") & target %in% c("NK", "Treg", "Tcyt_activated", "Tdp", "Teff_mem")) |
#         (target %in% c("Mel_nc-like", "Mel_trans", "Mel_trans_melan") & source %in% c("NK", "Treg", "Tcyt_activated", "Tdp", "Teff_mem")))
#     )
#   agg$sample=names(liana_aggs)[i]
#   agg
# })) %>% group_by(source, target, ligand.complex, receptor.complex) %>%
#   mutate(n=n()) %>%
#   ungroup() %>%
#   filter(n >= 1) %>%
#   unite(source, target, sep = " --> ", col = "celltypes") %>%
#   select(celltypes, ligand.complex, receptor.complex, n) %>%
#   group_by(ligand.complex, receptor.complex) %>%
#   slice_max(n, n = 1, with_ties = T) %>%
#   summarize(n = first(n), celltypes_tnk_mel = paste(unique(celltypes), collapse = ", ")) %>%
#   unite(ligand.complex, receptor.complex, sep = " -> ", col = "interaction") %>%
#   rename(hd_signif_max_TNK_mel = n)
#
# venn_joint = full_join(visium_venn, hd_venn, by = "interaction") %>%
#   left_join(hd_venn_TNK_Mel, by = "interaction") %>%
#   mutate(visium_signif_max = factor(replace_na(visium_signif_max, 0), levels = 0:4),
#          hd_signif_max = factor(replace_na(hd_signif_max, 0), levels = 0:8),
#          hd_signif_max_TNK_mel = factor(replace_na(hd_signif_max_TNK_mel, 0), levels = 0:8))
#
# table(venn_joint$visium_signif_max, venn_joint$hd_signif_max)
# table(venn_joint$hd_signif_max_TNK_mel, venn_joint$hd_signif_max)
#
# venn_plt_df = rbind(
#   venn_joint %>%
#     dplyr::count(visium_signif_max, hd_signif_max, .drop = F) %>%
#     mutate(type = "hd_all"),
#   venn_joint %>%
#     dplyr::count(visium_signif_max, hd_signif_max_TNK_mel, .drop = F) %>%
#     mutate(type = "hd_TNK_mel",
#            hd_signif_max = hd_signif_max_TNK_mel,
#            hd_signif_max_TNK_mel = NULL)
# )


# wrap_plots(a = ggplot(venn_plt_df, aes(hd_signif_max, visium_signif_max, fill = n, color = n)) +
#   geom_tile(data = ~ filter(.x, type == "hd_all")) +
#   geom_point(data = ~ filter(.x, type == "hd_TNK_mel"), size = 20) +
#   geom_tile(data = tibble(hd_signif_max=as.factor(0), visium_signif_max=as.factor(0), n=0), fill = "grey", linewidth = 0) +
#   scale_color_scico(palette = "lajolla", direction = 1,
#                     aesthetics = c("color", "fill"), limits = c(0, 250), oob = scales::oob_censor, na.value = "grey20") +
#   coord_equal() + theme(legend.position = "bottom"),
#   r = ggplot(venn_plt_df %>% group_by(visium_signif_max) %>% summarize(n = sum(n)),
#              aes(y = visium_signif_max, x = n)) +
#     geom_col(),
#   design = "aaar
#             aaar",
#   axes = "collect"
#   )
# p1 = ggplot(venn_plt_df, aes(hd_signif_max, visium_signif_max, fill = n, color = n)) +
#   geom_tile(data = ~ filter(.x, type == "hd_all")) +
#   geom_tile(data = ~ filter(.x, type == "hd_TNK_mel"), width = 0.7, height = 0.7) +
#   geom_tile(data = tibble(hd_signif_max=as.factor(0), visium_signif_max=as.factor(0), n=0), fill = "grey", linewidth = 0) +
#   scale_color_scico(palette = "lajolla", direction = 1, name = "Number of significant interactions",
#                     aesthetics = c("color", "fill"), limits = c(0, 250),
#                     oob = scales::oob_censor, na.value = "grey20") +
#   theme_cowplot() +
#   theme(legend.position = "bottom", plot.margin = margin(t = -8, r = -8, unit = "pt"),
#         legend.key.width = unit(2, "lines")) +
#   labs(x = "Interactions significant in # HD samples",
#        y = "Interactions significant in # Visium samples")
# p2 = ggplot(venn_plt_df %>% group_by(visium_signif_max) %>% summarize(n = sum(n)),
#             aes(y = visium_signif_max, x = n)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_cowplot() +
#   theme(axis.text.y = element_blank(), axis.title.y = element_blank())
# p3 = ggplot(venn_plt_df %>% group_by(hd_signif_max, type) %>% summarize(n = sum(n), .groups = "drop"),
#             aes(x = hd_signif_max, y = n, group = type, fill = type)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   theme_cowplot() +
#   theme(legend.position = "inside", legend.position.inside = c(0.7,0.7),
#         axis.text.x = element_blank(), axis.title.x = element_blank())
#
#
# (p3 + plot_spacer() + p1 + p2) + plot_layout(height = c(1,3), widths = c(4,1))
#

