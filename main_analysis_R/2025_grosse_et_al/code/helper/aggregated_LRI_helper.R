##
## TCGA for annotations
##

fetch_TCGA <- function(){
  library("recount")
  if (!(dir.exists("./data/TCGA_temp/"))) {
    dir.create("./data/TCGA_temp/", recursive = T)
  }
  cachepath = normalizePath("./data/TCGA_temp/")

  outpath = file.path(cachepath, "TCGAse.Rds")

  if (file.exists(outpath)) {
    rse_gene = readRDS(outpath)
  } else {
    if (file.exists(file.path(cachepath, "rse_gene.Rdata"))) {
      load(file.path(cachepath, "rse_gene.Rdata"))
    } else {
      options(timeout = max(300, getOption("timeout")))
      download_study('TCGA', type = 'rse-gene', outdir = cachepath)
      load(file.path(cachepath, "rse_gene.Rdata"))
    }

    # fixes missing slots issue when loading data with newer package versions
    rse_gene = updateObject(rse_gene)

    # Remove PAR locus genes (45 genes). This avoids EnemblID duplications after Suffix truncation
    rse_gene = rse_gene[!grepl("PAR_Y$", rownames(rse_gene)), ]

    assays(rse_gene)$tpm = recount::getTPM(rse_gene)

    saveRDS(rse_gene, file=outpath)
  }
  return(rse_gene)
}

make_TCGA.se <- function(cohort = "TCGA-SKCM", # SKCM == skin cutaneous melanoma
                         sites = c("Primary Tumor")){  # could also include Metastasis

  # make filename for caches
  cachepath = file.path("./data/TCGA_temp/", paste0("TCGA_cache_", hash(c(cohort, sites))))
  if (file.exists(cachepath)) {
    TCGA.se = readRDS(cachepath)
  } else {
    # Load data
    # confirmed that this is equivalent to michaels preprocessed TCGA download
    TCGA.se = fetch_TCGA()
    rownames(TCGA.se) = gsub("\\..+", "", rownames(TCGA.se))

    # filter to cohort
    TCGA.se = TCGA.se[, TCGA.se$gdc_cases.project.project_id == cohort]

    # filter to primary tumor (or metastases, as wanted)
    TCGA.se = TCGA.se[, TCGA.se$gdc_cases.samples.sample_type %in% sites]
    colData(TCGA.se) = droplevels(colData(TCGA.se))

    # cache
    if (!dir.exists(dirname(cachepath))) {dir.create(dirname(cachepath), recursive = T)}
    saveRDS(TCGA.se, file = cachepath)
  }
  return(TCGA.se)
}

geomean = function(x){exp(sum(log(x))/length(x))}

##
## function for filtering and joining results dataframes
##

make_consensus_df <- function(LR_results, # needs to be a named list of results
                              results.name = "border.markers.NeighborhoodToCell",
                              results.type = "logcounts" # logcounts or normalized
){
  sample.names = names(LR_results)
  LR_filtered = lapply(sample.names, function(sn){
    res = LR_results[[sn]]
    de.results = res[[results.type]][[results.name]]
    # keep only the interaction, mean_diff, log2FC, zscore,
    # pval and pval_adj with the sample name appended
    de.results = de.results %>%
      dplyr::select(mechanism = gene, mean_diff, log2FC, zscore,
                    pval_raw = pval, pval_adj) %>%
      dplyr::rename_with(function(colname){paste0(colname,"_",sn)},
                  all_of(c("mean_diff", "log2FC", "zscore",
                           "pval_raw", "pval_adj")))
    # undo the renaming of complexes from the NICHES patch
    de.results = de.results %>%
      mutate(mechanism = gsub(mechanism, pattern = "+",
                              replacement = "_", fixed = T))
    return(de.results)
  })
  # combine dataframes
  if (length(LR_filtered)<=1) {
    aggregated.results = LR_filtered[[1]]
  } else {
    x = LR_filtered[[1]]
    for (i in 2:length(LR_filtered)) {
      y = LR_filtered[[i]]
      x = full_join(x, y, by = "mechanism")
    }
    aggregated.results = x
    rm(x)
  }
  return(aggregated.results)
}

##
## function for counting in how many samples pvals are significant
##
annotate_df_with_significance_info = function(aggregated.df,
                                       pval.adj.threshold = 0.05,
                                       only.pos = T){
  # check significance of result in each sample
  if (only.pos) {
    signif.df = aggregated.df %>%
      dplyr::select(mechanism, starts_with("pval_adj"), starts_with("log2FC")) %>%
      mutate(across(starts_with("pval_adj"),
                    ~ (.x <= pval.adj.threshold)&(pick(starts_with("log2FC"))[[str_replace(cur_column(), pattern = "pval_adj", "log2FC")]] > 0) )) %>%
      dplyr::select(mechanism, starts_with("pval_adj")) %>%
      dplyr::rename_with(function(colname){
        gsub(colname, pattern = "pval_adj", replacement = "significant_in")})
  } else {
    signif.df = aggregated.df %>%
      dplyr::select(mechanism, starts_with("pval_adj")) %>%
      mutate(across(starts_with("pval_adj"),
                    function(pval){pval<=pval.adj.threshold})) %>%
      dplyr::rename_with(function(colname){
        gsub(colname, pattern = "pval_adj", replacement = "significant_in")})
  }
  aggregated.results = left_join(aggregated.df, signif.df, by = "mechanism")
  # count how often each mechanism is significant across samples
  aggregated.results = aggregated.results %>%
    mutate(count_signif =
             rowSums(across(starts_with("significant_in")), na.rm = T),
           summed_zscores =
             rowSums(across(starts_with("zscore")), na.rm = F)) %>%
    arrange(desc(count_signif))

  return(aggregated.results)
}

##
## function for adding annotations to dataframe
##

add.TCGA.annos <- function(aggregated.results,
                           TCGA_reference){
  library(biomaRt)

  mart = useMart("ENSEMBL_MART_ENSEMBL", host = "https://jul2022.archive.ensembl.org")
  mart = useDataset("hsapiens_gene_ensembl", mart=mart)
  attr = c('ensembl_gene_id', 'hgnc_symbol', 'wikigene_description', 'gene_biotype')

  # assemble list of all ligands and receptors to query
  res.df = aggregated.results %>% dplyr::select(mechanism)
  res.df = tidyr::separate(res.df, mechanism, sep = "—", into = c("source_genesymbol", "target_genesymbol"), remove = F)

  # extract relevant symbols
  # split complexes and rejoin later
  relevant.genesymbols = unique(c(res.df$source_genesymbol, res.df$target_genesymbol))
  complexes = grep(x = relevant.genesymbols, pattern = "_", fixed = T, value = T)
  relevant.genesymbols = setdiff(relevant.genesymbols, complexes)

  complex_components = str_split(complexes, pattern = "_", simplify = T)
  complex_components = setdiff(unique(c(complex_components)), "")

  relevant.genesymbols = unique(c(relevant.genesymbols, complex_components))
  # there's also 3 broken symbols: YARS and GPR1 (old) and SLURP2 (not in TCGA but is also an alias name for LYNX1, which is in TCGA)

  bm = getBM(attributes = attr, mart = mart, filters = "hgnc_symbol", values = relevant.genesymbols)
  bm = bm[bm$gene_biotype == "protein_coding", 1:3]

  tcga_tmpm_stats = data.frame(genesymbol = relevant.genesymbols)

  # match ENSEMBL IDs to gene symbols
  known.ensids = rownames(assays(TCGA_reference)$tpm)
  tcga_tmpm_stats$ensembl_gene_id = unlist(lapply(relevant.genesymbols, function(gs){
    ensids = unique(bm[bm$hgnc_symbol == gs, "ensembl_gene_id"])
    ensids = ensids[ensids %in% known.ensids]
    if (length(ensids)>1) {stop("Multiple ENSEMBLID matches for genesymbol: ", gs, " Check algorithm for correctnes.")}
    if (length(ensids)==0) {
      return("")
    } else {
      return(ensids)
    }
  }), use.names = F)

  # TPM mean and SD (display in log-scale later)
  tpmassay = assays(TCGA_reference)$tpm
  divisor = ncol(tpmassay)
  tcga_tmpm_stats$tpm_mean = sapply(tcga_tmpm_stats$ensembl_gene_id, function(gid){
    if (length(gid) == 0) {
      return(NA)
    } else {
      sum(tpmassay[known.ensids == gid,])/divisor
    }
  })
  tcga_tmpm_stats$tpm_sd = sapply(tcga_tmpm_stats$ensembl_gene_id, function(gid){
    if (length(gid) == 0) {
      return(NA)
    } else {
      sd(tpmassay[known.ensids == gid,])
    }
  })

  # robust z-score of log-expressions based on median and mad, because of non-normal distribution
  zscores = rowMeans(tpmassay)
  zscores = log(zscores)
  zscores[zscores == -Inf] = NA
  zscores = (zscores - median(zscores, na.rm = T)) / mad(zscores, na.rm = T)
  tcga_tmpm_stats$robust_z_score = sapply(tcga_tmpm_stats$ensembl_gene_id, function(gid){
    if (length(gid) == 0) {
      return(NA)
    } else {
      return(zscores[gid])
    }
  })

  # empirical quantiles in the distribution of mean TPM
  rms = rowMeans(tpmassay)
  distfun = ecdf(rms)
  tcga_tmpm_stats$quantiles = sapply(tcga_tmpm_stats$ensembl_gene_id, function(gid){
    if (length(gid) == 0) {
      return(NA)
    } else {
      return(distfun(rms[gid]))
    }
  })

  # deal with complexes
  rownames(tcga_tmpm_stats) = tcga_tmpm_stats$genesymbol
  complex_stats = data.frame(genesymbol = complexes)

  # ENSEMBL IDs from the stats df
  complex_stats$ensembl_gene_id = sapply(complex_stats$genesymbol, function(gs){
    parts = unlist(strsplit(gs, split = "_", fixed = T))
    gids = tcga_tmpm_stats[parts,"ensembl_gene_id"]
    return(paste(gids, collapse = "/"))
  })

  # geometric means of mean and sd of TPM
  complex_stats$tpm_mean = sapply(complex_stats$genesymbol, function(gs){
    parts = unlist(strsplit(gs, split = "_", fixed = T))
    mean_tpms = tcga_tmpm_stats[parts,"tpm_mean"]
    return(geomean(mean_tpms))
  })
  complex_stats$tpm_sd = sapply(complex_stats$genesymbol, function(gs){
    parts = unlist(strsplit(gs, split = "_", fixed = T))
    sd_tpms = tcga_tmpm_stats[parts,"tpm_sd"]
    return(geomean(sd_tpms))
  })

  # arithmetic mean of zscores
  complex_stats$robust_z_score = sapply(complex_stats$genesymbol, function(gs){
    parts = unlist(strsplit(gs, split = "_", fixed = T))
    zscores = tcga_tmpm_stats[parts,"robust_z_score"]
    return(mean(zscores))
  })

  # quantiles of geomean of expression in the original expression distribution
  complex_stats$quantiles = sapply(complex_stats$tpm_mean, function(gm){
    return(distfun(gm))
  })

  # combine complex and normal stats
  rownames(tcga_tmpm_stats) = NULL
  tcga_tmpm_stats = rbind(tcga_tmpm_stats, complex_stats)

  # join back
  res.df = left_join(x = res.df, by = "source_genesymbol", y = (tcga_tmpm_stats %>% dplyr::rename(source_genesymbol = genesymbol,
                                                                                           source_ensembl_gene_id = ensembl_gene_id,
                                                                                           source_TCGA_mean_tpm = tpm_mean,
                                                                                           source_TCGA_sd_tpm = tpm_sd,
                                                                                           source_TCGA_robust_z_score = robust_z_score,
                                                                                           source_TCGA_quantile = quantiles)))
  res.df = left_join(x = res.df, by = "target_genesymbol", y = (tcga_tmpm_stats %>% dplyr::rename(target_genesymbol = genesymbol,
                                                                                           target_ensembl_gene_id = ensembl_gene_id,
                                                                                           target_TCGA_mean_tpm = tpm_mean,
                                                                                           target_TCGA_sd_tpm = tpm_sd,
                                                                                           target_TCGA_robust_z_score = robust_z_score,
                                                                                           target_TCGA_quantile = quantiles)))


  if("source_genesymbol" %in% colnames(aggregated.results) & "target_genesymbol" %in% colnames(aggregated.results)){
    # remove "source_genesymbol" and "target_genesymbol" columns if already present
    res.df = mutate(res.df, source_genesymbol = NULL, target_genesymbol = NULL)
  }

  aggregated.results = left_join(aggregated.results, res.df, by = "mechanism")
  return(aggregated.results)
}

add.sc.annos <- function(aggregated.results,
                         expression.annos,
                         annos.name,
                         cell.type.filter){
  res.df = aggregated.results %>% dplyr::select(mechanism, source_genesymbol, target_genesymbol)

  # extract relevant symbols
  # split complexes and rejoin later
  relevant.genesymbols = unique(c(res.df$source_genesymbol, res.df$target_genesymbol))
  # no need to handle complexes separately, as they were already processed in the annotations file

  SC_stats = data.frame(genesymbol = relevant.genesymbols)

  expr.fraction = expression.annos$fraction.expressed[relevant.genesymbols,] %>%
    as.data.frame %>%
    dplyr::select(any_of(cell.type.filter)) %>%
    dplyr::rename_with(.fn = function(x)paste0(annos.name, "_expression_fraction_", x))
  expr.fraction$genesymbol = rownames(expr.fraction)

  log.avg.expr = expression.annos$log.avg.normalized.expression[relevant.genesymbols,] %>%
    as.data.frame %>%
    dplyr::select(any_of(cell.type.filter)) %>%
    dplyr::rename_with(.fn = function(x)paste0(annos.name, "_log_avg_expr_", x))
  log.avg.expr$genesymbol = rownames(log.avg.expr)

  norm.empirical.quantiles = expression.annos$across.normalized.celltypes.empirical.quantiles[relevant.genesymbols,] %>%
    as.data.frame %>%
    dplyr::select(any_of(cell.type.filter)) %>%
    dplyr::rename_with(.fn = function(x)paste0(annos.name, "_norm_empirical_quantiles_", x))
  norm.empirical.quantiles$genesymbol = rownames(norm.empirical.quantiles)

  unnorm.empirical.quantiles = expression.annos$per.celltype.empirical.quantile[relevant.genesymbols,] %>%
    as.data.frame %>%
    dplyr::select(any_of(cell.type.filter)) %>%
    dplyr::rename_with(.fn = function(x)paste0(annos.name, "_per_celltype_empirical_quantiles_", x))
  unnorm.empirical.quantiles$genesymbol = rownames(unnorm.empirical.quantiles)

  SC_stats = left_join(SC_stats, expr.fraction, by = "genesymbol")
  SC_stats = left_join(SC_stats, log.avg.expr, by = "genesymbol")
  SC_stats = left_join(SC_stats, norm.empirical.quantiles, by = "genesymbol")
  SC_stats = left_join(SC_stats, unnorm.empirical.quantiles, by = "genesymbol")

  # join back
  res.df = left_join(x = res.df, by = "source_genesymbol", y = (SC_stats %>% dplyr::rename_with(.fn = function(x)paste0("source_", x))))
  res.df = left_join(x = res.df, by = "target_genesymbol", y = (SC_stats %>% dplyr::rename_with(.fn = function(x)paste0("target_", x))))


  aggregated.results = left_join(aggregated.results, res.df, by = c("mechanism", "source_genesymbol", "target_genesymbol"))
  return(aggregated.results)
}


#####################################
#####################################
#####################################
zscore_limits = function(lims)c(-max(abs(lims), na.rm = T), max(abs(lims), na.rm = T))

clipped_mid_rescaler = function(x, min_clip = -4, max_clip = 4, ...){scales::rescale_mid(x = pmin(pmax(min_clip, x), max_clip), ...)}

plot_combined_df <- function(df,
                             exclude.samples = NULL,            # for subsetting samples
                             annos = sample_annos,              # to draw colored rectangles for sample grouping
                             anno.colmap = sample_anno_colmap,  # colormap for annotations
                             anno.height = 1.2,
                             order.by.hclust = T,               # whether to reorder rows by logFC clustering
                             order.signif.overwrite = F,        # if True, the number of samples a mechnaisms is significant in overwrites hclust order
                             hclust.on.samples = "all",         # enter vector of sample names to use subset for clustering
                             hclust.method = "ward.D",          # see hclust function
                             top_n = NULL,                      # plot all (NULL) rows or only the top n (before clustering)
                             plot_TCGA = c("none", "meansd", "quant", "both"),
                             plot_GTEx = c("none", "frac", "quant", "expr"),
                             plot_JA = c("none", "frac", "quant", "expr"),
                             plot_StubS = c("none", "frac", "quant", "expr"),
                             quant_anno_palette = scale_fill_gradientn(values = c(0., 0.5, 0.6, 0.7,0.8, 0.85, 0.9, 0.95, 1.),
                                                                       colors = c("grey", "#FFF7FB", "#D0D1E6", "#A6BDDB", "#64A3CB","#3075B3", "#02818A", "#016C59", "#014636"),
                                                                       rescaler = clipped_mid_rescaler),
                             frac_anno_palette = scale_fill_gradientn(values =  (0:8)/8,
                                                                      colors = scales::brewer_pal(palette = "PuBuGn")(9),
                                                                      rescaler = scales::rescale_none),
                             plot_additional_dfs_list = NULL,   # plot values in other dataframes (given as list) in separate panels, subset to main mechanisms
                             plot_signif = T,
                             legends.rel.heights = 0.1, legend_rows = 1,
                             use_geom_label = T,
                             anno.signifi.legend = c("plot", "none", "legend", "both"),   # none, plot, legend or both
                             signif_text_size = 3, width_lfc_extra = 4,
                             width_per_sample = 1, width_signif = 0.6,
                             width_TCGA = 5, with_GTEx = 8,
                             width_JA = 8, width_StubS = 10,
                             legend_barwidth = 8,
                             stubs_title = "z-transformed expression quantiles"
){
  plot_TCGA = match.arg(arg = plot_TCGA)
  plot_GTEx = match.arg(arg = plot_GTEx)
  plot_JA = match.arg(arg = plot_JA)
  plot_StubS = match.arg(arg = plot_StubS)
  anno.signifi.legend = match.arg(arg = anno.signifi.legend)
  if (!is.null(top_n)) {df = head(df, n = top_n)}
  if (nrow(df)==1) {order.by.hclust=F}
  if (!is.null(exclude.samples)) {df = dplyr::select(df, !contains(exclude.samples))}
  if (order.by.hclust){
    # draw and prepare log2FCs, replace NA with 0 and Inf with 1.1*max and -Inf with 1.1*min
    lfcs = df %>% dplyr::select(contains("log2FC"))
    if (hclust.on.samples != "all") {lfcs = lfcs %>% dplyr::select(contains(hclust.on.samples))}
    lfcs = as.matrix(lfcs)
    lfcs[is.na(lfcs)] = 0
    lfcs[lfcs == Inf] = NA
    lfcs[is.na(lfcs)] = 1.1 * max(lfcs, na.rm = T)
    lfcs[lfcs == -Inf] = NA
    lfcs[is.na(lfcs)] = 1.1 * min(lfcs, na.rm = T)
    rownames(lfcs) = df$mechanism
    dm = dist(lfcs)
    hclust.res = hclust(dm, method = hclust.method)

    df = df[hclust.res$order, ]
  }

  # significance
  # cant't really assume that the significant_in_stuff was already aggregated (may have excluded samples), so do that again
  signif.data = dplyr::select(df, mechanism, contains("significant_in_"))
  signif.data$significant_in_number = rowSums(dplyr::select(signif.data, !mechanism), na.rm = T)
  signif.data = melt(signif.data %>% dplyr::select(mechanism, significant_in_number), id.vars = "mechanism")
  signif.data$value = ordered(as.integer(signif.data$value))
  if (order.signif.overwrite) {
    signif.order = order(signif.data$value, decreasing = T)
    df = df[signif.order, ]
  }

  df$mechanism = factor(df$mechanism, levels = rev(unique(df$mechanism)))

  # plots
  ## log2FC
  log.fc.data = dplyr::select(df, mechanism, contains("log2FC"))
  colnames(log.fc.data) = gsub(colnames(log.fc.data), pattern = "log2FC_", replacement = "", fixed = T)
  log.fc.data = melt(log.fc.data, id.vars = "mechanism")
  log.fc.data$striptext = "Observed enrichment at tumor-immune border (log2FC)"

  # preparation: if there are additional dfs to plot beside, sync the log2FC color scale by calculating the maximum range
  realrange = range(log.fc.data$value[!(log.fc.data$value %in% c(-Inf, Inf))], na.rm = T)
  lfc_min_val = min(realrange)
  lfc_max_val = max(realrange)
  # update if plot_additional_dfs_list is not NULL
  for (rnge in lapply(plot_additional_dfs_list, function(add_df){
                      add_df = add_df %>% dplyr::filter(mechanism %in% levels(df$mechanism))
                      add_df$mechanism = factor(add_df$mechanism, levels = levels(df$mechanism))
                      add.log.fc.data = dplyr::select(add_df, mechanism, contains("log2FC"))
                      colnames(add.log.fc.data) = gsub(colnames(add.log.fc.data), pattern = "log2FC_", replacement = "", fixed = T)
                      add.log.fc.data = melt(add.log.fc.data, id.vars = "mechanism")
                      return(range(add.log.fc.data$value[!(add.log.fc.data$value %in% c(-Inf, Inf))], na.rm = T))})) {
    lfc_min_val = min(lfc_min_val, rnge)
    lfc_max_val = max(lfc_max_val, rnge)
  }

  p.lfc = ggplot(data = log.fc.data, aes(variable, mechanism)) +
    geom_tile(aes(fill = value), col = "white", linewidth = 0.5) + #coord_equal() +
    scale_fill_scico(palette = "broc", direction = -1, midpoint = 0., limits = c(lfc_min_val, lfc_max_val),
                     na.value = "white", oob = scales::squish_infinite,
                     guide = guide_colorbar(title.position = "top", title = "log2FC",
                                            barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"),
                                            direction = "horizontal", order = 1)) +
    ylab(NULL) + xlab(NULL) +
    facet_grid(cols = vars(striptext)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)),
          legend.position = "bottom", plot.margin = margin(r=2, l=2, unit = "pt"),
          strip.clip = "off")

  # ## pvals
  # pval.data = dplyr::select(df, mechanism, contains("pval_adj_"))
  # colnames(pval.data) = gsub(colnames(pval.data), pattern = "pval_adj_", replacement = "", fixed = T)
  # pval.data = melt(pval.data, id.vars = "mechanism")
  # pval.data$value = ifelse(pval.data$value > 0.05, NA, pval.data$value)
  #
  #
  # p.pvals = ggplot(data = pval.data, aes(variable, mechanism)) +
  #   geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
  #   scale_fill_scico(palette = "oslo", direction = 1, end = 0.9, begin = 0.2, na.value = "indianred1", trans = "log10", limits = c(0.00001, NA), oob = scales::squish) +
  #   ylab(NULL) + xlab(NULL) +
  #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom")
  #

  samples.in.use = unique(log.fc.data$variable)
  mapped.groups = annos[as.character(samples.in.use)]
  mapped.colors = anno.colmap[mapped.groups]
  groups.in.use = unique(mapped.groups)
  colors.in.use = unique(mapped.colors)
  anno.df <- data.frame (xmin = seq_along(samples.in.use) - 0.5,
                         xmax = seq_along(samples.in.use) + 0.5,
                         ymin = rep_along(samples.in.use, 0.45-anno.height),
                         ymax = rep_along(samples.in.use, 0.45),
                         fill = mapped.colors)
  p.lfc = p.lfc +
    new_scale_fill() +
    geom_rect(data = anno.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), inherit.aes = F) +
    scale_fill_identity(guide = guide_legend(title.position = "top", title = "Sample Type",
                                             order = 2, nrow = legend_rows),
                        breaks = colors.in.use, labels = groups.in.use) #+
    # guides(fill_new = guide_colorbar(title.position = "top", title = "Log2FC", available_aes = "fill_new",
                                     # barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), order = 1, direction = "horizontal"))


  ## Significance in Samples
  if (plot_signif) {
    signif.data = dplyr::select(df, mechanism, contains("significant_in_"))
    signif.data$significant_in_number = rowSums(dplyr::select(signif.data, !mechanism), na.rm = T)
    signif.data = melt(signif.data %>% dplyr::select(mechanism, significant_in_number), id.vars = "mechanism")
    signif.data$value = ordered(as.integer(signif.data$value))
    # clear names
    signif.data$variable = "Significant in number\nof samples"

    p.sig = ggplot(data = signif.data, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5)
    if (anno.signifi.legend %in% c("plot", "both")) {
      if (use_geom_label){
        p.sig = p.sig + geom_label(aes(label = value), size = signif_text_size, label.padding = unit(0.15, "lines"), label.size = 0.15)
      } else {
        p.sig = p.sig + geom_text(aes(label = value), size = signif_text_size)
      }
    }
    signif_guide = if (anno.signifi.legend %in% c("legend", "both")) {
      guide_legend(title.position = "top", title = "Number of\nSignificant Samples", barheight = unit(1, "lines"), barwidth = unit(1, "lines"))
      } else {"none"}
    p.sig = p.sig +
      scale_fill_scico_d(begin = 0.05, end = 0.65, palette = "lajolla", limits = function(lims)paste(1:max(as.numeric(lims)))) +
      guides(fill = signif_guide) +
      # guides(fill = guide_colorsteps(title.position = "top", title = "Number of\nSignificant Samples", show.limits = T, even.steps = F,
      #                                barheight = unit(.4, "lines"), barwidth = unit(6, "lines"), direction = "horizontal")) +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"))
    p.sig = p.sig + annotate("rect", xmin = 1, xmax = 1, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)
  }

  ## potentially plot additional dataframes
  if (!is.null(plot_additional_dfs_list)){
    ps.additional = lapply(plot_additional_dfs_list, function(add_df){
      # plots
      add_df = add_df %>% dplyr::filter(mechanism %in% levels(df$mechanism))
      add_df$mechanism = factor(add_df$mechanism, levels = levels(df$mechanism))
      ## log2FC
      add.log.fc.data = dplyr::select(add_df, mechanism, contains("log2FC"))
      colnames(add.log.fc.data) = gsub(colnames(add.log.fc.data), pattern = "log2FC_", replacement = "", fixed = T)
      add.log.fc.data = melt(add.log.fc.data, id.vars = "mechanism")

      p.add_lfc = ggplot(data = add.log.fc.data, aes(variable, mechanism)) +
        geom_tile(aes(fill = value), col = "white", size = 0.5) +
        scale_fill_scico(palette = "broc", direction = -1, midpoint = 0., limits = c(lfc_min_val, lfc_max_val),
                         na.value = "white", oob = scales::squish_infinite,
                         guide = "none") +
        ylab(NULL) + xlab(NULL) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
              axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"))

      samples.in.use = unique(add.log.fc.data$variable)
      mapped.groups = annos[as.character(samples.in.use)]
      mapped.colors = anno.colmap[mapped.groups]
      groups.in.use = unique(mapped.groups)
      colors.in.use = unique(mapped.colors)
      anno.df <- data.frame (xmin = seq_along(samples.in.use) - 0.5,
                             xmax = seq_along(samples.in.use) + 0.5,
                             ymin = rep_along(samples.in.use, 0.45-anno.height),
                             ymax = rep_along(samples.in.use, 0.45),
                             fill = mapped.colors)
      p.add_lfc = p.add_lfc +
        new_scale_fill() +
        geom_rect(data = anno.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), inherit.aes = F) +
        scale_fill_identity(guide = guide_legend(title.position = "top", title = "", nrow = legend_rows), breaks = colors.in.use, labels = groups.in.use)
      return(p.add_lfc)
    })
  }

  ## TCGA
  ##
  if (plot_TCGA != "none") {
    # mean/sd
    TCGA.data = dplyr::select(df, mechanism, contains("TCGA_mean")|contains("TCGA_sd"))
    TCGA.data = lapply(c("source", "target"), function(what){
      return(melt(TCGA.data %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(TCGA.data) = c("Ligand", "Receptor")
    TCGA.data = map_df(.x = (TCGA.data), .f = bind_rows, .id = "src")
    TCGA.data$variable = as.character(TCGA.data$variable)
    TCGA.data$variable = gsub(TCGA.data$variable, pattern = "source_TCGA_", replacement = "", fixed = T)
    TCGA.data$variable = gsub(TCGA.data$variable, pattern = "target_TCGA_", replacement = "", fixed = T)
    TCGA.data$variable = gsub(TCGA.data$variable, pattern = "mean_tpm", replacement = "TPM Mean", fixed = T)
    TCGA.data$variable = gsub(TCGA.data$variable, pattern = "sd_tpm", replacement = "TPM SD", fixed = T)
    TCGA.data$variable = as.factor(TCGA.data$variable)
    p.TCGA.meansd = ggplot(data = TCGA.data, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.clip = "off") +
      facet_grid(cols = vars(src)) +
      scale_fill_scico(palette = "hawaii", direction = -1, begin = 0.05, end = 0.9, na.value = "white", trans = "log10") +
      guides(fill = guide_colorbar(title.position = "top", title = "TCGA Expression", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # # zscores
    # TCGA.data = dplyr::select(df, mechanism, contains("TCGA_robust_z_score"))
    # TCGA.data = lapply(c("source", "target"), function(what){
    #   return(melt(TCGA.data %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    # })
    # names(TCGA.data) = c("Ligand", "Receptor")
    # TCGA.data = map_df(.x = (TCGA.data), .f = bind_rows, .id = "src")
    # TCGA.data$variable = as.character(TCGA.data$variable)
    # TCGA.data$variable = gsub(TCGA.data$variable, pattern = "source_TCGA_robust_z_score", replacement = "Robust z-Score", fixed = T)
    # TCGA.data$variable = gsub(TCGA.data$variable, pattern = "target_TCGA_robust_z_score", replacement = "Robust z-Score", fixed = T)
    # TCGA.data$variable = as.factor(TCGA.data$variable)
    # p.TCGA.zscores = ggplot(data = TCGA.data, aes(variable, mechanism)) +
    #   geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
    #   ylab(NULL) + xlab(NULL) +
    #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
    #         axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    #   facet_grid(cols = vars(src)) +
    #   scale_fill_scico(palette = "broc", direction = -1, midpoint = 0, begin = 0.1, end = 0.9, na.value = "white", limits = zscore_limits) +
    #   guides(fill = guide_colorbar(title.position = "top", title = "TCGA Expression", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
    #   annotate("rect", xmin = 1, xmax = 1, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # quantiles
    TCGA.data = dplyr::select(df, mechanism, contains("TCGA_quantile"))
    TCGA.data = lapply(c("source", "target"), function(what){
      return(melt(TCGA.data %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(TCGA.data) = c("Ligand", "Receptor")
    TCGA.data = map_df(.x = (TCGA.data), .f = bind_rows, .id = "src")
    TCGA.data$variable = as.character(TCGA.data$variable)
    TCGA.data$variable = gsub(TCGA.data$variable, pattern = "source_TCGA_quantile", replacement = "TCGA SKCM", fixed = T)
    TCGA.data$variable = gsub(TCGA.data$variable, pattern = "target_TCGA_quantile", replacement = "TCGA SKCM", fixed = T)
    TCGA.data$variable = as.factor(TCGA.data$variable)
    TCGA.data$value = qnorm(TCGA.data$value)
    p.TCGA.quantiles = ggplot(data = TCGA.data, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src)) +
      quant_anno_palette +
      guides(fill = guide_colorbar(title.position = "top", title = "TCGA Empirical Quantiles", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 1, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)
  }

  ### Expression in GTEx
  if (plot_GTEx != "none") {
    GTEx.data = dplyr::select(df, mechanism, contains("GTEx"))

    # avg expr
    GTEx.plt.df = dplyr::select(GTEx.data, mechanism, contains("log_avg_expr"))
    GTEx.plt.df = lapply(c("source", "target"), function(what){
      return(melt(GTEx.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(GTEx.plt.df) = c("Ligand", "Receptor")
    GTEx.plt.df = map_df(.x = (GTEx.plt.df), .f = bind_rows, .id = "src")
    GTEx.plt.df$variable = as.character(GTEx.plt.df$variable)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "source_GTEx_", replacement = "", fixed = T)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "target_GTEx_", replacement = "", fixed = T)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "log_avg_expr_", replacement = "", fixed = T)
    GTEx.plt.df$variable = factor(x = GTEx.plt.df$variable, levels = rev(unique(GTEx.plt.df$variable)))

    GTEx.plt.df$topnote = "Log Average Normalized Expr."
    GTEx.plt.df$value = exp(GTEx.plt.df$value)

    p.GTEx.avgexp = ggplot(data = GTEx.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      scale_fill_scico(palette = "hawaii", direction = -1, begin = 0.05, end = 0.9, na.value = "white", trans = "log10") +
      guides(fill = guide_colorbar(title.position = "top", title = "GTEx Expression", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # fraction expressed
    GTEx.plt.df = dplyr::select(GTEx.data, mechanism, contains("expression_fraction"))
    GTEx.plt.df = lapply(c("source", "target"), function(what){
      return(melt(GTEx.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(GTEx.plt.df) = c("Ligand", "Receptor")
    GTEx.plt.df = map_df(.x = (GTEx.plt.df), .f = bind_rows, .id = "src")
    GTEx.plt.df$variable = as.character(GTEx.plt.df$variable)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "source_GTEx_", replacement = "", fixed = T)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "target_GTEx_", replacement = "", fixed = T)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "expression_fraction_", replacement = "", fixed = T)
    GTEx.plt.df$variable = factor(x = GTEx.plt.df$variable, levels = rev(unique(GTEx.plt.df$variable)))

    GTEx.plt.df$topnote = "Fraction of Cells w/ Detectable Expr."

    p.GTEx.expfrac = ggplot(data = GTEx.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      frac_anno_palette +
      guides(fill = guide_colorbar(title.position = "top", title = "GTEx Expression Fraction", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # Expression quantiles
    # normalized
    GTEx.plt.df = dplyr::select(GTEx.data, mechanism, contains("norm_empirical_quantiles"))
    GTEx.plt.df = lapply(c("source", "target"), function(what){
      return(melt(GTEx.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(GTEx.plt.df) = c("Ligand", "Receptor")
    GTEx.plt.df = map_df(.x = (GTEx.plt.df), .f = bind_rows, .id = "src")
    GTEx.plt.df$variable = as.character(GTEx.plt.df$variable)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "source_GTEx_", replacement = "", fixed = T)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "target_GTEx_", replacement = "", fixed = T)
    GTEx.plt.df$variable = gsub(GTEx.plt.df$variable, pattern = "norm_empirical_quantiles_", replacement = "", fixed = T)
    GTEx.plt.df$variable = factor(x = GTEx.plt.df$variable, levels = rev(unique(GTEx.plt.df$variable)))

    GTEx.plt.df$topnote = "Empirical Quantiles of Norm. Expr."
    GTEx.plt.df$value = qnorm(GTEx.plt.df$value)

    p.GTEx.normq = ggplot(data = GTEx.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      quant_anno_palette +
      guides(fill = guide_colorbar(title.position = "top", title = "GTEx Empirical Quantiles", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)
  }

  ### Expression in Jerby-Arnon
  ###
  if (plot_JA != "none") {
    JA.data = dplyr::select(df, mechanism, contains("Jerby_Arnon"))

    # avg expr
    JA.plt.df = dplyr::select(JA.data, mechanism, contains("log_avg_expr"))
    JA.plt.df = lapply(c("source", "target"), function(what){
      return(melt(JA.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(JA.plt.df) = c("Ligand", "Receptor")
    JA.plt.df = map_df(.x = (JA.plt.df), .f = bind_rows, .id = "src")
    JA.plt.df$variable = as.character(JA.plt.df$variable)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "source_Jerby_Arnon_", replacement = "", fixed = T)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "target_Jerby_Arnon_", replacement = "", fixed = T)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "log_avg_expr_", replacement = "", fixed = T)
    JA.plt.df$variable = factor(x = JA.plt.df$variable, levels = rev(unique(JA.plt.df$variable)))

    JA.plt.df$topnote = "Log Average Normalized Expr."
    JA.plt.df$value = exp(JA.plt.df$value)

    p.JA.avgexp = ggplot(data = JA.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      scale_fill_scico(palette = "hawaii", direction = -1, begin = 0.05, end = 0.9, na.value = "white", trans = "log10") +
      guides(fill = guide_colorbar(title.position = "top", title = "Jerby-Arnon Expression", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # fraction expressed
    JA.plt.df = dplyr::select(JA.data, mechanism, contains("expression_fraction"))
    JA.plt.df = lapply(c("source", "target"), function(what){
      return(melt(JA.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(JA.plt.df) = c("Ligand", "Receptor")
    JA.plt.df = map_df(.x = (JA.plt.df), .f = bind_rows, .id = "src")
    JA.plt.df$variable = as.character(JA.plt.df$variable)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "source_Jerby_Arnon_", replacement = "", fixed = T)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "target_Jerby_Arnon_", replacement = "", fixed = T)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "expression_fraction_", replacement = "", fixed = T)
    JA.plt.df$variable = factor(x = JA.plt.df$variable, levels = rev(unique(JA.plt.df$variable)))

    JA.plt.df$topnote = "Fraction of Cells w/ Detectable Expr."

    p.JA.expfrac = ggplot(data = JA.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      frac_anno_palette +
      guides(fill = guide_colorbar(title.position = "top", title = "Jerby-Arnon Expression Fraction", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # Expression quantiles
    # normalized
    JA.plt.df = dplyr::select(JA.data, mechanism, contains("norm_empirical_quantiles"))
    JA.plt.df = lapply(c("source", "target"), function(what){
      return(melt(JA.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(JA.plt.df) = c("Ligand", "Receptor")
    JA.plt.df = map_df(.x = (JA.plt.df), .f = bind_rows, .id = "src")
    JA.plt.df$variable = as.character(JA.plt.df$variable)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "source_Jerby_Arnon_", replacement = "", fixed = T)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "target_Jerby_Arnon_", replacement = "", fixed = T)
    JA.plt.df$variable = gsub(JA.plt.df$variable, pattern = "norm_empirical_quantiles_", replacement = "", fixed = T)
    JA.plt.df$variable = factor(x = JA.plt.df$variable, levels = rev(unique(JA.plt.df$variable)))

    JA.plt.df$topnote = "Empirical Quantiles of Norm. Expr."
    JA.plt.df$value = qnorm(JA.plt.df$value)

    p.JA.normq = ggplot(data = JA.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      quant_anno_palette +
      guides(fill = guide_colorbar(title.position = "top", title = "Jerby-Arnon Empirical Quantiles", barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)
  }

  ### Expression in Stubenvoll, Schmidt et al.
  ###
  if (plot_StubS != "none") {
    StubS.data = dplyr::select(df, mechanism, contains("Stubenvoll_Schmidt"))

    # avg expr
    StubS.plt.df = dplyr::select(StubS.data, mechanism, contains("log_avg_expr"))
    StubS.plt.df = lapply(c("source", "target"), function(what){
      return(melt(StubS.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(StubS.plt.df) = c("Ligand", "Receptor")
    StubS.plt.df = map_df(.x = (StubS.plt.df), .f = bind_rows, .id = "src")
    StubS.plt.df$variable = as.character(StubS.plt.df$variable)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "source_Stubenvoll_Schmidt_", replacement = "", fixed = T)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "target_Stubenvoll_Schmidt_", replacement = "", fixed = T)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "log_avg_expr_", replacement = "", fixed = T)
    StubS.plt.df$variable = factor(x = StubS.plt.df$variable, levels = rev(unique(StubS.plt.df$variable)))

    StubS.plt.df$topnote = "Log Average Normalized Expr."
    StubS.plt.df$value = exp(StubS.plt.df$value)

    p.StubS.avgexp = ggplot(data = StubS.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      scale_fill_scico(palette = "hawaii", direction = -1, begin = 0.05, end = 0.9, na.value = "white", trans = "log10") +
      guides(fill = guide_colorbar(title.position = "top", title = stubs_title, barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # fraction expressed
    StubS.plt.df = dplyr::select(StubS.data, mechanism, contains("expression_fraction"))
    StubS.plt.df = lapply(c("source", "target"), function(what){
      return(melt(StubS.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(StubS.plt.df) = c("Ligand", "Receptor")
    StubS.plt.df = map_df(.x = (StubS.plt.df), .f = bind_rows, .id = "src")
    StubS.plt.df$variable = as.character(StubS.plt.df$variable)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "source_Stubenvoll_Schmidt_", replacement = "", fixed = T)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "target_Stubenvoll_Schmidt_", replacement = "", fixed = T)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "expression_fraction_", replacement = "", fixed = T)
    StubS.plt.df$variable = factor(x = StubS.plt.df$variable, levels = rev(unique(StubS.plt.df$variable)))

    StubS.plt.df$topnote = "Fraction of Cells w/ Detectable Expr."

    p.StubS.expfrac = ggplot(data = StubS.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      frac_anno_palette +
      guides(fill = guide_colorbar(title.position = "top", title = stubs_title, barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)

    # Expression quantiles
    # normalized
    StubS.plt.df = dplyr::select(StubS.data, mechanism, contains("norm_empirical_quantiles"))
    StubS.plt.df = lapply(c("source", "target"), function(what){
      return(melt(StubS.plt.df %>% dplyr::select(mechanism, contains(what)), id.vars = "mechanism"))
    })
    names(StubS.plt.df) = c("Ligand", "Receptor")
    StubS.plt.df = map_df(.x = (StubS.plt.df), .f = bind_rows, .id = "src")
    StubS.plt.df$variable = as.character(StubS.plt.df$variable)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "source_Stubenvoll_Schmidt_", replacement = "", fixed = T)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "target_Stubenvoll_Schmidt_", replacement = "", fixed = T)
    StubS.plt.df$variable = gsub(StubS.plt.df$variable, pattern = "norm_empirical_quantiles_", replacement = "", fixed = T)
    StubS.plt.df$variable = factor(x = StubS.plt.df$variable, levels = rev(unique(StubS.plt.df$variable)))

    StubS.plt.df$topnote = "Empirical Quantiles of Norm. Expr."
    StubS.plt.df$value = qnorm(StubS.plt.df$value)

    p.StubS.normq = ggplot(data = StubS.plt.df, aes(variable, mechanism)) +
      geom_tile(aes(fill = value), col = "white", size = 0.5) + #coord_equal() +
      ylab(NULL) + xlab(NULL) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = rel(1)), legend.position = "bottom",
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(r=2, l=2, unit = "pt"),
            strip.clip = "off") +
      facet_grid(cols = vars(src, topnote)) +
      quant_anno_palette +
      guides(fill = guide_colorbar(title.position = "top", title = stubs_title, barheight = unit(.4, "lines"), barwidth = unit(legend_barwidth, "lines"), direction = "horizontal")) +
      annotate("rect", xmin = 1, xmax = 2, ymin = 0.45-anno.height, ymax = 0.45, alpha = 0.)
  }

  # adjust as appropriate
  nsamples = ncol(df %>% dplyr::select(contains("log2FC")))
  width.lfc = length(nsamples)*width_per_sample + width_lfc_extra # mechanism names
  width.sig = width_signif
  if (!(is.null(plot_additional_dfs_list))) {# additional dfs can be a bit squished
   width.additionals = lapply(plot_additional_dfs_list, function(add_df){width_per_sample*ncol(add_df %>% dplyr::select(contains("log2FC")))})
  }
  width.TCGA_meansd = width_TCGA
  width.TCGA_quant = width_TCGA
  width.GTEx = with_GTEx
  width.JA = width_JA
  width.StubS = width_StubS

  p_widths = c(width.lfc)
  pps = list(lfc = p.lfc)
  if (plot_signif) {
    p_widths = c(p_widths, width.sig)
    pps[["sig"]] = p.sig
  }
  if (!(is.null(plot_additional_dfs_list))){
    for (i in seq_along(plot_additional_dfs_list)) {
      p_widths = c(p_widths, width.additionals[[i]])
      pps[[paste0("additional_",i)]] = ps.additional[[i]]
    }
  }

    if (plot_TCGA == "meansd"){
    p_widths = c(p_widths, width.TCGA_meansd)
    pps[["TCGAmeansd"]] = p.TCGA.meansd
  }
  if (plot_TCGA == "quant"){
    p_widths = c(p_widths, width.TCGA_quant)
    pps[["TCGAquant"]] = p.TCGA.quantiles
  }
  if (plot_TCGA == "both"){
    p_widths = c(p_widths, width.TCGA_meansd, width.TCGA_quant)
    pps[["TCGAmeansd"]] = p.TCGA.meansd
    pps[["TCGAquant"]] = p.TCGA.quantiles
  }

  if (plot_GTEx == "frac"){
    p_widths = c(p_widths, width.GTEx)
    pps[["GTExfrac"]] = p.GTEx.expfrac
  }
  if (plot_GTEx == "quant"){
    p_widths = c(p_widths, width.GTEx)
    pps[["GTExquant"]] = p.GTEx.normq
  }
  if (plot_GTEx == "expr"){
    p_widths = c(p_widths, width.GTEx)
    pps[["GTExexpr"]] = p.GTEx.avgexp
  }

  if (plot_JA == "frac"){
    p_widths = c(p_widths, width.JA)
    pps[["JAfrac"]] = p.JA.expfrac
  }
  if (plot_JA == "quant"){
    p_widths = c(p_widths, width.JA)
    pps[["JAquant"]] = p.JA.normq
  }
  if (plot_JA == "expr"){
    p_widths = c(p_widths, width.JA)
    pps[["JAexpr"]] = p.JA.avgexp
  }

  if (plot_StubS == "frac"){
    p_widths = c(p_widths, width.StubS)
    pps[["StubSfrac"]] = p.StubS.expfrac
  }
  if (plot_StubS == "quant"){
    p_widths = c(p_widths, width.StubS)
    pps[["StubSquant"]] = p.StubS.normq
  }
  if (plot_StubS == "expr"){
    p_widths = c(p_widths, width.StubS)
    pps[["StubSexpr"]] = p.StubS.avgexp
  }

  # legs = list(NULL)
  # legsfrominds = seq(length(pps))
  # if (plot_signif && (anno.signifi.legend %in% c("none", "plot"))) {
  #   legsfrominds = setdiff(legsfrominds, which(names(pps)=="sig"))
  # }
  # i = 2 # leaves 1 legend "box" empty, such that legends are more centered
  # for (ppsi in legsfrominds) {
  #   legs[[i]] = cowplot::get_legend(pps[[ppsi]])
  #   i = i+1
  # }
  # pps = lapply(pps, function(p){p+theme(legend.position = "none")})
  #
  # ret = cowplot::plot_grid(
  #   cowplot::plot_grid(plotlist = pps, rel_widths = p_widths, nrow = 1, axis = "tb", align = "hv"),
  #   cowplot::plot_grid(plotlist = legs, rel_widths = c(0.3, rep(1, length(p_widths))), nrow = 1, axis = "tb", align = "hv"),
  #   ncol = 1, rel_heights = c(1, legends.rel.heights)
  # )

  ret = patchwork::wrap_plots(pps, nrow=1, guides = "collect", widths = p_widths) &
    theme(legend.position = "bottom", strip.text.x.top = element_text(size = rel(1.1)))

  return(ret)
}


## Circle Plot for results

circle.vis.LRI <- function(x,
                           fontsize = 0.85, # fontsize, duh..
                           text.offset = 0.05, # text offset
                           pos = NULL, colors = NULL,
                           small.gap = 1.25, big.gap = 20,
                           h.ratio = 0.666, # curvature, 1. = straight lines, 0. = curves almost meeting in middle
                           text.height.mult = fontsize, # multiplier for the text track
                           ...){
  require(circlize)
  circ.df = tidyr::separate(x["mechanism"], mechanism, sep = "—", into = c("from", "to"), remove = T)
  circ.df$value = 1
  #circ.df = dplyr::arrange(circ.df, from)
  #group ligands, receptors and "either"
  all.genes = unique(c(circ.df$from, circ.df$to))
  grp = setNames(nm = all.genes, object =
                   factor(ifelse(all.genes %in% circ.df$from,
                                 ifelse(all.genes %in% circ.df$to,
                                        "either",
                                        "source"),
                                 "target"), levels = c("source", "either", "target")))
  # for better readability, the order of genes is important:
  # to reduce crossings, start (going clockwise) by ligands
  # whose targets are only receptors, then ligands with
  # "either" receptors, then the "either receptor or ligand"s
  # then receptors for "either" ligands and then receptors that
  # only connect to pure ligands (reversed)
  #
  g.either = all.genes[(all.genes %in% circ.df$from) & (all.genes %in% circ.df$to)]
  g.ligs = setdiff(all.genes[all.genes %in% circ.df$from], g.either)
  g.ligs.mixed = unique(dplyr::filter(circ.df, from %in% g.ligs, to %in% g.either)$from)
  g.ligs.pure = setdiff(g.ligs, g.ligs.mixed)
  g.recs = setdiff(all.genes[all.genes %in% circ.df$to], g.either)
  g.recs.mixed = unique(dplyr::filter(circ.df, to %in% g.recs, from %in% g.either)$to)
  g.recs.pure = setdiff(g.recs, g.recs.mixed)

  g.pure.pure = dplyr::filter(circ.df, from %in% g.ligs.pure, to %in% g.recs.pure)
  g.pure.mixed = dplyr::filter(circ.df, from %in% g.ligs.pure, to %in% g.recs.mixed)
  g.mixed.pure = dplyr::filter(circ.df, from %in% g.ligs.mixed, to %in% g.recs.pure)

  g.ligs.pure = rev(unique(c(rev(g.pure.mixed$from), (g.pure.pure$from))))
  g.recs.pure = unique(c(rev(g.mixed.pure$to), (g.pure.pure$to)))

  g.order = c(g.ligs.pure, g.ligs.mixed, g.either, g.recs.mixed, g.recs.pure)

  grp = grp[g.order]
  circos.par(...)
  chordDiagram(circ.df, directional = 1, group = droplevels(grp),
               direction.type = c("diffHeight+arrows"), diffHeight = convert_height(-1, "mm"),
               link.arr.type = "big.arrow", annotationTrack = "grid",
               small.gap = small.gap, big.gap = big.gap, grid.col = colors,
               preAllocateTracks = list(track.height = text.height.mult*max(strwidth(c(circ.df$from, circ.df$to)))), h.ratio = h.ratio)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + text.offset, CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0., 0.5),
                cex = fontsize)
  }, bg.border = NA)
  circos.clear()
}
