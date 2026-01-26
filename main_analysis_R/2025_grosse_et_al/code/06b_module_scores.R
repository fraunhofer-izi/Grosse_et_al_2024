# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# INFO:
# Calculate all sorts of module scores for each spot in each sample
# Take:
#   04_seurat.rds
# Store:
#   06_seurat.rds
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Packages
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
library("tidyverse")
library("Seurat")
library("yaml")
library("Matrix")
library("fgsea")
library("readxl")
library("stringi")
library("robustbase")

# bioc
library("GSEABase")
library("BiocParallel")


# load a fixed version of AddModuleScore
source("code/helper/FG_modulescores_helper.R")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

manifest = yaml.load_file("assets/manifest.yaml")
all_samples = c(manifest$samples, manifest$external_samples)
random_seed = 1234
set.seed(random_seed)

all_gene_symbols_in_ses = character()

ses = lapply(all_samples, function(sample_id){
  print(paste("Reading input for sample", sample_id))
  sample_manifest = manifest[[sample_id]]
  work.dir = file.path(manifest$workdata_parent, sample_manifest$workdata)
  se = readRDS(file.path(work.dir, "04_seurat.rds"))
  DefaultAssay(se) = "RNA_logcounts"
  return(list(sid = sample_id, se = se, workdir = work.dir))
})

for (obj in ses) {
  all_gene_symbols_in_ses = unique(c(all_gene_symbols_in_ses, rownames(obj$se@assays$RNA_logcounts)))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# all precomputed marker lists
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


obj = readRDS(file = "./data/genesets/collected_markers.rds")
panglao_gs = obj$panglao_gs
# additionally, make double sure that nothing in R complains about wrong symbols (-/ ) in the names
names(panglao_gs) = make.names(names(panglao_gs))
jerby_arnon_sole_boldo_markers = obj$final.markers
names(jerby_arnon_sole_boldo_markers) = make.names(names(jerby_arnon_sole_boldo_markers))

# all filtering has been changed to happen in the 05_generate_marker_lists script

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# additional genesets
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Hypoxia scores from Buffa et al. 2010 (PMID: 20087356)
buffa_hypoxia_full = c("GPI",    "UTP11",  "ENO1",   "P4HA1",   "PFKP",   "CA9",   "PGAM1",  "PNP",    "HILPDA",   "MRGBP",  "AK4",   "ALDOA",   "ADM",   "SEC61G",
                       "ACOT7",  "PSMA7",  "LDHA",   "HK2",     "NDRG1",  "TPI1",  "SLC2A1", "MRPL15", "SLC25A32", "TUBB6",  "DDIT4", "CDKN3",   "VEGFA", "MRPS17",
                       "PGK1",   "BNIP3",  "CORO1C", "ANKRD37", "MAP7D1", "MIF",   "MCTS1",  "MAD2L2", "MRPL13",   "SHCBP1", "GAPDH", "SLC16A1", "YKT6",  "ESRP1",
                       "KIF20A", "TUBA1B", "TUBA1C", "CHCHD2",  "ANLN",   "PSRC1", "KIF4A",  "CTSV",   "LRRC42")
buffa_hypoxia_top15 = c("VEGFA", "SLC2A1", "PGAM1", "ENO1", "LDHA", "TPI1", "P4HA1", "MRPS17", "CDKN3", "ADM", "NDRG1", "TUBB6", "ALDOA", "MIF", "ACOT7")

other_scores = list(buffa_hypoxia_full = buffa_hypoxia_full, buffa_hypoxia_top15 = buffa_hypoxia_top15)

panglao_gs = lapply(panglao_gs, function(x)update_genes(x, known.genes = all_gene_symbols_in_ses))
jerby_arnon_sole_boldo_markers = lapply(jerby_arnon_sole_boldo_markers, function(x)update_genes(x, known.genes = all_gene_symbols_in_ses))
other_scores = lapply(other_scores, function(x)update_genes(x, known.genes = all_gene_symbols_in_ses))

# Add module scores to samples and save

for (se.i in seq_along(ses)) {
  entry = ses[[se.i]]
  sample_id = entry$sid
  se = entry$se
  work.dir = entry$workdir
  ms <- list()

  print(paste("Calculating module scores for sample:", sample_id))

  ### Non-normalized module scores
  DefaultAssay(se) = "RNA_logcounts"

  # panglao logcounts
  nm = paste0(names(panglao_gs), "_logcounts")
  se = add_ms_with_checks(se = se, features = panglao_gs, names = nm)
  ms$panglao_logcounts = nm

  # jerby-arnon/sole-boldo logcounts
  nm = paste0(names(jerby_arnon_sole_boldo_markers), "_logcounts")
  se = add_ms_with_checks(se = se, features = jerby_arnon_sole_boldo_markers, names = nm)
  ms$jerby_logcounts = nm

  # other logcounts
  nm = paste0(names(other_scores), "_logcounts")
  se = add_ms_with_checks(se = se, features = other_scores, names = nm)
  ms$other_scores_logcounts = nm

  ### Normalized module scores
  DefaultAssay(se) = "Spatial"

  # panglao nprmalized
  nm = paste0(names(panglao_gs), "_normalized")
  se = add_ms_with_checks(se = se, features = panglao_gs, names = nm)
  ms$panglao_normalized = nm

  # jerby-arnon/sole-boldo normalized
  nm = paste0(names(jerby_arnon_sole_boldo_markers), "_normalized")
  se = add_ms_with_checks(se = se, features = jerby_arnon_sole_boldo_markers, names = nm)
  ms$jerby_normalized = nm

  # other normalized
  nm = paste0(names(other_scores), "_normalized")
  se = add_ms_with_checks(se = se, features = other_scores, names = nm)
  ms$other_scores_normalized = nm

  Misc(se, "module_scores") <- ms

  # also record the size of each applied geneset for later
  Misc(se, "applied_genesets") <- c(panglao_gs, jerby_arnon_sole_boldo_markers, other_scores)
  Misc(se, "module_score_sizes") <- lapply(Misc(se, "applied_genesets"), length)

  # Save and quit
  saveRDS(se, file = file.path(work.dir, "06_seurat.rds"))
}
print("FIN")
