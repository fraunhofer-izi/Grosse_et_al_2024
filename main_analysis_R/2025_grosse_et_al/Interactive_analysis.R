# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# This file is supposed to guide through the analysis in case that somebody
# wants to reproduce it all. If you don't intend to set cutoffs and cluster
# numbers and the like yourself again, feel free to keep the values in the
# manifest file. This is were I saved all sample-specific parameters. 
# 
# If you do however want to redo the whole analysis, maybe even on new samples,
# then I suggest deleting all the values in the manifest that are not paths.
# With that said, I do warn the reader that I initially had no clue how to properly
# set up such a workflowr repo for the repeated analyis of multiple samples,
# so forgive me that I did hack the part that renders the same markdown document
# for different samples together. It should work relatively decently though.
#
# Since the analysis is somewhat of an interactive process where you have
# to control and possibly adjust parameters based on the specific sample
# at various steps, the idea here is to run the analysis up to each such
# point, build and view the corresponding output .html page, possibly adjust
# some settings (QC cutoffs, appropriate clustering resolution, annotations,
# deconvolution references, filtering settings, etc.) and move on to the next
# slide until all are done with a specific step.
# The .html output of workflowr is finally moved to the proper
# output directory of that sample.
# One last word: A .tmp folder and file is created to keep track of the active 
# sample between wflow_build()'s. This is somewhat not-pretty, I know.
#
# Addendum from in between the work:
# The initial parts are quite fast and best used interactively
# (definition of functions, setup of indices, spot QC).
# Starting from ~ the transfer of annotations onwards, it's more comfortable to
# run the scripts as jobs. This is easy when working in rstudio ("run as local job"),
# with the r_bg function from callr (use instead of rscript) or your favorite
# cluster manager (at least for all the script (everything in ./code but not ./analysis)).
# The workflowr knitting currently has to be done in serial though, unfortunately. 
#
# Final words: if you only intend to reproduce the images from intermediate
# objects, have a look at the Reproduction_of_paper_figures.R in this folder.
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Setup and stuff
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Before I incorporate a huge installation chunk here, just make sure that the 
# libraries are there. Ideally restore the correct versions using renv first.
# Specific versions are borderline irrelevant (only!) here,
# since these functions never directly touch the sample data.
library(workflowr)
library(rmarkdown)
library(yaml)
library(callr)

# please make sure that the manifest is in proper shape and all samples have at least some entry
manifest = yaml.load_file("./assets/manifest.yaml")

# set flags in workdir to build workflowr pages for specific samples
if(!dir.exists(file.path(manifest$workdata_parent, ".tmp"))){
  dir.create(file.path(manifest$workdata_parent, ".tmp"))
}
tmp.file <- file.path(file.path(manifest$workdata_parent, ".tmp/active_sample"))
if(file.exists(tmp.file)){
  print("Warning: There's already an active_sample temporary file. You may want to look at it if you stopped mid-analysis last time")
} else {
  file.create(tmp.file)
  cat("TEST", file = tmp.file)
  test.content <- scan(file.path(manifest$workdata_parent, ".tmp/active_sample"), what = "character", quiet = T)
  if (test.content == "TEST") {
    print("Test of temporary file successful.")
  } else {
    stop("Temp file setup doesn't seem to work!")
  }
}
# makes setting the active sample somewhat prettier
active_sample <- function(sample) {
  print(paste("Active sample is now:", sample))
  cat(sample, file = tmp.file)
}

# small helper function to copy all resuls to a new location in a way such that the
# report pages are viewable in the browser (i.e. links point to the right location)
# currently assumes that all sample's results folders (manifest[[sample_id]]$tertiary)
# are named by their sample ID, because that is what the global_index.html file
# currently assumes as relative path
#
# If you do not deploy to GitLab pages in the end, use this function to copy the 
# final report pages to their final location in an accessible format.
#
copy_results <- function(destination, manifest){
  if (!dir.exists(destination)){
    dir.create(destination, recursive = T)
  }
  file.copy(from = file.path(manifest$results_directory_parent, manifest$global_index_dir, "."),
            to = destination, recursive = T)
  for (s in c(manifest$samples, manifest$external_samples)){
    file.copy(from = file.path(manifest$results_directory_parent, manifest[[s]]$tertiary),
              to = destination, recursive = T)
  }
}

# So, this project was originally intended to be viewed via GitLab pages. 
# Each sample has it's own index file and each index refers to a shared global 
# index. We let workflowr handle the sample indices for now and start by setting 
# up the global index ourselves. We have slightly different _site.yml files for
# normal pages, index pages and the global index, so we kinda need to switch
# between those files a few times.

# This function is the main workhorse here. It renders each analysis file with
# the appropriate base template (since e.g. the analysis files should refer back
# to either the global or the appropriate sample index, and each sample index
# should refer back to the global index etc.) which gets augmented by workflowr
# with the reproducibility section. The generated files are then moved to the
# specified destination and the previous state of the folder restored at the end.

render_with_template <- function(rmd.file, rename.to.before.render = NULL,
                                 destination.dir, overwrite.results = T,
                                 file.and.template.folder = "./analysis", template = "_site.yml",
                                 use.rmarkdown.render = F, unlink.after = T){
  if (!(dir.exists(destination.dir))) {dir.create(destination.dir, recursive = T)}
  # clean up public in case there are files left from something else
  unlink("./public/*", recursive = T, expand = T)
  maybe_unlink <- function(unlink.after){
    if (unlink.after) {
      unlink("./public/*", recursive = T, expand = T)
      message("Successfully cleaned public/ folder.")
    } else {
      message("No cleaning of public/ folder was performed.")
    }
  }
  tryCatch({
    # main part here
    templatebackup = NULL
    maybe_recover_template = function(from, basepath){
      if (!(is.null(from))) {
        file.remove(file.path(basepath, "_site.yml"))
        file.rename(from = from, file.path(basepath, "_site.yml"))
      }
    }
    if (template != "_site.yml") {
      # switch template only if necessary
      templatebackup = tempfile(pattern = "template_backup",
                                tmpdir = file.and.template.folder,
                                fileext = ".yml")
      file.rename(file.path(file.and.template.folder, "_site.yml"), templatebackup)
      # copy is safer than remove
      file.copy(file.path(file.and.template.folder, template),
                file.path(file.and.template.folder, "_site.yml"))
    }
    # template is correct now, poke renderer and make sure to restore template in any case
    tryCatch({
      renderfile = file.path(file.and.template.folder, rmd.file)
      originalname = NULL
      file.backup = NULL
      maybe_recover_rmd = function(renderfile, originalname, file.backup){
        if (!(is.null(originalname))) { # if file was renamed/copied
          file.remove(renderfile)
        }
        if (!(is.null(file.backup))) {
          file.rename(file.backup, renderfile)
        }
      }
      if (!(is.null(rename.to.before.render))){
        originalname = renderfile
        renderfile = file.path(file.and.template.folder, rename.to.before.render)
        # if such a file already exists
        if (file.exists(renderfile)) {
          file.backup = tempfile(pattern = "file_backup",
                                 tmpdir = file.and.template.folder,
                                 fileext = ".rmd")
          file.rename(renderfile, file.backup)
          # copy is safer than remove
          file.copy(originalname, renderfile)
        }
      }
      tryCatch({
        if (use.rmarkdown.render) {
          render_site(input = renderfile)
        } else {
          wflow_build(files = renderfile)
        }
      },
      finally = maybe_recover_rmd(renderfile = renderfile, originalname = originalname, file.backup = file.backup))
      file.copy("./public/.", destination.dir, recursive = T, overwrite = overwrite.results)
    },
    finally = maybe_recover_template(templatebackup, basepath = file.and.template.folder))
  },
  finally = maybe_unlink(unlink.after = unlink.after))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Actual analysis
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

all_samples = c(manifest$samples, manifest$external_samples)
good_samples = setdiff(c(manifest$samples, manifest$external_samples), manifest$rejected_samples)




# render global index site with appropriate template
render_with_template(rmd.file = "global_index.Rmd", rename.to.before.render = "index.Rmd",
                     template = "_global_index_site.yml", use.rmarkdown.render = T,
                     destination.dir = file.path(manifest$results_directory_parent, manifest$global_index_dir))

# now we set up the index pages for all samples to populate our folder structure
# a little remember that we switch to _index_site.yml (this one has the home
# button in the top left corner replaced with a button "Switch Samples" linking
# to the global index page)
for (sample in all_samples){
  # remember to store the active sample in the temporary file
  active_sample(sample)
  destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
  render_with_template(rmd.file = "index.Rmd", template = "_index_site.yml",
                       destination.dir = destination)
}

# We also need to copy the logo file to the results directories,
# so that it appears in the page headers
logo.destinations = normalizePath(file.path(manifest$results_directory_parent, c(manifest$global_index_dir, sapply(manifest$samples, function(s){manifest[[s]]$tertiary}))))
for (dest in file.path(logo.destinations, "assets/fhglogoizi.gif")){
  if (!(dir.exists(dirname(dest)))) {dir.create(dirname(dest))}
  file.copy(from = "./assets/fhglogoizi.gif", to = dest)
}

# Now comes the real analysis part. This will consist of steps alternating
# between R scripts to do the heavy working, possibly in parallel and building
# the workflowr sites for analysis

# The first step is sample preparation, loading pictures into Staffli objects,
# normalization and scaling and the like. The script is pretty self-contained
# and processes all samples in the manifest sequentially. There's not much to
# look at or adjust except from the cropping maybe. If you need to adjust
# the image cropping, you can set the sample manually and follow the steps
# in the script up to the cropping step and play around with it until you're
# happy with the result. The final cropping settings are saved in the manifest

# If everything is alright with script and manifest, simply run this:
rscript("./code/01_pipeline.R")
# You have to wait until this is finished of course.

# The results of these basic steps should reside in the working directories now
# (see manifest). We can now build the basic QC overview for each sample
for (sample in all_samples){
  active_sample(sample)
  destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
  render_with_template(rmd.file = "01_qc.Rmd", destination.dir = destination)
}

# The next step is adding annotations to the datasets (removing erroneous tissue spots)
# then do QC and build the page then build a page where different clustering
# resolutions can be compared to the expert annotations

  # Make sure that all annotations you have (as .csv files) are in the corresponding
  # folders for each sample (see manifest). The annotations file should be called
  # annotations.csv, and there should be a .png file in the annotations folder:
  #   annotations.png       with the expert annotation to compare to the resulting per-spot annotation
  rscript("./code/02_Kunz_annotations.R")

  ## Spot QC and Subsetting

  # We are going to discard spots for the coming analysis that are low-quality,
  # meaning that either the amount of detected genes or reads is unusually low
  # (this makes the expression patterns extremely noisy) or the number of reads
  # mapping to mitochondrial genes is very high (pointing to popped cells in the
  # sample preparation) or both

  # To determine appropriate measures and cutoffs, all samples are pooled first
  # and jointly clustered (after batch correction probably). Then the quality
  # metric distribution is judged across the combined sample clusters to see if
  # a cluster attracts low-quality spots or not.
  # In this analysis, no low-Q cluster was found, so we went on with determining
  # a threshold for low-quality spots based on the count/gene distributions
  # across samples.

  # This step uses the settings from the manifest under
  # manifest > QC_cutoffs > [genes|counts|mito]
  # If those are not already set, enter anything to get it running,
  # then run the script and .html build once, look at the results and determine
  # appropriate cutoffs from the histograms/prio knowledge or best practices
  # (enter in manifest!) and build again to look at the results.
  # If everything worked as expected, run the second script to do the final
  # subsetting and continue below.
  # Alternatively: run the first script and manually plot the global statistics from
  # analysis/qc_select.Rmd, then set the cutoffs and proceed with second script and builds

  # This script is self-contained like the first and processes all samples
  rscript("./code/03_global_QC.R")

  # This may be the first step that warrants looking at the whole thing,
  # so either copy locally using the function defined above, or trigger the
  # Git[Lab|Hub] CI to publish results.
  #
  # IZI Internal note:
  # If something fails with the CI, it is often the case that some of the
  # generated files lack group read permissions. You may need to set
  # those in that case.
  # Usually to do so go to: /mnt/.../2021-uccl-melanom-spacialseq-tertiary
  # And issue: chmod -c -R g+rx .
  # To make everything accessible to the CI runner
  #
  # Now build the pages.

  for (sample in all_samples){
    active_sample(sample)
    destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
    render_with_template(rmd.file = "03_qc_select.Rmd", destination.dir = destination)
  }

  # After cutoffs have been determined, run this to subset each sample accordingly
  rscript("./code/04_QC_selection.R")

  # And finally, visualize histological annotations and generated clusterings

  for (sample in all_samples){
    active_sample(sample)
    destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
    render_with_template(rmd.file = "04_annotations.Rmd", destination.dir = destination)
  }


## Integrated clustering worked somewhat fine but this direction was not
## pursued any further, so this step is currently skipped

# # we can also now do the umap plots for groups of samples
# rscript("./code/05_global_01_UMAPS.R")
#
# destination = file.path(manifest$results_directory_parent, manifest$global_index_dir)
# render_with_template(rmd.file = "05_global_01_umaps.Rmd", destination.dir = destination)

# Prepare the Stubenvoll, Schmidt et al. 2025 single cell dataset as deconvolution reference
rscript("./code/05a_prep_sc_ref.R")

# Run deconvolution on multiple levels
rscript("./code/05b_run_CARD_RCTD.R")

# Report deconvolution results
for (sample in all_samples){
  active_sample(sample)
  destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
  render_with_template(rmd.file = "05_deconv_results.Rmd", destination.dir = destination)
}

# this script generates auxillary marker tables from
# the jerby-arnon publication data and panglaoDB
# it will cache results and be pretty much instantaneous after
# the first run, but that first run will take a while
# IF: NCBI acts up again and won't let you download certain files from the script,
# you can go through the code, set up all directories and download the input files
# via the links provided or by searching the GEO entry with the corresponding ID
# and downloading it yourself. This script will not download already present files
rscript("./code/06a_generate_marker_lists.R")

# this script takes the consolidated objects and adds the module scores
rscript("./code/06b_module_scores.R")

# Now build reports pages for all module score analyses
for (sample in all_samples){
  active_sample(sample)
  destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
  render_with_template(rmd.file = "06_module_scores.Rmd", destination.dir = destination)
}

# make sure to take a look at the pages, specifically module scores and the generated
# clusterings and then select the best possible clustering resolution
# for each sample and run the following to generate a deeper analysis of this clustering
rscript("./code/07_clusterMarker_and_ora.R")

for (sample in all_samples){
  active_sample(sample)
  destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
  render_with_template(rmd.file = "07_clustering.Rmd", destination.dir = destination)
}

# This was the first attempt at rendering plots necessary for the paper in a structured
# form. It's no longer directly used but may be nice to look at nontheless.
for (sample in all_samples){
  active_sample(sample)
  destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
  render_with_template(rmd.file = "08_paper_plots.Rmd", destination.dir = destination)
}

###
# From here on, samples that have been determined to be of generally low quality
# or otherwise unsuitable (e.g. due to laying mostly outside of capture area or
# being damaged, infiltration too low for LRI) can be ignored when building reports
# (uses good_samples instead of all_samples)
###

# ligand-receptor interaction analyses for each sample separately
rscript("./code/09_combined_LRI_test.R")

# per sample reports
for (sample in good_samples){
  active_sample(sample)
  destination = file.path(manifest$results_directory_parent, manifest[[sample]]$tertiary)
  render_with_template(rmd.file = "09_LRI.Rmd", destination.dir = destination)
}

# prepare GTEx, Jerby-Arnon and Stubenvoll Schmidt datasets for heatmap annotations
rscript("./code/10_prep_SC_annos.R")

# using a combined approach, aggregate individual ligand-receptor analyses
render_with_template(rmd.file = "10_global_02_aggregated_LRI.Rmd", destination.dir = file.path(manifest$results_directory_parent, manifest$global_index_dir))

# basic stats for selected samples
render_with_template(rmd.file = "11_global_03_basic_spot_stats.Rmd", destination.dir = file.path(manifest$results_directory_parent, manifest$global_index_dir))

# steps 12 and 13 used to be analysis of interferone pathways and another collection of
# plots for each sample, which were discarded eventually

# prepare survival data for analysis
#
# TODO: this file needs some manual work from you, 
# please take a look at it, download and process external files for the reproduction
# 
rscript("./code/14_prep_survival_cohorts.R")

# Initial survival analysis, contains a lot of exploratory stuff
render_with_template(rmd.file = "14_global_04_TCGA_survival.Rmd", destination.dir = file.path(manifest$results_directory_parent, manifest$global_index_dir))


### Compiled plots for the paper:
# now finally, we compiled two pages that contain all plots needed for the paper
# (split into non-HD + survival, which uses samples analyzed in this repo +
# HD samples, which were mostly preprocessed with python (scikit-image/stardist/scanpy/TACCO)
# and handed over as anndata .h5ad files )

render_with_template(rmd.file = "15_paper_plots_merged.Rmd", destination.dir = file.path(manifest$results_directory_parent, manifest$global_index_dir))
render_with_template(rmd.file = "16_HD_plots.Rmd", destination.dir = file.path(manifest$results_directory_parent, manifest$global_index_dir))

# either use this to copy the results to a destination in proper 
# directory structure (such that you only need to open the topmost 
# index.html and everything is linked from there), or use GitLab pages
# with a similar copy script to deploy results in the web.
# copy_results(manifest=manifest, ~some accessible path to save to~)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# (╮°-°)╮┳━━┳                        BYE                         ( ╯°□°)╯ ┻━━┻
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
