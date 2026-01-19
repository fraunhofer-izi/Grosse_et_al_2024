# Grosse etal. 2026

This repository contains code used to produce the results in:

> "Florian Große et al. Spatial transcriptomics of melanoma tissues reveals new tumor– immune cell interactions in humans. (not released yet)"

# Singularity

All scripts were developed in an Rstudio server running in a Singularity container. Package versions are controlled via renv which should allow anyone to replicate the exact state of package versions used to produce the results.

# Analysis

1. Stain deconvolution (using [`stain_deconvolution.py`](scripts/stain_deconvolution/stain_deconvolution.py))
2. Segmentation (using [`bin2cell-workflow.py`](scripts/bin2cell/bin2cell-workflow.py))
3. Generate single cell reference (using [`generate_sc_reference.py`](scripts/generate_sc_reference/generate_sc_reference.py))
3. Cell type annotation (using [`tacco-annotate.py`](scripts/tacco/tacco-annotate.py))
4. Rotate image data (using [`rotate-sample.py`](scripts/rotate_sample/rotate-sample.py))


# Content

# Reproduction
