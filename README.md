# Grosse etal. 2026

This repository contains code used to produce the results in:

> "Florian Große et al. Spatial transcriptomics of melanoma tissues reveals new tumor– immune cell interactions in humans. (not released yet)"


## Analysis

This repository contains the code to analyse spatial transcriptomics data generated using 10x Genomics Visium and Visium HD platforms.

### Visium Data Analysis

The analysis of the lower resolution Visium samples was conducted in R.
Please refer to the [specific README file for the Visium analysis](main_analysis_R/README.md) for more information and usage instructions.

### Visium HD Data Analysis

The analysis of the higher resolution Visium HD samples was conducted both in Python and R.


1. Stain deconvolution (using [`stain_deconvolution.py`](scripts/stain_deconvolution/stain_deconvolution.py))
2. Segmentation (using [`bin2cell-workflow.py`](scripts/bin2cell/bin2cell-workflow.py))
3. Generate single cell reference (using [`generate_sc_reference.py`](scripts/generate_sc_reference/generate_sc_reference.py))
3. Cell type annotation (using [`tacco-annotate.py`](scripts/tacco/tacco-annotate.py))
4. Rotate image data (using [`rotate-sample.py`](scripts/rotate_sample/rotate-sample.py))

## Coding Guidelines

The Python code was developed without using a particular coding guideline.

## Linter

We rely on the [`black`](https://github.com/psf/black) linter to format code in a consistent manner.

## Reproduction

We use `conda` for managing software dependencies.
Thus you need to install [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) to run the Visium HD analysis.
The repository comes with the defintion of a `conda` environment.
This environemnt can be restored by:

```bash
$ make conda_env_create
```
