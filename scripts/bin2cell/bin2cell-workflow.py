#!/usr/bin/env python3

import os

import argparse
import json
import logging

logger = logging.getLogger(__name__)
from pathlib import Path
import sys

import anndata
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

import bin2cell as b2c

from spatial_transcriptomics_analysis.io import (
    load_scalefactors_json,
    load_visium_hd,
)
from spatial_transcriptomics_analysis.bin2cell_wrapper import (
    scale_he_image,
    destripe_anndata,
    apply_stardist,
    generate_gex_image,
    combine_labels,
    write_h5ad,
    bin_to_cell,
)


def filter_data(adata, min_cells_with_gene, min_counts_per_cell):
    logging.info("Filter Data")
    sc.pp.filter_genes(adata, min_cells=min_cells_with_gene)
    sc.pp.filter_cells(adata, min_counts=min_counts_per_cell)
    return adata


def main():
    logging.basicConfig(format="%(levelname)s:%(message)s", level=logging.INFO)
    # Define argument parser
    parser = argparse.ArgumentParser()
    # Tell them we started
    logger.info(f"Started {parser.prog}")
    parser.add_argument(
        "--json",
        help="Path to json file containing configuration parameter",
        required=True,
    )
    parser.add_argument(
        "--verbose", help="increase output verbosity", action="store_true"
    )
    args = parser.parse_args()
    # args = parser.parse_args(["--json", "assets/manifest.json"])
    # Read params from JSON file
    json_path = Path(args.json)
    params = None
    if json_path.is_file():
        with open(json_path) as json_file:
            params = json.load(json_file)
    else:
        logging.error(f"File {json_path} does not exist!")
        sys.exit()

    # Iterate over all samples in JSON file
    for sample in params.keys():
        logging.info(f"Start 'bin2cell' workflow on sample: {sample}")
        # Get parameters for each sample
        visium_params = params[sample]["visium_hd"]
        microscope_image_path = Path(visium_params["microscope_image_path"])
        spaceranger_image_path = Path(visium_params["spaceranger_image_path"])
        square_002um_path = Path(visium_params["square_002um_path"])

        output_dir = visium_params["output_dir"]
        stardist = f"{output_dir}/{sample}/stardist"
        h5_out_dir = f"{output_dir}/{sample}/h5-files"
        os.makedirs(stardist, exist_ok=True)
        os.makedirs(h5_out_dir, exist_ok=True)
        # Params for bin2cell.scaled_he_image()
        stardist_he_tiff = f"{stardist}/HE.tiff"
        stardist_dapi_tiff = f"{stardist}/DAPI.tiff"

        he_labels_npz = f"{stardist}/HE.npz"
        dapi_labels_npz = f"{stardist}/DAPI.npz"

        # Params for bin2cell.scaled_
        gex_pre_tiff = f"{stardist}/GEX.pre_cac.tiff"
        gex_labels_pre_npz = f"{stardist}/GEX.pre_cac.npz"
        gex_tiff = f"{stardist}/GEX.tiff"
        gex_labels_npz = f"{stardist}/GEX.npz"

        # Work through bin2cell workflow

        # Load Visium HD data to AnnData
        adata = load_visium_hd(
            microscope_image_path=microscope_image_path,
            square_002um_path=square_002um_path,
            spaceranger_image_path=spaceranger_image_path,
        )

        # Store raw counts in adata.raw.X
        adata.raw = adata.copy()

        # Filter AnnData
        adata = filter_data(adata=adata, min_cells_with_gene=3, min_counts_per_cell=1)

        # Let's try if bin2cell works without image scaling

        scalefactors = load_scalefactors_json(path_to_search=square_002um_path)
        mpp = scalefactors["microns_per_pixel"]
        spatial_key = "spatial"
        if params[sample]["scale_image"]:
            logger.info("Scale HE image")
            scale_he_image(adata=adata, mpp=mpp, save_path=stardist_he_tiff)
            # if we scale the image:
            # - the scaled image is used down-stream
            microscope_image_path = stardist_he_tiff
            # - adata.obsm contains the key "spatial_cropped_150_buffer"
            spatial_key = 'spatial_cropped_150_buffer'
        else:
            logger.info("Not scaling HE image")

        adata = destripe_anndata(adata=adata, h5_out_dir=h5_out_dir)
        # Apply stardist to microscopy image
        adata, exp_lbl_he = apply_stardist(
            adata=adata,
            image_path=microscope_image_path,
            labels_npz_path=he_labels_npz,
            labels_key="labels_he",
            mpp=mpp,
            stardist_model="2D_versatile_he",
            spatial_key=spatial_key,
        )
        generate_gex_image(adata=adata, mpp=mpp, save_path=gex_tiff)
        # Apply stardist to 'fluorescence' image created from expression data
        adata, exp_lbl_gex = apply_stardist(
            adata=adata,
            image_path=gex_tiff,
            labels_npz_path=gex_labels_npz,
            labels_key="labels_gex",
            mpp=mpp,
            stardist_model="2D_versatile_fluo",
            spatial_key=spatial_key
        )

        adata = combine_labels(
            adata=adata,
            primary_labels=exp_lbl_he,
            secondary_labels=exp_lbl_gex,
            joint_labels="labels_joint",
        )
        write_h5ad(adata=adata, save_path=f"{h5_out_dir}/adata-pre-b2c.h5ad")
        adata = bin_to_cell(adata=adata, labels_key="labels_joint")
        write_h5ad(adata=adata, save_path=f"{h5_out_dir}/adata-post-b2c.h5ad")
    logger.info("Finished")


if __name__ == "__main__":
    main()
