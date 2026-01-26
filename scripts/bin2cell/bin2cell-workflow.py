#!/usr/bin/env python3

import os

# os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(pow(2,40))
# os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import argparse
import json
import logging

logger = logging.getLogger(__name__)
from pathlib import Path
import sys

import scanpy as sc
import scipy
import bin2cell as b2c


def read_visisum_hd(microscope_image_path, square_002um_path, spaceranger_image_path):
    logging.info("Read Visium HD Data")
    adata = b2c.read_visium(
        square_002um_path,
        source_image_path=microscope_image_path,
        spaceranger_image_path=spaceranger_image_path,
    )
    adata.var_names_make_unique()
    return adata


def filter_data(adata, min_cells, min_counts):
    logging.info("Filter Data")
    sc.pp.filter_genes(adata, min_cells=min_cells)
    sc.pp.filter_cells(adata, min_counts=min_counts)
    return adata


def load_scalefactors_json(square_002um_path):
    # We expect to know where the 'scalefactors_json.json' file
    # is located relative to the square_002um_path
    scalefactors_path = Path(f"{square_002um_path}spatial/scalefactors_json.json")
    scalefactors = None
    if scalefactors_path.is_file():
        logging.info(f"File '{scalefactors_path}' exists")
        # Read scalefactors
        with open(scalefactors_path) as scalefactors_file:
            scalefactors = json.load(scalefactors_file)
            return scalefactors
    else:
        logging.error(f"File '{scalefactors_path}' does NOT exist")
        sys.exit()


def scale_he_image(adata, mpp, save_path):
    logging.info("Scale HE Image")
    b2c.scaled_he_image(adata, mpp=mpp, save_path=save_path)


def destripe_anndata(adata, h5_out_dir):
    logging.info("Destripe AnnData")
    b2c.destripe(adata)
    adata.write_h5ad(f"{h5_out_dir}/adata-destripe.h5ad")
    return adata


def apply_stardist(image_path, labels_npz_path, stardist_model):
    # Apply Stardist on HE image
    # 2025-03-05 Changed to use "2D_versatile_fluo" model as test!!!
    if not stardist_model in ["2D_versatile_fluo", "2D_versatile_he"]:
        logging.error(f"Unknown stardist model: {stardist_model}")
        sys.exit()
    logging.info(f"Apply Stardist to Image: {image_path}")
    b2c.stardist(
        image_path=image_path,
        labels_npz_path=labels_npz_path,
        stardist_model=stardist_model,
    )


def insert_labels_from_npz(adata, labels_key, labels_npz_path, mpp):
    logging.info("Insert Stardist Labels")
    b2c.insert_labels(
        adata,
        labels_npz_path=labels_npz_path,
        basis="spatial",
        spatial_key="spatial_cropped_150_buffer",
        mpp=mpp,
        labels_key=labels_key,
    )


def expand_labels(adata, labels_key):
    logging.info("Expand Stardist Labels")
    expanded_labels_key = f"{labels_key}_expanded"
    b2c.expand_labels(
        adata, labels_key=labels_key, expanded_labels_key=f"{labels_key}_expanded"
    )
    return adata, expanded_labels_key


def generate_gex_image(adata, mpp, save_path):
    logging.info("Check Array Coordinates")
    b2c.check_array_coordinates(adata)
    logging.info("Generate Image from Expression Data Pre check_array_coordinates")
    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, save_path=save_path)


def combine_labels(adata, primary_labels, secondary_labels, joint_labels):
    ## Combine labels "labels_he_expanded" and "labels_gex_bdata"
    logging.info("Combine Labels from Microscopy and Expression")
    b2c.salvage_secondary_labels(
        adata,
        primary_label=primary_labels,
        secondary_label=secondary_labels,
        labels_key=joint_labels,
    )
    return adata


def write_h5ad(adata, save_path):
    logging.info("Save AnnData object to h5ad file")
    adata.write_h5ad(save_path)


def bin_to_cell(adata, labels_key):
    logging.info(f"Construct Binned Cells based on Label: {labels_key}")
    b2c_adata = b2c.bin_to_cell(
        adata,
        labels_key="labels_joint",
        spatial_keys=["spatial"],
        diameter_scale_factor=None,
    )
    return b2c_adata


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
        microscope_image_path = visium_params["microscope_image_path"]
        spaceranger_image_path = visium_params["spaceranger_image_path"]
        square_002um_path = visium_params["square_002um_path"]

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
        adata = read_visisum_hd(
            microscope_image_path=microscope_image_path,
            square_002um_path=square_002um_path,
            spaceranger_image_path=spaceranger_image_path,
        )
        adata = filter_data(adata=adata, min_cells=3, min_counts=1)
        scalefactors = load_scalefactors_json(square_002um_path=square_002um_path)
        mpp = scalefactors["microns_per_pixel"]
        scale_he_image(adata=adata, mpp=mpp, save_path=stardist_he_tiff)
        adata = destripe_anndata(adata=adata, h5_out_dir=h5_out_dir)
        # Apply stardist to microscopy image
        #adata, exp_lbl_he = apply_stardist(
        apply_stardist(
            image_path=stardist_he_tiff,
            labels_npz_path=he_labels_npz,
            stardist_model="2D_versatile_he",
        )
        exp_lbl_he = scipy.sparse.load_npz(he_labels_npz)
        # Insert segmentation labels created by stardist
        he_label_key = "labels_he"
        insert_labels_from_npz(
            adata=adata, labels_key=he_label_key, labels_npz_path=he_labels_npz, mpp=mpp
        )
        # Expand segmentgation labels to bins surrounding the initial segmentation
        # Goal is to expand nuclei to "cells"
        expand_labels(adata=adata, labels_key=he_label_key)
        # Generate a TIFF image from the expresssion data
        generate_gex_image(adata=adata, mpp=mpp, save_path=gex_tiff)
        # Apply stardist to 'fluorescence' image created from expression data
        adata, exp_lbl_gex = apply_stardist(
            image_path=gex_tiff,
            labels_npz_path=gex_labels_npz,
            stardist_model="2D_versatile_fluo",
        )

        gex_label_key = "labels_gex"
        insert_labels_from_npz(
            adata=adata, labels_key=gex_label_key, labels_npz_path=gex_labels_npz
        )
        expand_labels(adata=adata, labels_key=gex_label_key)

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
