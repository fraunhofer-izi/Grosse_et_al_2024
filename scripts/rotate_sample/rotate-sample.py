#!/usr/bin/env python3

# Import Required Libraries/Functions

import argparse
import json
import logging

logger = logging.getLogger(__name__)
import numpy as np
from pathlib import Path
import sys
import anndata as ad
import tacco as tc

import matplotlib.patches as patches
from matplotlib import pyplot as plt
import matplotlib as mpl

import pandas as pd
import scanpy as sc
from scipy.ndimage import rotate
from skimage.color import rgb2gray
import spatialdata as sd
import spatialdata_io as sd_io
import spatialdata_plot
import squidpy as sq

# Functions


def rotate_images_in_adata(adata, sample_id, angle=0):

    # Rotate the spatial coordinates by a certain angle (e.g., 45 degrees)
    theta = np.radians(angle)  # angle in degrees
    rotation_matrix = np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
    )

    ## Rotate the spatial coordinates
    coords = adata.obsm["spatial"]
    rotated_coords = coords @ rotation_matrix.T

    # Save original coordinates if needed
    adata.obsm["spatial_original"] = adata.obsm["spatial"].copy()

    # Replace spatial coordinates with rotated version
    adata.obsm["spatial"] = rotated_coords

    # Save original image
    adata.obsm["spatial_original"] = adata.obsm["spatial"].copy()

    ## Rotate the images
    for res in ["hires", "lowres"]:

        # Extract image
        img = adata.uns["spatial"][sample_id]["images"][res]

        # Rotate image 90 degrees counter-clockwise
        rotated_img = rotate(
            img, angle=-angle, reshape=True
        )  # reshape=True to keep full image

        # If needed, overwrite or store separately
        adata.uns["spatial"][sample_id]["images"][res] = rotated_img

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

    # Read params from JSON file

    json_path = Path(args.json)
    params = None
    if json_path.is_file():
        with open(json_path) as json_file:
            params = json.load(json_file)
        logging.info(f"Read params from JSON file: {json_path}")
    else:
        logging.error(f"File {json_path} does not exist!")
        sys.exit()

    # iterate over samples in JSON file
    for sample in params.keys():

        # Make sure out path exists
        output_dir = Path(params[sample]["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)

        for ad_name, ad_path in params[sample]["anndata_files"].items():
            #
            ad_path = Path(ad_path)

            # Load AnnData with cells/rows to be annotated/categorized
            logging.info(f"Load AnnData from {ad_path}")
            spatial_adata = ad.read_h5ad(ad_path)

            rotate_adata = rotate_images_in_adata(
                adata=spatial_adata,
                sample_id=params[sample]["image_id"],
                angle=int(params[sample]["rotation_angle"]),
            )

            # Assemble output file name
            rotate_file_name = f"{ad_path.stem}-rotate{ad_path.suffix}"
            rotate_adata_path = output_dir.joinpath(rotate_file_name)

            logging.info(f"""Write rotated AnnData to file: {rotate_adata_path}""")

            rotate_adata.write_h5ad(rotate_adata_path)


if __name__ == "__main__":
    main()
