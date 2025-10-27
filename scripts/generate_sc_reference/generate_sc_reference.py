#!/usr/bin/env python3

# Import Required Libraries/Functions

import argparse
import json
import logging

logger = logging.getLogger(__name__)
from pathlib import Path
import sys

import anndata as ad
import pandas as pd
import scanpy as sc

# Functions

def create_adata_per_sample(path_to_sc_data, gsm, geo_sample_name, cell_anno, cell_type_state):
    logger.info(f"Try loading GSM {gsm} with sample name {geo_sample_name}.")
    # load scRNA data for samples: MM01, MM02, etc.
    gsm_adata = sc.read_10x_mtx(
        path=path_to_sc_data,
        prefix=f"{gsm}_{geo_sample_name}_",
        )
    # filter cell_anno for "GEO_samplename" == "MM01", "MM02", etc.
    geo_sample_anno = cell_anno[cell_anno["GEO_samplename"] == geo_sample_name]

    # Add information from cell_type_state to geo_sample_anno
    geo_sample_anno = geo_sample_anno.merge(cell_type_state, on='celltype', how='left')

    # Add index
    geo_sample_anno.index = geo_sample_anno["in_sample_barcode"]

    # Filter for cells in sample-specific AnnData object, that have cell type annotation
    gsm_adata = gsm_adata[gsm_adata.obs_names.isin(geo_sample_anno.index)].copy()

    # Filter for cells in sample-specific annotation data, that exist in scRNA data/AnnData object
    geo_sample_anno = geo_sample_anno[geo_sample_anno.index.isin(gsm_adata.obs_names)]

    # Reindex GEO sample annotation data ...
    geo_sample_anno.reindex(gsm_adata.obs_names)

    # ... so we can add it to the AnnData object
    geo_sample_anno_cat = geo_sample_anno[[
        "celltype",
        "celltype_medium",
        "celltype_broad",
        "celltype_heatmap_anno",
        "sum_type"]]
    gsm_adata.obs = pd.concat([gsm_adata.obs, geo_sample_anno_cat], axis=1)
    gsm_adata.obs = gsm_adata.obs.astype('category')
    return gsm_adata

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
    # args = parser.parse_args(["--json", "assets/segmentation2cells/LE-585-HM.json"])
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

    # Let's load the cell annotation information ...
    cell_anno = pd.read_csv(
        filepath_or_buffer=params["celltype_annotation"],
        sep="\t")

    # ... and the information about cell type states
    cell_type_state = pd.read_csv(
        filepath_or_buffer=params["cell_type_state"],
        sep="\t")

    # Path to directory containing the scRNA for all samples
    path_to_sc_data = params["path_to_sc_data"]

    # Map GSM accession to sample name 
    sample_name_map = params["sample_name_map"]

    # Initialize empty hash that will hold all AnnData objects of all scRNA samples
    adatas = {}
    # Iterate through sample_name_map and fill adatas with AnnData objects for each scRNA sample
    for k, v in sample_name_map.items():
        adatas[v] = create_adata_per_sample(
            path_to_sc_data=path_to_sc_data,
            gsm=k,
            geo_sample_name=v,
            cell_anno=cell_anno,
            cell_type_state=cell_type_state)

    sc_ref = ad.concat(adatas)
    sc_ref.obs_names_make_unique()

    # Make sure out path exists
    output_dir = Path(params["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Assemble path to h5ad file
    h5ad_file = output_dir.joinpath(params["h5ad_file"])
    # Save single cell reference to h5ad file
    sc_ref.write_h5ad(filename=h5ad_file)

if __name__ == "__main__":
    main()
