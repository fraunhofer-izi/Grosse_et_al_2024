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

# Functions

def preprocess_reference(sc_ref: ad.AnnData, anno_key: str):
    # Construct reference profiles from categorical annotation
    logging.info(f"Construct reference profiles for annotation {anno_key}")
    sc_ref = tc.pp.construct_reference_profiles(
        sc_ref,
        annotation_key=anno_key,
        counts_location="X",
        inplace=True,
        normalize=True,
        target_sum=10000,
        trafo="log1p",
    )

    # Filter reference genes
    logging.info(f"Filter reference genes given annotation {anno_key}")
    sc_ref = tc.pp.filter_reference_genes(
        sc_ref,
        annotation_key=anno_key,
        min_log2foldchange=1,
        min_expression=1e-05,
        remove_mito=True,
        n_hvg=None,
        return_mask=False,
        return_view=True,
    )
    return(sc_ref)

def run_tacco_OT(spatial_adata: ad.AnnData, sc_ref: ad.AnnData, anno_key: str):

    # Preprocess annotation reference
    sc_ref = preprocess_reference(sc_ref, anno_key)
    
    # Annotate spatial_data using OT method
    result_key = f"tacco_OT_{anno_key}"
    logging.info("Annotate AnnData using OT method")
    spatial_adata = tc.tl.annotate(
        spatial_adata,
        sc_ref,
        method="OT",
        result_key=result_key,
        assume_valid_counts=True,
        annotation_key=anno_key,
        counts_location=("X"),
    )

    spatial_adata = tc.utils.get_maximum_annotation(
        spatial_adata,
        obsm_key=result_key,
        result_key=result_key)

    return(spatial_adata)

def run_tacco_RCTD(spatial_adata: ad.AnnData,
                   sc_ref: ad.AnnData,
                   anno_key: str,
                   conda_env: Path,
                   n_cores: int = 20):

    # Preprocess annotation reference
    sc_ref = preprocess_reference(sc_ref, anno_key)

    # Round counts (we do it here because it fails in the tacco code for me)
    spatial_adata.X.data = np.round(spatial_adata.X.data)
    sc_ref.X.data = np.round(sc_ref.X.data)
    sc_ref.X.data = sc_ref.X.data.astype(np.int32)

    # Annotate spatial_data using RCTD method
    result_key = f"tacco_RCTD_{anno_key}"
    logging.info("Annotate AnnData using RCTD method")

    # Tell tacco to run RCTD
    spatial_adata = tc.tl.annotate_RCTD(
        spatial_adata,
        sc_ref,
        annotation_key=anno_key,
        counts_location="X",
        conda_env=conda_env,
        x_coord_name="array_col",
        y_coord_name="array_row",
        doublet=True,
        n_cores=n_cores,
        verbose=True,
        working_directory="."
    )

    spatial_adata = tc.utils.get_maximum_annotation(
        spatial_adata,
        obsm_key=result_key,
        result_key=result_key)

    return(spatial_adata)


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

        # Read single cell reference data (unfiltered, raw counts)
        logging.info("Read single cell reference data")
        sc_ref_path = Path(params[sample]["single_cell_reference"]["h5ad_file"])
        sc_ref = ad.read_h5ad(sc_ref_path)

        # # Check count matrix for validity
        # logging.info(f"Check count matrix from {sc_ref_path} for validity")
        # tc.pp.check_counts_validity(sc_ref.X)

        # Preprocess single cell reference data as input for tacco.tools.annotate*
        logging.info("Preprocess single cell data (reference)")
        tc.utils.preprocess_single_cell_data(
            sc_ref,
            hvg=False,
            scale=False,
            pca=True,
            inplace=True,
            min_cells=10,
            min_genes=20,
            verbose=2,
        )

        for ad_name, ad_path in params[sample]["anndata_files"].items():
            #
            ad_path = Path(ad_path)

            # Load AnnData with cells/rows to be annotated/categorized
            logging.info(f"Load AnnData from {ad_path}")
            spatial_adata = ad.read_h5ad(ad_path)

            # # Check count matrix for validity
            # logging.info(f"Check count matrix from {ad_path} for validity")
            # tc.pp.check_counts_validity(spatial_adata.X)

            for method in params[sample]["tacco"]["method"]:
                logging.info(f"Apply annotation method {method} to {ad_path.name}")

                for anno_key in params[sample]["single_cell_reference"]["annotation_key"]:
                    logging.info("")
                    if method == "OT":
                        spatial_adata = run_tacco_OT(spatial_adata.copy(), sc_ref, anno_key)
                    elif method == "RCTD":
                        adata = spatial_adata.copy()
                        # We need to do the rounding here, because it does not work in tacco itself:
                        # - https://github.com/simonwm/tacco/blob/87252c6ffcf5f616ffbf167a8d4d34bb8d86d276/tacco/tools/_RCTD.py#L17
                        #   does not alter the values in adata.X.adata (I do not understand why)
                        
                        run_tacco_RCTD(
                            spatial_adata,
                            sc_ref,
                            anno_key,
                            params[sample]["conda_env"],
                            n_cores = 20)

                # Assemble output file name
                spatial_adata_path = output_dir.joinpath(f"tacco-annotate-{method}-{sample}-{ad_name}.h5ad")            

                logging.info(f"""Write annotated AnnData for method {method}
                                using annotaton "{anno_key}" to file 
                                {spatial_adata_path}""")

                spatial_adata.write_h5ad(spatial_adata_path)

if __name__ == "__main__":
    main()
