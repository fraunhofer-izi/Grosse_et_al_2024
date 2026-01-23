#!/usr/bin/env python3

import argparse
import gc
import json
import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from PIL import Image
import scanpy as sc
import skimage as skim
import sys
import tifffile
import tiffslide

# Initialize logger
logger = logging.getLogger(__name__)

## Color Deconvolution Functions

# Functions for stain matrices


def join_vertically(*args):
    """Joins many PIL images of the same dimensions vertically"""
    w, h = args[0].size
    n = len(args)
    joined = Image.new("RGB", (w, n * h))
    for y_off, img in zip(range(0, n * h, h), args):
        joined.paste(img, (0, y_off))
    return joined


def join_horizontally(*args):
    """Joins many PIL images of the same dimensions horizontally"""
    w, h = args[0].size
    n = len(args)
    joined = Image.new("RGB", (n * w, h))
    for x_off, img in zip(range(0, n * w, w), args):
        joined.paste(img, (x_off, 0))
    return joined


def create_matrix_image(matrix, inverse=False):
    """Vertical visualisation of stain-matrix.

    Args:
        matrix (list): List of length >= 3; first three entries must have values between 0 and 1
        inverse (bool, optional): If True the inverse of the color values are used. Defaults to False.
    """

    if inverse == False:
        pil_matrix = join_vertically(
            Image.new("RGB", (100, 100), tuple([int(255 * c) for c in matrix[0]])),
            Image.new("RGB", (100, 100), tuple([int(255 * c) for c in matrix[1]])),
            Image.new("RGB", (100, 100), tuple([int(255 * c) for c in matrix[2]])),
        )
    else:
        pil_matrix = join_vertically(
            Image.new(
                "RGB", (100, 100), tuple([int(255 * (1 - c)) for c in matrix[0]])
            ),
            Image.new(
                "RGB", (100, 100), tuple([int(255 * (1 - c)) for c in matrix[1]])
            ),
            Image.new(
                "RGB", (100, 100), tuple([int(255 * (1 - c)) for c in matrix[2]])
            ),
        )

    return pil_matrix


def display_stain_matrix(matrix, inverse=False):
    """Horizontal visualisation of stain-matrix.

    Args:
        matrix (list): List of length >= 3; first three entries must have values between 0 and 1
        inverse (bool, optional): If True the inverse of the color values are used. Defaults to False.
    """
    plt.axis("off")
    if inverse == False:
        plt.imshow(
            join_horizontally(
                Image.new("RGB", (100, 100), tuple([int(255 * c) for c in matrix[0]])),
                Image.new("RGB", (100, 100), tuple([int(255 * c) for c in matrix[1]])),
                Image.new("RGB", (100, 100), tuple([int(255 * c) for c in matrix[2]])),
            )
        )
    else:
        plt.imshow(
            join_horizontally(
                Image.new(
                    "RGB", (100, 100), tuple([int(255 * (1 - c)) for c in matrix[0]])
                ),
                Image.new(
                    "RGB", (100, 100), tuple([int(255 * (1 - c)) for c in matrix[1]])
                ),
                Image.new(
                    "RGB", (100, 100), tuple([int(255 * (1 - c)) for c in matrix[2]])
                ),
            )
        )


def create_hed_stains(
    img_array,
    stain_matrix,
    clipping_function=None,
    to_clip=[],
    clip_args={},
    contrast_fct=None,
    to_adjust=[],
    contrast_args=[],
):
    """Computes stain deconvolution of image in 'img_array' based on 'stain_matrix'.
    Usually used for HED color separation (H (hematoxylin), E (Eosin) and D (DAB)).
    Additionally

    Args:
        img_array (np.array): Image to deconvolute
        stain_matrix (np.array): Stain matrix to be used for stain deconvolution
        clipping_function(function, optional): Function used for clipping
        to_clip (list, optional): List of channels to clip (H, E or D). If the list is empty, nothing will be clipped.
        clip_args(dict, optional): Dictionary with arguments required for function given to clipping_function()
        contrast_fct (function, optional): Function used for contrast adjustment
        to_adjust (list, optional): List of channels to adjust contrast for
        contrast_args (list, optional): List of dictionaries, each dictionary contains arguments for each channel to adjust contrast

    Returns:
        {
        'my_h' (np.array): Hematoxylin stain, for visualisation
        'my_e' (np.array): Eosin stain, for visualisation
        'my_d' (np.array): DAB stain, for visualisation

        'my_hed' (np.array): Reassembled image in HED color space, e. g. to calculate histograms of pixel intensities
        'my_ihc' (np.array): Reassembled image in RGB color space.
        }

    """

    # stain matrix has to be normalised and inverted
    normalised_stain_matrix = stain_matrix
    normalised_stain_matrix = normalised_stain_matrix / np.reshape(
        np.sum(normalised_stain_matrix**2, axis=1) ** (1 / 2), (-1, 1)
    )
    # L2-normalisation of matrix
    # every value will be squared and row-wise summed up; square root of this sum == euclidian distance
    # reshape column vector
    # afterwards each vector (= rows of the matrix) is divided by its length === each vector has length of 1

    inverse_norm_stain_mat = np.linalg.inv(normalised_stain_matrix)
    # inverse of normalised stain matrix

    channel_axis = -1
    logger.info(f"Image array shape: {img_array.shape[channel_axis]}")
    if img_array.shape[channel_axis] != 3:
        msg = (
            f"the input array must have size 3 along `channel_axis`, "
            f"got {img_array.shape}"
        )
        raise ValueError(msg)
    my_hed = skim.color.separate_stains(img_array, inverse_norm_stain_mat)
    # single h,e,d stains; should be black, if displayed

    use_h = my_hed[:, :, 0]
    use_e = my_hed[:, :, 1]
    use_d = my_hed[:, :, 2]

    # Clipping:
    if "H" in to_clip:
        use_h = clipping_function(my_hed[:, :, 0], **clip_args)
        print("Clipped H with: " + str(clip_args))

    if "E" in to_clip:
        use_e = clipping_function(my_hed[:, :, 1], **clip_args)
        print("Clipped E with: " + str(clip_args))

    if "D" in to_clip:
        use_d = clipping_function(my_hed[:, :, 2], **clip_args)
        print("Clipped D with: " + str(clip_args))

    # Adjust contrast:
    if "H" in to_adjust:
        use_h = contrast_fct(use_h, **contrast_args[0])
        print("Adjusted H with: " + str(contrast_args[0]))
        # print("Adjusted H: ", use_h)

    if "E" in to_adjust:
        use_e = contrast_fct(use_e, **contrast_args[1])
        print("Adjusted E with: " + str(contrast_args[1]))
        # print("Adjusted E: ", use_e)

    if "D" in to_adjust:
        use_d = contrast_fct(use_d, **contrast_args[2])
        print("Adjusted D with: " + str(contrast_args[2]))

    # contrast adjustment is also possible on clipped channels!

    my_h = skim.color.combine_stains(
        np.stack(
            (
                use_h,
                np.zeros_like(my_hed[:, :, 0]),
                np.zeros_like(my_hed[:, :, 0]),
            ),
            axis=-1,
        ),
        normalised_stain_matrix,
    )

    my_e = skim.color.combine_stains(
        np.stack(
            (
                np.zeros_like(my_hed[:, :, 0]),
                use_e,
                np.zeros_like(my_hed[:, :, 0]),
            ),
            axis=-1,
        ),
        normalised_stain_matrix,
    )

    my_d = skim.color.combine_stains(
        np.stack(
            (
                np.zeros_like(my_hed[:, :, 0]),
                np.zeros_like(my_hed[:, :, 0]),
                use_d,
            ),
            axis=-1,
        ),
        normalised_stain_matrix,
    )

    # stack h, e, d channels together. Calculate histograms based on this
    my_hed = np.stack((use_h, use_e, use_d), axis=-1)

    # Back transformation of Rücktransformation von HED zu RGB
    my_ihc = skim.color.combine_stains(my_hed, normalised_stain_matrix)
    # = zusammengebautes bild aus h,e,d in rgb --> für Visualisierung

    return my_h, my_e, my_d, my_hed, my_ihc


def stain_deconvolution(img_array, stain_matrix) -> np.array:
    img_array_deconv = skim.color.separate_stains(img_array, stain_matrix)
    return img_array_deconv


def imwrite(path: Path, img_arr: np.array, **kwargs):
    """Function to wrap tifffile.imwrite

    Args:
        path (Path): Path to write image to
        img_arr (np.array): Numpy array with image data
    """
    # Test if array to save hast data type float
    if np.issubdtype(img_arr.dtype, np.floating):
        logger.info("Convert float to int array")
        img_arr = (img_arr * 255).astype(np.uint8)
    if not np.issubdtype(img_arr.dtype, np.integer):
        logger.error("'img_arr' has not integer data type")
        sys.exit(1)
    if img_arr.min() < 0:
        logger.error("'img_arr' has min value < 0")
        sys.exit(1)
    if img_arr.max() > 255:
        logger.error("'img_arr' has max value > 255")
        sys.exit(1)
    # Compression type "zstd" seems to work on HPC and
    # VS Code can display the resulting files
    tifffile.imwrite(path, img_arr, **kwargs)


def save_image_array(img_arr, output_dir, sample: str, matrix: str, stain: str):
    """Save numpy array to image file on disc.

    Args:
        img_arr (numpy.array): Numpy array with image data
        output_dir (str): Path to directory to write image to
        sample (str): Name of sample
        matrix (str): Name of used stain matrix
        stain (str): Name of stain
    """
    logger.info(f"img_arr for {stain} looks like:\nShape: {img_arr.shape}\n")
    out_path = Path(output_dir, sample)
    out_path.mkdir(parents=True, exist_ok=True)
    tiff_file = out_path.joinpath(f"{matrix}_{stain}.tiff")
    logger.info(f"Write {stain} TIFF file:\n{tiff_file}")
    imwrite(tiff_file, img_arr, photometric="rgb", compression="jpeg")

    # Scale image for display
    logger.info(f"Scale {stain} image")
    tiff_scaled_file = out_path.joinpath(f"{matrix}_{stain}_scale_10-percent.tiff")
    img_arr_rescaled = skim.transform.rescale(img_arr, 0.1, channel_axis=-1)
    logger.info(f"Write {stain} TIFF file:\n{tiff_scaled_file}")
    imwrite(tiff_scaled_file, img_arr_rescaled, photometric="rgb", compression="jpeg")


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

    # Run stain deconvolution for each sample in config
    for sample in params.keys():
        logging.info(f"Start 'stain deconvolution' on sample: {sample}")
        sd_params = params[sample]

        # Read original image (using package tiffslide)
        full_image_slide = tiffslide.TiffSlide(sd_params["microscope_image_path"])
        width, height = full_image_slide.dimensions[0], full_image_slide.dimensions[1]
        full_image_array = full_image_slide.read_region(
            (0, 0), 0, (width, height), as_array=True
        )
        logger.info(f"Full image looks like:\n{full_image_array}")

        # Perform stain deconvolution for each stain matrix in config
        for matrix in sd_params["stain_matrices"]:
            # img = stain_deconvolution(arr_crop1)
            stain_matrix = np.array(sd_params["stain_matrices"][matrix])
            logger.info(stain_matrix)
            logger.info(type(stain_matrix))
            hema_arr, eosin_arr, dab_arr, hed_arr, rgb_arr = create_hed_stains(
                full_image_array, stain_matrix
            )

            out_dir = sd_params["output_dir"]
            # save eosin image
            # save_image_array(eosin_arr, out_dir, sample=sample, matrix=matrix, stain="eosin")

            # save hematoxylin image
            save_image_array(
                hema_arr, out_dir, sample=sample, matrix=matrix, stain="hema"
            )

            # Clean up unused objects

            del hema_arr, eosin_arr, dab_arr, hed_arr, rgb_arr
            gc.collect()


if __name__ == "__main__":
    main()
