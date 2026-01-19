# Set SHELL to bash

SHELL:=/bin/bash

# Variables

RSCRIPT=Rscript -e

download_data: download_single_cell_data download_spatial_seq_data

download_single_cell_data:
	./download_single_cell_data.sh

download_spatial_seq_data:
	./download_spatial_seq_data.sh