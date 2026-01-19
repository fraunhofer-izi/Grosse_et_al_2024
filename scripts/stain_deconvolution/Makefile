# Set SHELL to bash

SHELL:=/bin/bash

# Variables

tertiary_analysis = ../../data/output/
RSCRIPT=Rscript -e

# Use sbatch to run stain_deconvolution

sbatch:
	sbatch stain_deconvolution.sh

# Run stain_deconvolution workflow

stain-deconvolution: stain-deconvolution-HD_PM1 stain-deconvolution-HD_PM2 stain-deconvolution-HD_PM3 stain-deconvolution-HD_PM4 stain-deconvolution-HD_MM1 stain-deconvolution-HD_MM2 stain-deconvolution-HD_MM3

stain-deconvolution-HD_PM1: $(tertiary_analysis)stain_deconvolution/HD_PM1/*.tiff

stain-deconvolution-HD_PM2: $(tertiary_analysis)stain_deconvolution/HD_PM2/*.tiff

stain-deconvolution-HD_PM3: $(tertiary_analysis)stain_deconvolution/HD_PM3/*.tiff

stain-deconvolution-HD_PM4: $(tertiary_analysis)stain_deconvolution/HD_PM4/*.tiff

stain-deconvolution-HD_MM1: $(tertiary_analysis)stain_deconvolution/HD_MM1/*.tiff

stain-deconvolution-HD_MM2: $(tertiary_analysis)stain_deconvolution/HD_MM2/*.tiff

stain-deconvolution-HD_MM3: $(tertiary_analysis)stain_deconvolution/HD_MM3/*.tiff

$(tertiary_analysis)stain_deconvolution/HD_PM1/*.tiff: ../../configs/stain_deconvolution/HD_PM1.json ./stain_deconvolution.py
	conda activate Grosse_et_al_2026 ; \
	time ./stain_deconvolution.py --json $<

$(tertiary_analysis)stain_deconvolution/HD_PM2/*.tiff: ../../configs/stain_deconvolution/HD_PM2.json ./stain_deconvolution.py
	conda activate Grosse_et_al_2026 ; \
	time ./stain_deconvolution.py --json $<

$(tertiary_analysis)stain_deconvolution/HD_PM3/*.tiff: ../../configs/stain_deconvolution/HD_PM3.json ./stain_deconvolution.py
	conda activate Grosse_et_al_2026 ; \
	time ./stain_deconvolution.py --json $<

$(tertiary_analysis)stain_deconvolution/HD_PM4/*.tiff: ../../configs/stain_deconvolution/HD_PM4.json ./stain_deconvolution.py
	conda activate Grosse_et_al_2026 ; \
	time ./stain_deconvolution.py --json $<

$(tertiary_analysis)stain_deconvolution/HD_MM1/*.tiff: ../../configs/stain_deconvolution/HD_MM1.json ./stain_deconvolution.py
	conda activate Grosse_et_al_2026 ; \
	time ./stain_deconvolution.py --json $<

$(tertiary_analysis)stain_deconvolution/HD_MM2/*.tiff: ../../configs/stain_deconvolution/HD_MM2.json ./stain_deconvolution.py
	conda activate Grosse_et_al_2026 ; \
	time ./stain_deconvolution.py --json $<

$(tertiary_analysis)stain_deconvolution/HD_MM3/*.tiff: ../../configs/stain_deconvolution/HD_MM3.json ./stain_deconvolution.py
	conda activate Grosse_et_al_2026 ; \
	time ./stain_deconvolution.py --json $<

.PHONY: sbatch stain-deconvolution stain-deconvolution-HD_PM1 stain-deconvolution-HD_PM2 stain-deconvolution-HD_PM3 stain-deconvolution-HD_PM4 stain-deconvolution-HD_MM1 stain-deconvolution-HD_MM2 stain-deconvolution-HD_MM3