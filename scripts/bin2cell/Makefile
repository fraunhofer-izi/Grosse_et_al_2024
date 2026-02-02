# Set SHELL to bash

SHELL:=/bin/bash

# Variables

tertiary_analysis = ../../data/output/
RSCRIPT=Rscript -e

# Run bin2cell workflow

bin2cell: bin2cell-HD_PM1 bin2cell-HD_PM2 bin2cell-HD_PM3 bin2cell-HD_PM4 bin2cell-HD_MM1 bin2cell-HD_MM2 bin2cell-HD_MM3

# bin2cell using images manipulated using stain-deconvolution.py

bin2cell-HD_PM1: $(tertiary_analysis)bin2cell/HD_PM1/h5-files/adata-post-b2c.h5ad

bin2cell-HD_PM2: $(tertiary_analysis)bin2cell/HD_PM2/h5-files/adata-post-b2c.h5ad

bin2cell-HD_PM3: $(tertiary_analysis)bin2cell/HD_PM3/h5-files/adata-post-b2c.h5ad

bin2cell-HD_PM4: $(tertiary_analysis)bin2cell/HD_PM4/h5-files/adata-post-b2c.h5ad

bin2cell-HD_MM1: $(tertiary_analysis)bin2cell/HD_MM1/h5-files/adata-post-b2c.h5ad

bin2cell-HD_MM2: $(tertiary_analysis)bin2cell/HD_MM2/h5-files/adata-post-b2c.h5ad

bin2cell-HD_MM3: $(tertiary_analysis)bin2cell/HD_MM3/h5-files/adata-post-b2c.h5ad

$(tertiary_analysis)bin2cell/HD_PM1/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/HD_PM1.json $(tertiary_analysis)stain_deconvolution/HD_PM1/scikitimage_rgb_from_hed_hema.tiff
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/HD_PM2/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/HD_PM2.json $(tertiary_analysis)stain_deconvolution/HD_PM2/scikitimage_rgb_from_hed_hema.tiff
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/HD_PM3/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/HD_PM3.json $(tertiary_analysis)stain_deconvolution/HD_PM3/scikitimage_rgb_from_hed_hema.tiff
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/HD_PM4/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/HD_PM4.json $(tertiary_analysis)stain_deconvolution/HD_PM4/scikitimage_rgb_from_hed_hema.tiff
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/HD_MM1/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/HD_MM1.json $(tertiary_analysis)stain_deconvolution/HD_MM1/scikitimage_rgb_from_hed_hema.tiff
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/HD_MM2/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/HD_MM2.json $(tertiary_analysis)stain_deconvolution/HD_MM2/scikitimage_rgb_from_hed_hema.tiff
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/HD_MM3/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/HD_MM3.json $(tertiary_analysis)stain_deconvolution/HD_MM3/scikitimage_rgb_from_hed_hema.tiff
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

# bin2cell using the original microscopy images

.PHONY: bin2cell-HD_PM1 bin2cell-HD_PM2 bin2cell-HD_PM3 bin2cell-HD_PM4 bin2cell-HD_MM1  bin2cell-HD_MM2 bin2cell-HD_MM3