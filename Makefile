PHONY:

conda_env_create: environment.yml
	conda env create -f environment.yml

conda_env_export:
	conda env export > environment.yml

00_download_data: scripts/download_data/Makefile
	$(MAKE) -C scripts/download_data download_data

01_stain_deconvolution: scripts/stain_deconvolution/Makefile
	$(MAKE) -C scripts/stain_deconvolution stain-deconvolution

02_bin2cell: scripts/bin2cell/Makefile
	$(MAKE) -C scripts/bin2cell/