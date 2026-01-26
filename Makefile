PHONY: conda_env_create conda_env_export 00_download_data 00_generate_sc_reference 01_stain_deconvolution 02_bin2cell 03_tacco 04_rotate_sample

conda_env_create: environment.yml
	conda env create -f environment.yml

conda_env_export:
	conda env export > environment.yml

00_download_data: scripts/download_data/Makefile
	$(MAKE) -C scripts/download_data download_data

00_generate_sc_reference: scripts/generate_sc_reference/Makefile
	$(MAKE) -C scripts/generate_sc_reference generate-sc-reference

01_stain_deconvolution: scripts/stain_deconvolution/Makefile
	$(MAKE) -C scripts/stain_deconvolution stain-deconvolution

02_bin2cell: scripts/bin2cell/Makefile
	$(MAKE) -C scripts/bin2cell bin2cell

03_tacco: scripts/tacco/Makefile
	$(MAKE) -C scripts/tacco tacco

04_rotate_sample: scripts/rotate_sample/Makefile
	$(MAKE) -C scripts/rotate_sample rotate-sample
