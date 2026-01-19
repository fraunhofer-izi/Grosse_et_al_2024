# Set SHELL to bash

SHELL:=/bin/bash

# Variables

tertiary_analysis = ../../data/output/
RSCRIPT=Rscript -e

# Use sbatch to run anndata2csv4rctd

sbatch: sbatch-eosin sbatch-hema sbatch-nuclei sbatch-original

sbatch-eosin:
	sbatch bin2cell-workflow-eosin.sh

sbatch-hema:
	sbatch bin2cell-workflow-hema.sh

sbatch-nuclei:
	sbatch bin2cell-workflow-nuclei.sh

sbatch-original:
	sbatch bin2cell-workflow-original.sh

# Run bin2cell workflow

bin2cell-workflow-nuclei: bin2cell-LE-589-BE-nuclei bin2cell-LE-595-SV-nuclei bin2cell-LE-597-EG-nuclei

# bin2cell using images manipulated using imagej
# The following target failed
#bin2cell-LE-585-HM-nuclei: $(tertiary_analysis)bin2cell/LE-585-HM-nuclei/h5-files/adata-post-b2c.h5ad

bin2cell-LE-589-BE-nuclei: $(tertiary_analysis)bin2cell/LE-589-BE-nuclei/h5-files/adata-post-b2c.h5ad

bin2cell-LE-595-SV-nuclei: $(tertiary_analysis)bin2cell/LE-595-SV-nuclei/h5-files/adata-post-b2c.h5ad

bin2cell-LE-597-EG-nuclei: $(tertiary_analysis)bin2cell/LE-597-EG-nuclei/h5-files/adata-post-b2c.h5ad

$(tertiary_analysis)bin2cell/LE-585-HM-nuclei/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-585-HM-nuclei.json
	conda activate Grosse_et_al_2026 ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-589-BE-nuclei/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-589-BE-nuclei.json
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-595-SV-nuclei/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-595-SV-nuclei.json
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-597-EG-nuclei/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-597-EG-nuclei.json
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

# bin2cell using images manipulated using stain-deconvolution.py
bin2cell-workflow-eosin: bin2cell-LE-585-HM-eosin bin2cell-LE-589-BE-eosin bin2cell-LE-595-SV-eosin bin2cell-LE-597-EG-eosin bin2cell-LE-713-SI-eosin bin2cell-LE-773-SK-eosin bin2cell-LE-781-KHJ-eosin

bin2cell-LE-585-HM-eosin: $(tertiary_analysis)bin2cell/LE-585-HM-eosin/h5-files/adata-post-b2c.h5ad

bin2cell-LE-589-BE-eosin: $(tertiary_analysis)bin2cell/LE-589-BE-eosin/h5-files/adata-post-b2c.h5ad

bin2cell-LE-595-SV-eosin: $(tertiary_analysis)bin2cell/LE-595-SV-eosin/h5-files/adata-post-b2c.h5ad

bin2cell-LE-597-EG-eosin: $(tertiary_analysis)bin2cell/LE-597-EG-eosin/h5-files/adata-post-b2c.h5ad

# The following target failed (bin2cell-983-1-eosin)
#bin2cell-983-1-eosin: $(tertiary_analysis)bin2cell/983-1-eosin/h5-files/adata-post-b2c.h5ad

bin2cell-LE-713-SI-eosin: $(tertiary_analysis)bin2cell/LE-713-SI-eosin/h5-files/adata-post-b2c.h5ad

bin2cell-LE-773-SK-eosin: $(tertiary_analysis)bin2cell/LE-773-SK-eosin/h5-files/adata-post-b2c.h5ad

bin2cell-LE-781-KHJ-eosin: $(tertiary_analysis)bin2cell/LE-781-KHJ-eosin/h5-files/adata-post-b2c.h5ad

$(tertiary_analysis)bin2cell/LE-585-HM-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-585-HM-eosin.json $(tertiary_analysis)stain_deconvolution/LE-585-HM/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-589-BE-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-589-BE-eosin.json $(tertiary_analysis)stain_deconvolution/LE-589-BE/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-595-SV-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-595-SV-eosin.json $(tertiary_analysis)stain_deconvolution/LE-595-SV/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-597-EG-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-597-EG-eosin.json $(tertiary_analysis)stain_deconvolution/LE-597-EG/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/983-1-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/983-1-eosin.json $(tertiary_analysis)stain_deconvolution/983-1/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-713-SI-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-713-SI-eosin.json $(tertiary_analysis)stain_deconvolution/LE-713-SI/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-773-SK-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-773-SK-eosin.json $(tertiary_analysis)stain_deconvolution/LE-773-SK/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-781-KHJ-eosin/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-781-KHJ-eosin.json $(tertiary_analysis)stain_deconvolution/LE-781-KHJ/scikitimage_rgb_from_hed_eosin.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<


bin2cell-workflow-hema: bin2cell-LE-585-HM-hema bin2cell-LE-589-BE-hema bin2cell-LE-595-SV-hema bin2cell-LE-597-EG-hema bin2cell-983-1-hema bin2cell-LE-713-SI-hema bin2cell-LE-773-SK-hema bin2cell-LE-781-KHJ-hema

bin2cell-LE-585-HM-hema: $(tertiary_analysis)bin2cell/LE-585-HM-hema/h5-files/adata-post-b2c.h5ad

bin2cell-LE-589-BE-hema: $(tertiary_analysis)bin2cell/LE-589-BE-hema/h5-files/adata-post-b2c.h5ad

bin2cell-LE-595-SV-hema: $(tertiary_analysis)bin2cell/LE-595-SV-hema/h5-files/adata-post-b2c.h5ad

bin2cell-LE-597-EG-hema: $(tertiary_analysis)bin2cell/LE-597-EG-hema/h5-files/adata-post-b2c.h5ad

bin2cell-983-1-hema: $(tertiary_analysis)bin2cell/983-1-hema/h5-files/adata-post-b2c.h5ad

bin2cell-LE-713-SI-hema: $(tertiary_analysis)bin2cell/LE-713-SI-hema/h5-files/adata-post-b2c.h5ad

bin2cell-LE-773-SK-hema: $(tertiary_analysis)bin2cell/LE-773-SK-hema/h5-files/adata-post-b2c.h5ad

bin2cell-LE-781-KHJ-hema: $(tertiary_analysis)bin2cell/LE-781-KHJ-hema/h5-files/adata-post-b2c.h5ad

$(tertiary_analysis)bin2cell/LE-585-HM-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-585-HM-hema.json $(tertiary_analysis)stain_deconvolution/LE-585-HM/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-589-BE-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-589-BE-hema.json $(tertiary_analysis)stain_deconvolution/LE-589-BE/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-595-SV-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-595-SV-hema.json $(tertiary_analysis)stain_deconvolution/LE-595-SV/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-597-EG-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-597-EG-hema.json $(tertiary_analysis)stain_deconvolution/LE-597-EG/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/983-1-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/983-1-hema.json $(tertiary_analysis)stain_deconvolution/983-1/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-713-SI-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-713-SI-hema.json $(tertiary_analysis)stain_deconvolution/LE-713-SI/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-773-SK-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-773-SK-hema.json $(tertiary_analysis)stain_deconvolution/LE-773-SK/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-781-KHJ-hema/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-781-KHJ-hema.json $(tertiary_analysis)stain_deconvolution/LE-781-KHJ/scikitimage_rgb_from_hed_hema.tiff
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

# bin2cell using the original microscopy images
bin2cell-workflow-original: bin2cell-LE-585-HM-original bin2cell-LE-589-BE-original bin2cell-LE-595-SV-original bin2cell-LE-597-EG-original bin2cell-983-1-original bin2cell-LE-713-SI-original bin2cell-LE-773-SK-original bin2cell-LE-781-KHJ-original

bin2cell-LE-585-HM-original: $(tertiary_analysis)bin2cell/LE-585-HM-original/h5-files/adata-post-b2c.h5ad

bin2cell-LE-589-BE-original: $(tertiary_analysis)bin2cell/LE-589-BE-original/h5-files/adata-post-b2c.h5ad

bin2cell-LE-595-SV-original: $(tertiary_analysis)bin2cell/LE-595-SV-original/h5-files/adata-post-b2c.h5ad

bin2cell-LE-597-EG-original: $(tertiary_analysis)bin2cell/LE-597-EG-original/h5-files/adata-post-b2c.h5ad

bin2cell-983-1-original: $(tertiary_analysis)bin2cell/983-1-original/h5-files/adata-post-b2c.h5ad

bin2cell-LE-713-SI-original: $(tertiary_analysis)bin2cell/LE-713-SI-original/h5-files/adata-post-b2c.h5ad

bin2cell-LE-773-SK-original: $(tertiary_analysis)bin2cell/LE-773-SK-original/h5-files/adata-post-b2c.h5ad

bin2cell-LE-781-KHJ-original: $(tertiary_analysis)bin2cell/LE-781-KHJ-original/h5-files/adata-post-b2c.h5ad

$(tertiary_analysis)bin2cell/LE-585-HM-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-585-HM-original.json /mnt/ribolution/secondary_analysis/2024_SpatSeq_VisiumHD/240726_VH00421_54_2223G3WNX/microscope_images/Hautklinik_Kunz/20240705_A1-538_4.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-589-BE-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-589-BE-original.json /mnt/ribolution/secondary_analysis/2024_SpatSeq_VisiumHD/240726_VH00421_54_2223G3WNX/microscope_images/Hautklinik_Kunz/20240704_D1_846-1.6.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-595-SV-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-595-SV-original.json /mnt/ribolution/secondary_analysis/2024_SpatSeq_VisiumHD/240826_VH00421_59_AACGNK5HV/microscope_images/Hautklinik_Kunz/20240704_D1_941-4.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-597-EG-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-597-EG-original.json /mnt/ribolution/secondary_analysis/2024_SpatSeq_VisiumHD/240826_VH00421_59_AACGNK5HV/microscope_images/Hautklinik_Kunz/20240704_A1_942-4.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/983-1-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/983-1-original.json /mnt/ribolution/secondary_analysis/2025_SpatSeq_VisiumHD/250627_VH00421_91_222CCLVNX/microscope_images/983-1-25_20x.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-713-SI-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-713-SI-original.json /mnt/ribolution/secondary_analysis/2025_SpatSeq_VisiumHD/250627_VH00421_91_222CCLVNX/microscope_images/LE-713-SI_croop.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-773-SK-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-773-SK-original.json /mnt/ribolution/secondary_analysis/2025_SpatSeq_VisiumHD/250627_VH00421_91_222CCLVNX/microscope_images/LE-773-SK_20x.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<

$(tertiary_analysis)bin2cell/LE-781-KHJ-original/h5-files/adata-post-b2c.h5ad: ../../configs/bin2cell/LE-781-KHJ-original.json /mnt/ribolution/secondary_analysis/2025_SpatSeq_VisiumHD/250627_VH00421_91_222CCLVNX/microscope_images/LE-781-KHJ_20x.tif
	source ../../virtualenv/bin/activate ; \
	source /etc/profile.d/modules.sh ; \
	module load CUDA ; \
	module load cuDNN ; \
	time ./bin2cell-workflow.py --json $<


.PHONY: sbatch-eosin sbatch-nuclei sbatch-original bin2cell-workflow-nuclei bin2cell-LE-585-HM-nuclei bin2cell-LE-589-BE-nuclei bin2cell-LE-595-SV-nuclei bin2cell-LE-597-EG-nuclei bin2cell-workflow-eosin bin2cell-LE-585-HM-eosin bin2cell-LE-589-BE-eosin bin2cell-LE-595-SV-eosin bin2cell-LE-597-EG-eosin bin2cell-workflow-hema bin2cell-LE-585-HM-hema bin2cell-LE-589-BE-hema bin2cell-LE-595-SV-hema bin2cell-LE-597-EG-hema bin2cell-workflow-original bin2cell-LE-585-HM-original bin2cell-LE-589-BE-original bin2cell-LE-595-SV-original bin2cell-LE-597-EG-original