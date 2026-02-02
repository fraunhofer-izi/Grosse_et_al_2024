# Set SHELL to bash

SHELL:=/bin/bash

# Variables

tertiary_analysis = ../../data/output/
sc_reference = ../../data/output/single_cell_reference/Stubenvoll_2025_melanoma_single_cell_reference.h5ad

# Run tacco annotate_OT

tacco-annotate-OT: tacco-annotate-OT-HD_PM1 tacco-annotate-OT-HD_PM2 tacco-annotate-OT-HD_PM3 tacco-annotate-OT-HD_PM4 tacco-annotate-OT-HD_MM1 tacco-annotate-OT-HD_MM2 tacco-annotate-OT-HD_MM3

tacco-annotate-OT-HD_PM1: $(tertiary_analysis)tacco/HD_PM1/tacco-annotate-OT-HD_PM1-hema.h5ad
tacco-annotate-OT-HD_PM2: $(tertiary_analysis)tacco/HD_PM2/tacco-annotate-OT-HD_PM2-hema.h5ad
tacco-annotate-OT-HD_PM3: $(tertiary_analysis)tacco/HD_PM3/tacco-annotate-OT-HD_PM3-hema.h5ad
tacco-annotate-OT-HD_PM4: $(tertiary_analysis)tacco/HD_PM4/tacco-annotate-OT-HD_PM4-hema.h5ad
tacco-annotate-OT-HD_MM1: $(tertiary_analysis)tacco/HD_MM1/tacco-annotate-OT-HD_MM1-hema.h5ad
tacco-annotate-OT-HD_MM2: $(tertiary_analysis)tacco/HD_MM2/tacco-annotate-OT-HD_MM2-hema.h5ad
tacco-annotate-OT-HD_MM3: $(tertiary_analysis)tacco/HD_MM3/tacco-annotate-OT-HD_MM3-hema.h5ad

$(tertiary_analysis)tacco/HD_PM1/tacco-annotate-OT-HD_PM1-hema.h5ad: ../../configs/tacco/HD_PM1-OT.json $(sc_reference)
	conda activate Grosse_et_al_2026 ; \
	time ./tacco-annotate.py --json $<

$(tertiary_analysis)tacco/HD_PM2/tacco-annotate-OT-HD_PM2-hema.h5ad: ../../configs/tacco/HD_PM2-OT.json $(sc_reference)
	conda activate Grosse_et_al_2026 ; \
	time ./tacco-annotate.py --json $<

$(tertiary_analysis)tacco/HD_PM3/tacco-annotate-OT-HD_PM3-hema.h5ad: ../../configs/tacco/HD_PM3-OT.json $(sc_reference)
	conda activate Grosse_et_al_2026 ; \
	time ./tacco-annotate.py --json $<

$(tertiary_analysis)tacco/HD_PM4/tacco-annotate-OT-HD_PM4-hema.h5ad: ../../configs/tacco/HD_PM4-OT.json $(sc_reference)
	conda activate Grosse_et_al_2026 ; \
	time ./tacco-annotate.py --json $<

$(tertiary_analysis)tacco/HD_MM1/tacco-annotate-OT-HD_MM1-hema.h5ad: ../../configs/tacco/HD_MM1-OT.json $(sc_reference)
	conda activate Grosse_et_al_2026 ; \
	time ./tacco-annotate.py --json $<

$(tertiary_analysis)tacco/HD_MM2/tacco-annotate-OT-HD_MM2-hema.h5ad: ../../configs/tacco/HD_MM2-OT.json $(sc_reference)
	conda activate Grosse_et_al_2026 ; \
	time ./tacco-annotate.py --json $<

$(tertiary_analysis)tacco/HD_MM3/tacco-annotate-OT-HD_MM3-hema.h5ad: ../../configs/tacco/HD_MM3-OT.json $(sc_reference)
	conda activate Grosse_et_al_2026 ; \
	time ./tacco-annotate.py --json $<
.PHONY: tacco-annotate-OT tacco-annotate-OT-HD_PM1 tacco-annotate-OT-HD_PM2 tacco-annotate-OT-HD_PM3 tacco-annotate-OT-HD_PM4 tacco-annotate-OT-HD_MM1 tacco-annotate-OT-HD_MM2 tacco-annotate-OT-HD_MM3
