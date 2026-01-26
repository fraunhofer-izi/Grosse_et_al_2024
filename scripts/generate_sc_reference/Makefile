# Set SHELL to bash

SHELL:=/bin/bash

# Variables

tertiary_analysis = ../../data/output/

# Run tacco annotate_OT

generate-sc-reference: ../../configs/single_cell_reference.json
	conda activate Grosse_et_al_2026 ; \
	time ./generate_sc_reference.py --json $<

.PHONY: generate-sc-reference
