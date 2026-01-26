# Set SHELL to bash

SHELL:=/bin/bash

# Variables

tertiary_analysis = ../../data/output/

# Run tacco annotate_OT

rotate-sample: ../../configs/rotate_sample/rotate_tacco_OT_hema.json
	conda activate Grosse_et_al_2026 ; \
	time ./rotate-sample.py --json $<

.PHONY: rotate-sample
