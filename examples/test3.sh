#!/usr/bin/env bash
panopticon multireference-dna-correspondence --output mdc.pdf --seg Seg2 --seg Seg4 --query "cell_type == tumor" --query "patient_ID == 43" example_data/example.loom
