#!/bin/bash

CONDITION="primary_tumors/centers/node_positive"

DATA_DIRECTORY="$HOME/Desktop/plevritis_analysis/data"
CELESTA_DIRECTORY="$HOME/Desktop/plevritis_analysis/celesta"
SCRIPT_DIRECTORY="$HOME/Desktop/EPOCHS/topology/source"

COLORMAP="$HOME/Desktop/EPOCHS/celesta/source/colormaps/plevritis/colormap_primary_plevritis.json"

SIGNATURE_MATRIX="$CELESTA_DIRECTORY/signature_matrices/$CONDITION/center_positive_signature_matrix.csv"

IMAGE_DIRECTORY="$DATA_DIRECTORY/$CONDITION"
ASSIGNMENTS_DIRECTORY="$IMAGE_DIRECTORY/assignments"

TOPOLOGY_SCRIPT="$SCRIPT_DIRECTORY/apply_persistent_homology.py"

python3 "$TOPOLOGY_SCRIPT" \
    --data_directory "$ASSIGNMENTS_DIRECTORY" \
    --save_path "$IMAGE_DIRECTORY" \
    --signature_matrix "$SIGNATURE_MATRIX"