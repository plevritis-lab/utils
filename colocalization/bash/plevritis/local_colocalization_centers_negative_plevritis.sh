#!/bin/bash

CONDITION="primary_tumors/centers/node_negative"

DATA_DIRECTORY="$HOME/Desktop/plevritis_analysis/data"
CELESTA_DIRECTORY="$HOME/Desktop/plevritis_analysis/celesta"
SCRIPT_DIRECTORY="$HOME/Desktop/EPOCHS/colocalization/source/local_colocalization"

COLORMAP="$HOME/Desktop/EPOCHS/celesta/source/colormaps/plevritis/colormap_primary_plevritis.json"

SIGNATURE_MATRIX="$CELESTA_DIRECTORY/signature_matrices/$CONDITION/center_negative_signature_matrix.csv"

IMAGE_DIRECTORY="$DATA_DIRECTORY/$CONDITION"
ASSIGNMENTS_DIRECTORY="$IMAGE_DIRECTORY/assignments"

LOCAL_COLOCALIZATION_SCRIPT="$SCRIPT_DIRECTORY/generate_spatial_embeddings.py"

python3 "$LOCAL_COLOCALIZATION_SCRIPT" \
    --colormap_path "$COLORMAP" \
    --data_directory "$ASSIGNMENTS_DIRECTORY" \
    --save_path "$IMAGE_DIRECTORY" \
    --signature_matrix "$SIGNATURE_MATRIX"