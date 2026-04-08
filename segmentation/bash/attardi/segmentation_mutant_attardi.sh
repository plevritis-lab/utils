#!/bin/bash

SEGMENTATION_METHOD="cellpose"

IMAGE_DIRECTORY="$HOME/Desktop/attardi_analysis/data/mutant"
QUANTIFICATION_SCRIPT="$HOME/Desktop/utils/segmentation/source/quantify_expression.py"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")
    SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.qptiff"
    
    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" ]]; then
        python3 "$SEGMENTATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --panel_path $PANEL_PATH \
            --nuclear_channel "NUCLEAR_ONE" \
            --segment_channel "E_CADHERIN" \
            --apply_$SEGMENTATION_METHOD
    fi
done