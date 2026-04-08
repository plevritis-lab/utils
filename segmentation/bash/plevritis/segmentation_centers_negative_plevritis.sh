#!/bin/bash

SEGMENTATION_METHOD="cellpose"

PANEL_PATH="$HOME/Desktop/plevritis_analysis/data/channel_names.txt"
IMAGE_DIRECTORY="$HOME/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative"
QUANTIFICATION_SCRIPT="$HOME/Desktop/utils/segmentation/source/quantify_expression.py"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")
    SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.tif"
    
    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" ]]; then
        python3 "$SEGMENTATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --panel_path $PANEL_PATH \
            --nuclear_channel "DRAQ5" \
            --segment_channel "PDGFRB,CD3,CYTOKERATIN" \
            --apply_$SEGMENTATION_METHOD
    fi
done