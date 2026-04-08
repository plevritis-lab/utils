#!/bin/bash

SEGMENTATION_METHOD="mesmer"

PANEL_PATH="$HOME/Desktop/plevritis_analysis/data/channel_names.txt"
IMAGE_DIRECTORY="$HOME/Desktop/plevritis_analysis/data/primary_tumors/centers/node_negative"
QUANTIFICATION_SCRIPT="$HOME/Desktop/utils/segmentation/source/quantify_expression.py"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")
    SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.tif"
    SAMPLE_SEGMENTATION="$SAMPLE_DIRECTORY/full/$SEGMENTATION_METHOD/PDGFRB + CD45 + CYTOKERATIN/image_1_seg.npy"

    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" && "$SAMPLE_NAME" != "assignments" ]]; then
        python3 "$QUANTIFICATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --mask_path "$SAMPLE_SEGMENTATION" \
            --panel_path "$PANEL_PATH" \
            --save_path "$IMAGE_DIRECTORY" \
            --apply_$SEGMENTATION_METHOD
    fi
done