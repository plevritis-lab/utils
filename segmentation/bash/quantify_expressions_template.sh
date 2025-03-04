#!/bin/bash

PANEL_PATH="<TODO>"
IMAGE_DIRECTORY="<TODO>"
QUANTIFICATION_SCRIPT="<TODO>"

for SAMPLE_DIRECTORY in "$IMAGE_DIRECTORY"/*; do
    SAMPLE_NAME=$(basename "$SAMPLE_DIRECTORY")

    if [[ "$SAMPLE_NAME" != "histology" && "$SAMPLE_NAME" != "quantifications" ]]; then
        SAMPLE_DATA="$SAMPLE_DIRECTORY/data/${SAMPLE_NAME}.tif"
        SAMPLE_SEGMENTATION="$SAMPLE_DIRECTORY/full/<TODO>/<TODO>/image_1_seg.npy"
        
        python3 "$QUANTIFICATION_SCRIPT" \
            --image_path "$SAMPLE_DATA" \
            --mask_path "$SAMPLE_SEGMENTATION" \
            --panel_path "$PANEL_PATH" \
            --save_path "$IMAGE_DIRECTORY" \
            --<TODO>
    fi
done